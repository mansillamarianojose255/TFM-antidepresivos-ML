# ==============================================================================
# TFM — PIPELINE REPRODUCIBLE
# Análisis predictivo del uso de antidepresivos mediante modelos de Machine
# Learning utilizando datos genéticos públicos de estudios GWAS
#
# Autor: Mariano José Mansilla
#
# Este script está ordenado para seguir la secuencia de Materiales y Métodos:
# 1) procesamiento del GWAS;
# 2) exploración gráfica;
# 3) distribución cromosómica;
# 4) LD/haplotipos (paso externo documentado);
# 5) anotación funcional y eQTL;
# 6) construcción del pool genotípico;
# 7) construcción de E1, E2 y E3;
# 8) preparación de datos para Machine Learning;
# 9) modelos;
# 10) métricas, tuning, evaluación final, comparaciones e importancia.
#
# IMPORTANTE
# - El VCF original NO se incluye en GitHub por su tamaño.
# - Colocar finn-b-ANTIDEPRESSANTS.vcf.gz en datos/raw/.
# - El análisis de LD/LDhap se realizó en la web de LDlink y se documenta aquí,
#   pero no se automatiza porque no formó parte del código R original.
# ==============================================================================

options(stringsAsFactors = FALSE)
options(warn = 1)

# ==============================================================================
# 0. CONFIGURACIÓN GENERAL
# ==============================================================================

CONFIG <- list(
  archivo_vcf = file.path("datos", "raw", "finn-b-ANTIDEPRESSANTS.vcf.gz"),
  umbral_lp_inicial = 4,
  umbral_genomico = -log10(5e-8),
  tam_chunk = 250000L,

  n_pool = 10000L,
  n_mujeres_pool = 5000L,
  n_varones_pool = 5000L,

  n_casos = 1000L,
  n_controles = 1000L,
  n_casos_por_sexo = 500L,
  n_controles_por_sexo = 500L,

  lambda = log(1.5),

  prop_train = 0.80,
  cv_folds = 5L,
  cv_repeticiones = 5L,

  semilla_pool = 123L,
  semilla_e1 = 1101L,
  semilla_e2 = 1202L,
  semilla_e3 = 1303L,
  semilla_split_base = 20260816L,
  semilla_cv_base = 20260817L,

  clase_negativa = "control",
  clase_positiva = "caso"
)

DIRS <- c(
  "datos",
  file.path("datos", "raw"),
  file.path("datos", "ldlink"),
  "resultados",
  "figuras",
  "modelos"
)

invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

# ------------------------------------------------------------------------------
# Paquetes
# ------------------------------------------------------------------------------

PAQUETES_CRAN <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "ggplot2",
  "ggrepel",
  "scales",
  "purrr",
  "tibble",
  "stringr",
  "httr2",
  "jsonlite",
  "gtexr",
  "pROC",
  "PRROC",
  "ranger",
  "xgboost",
  "lightgbm",
  "nnet"
)

PAQUETES_BIOC <- c(
  "GenomicRanges",
  "IRanges",
  "S4Vectors",
  "karyoploteR"
)

faltan_cran <- PAQUETES_CRAN[
  !vapply(PAQUETES_CRAN, requireNamespace, logical(1), quietly = TRUE)
]

faltan_bioc <- PAQUETES_BIOC[
  !vapply(PAQUETES_BIOC, requireNamespace, logical(1), quietly = TRUE)
]

if (length(faltan_cran) > 0L || length(faltan_bioc) > 0L) {
  stop(
    paste0(
      "Faltan paquetes.\n\nCRAN: ",
      ifelse(length(faltan_cran) == 0L, "ninguno", paste(faltan_cran, collapse = ", ")),
      "\nBioconductor: ",
      ifelse(length(faltan_bioc) == 0L, "ninguno", paste(faltan_bioc, collapse = ", ")),
      "\n\nConsulte README.md para la instalación."
    )
  )
}

# ------------------------------------------------------------------------------
# Funciones utilitarias generales
# ------------------------------------------------------------------------------

normalizar_chr <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y <- gsub("^CHR", "", y)
  y[y == "23"] <- "X"
  y
}

es_chr_x <- function(x) {
  normalizar_chr(x) == "X"
}

guardar_csv <- function(x, nombre, carpeta = "resultados") {
  utils::write.csv(
    x,
    file = file.path(carpeta, nombre),
    row.names = FALSE,
    na = ""
  )
}

clip_prob <- function(x, eps = 1e-15) {
  pmin(pmax(as.numeric(x), eps), 1 - eps)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}


# ==============================================================================
# 1. DESCARGA, EXTRACCIÓN Y PROCESAMIENTO DE LAS ESTADÍSTICAS RESUMIDAS
# ==============================================================================

# El archivo VCF se descargó previamente desde IEU OpenGWAS.
# El script parte del archivo local:
#
# datos/raw/finn-b-ANTIDEPRESSANTS.vcf.gz

if (!file.exists(CONFIG$archivo_vcf)) {
  stop(
    paste0(
      "No se encontró el VCF original:\n",
      CONFIG$archivo_vcf,
      "\n\nColóquelo en datos/raw/ antes de ejecutar el pipeline."
    )
  )
}

# ------------------------------------------------------------------------------
# 1.1. Lectura de la cabecera del VCF
# ------------------------------------------------------------------------------

leer_cabecera_vcf <- function(archivo_vcf) {
  con <- gzfile(archivo_vcf, open = "rt")
  on.exit(close(con), add = TRUE)

  cabecera <- character()

  repeat {
    linea <- readLines(con, n = 1L, warn = FALSE)

    if (length(linea) == 0L) {
      break
    }

    if (!startsWith(linea, "#")) {
      break
    }

    cabecera <- c(cabecera, linea)

    if (startsWith(linea, "#CHROM")) {
      break
    }
  }

  cabecera
}

cabecera_vcf <- leer_cabecera_vcf(CONFIG$archivo_vcf)

linea_columnas <- grep("^#CHROM", cabecera_vcf, value = TRUE)

if (length(linea_columnas) != 1L) {
  stop("No se pudo identificar de forma inequívoca la línea #CHROM del VCF.")
}

linea_lp <- grep("^##FORMAT=<ID=LP,", cabecera_vcf, value = TRUE)

if (length(linea_lp) == 0L) {
  warning("No se encontró la definición de LP en la cabecera del VCF.")
}

# Extraer longitudes cromosómicas declaradas en el VCF.
lineas_contig <- grep("^##contig=", cabecera_vcf, value = TRUE)

extraer_contig <- function(linea) {
  id <- sub(".*ID=([^,>]+).*", "\\1", linea)
  longitud <- suppressWarnings(
    as.numeric(sub(".*length=([0-9]+).*", "\\1", linea))
  )

  data.frame(
    CHR = normalizar_chr(id),
    LONGITUD = longitud,
    stringsAsFactors = FALSE
  )
}

tabla_contigs <- dplyr::bind_rows(lapply(lineas_contig, extraer_contig)) |>
  dplyr::filter(!is.na(.data$LONGITUD)) |>
  dplyr::distinct(.data$CHR, .keep_all = TRUE)

guardar_csv(tabla_contigs, "01_longitudes_cromosomicas_GRCh37.csv", "datos")

# ------------------------------------------------------------------------------
# 1.2. Filtrado por LP > 4 mediante lectura por fragmentos
# ------------------------------------------------------------------------------

filtrar_vcf_lp <- function(
    archivo_entrada,
    umbral = 4,
    tam_chunk = 250000L
) {

  con <- gzfile(archivo_entrada, open = "rt")
  on.exit(close(con), add = TRUE)

  # Avanzar hasta #CHROM.
  repeat {
    linea <- readLines(con, n = 1L, warn = FALSE)

    if (length(linea) == 0L) {
      stop("Se alcanzó el final del VCF antes de encontrar #CHROM.")
    }

    if (startsWith(linea, "#CHROM")) {
      break
    }
  }

  salida <- list()
  k <- 1L
  n_leidas <- 0L

  repeat {
    lineas <- readLines(
      con,
      n = tam_chunk,
      warn = FALSE
    )

    if (length(lineas) == 0L) {
      break
    }

    n_leidas <- n_leidas + length(lineas)

    dt <- data.table::fread(
      text = paste(lineas, collapse = "\n"),
      sep = "\t",
      header = FALSE,
      select = c(1L, 2L, 3L, 4L, 5L, 9L, 10L),
      showProgress = FALSE
    )

    data.table::setnames(
      dt,
      c("CHROM", "POS", "ID", "REF", "ALT", "FORMAT", "SAMPLE")
    )

    formatos_permitidos <- c(
      "ES:SE:LP:AF",
      "ES:SE:LP:AF:ID"
    )

    formatos_inesperados <- setdiff(
      unique(dt$FORMAT),
      formatos_permitidos
    )

    if (length(formatos_inesperados) > 0L) {
      stop(
        paste0(
          "Se encontró un FORMAT no esperado: ",
          paste(formatos_inesperados, collapse = ", ")
        )
      )
    }

    campos <- data.table::tstrsplit(
      dt$SAMPLE,
      ":",
      fixed = TRUE,
      keep = 1:4
    )

    dt[, ES := suppressWarnings(as.numeric(campos[[1]]))]
    dt[, SE := suppressWarnings(as.numeric(campos[[2]]))]
    dt[, LP := suppressWarnings(as.numeric(campos[[3]]))]
    dt[, AF := suppressWarnings(as.numeric(campos[[4]]))]

    seleccion <- dt[
      !is.na(LP) & LP > umbral,
      .(CHROM, POS, ID, REF, ALT, ES, SE, LP, AF)
    ]

    if (nrow(seleccion) > 0L) {
      salida[[k]] <- seleccion
      k <- k + 1L
    }

    rm(dt, seleccion, campos, lineas)
  }

  gwas_lp4 <- data.table::rbindlist(salida, use.names = TRUE)

  list(
    datos = gwas_lp4,
    n_variantes_vcf = n_leidas
  )
}

filtrado <- filtrar_vcf_lp(
  archivo_entrada = CONFIG$archivo_vcf,
  umbral = CONFIG$umbral_lp_inicial,
  tam_chunk = CONFIG$tam_chunk
)

gwas_lp4 <- filtrado$datos
n_variantes_vcf <- filtrado$n_variantes_vcf

gwas_lp4[, CHROM := normalizar_chr(CHROM)]
gwas_lp4[, POS := as.integer(POS)]
gwas_lp4[, ID := as.character(ID)]
gwas_lp4[, REF := toupper(as.character(REF))]
gwas_lp4[, ALT := toupper(as.character(ALT))]
gwas_lp4[, ES := as.numeric(ES)]
gwas_lp4[, SE := as.numeric(SE)]
gwas_lp4[, LP := as.numeric(LP)]
gwas_lp4[, AF := as.numeric(AF)]
gwas_lp4[, P := 10^(-LP)]

if (nrow(gwas_lp4) != 10234L) {
  warning(
    paste0(
      "El archivo LP > 4 contiene ",
      nrow(gwas_lp4),
      " variantes; en el TFM se documentaron 10.234."
    )
  )
}

saveRDS(
  gwas_lp4,
  file.path("datos", "antidepressants_LP_gt4.rds"),
  compress = TRUE
)

data.table::fwrite(
  gwas_lp4,
  file.path("datos", "antidepressants_LP_gt4.tsv"),
  sep = "\t"
)

resumen_filtrado <- data.frame(
  N_VARIANTES_VCF = n_variantes_vcf,
  UMBRAL_LP_INICIAL = CONFIG$umbral_lp_inicial,
  N_LP_GT4 = nrow(gwas_lp4),
  UMBRAL_GENOMICO = CONFIG$umbral_genomico,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 1.3. Selección de variantes con significación genómica
# ------------------------------------------------------------------------------

variantes_significativas <- gwas_lp4[
  !is.na(LP) & LP > CONFIG$umbral_genomico
]

data.table::setorder(
  variantes_significativas,
  CHROM,
  POS
)

if (nrow(variantes_significativas) != 11L) {
  stop(
    paste0(
      "Se esperaban 11 variantes con significación genómica y se obtuvieron ",
      nrow(variantes_significativas),
      "."
    )
  )
}

variantes <- as.data.frame(variantes_significativas) |>
  dplyr::transmute(
    CHR = normalizar_chr(.data$CHROM),
    POS = as.numeric(.data$POS),
    SNP = as.character(.data$ID),
    REF = toupper(as.character(.data$REF)),
    ALT = toupper(as.character(.data$ALT)),
    ES = as.numeric(.data$ES),
    SE = as.numeric(.data$SE),
    LP = as.numeric(.data$LP),
    AF_ALT = as.numeric(.data$AF),
    P = 10^(-.data$LP),
    AF_REF = 1 - .data$AF_ALT,
    MAF = pmin(.data$AF_ALT, 1 - .data$AF_ALT),
    MINOR_IS_ALT = .data$AF_ALT <= 0.5,
    ALELO_MINORITARIO = ifelse(.data$MINOR_IS_ALT, .data$ALT, .data$REF),
    TIPO_VARIANTE = dplyr::case_when(
      nchar(.data$REF) == 1L & nchar(.data$ALT) == 1L ~ "SNV",
      nchar(.data$REF) < nchar(.data$ALT) ~ "INS",
      nchar(.data$REF) > nchar(.data$ALT) ~ "DEL",
      TRUE ~ "COMPLEJA"
    )
  )

guardar_csv(
  variantes,
  "02_variantes_significativas.csv",
  "datos"
)

saveRDS(
  variantes,
  file.path("datos", "02_variantes_significativas.rds"),
  compress = TRUE
)

resumen_filtrado$N_GENOME_WIDE <- nrow(variantes)
guardar_csv(resumen_filtrado, "00_resumen_filtrado.csv", "resultados")


# ==============================================================================
# 2. EXPLORACIÓN GRÁFICA
#    2.1. MANHATTAN
#    2.2. VOLCANO
# ==============================================================================

# ------------------------------------------------------------------------------
# 2.1. Preparación del eje acumulado para Manhattan
# ------------------------------------------------------------------------------

orden_chr <- c(as.character(1:22), "X")

gwas_plot <- data.table::copy(gwas_lp4)
gwas_plot[, CHROM_plot := normalizar_chr(CHROM)]

gwas_plot <- gwas_plot[
  CHROM_plot %in% orden_chr &
    !is.na(POS) &
    !is.na(LP)
]

layout_chr <- data.table::data.table(
  CHROM_plot = orden_chr
)

layout_chr[, CHROM_num := seq_len(.N)]

map_longitud <- setNames(
  tabla_contigs$LONGITUD,
  tabla_contigs$CHR
)

layout_chr[, longitud := as.numeric(map_longitud[CHROM_plot])]

if (any(is.na(layout_chr$longitud))) {
  faltan <- layout_chr$CHROM_plot[is.na(layout_chr$longitud)]
  stop(
    paste0(
      "Faltan longitudes cromosómicas para: ",
      paste(faltan, collapse = ", ")
    )
  )
}

layout_chr[, offset_chr := data.table::shift(
  cumsum(longitud),
  fill = 0
)]

layout_chr[, centro_chr := offset_chr + longitud / 2]

map_chr_num <- setNames(
  layout_chr$CHROM_num,
  layout_chr$CHROM_plot
)

map_offset <- setNames(
  layout_chr$offset_chr,
  layout_chr$CHROM_plot
)

gwas_plot[, CHROM_num := as.integer(map_chr_num[CHROM_plot])]
gwas_plot[, offset_chr := as.numeric(map_offset[CHROM_plot])]
gwas_plot[, POS_acumulada := as.numeric(POS) + offset_chr]
gwas_plot[, grupo_color := factor(CHROM_num %% 2L)]

data.table::setorder(gwas_plot, CHROM_num, POS)

lp_max <- max(gwas_plot$LP, na.rm = TRUE)

p_manhattan <- ggplot2::ggplot(
  gwas_plot,
  ggplot2::aes(
    x = POS_acumulada,
    y = LP,
    color = grupo_color
  )
) +
  ggplot2::geom_point(
    size = 0.72,
    alpha = 0.72
  ) +
  ggplot2::geom_hline(
    yintercept = CONFIG$umbral_lp_inicial,
    linetype = "dotted",
    linewidth = 0.75,
    color = "grey25"
  ) +
  ggplot2::geom_hline(
    yintercept = CONFIG$umbral_genomico,
    linetype = "dashed",
    linewidth = 0.75,
    color = "black"
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "0" = "#8C8C8C",
      "1" = "#303030"
    ),
    guide = "none"
  ) +
  ggplot2::scale_x_continuous(
    breaks = layout_chr$centro_chr,
    labels = layout_chr$CHROM_plot,
    expand = ggplot2::expansion(mult = c(0.004, 0.004)),
    guide = ggplot2::guide_axis(n.dodge = 2)
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(
      CONFIG$umbral_lp_inicial,
      ceiling(lp_max + 0.4),
      by = 1
    ),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::coord_cartesian(
    ylim = c(3.75, ceiling(lp_max + 0.4)),
    clip = "off"
  ) +
  ggplot2::labs(
    x = "Cromosoma",
    y = expression(-log[10](italic(p)))
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(size = 12, margin = ggplot2::margin(t = 10)),
    axis.title.y = ggplot2::element_text(size = 12, margin = ggplot2::margin(r = 9)),
    axis.text.x = ggplot2::element_text(size = 9.5),
    axis.text.y = ggplot2::element_text(size = 9.5),
    axis.ticks.x = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(linewidth = 0.45),
    plot.margin = ggplot2::margin(10, 15, 10, 10)
  )

ggplot2::ggsave(
  filename = file.path("figuras", "Figura_Manhattan_LP_gt4.png"),
  plot = p_manhattan,
  width = 18,
  height = 11.5,
  units = "cm",
  dpi = 600,
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path("figuras", "Figura_Manhattan_LP_gt4.pdf"),
  plot = p_manhattan,
  width = 18,
  height = 11.5,
  units = "cm",
  device = "pdf",
  bg = "white"
)

# ------------------------------------------------------------------------------
# 2.2. Volcano plot
# ------------------------------------------------------------------------------

datos_volcano <- as.data.frame(gwas_lp4)

datos_volcano$grupo_volcano <- ifelse(
  datos_volcano$LP > CONFIG$umbral_genomico,
  "Significación genómica",
  "LP > 4"
)

datos_volcano$grupo_volcano <- factor(
  datos_volcano$grupo_volcano,
  levels = c("LP > 4", "Significación genómica")
)

datos_significativos <- datos_volcano |>
  dplyr::filter(.data$LP > CONFIG$umbral_genomico)

significativos_negativos <- datos_significativos |>
  dplyr::filter(.data$ES < 0)

significativos_positivos <- datos_significativos |>
  dplyr::filter(.data$ES >= 0)

p_volcano <- ggplot2::ggplot(
  datos_volcano,
  ggplot2::aes(x = ES, y = LP)
) +
  ggplot2::geom_vline(
    xintercept = 0,
    color = "grey55",
    linewidth = 0.4
  ) +
  ggplot2::geom_hline(
    yintercept = CONFIG$umbral_lp_inicial,
    linetype = "dotted",
    color = "grey55",
    linewidth = 0.5
  ) +
  ggplot2::geom_hline(
    yintercept = CONFIG$umbral_genomico,
    linetype = "dashed",
    color = "grey20",
    linewidth = 0.7
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      color = grupo_volcano,
      alpha = grupo_volcano,
      size = grupo_volcano
    )
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "LP > 4" = "grey65",
      "Significación genómica" = "#C62828"
    ),
    name = NULL
  ) +
  ggplot2::scale_alpha_manual(
    values = c(
      "LP > 4" = 0.50,
      "Significación genómica" = 0.95
    ),
    guide = "none"
  ) +
  ggplot2::scale_size_manual(
    values = c(
      "LP > 4" = 1.15,
      "Significación genómica" = 2.25
    ),
    guide = "none"
  ) +
  ggrepel::geom_text_repel(
    data = significativos_negativos,
    ggplot2::aes(x = ES, y = LP, label = ID),
    inherit.aes = FALSE,
    nudge_x = -0.30,
    direction = "y",
    hjust = 1,
    ylim = c(CONFIG$umbral_genomico + 0.12, Inf),
    size = 3.2,
    box.padding = 0.40,
    point.padding = 0.30,
    segment.color = "grey45",
    min.segment.length = 0,
    force = 2.5,
    max.overlaps = Inf,
    seed = 123
  ) +
  ggrepel::geom_text_repel(
    data = significativos_positivos,
    ggplot2::aes(x = ES, y = LP, label = ID),
    inherit.aes = FALSE,
    nudge_x = 0.30,
    direction = "y",
    hjust = 0,
    ylim = c(CONFIG$umbral_genomico + 0.12, Inf),
    size = 3.2,
    box.padding = 0.40,
    point.padding = 0.30,
    segment.color = "grey45",
    min.segment.length = 0,
    force = 2.5,
    max.overlaps = Inf,
    seed = 123
  ) +
  ggplot2::coord_cartesian(
    ylim = c(
      3.8,
      max(datos_volcano$LP, na.rm = TRUE) + 0.5
    ),
    clip = "off"
  ) +
  ggplot2::labs(
    x = "Tamaño del efecto (ES)",
    y = expression(-log[10](italic(p)))
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(size = 12, margin = ggplot2::margin(t = 8)),
    axis.title.y = ggplot2::element_text(size = 12, margin = ggplot2::margin(r = 8)),
    axis.text = ggplot2::element_text(size = 10, color = "black"),
    legend.position = "top",
    plot.margin = ggplot2::margin(10, 50, 10, 50)
  )

ggplot2::ggsave(
  filename = file.path("figuras", "Figura_Volcano.png"),
  plot = p_volcano,
  width = 18,
  height = 12,
  units = "cm",
  dpi = 600,
  bg = "white"
)


# ==============================================================================
# 3. DISTRIBUCIÓN CROMOSÓMICA DE LAS VARIANTES SIGNIFICATIVAS
# ==============================================================================

snps_chr6 <- variantes |>
  dplyr::filter(.data$CHR == "6") |>
  dplyr::arrange(.data$POS)

if (nrow(snps_chr6) != 10L) {
  stop("Se esperaban 10 variantes significativas en el cromosoma 6.")
}

# ------------------------------------------------------------------------------
# 3.1. Ideograma del cromosoma 6 con karyoploteR
# ------------------------------------------------------------------------------

gr_chr6 <- GenomicRanges::GRanges(
  seqnames = rep("chr6", nrow(snps_chr6)),
  ranges = IRanges::IRanges(
    start = snps_chr6$POS,
    width = 1
  )
)

S4Vectors::mcols(gr_chr6)$SNP <- snps_chr6$SNP

png(
  filename = file.path("figuras", "Figura_chr6_ideograma.png"),
  width = 1600,
  height = 1000,
  res = 180
)

kp <- karyoploteR::plotKaryotype(
  genome = "hg19",
  chromosomes = "chr6"
)

karyoploteR::kpPlotMarkers(
  kp,
  data = gr_chr6,
  labels = S4Vectors::mcols(gr_chr6)$SNP,
  text.orientation = "vertical"
)

dev.off()

# ------------------------------------------------------------------------------
# 3.2. Ampliación 29,0–33,5 Mb
# ------------------------------------------------------------------------------

snps_chr6_plot <- snps_chr6 |>
  dplyr::mutate(
    POS_Mb = .data$POS / 1e6,
    region = dplyr::case_when(
      .data$POS_Mb < 31 ~ "Región 29–31 Mb",
      .data$POS_Mb < 32.5 ~ "Región 31–32,5 Mb",
      TRUE ~ "Región 32,5–33 Mb"
    ),
    SNP = factor(.data$SNP, levels = rev(.data$SNP))
  )

p_chr6_zoom <- ggplot2::ggplot(
  snps_chr6_plot,
  ggplot2::aes(
    x = POS_Mb,
    y = SNP,
    color = region
  )
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = 29,
      xend = POS_Mb,
      y = SNP,
      yend = SNP
    ),
    linewidth = 0.4,
    color = "grey75"
  ) +
  ggplot2::geom_point(size = 3) +
  ggplot2::scale_x_continuous(
    limits = c(29, 33.5),
    breaks = seq(29, 33.5, by = 0.5),
    labels = scales::label_number(decimal.mark = ",")
  ) +
  ggplot2::labs(
    title = "Distribución de SNPs significativos en cromosoma 6",
    subtitle = "Región 29,0–33,5 Mb | GRCh37/hg19",
    x = "Posición genómica en cromosoma 6 (Mb)",
    y = "SNP",
    color = "Subregión"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank()
  )

ggplot2::ggsave(
  filename = file.path("figuras", "Figura_chr6_region_29_33_5Mb.png"),
  plot = p_chr6_zoom,
  width = 10,
  height = 6,
  dpi = 300
)


# ==============================================================================
# 4. ANÁLISIS DE DESEQUILIBRIO DE LIGAMIENTO Y ESTRUCTURA HAPLOTÍPICA
# ==============================================================================

# Este paso se realizó EXTERNAMENTE mediante la plataforma web LDlink.
# No se inventa aquí una automatización que no haya formado parte del análisis.
#
# Parámetros documentados en Materiales y Métodos:
#
# LDmatrix
# - Población: FIN (Finnish in Finland)
# - Ensamblado: GRCh37
# - Variantes evaluadas:
#   rs1612063
#   rs2857702
#   rs140872668
#   rs2858877
#   rs2071536
#   rs115493740
#   rs3749981
#   rs241407
#   rs114463908
# - Métricas: r² y D'
#
# Grupo examinado con mayor detalle:
#   rs2071536
#   rs115493740
#   rs3749981
#   rs241407
#   rs114463908
#
# LDhap
# - Misma población FIN
# - Mismas cinco variantes del grupo anterior
#
# Guardar las exportaciones originales de LDlink, si se desea conservarlas
# junto con el proyecto, dentro de datos/ldlink/.
#
# IMPORTANTE:
# estos datos se utilizaron como aproximación poblacional basada en 1000 Genomes
# y NO como reconstrucción del LD de la cohorte FinnGen.


# ==============================================================================
# 5. ANOTACIÓN FUNCIONAL DE LAS VARIANTES PRIORIZADAS
# ==============================================================================

# ------------------------------------------------------------------------------
# 5.1. Ensembl Variant Effect Predictor (VEP), GRCh37
# ------------------------------------------------------------------------------

snps_stats <- variantes

consultar_vep <- function(ids) {

  req <- httr2::request("https://grch37.rest.ensembl.org") |>
    httr2::req_url_path_append("vep", "human", "id") |>
    httr2::req_url_query(
      pick = "1",
      canonical = "1",
      symbol = "1",
      hgvs = "1"
    ) |>
    httr2::req_headers(
      "Content-Type" = "application/json",
      "Accept" = "application/json"
    ) |>
    httr2::req_body_json(list(ids = ids))

  resp <- httr2::req_perform(req)

  httr2::resp_body_json(
    resp,
    simplifyVector = FALSE
  )
}

parsear_vep <- function(v) {

  tc <- v$transcript_consequences

  if (!is.null(tc) && length(tc) > 0L) {

    pick_idx <- which(
      purrr::map_lgl(
        tc,
        function(x) {
          isTRUE(x$pick == 1) || isTRUE(x$pick == TRUE)
        }
      )
    )

    elegido <- if (length(pick_idx) > 0L) {
      tc[[pick_idx[[1]]]]
    } else {
      tc[[1]]
    }

    tibble::tibble(
      SNP = v$id %||% NA_character_,
      input = v$input %||% NA_character_,
      vep_chr = v$seq_region_name %||% NA_character_,
      vep_start = v$start %||% NA_real_,
      vep_end = v$end %||% NA_real_,
      allele_string = v$allele_string %||% NA_character_,
      most_severe_consequence = v$most_severe_consequence %||% NA_character_,
      selected_consequence = paste(
        unlist(elegido$consequence_terms %||% NA),
        collapse = ";"
      ),
      impact = elegido$impact %||% NA_character_,
      gene_symbol = elegido$gene_symbol %||% NA_character_,
      gene_id = elegido$gene_id %||% NA_character_,
      transcript_id = elegido$transcript_id %||% NA_character_,
      biotype = elegido$biotype %||% NA_character_,
      canonical = elegido$canonical %||% NA
    )

  } else {

    tibble::tibble(
      SNP = v$id %||% NA_character_,
      input = v$input %||% NA_character_,
      vep_chr = v$seq_region_name %||% NA_character_,
      vep_start = v$start %||% NA_real_,
      vep_end = v$end %||% NA_real_,
      allele_string = v$allele_string %||% NA_character_,
      most_severe_consequence = v$most_severe_consequence %||% NA_character_,
      selected_consequence = NA_character_,
      impact = NA_character_,
      gene_symbol = NA_character_,
      gene_id = NA_character_,
      transcript_id = NA_character_,
      biotype = NA_character_,
      canonical = NA
    )
  }
}

ids_vep <- snps_stats$SNP

vep_raw_list <- purrr::map(
  split(ids_vep, ceiling(seq_along(ids_vep) / 50)),
  consultar_vep
)

vep_raw <- unlist(
  vep_raw_list,
  recursive = FALSE
)

vep_table <- purrr::map_dfr(
  vep_raw,
  parsear_vep
)

vep_con_gwas <- vep_table |>
  dplyr::left_join(
    snps_stats |>
      dplyr::select(
        .data$SNP,
        .data$CHR,
        .data$POS,
        .data$REF,
        .data$ALT,
        .data$ES,
        .data$SE,
        .data$LP,
        .data$AF_ALT,
        .data$P
      ),
    by = "SNP"
  )

vep_primary <- vep_con_gwas |>
  dplyr::filter(
    normalizar_chr(.data$vep_chr) == normalizar_chr(.data$CHR)
  )

priorizar_consecuencia <- function(consequence) {
  dplyr::case_when(
    stringr::str_detect(
      consequence,
      "stop_gained|frameshift_variant|splice_acceptor_variant|splice_donor_variant"
    ) ~ 1,
    stringr::str_detect(consequence, "missense_variant") ~ 2,
    stringr::str_detect(
      consequence,
      "5_prime_UTR_variant|3_prime_UTR_variant"
    ) ~ 3,
    stringr::str_detect(
      consequence,
      "non_coding_transcript_variant"
    ) ~ 4,
    stringr::str_detect(consequence, "intron_variant") ~ 5,
    stringr::str_detect(consequence, "upstream_gene_variant") ~ 6,
    stringr::str_detect(consequence, "downstream_gene_variant") ~ 7,
    stringr::str_detect(consequence, "intergenic_variant") ~ 8,
    TRUE ~ 99
  )
}

vep_resumen <- vep_primary |>
  dplyr::mutate(
    prioridad = priorizar_consecuencia(.data$selected_consequence),
    canonical_num = dplyr::case_when(
      as.character(.data$canonical) %in% c("1", "TRUE") ~ 1,
      TRUE ~ 0
    )
  ) |>
  dplyr::arrange(
    .data$SNP,
    .data$prioridad,
    dplyr::desc(.data$canonical_num)
  ) |>
  dplyr::group_by(.data$SNP) |>
  dplyr::slice(1) |>
  dplyr::ungroup()

guardar_csv(
  vep_table,
  "03_VEP_completo.csv",
  "resultados"
)

guardar_csv(
  vep_resumen,
  "04_VEP_resumen_variantes.csv",
  "resultados"
)

# ------------------------------------------------------------------------------
# 5.2. GTEx v8: eQTL cerebrales
# ------------------------------------------------------------------------------

options(gtexr.itemsPerPage = 100000)
options(gtexr.verbose = FALSE)

brain_tissues <- c(
  "Brain_Cortex",
  "Brain_Frontal_Cortex_BA9",
  "Brain_Hippocampus",
  "Brain_Amygdala",
  "Brain_Anterior_cingulate_cortex_BA24",
  "Brain_Caudate_basal_ganglia",
  "Brain_Nucleus_accumbens_basal_ganglia",
  "Brain_Putamen_basal_ganglia",
  "Brain_Hypothalamus",
  "Brain_Cerebellum",
  "Brain_Cerebellar_Hemisphere"
)

variant_map <- purrr::map_dfr(
  snps_stats$SNP,
  function(rs) {
    tryCatch(
      {
        gtexr::get_variant(snpId = rs) |>
          dplyr::mutate(SNP = rs)
      },
      error = function(e) {
        tibble::tibble(
          SNP = rs,
          snpId = rs,
          variantId = NA_character_,
          b37VariantId = NA_character_
        )
      }
    )
  }
)

variant_map_clean <- variant_map |>
  dplyr::filter(!is.na(.data$variantId)) |>
  dplyr::distinct(
    .data$SNP,
    .data$variantId,
    .keep_all = TRUE
  )

gtex_eqtl_brain_raw <- tryCatch(
  {
    gtexr::get_significant_single_tissue_eqtls(
      variantIds = variant_map_clean$variantId,
      tissueSiteDetailIds = brain_tissues,
      datasetId = "gtex_v8"
    )
  },
  error = function(e) {
    warning(
      paste0(
        "No se pudo completar la consulta de GTEx: ",
        conditionMessage(e)
      )
    )
    tibble::tibble()
  }
)

if (nrow(gtex_eqtl_brain_raw) > 0L) {

  gtex_eqtl_brain <- gtex_eqtl_brain_raw |>
    dplyr::left_join(
      variant_map_clean |>
        dplyr::select(.data$SNP, .data$variantId),
      by = "variantId"
    )

  gtex_summary <- gtex_eqtl_brain |>
    dplyr::group_by(.data$SNP) |>
    dplyr::summarise(
      n_brain_eqtl = dplyr::n(),
      eQTL_genes = paste(
        unique(.data$geneSymbol),
        collapse = "; "
      ),
      eQTL_tissues = paste(
        unique(.data$tissueSiteDetailId),
        collapse = "; "
      ),
      min_pvalue_gtex = min(.data$pValue, na.rm = TRUE),
      .groups = "drop"
    )

} else {

  gtex_eqtl_brain <- tibble::tibble()

  gtex_summary <- tibble::tibble(
    SNP = snps_stats$SNP,
    n_brain_eqtl = 0L,
    eQTL_genes = NA_character_,
    eQTL_tissues = NA_character_,
    min_pvalue_gtex = NA_real_
  )
}

guardar_csv(
  variant_map,
  "05_GTEx_mapeo_variantes.csv",
  "resultados"
)

guardar_csv(
  gtex_eqtl_brain,
  "06_GTEx_eQTL_cerebro_detalle.csv",
  "resultados"
)

guardar_csv(
  gtex_summary,
  "07_GTEx_eQTL_cerebro_resumen.csv",
  "resultados"
)


# ==============================================================================
# 6. SIMULACIÓN DE DATOS GENOTÍPICOS INDIVIDUALES
#    CONSTRUCCIÓN DEL POOL GENOTÍPICO
# ==============================================================================

set.seed(CONFIG$semilla_pool)

sexo_pool <- sample(
  c(
    rep("F", CONFIG$n_mujeres_pool),
    rep("M", CONFIG$n_varones_pool)
  ),
  size = CONFIG$n_pool,
  replace = FALSE
)

pool <- data.frame(
  ID = sprintf("IND_%05d", seq_len(CONFIG$n_pool)),
  SEXO = sexo_pool,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 6.1. Generación genotípica bajo HWE
#
# Autosomas:
# rbinom(n, size = 2, prob = q)
# -> 0, 1 o 2 copias ALT.
#
# Cromosoma X:
# mujeres: size = 2
# varones: size = 1
# ------------------------------------------------------------------------------

for (i in seq_len(nrow(variantes))) {

  snp <- variantes$SNP[[i]]
  chr <- variantes$CHR[[i]]
  q <- variantes$AF_ALT[[i]]

  if (!es_chr_x(chr)) {

    pool[[snp]] <- stats::rbinom(
      n = CONFIG$n_pool,
      size = 2L,
      prob = q
    )

  } else {

    g <- integer(CONFIG$n_pool)

    idx_f <- which(pool$SEXO == "F")
    idx_m <- which(pool$SEXO == "M")

    g[idx_f] <- stats::rbinom(
      n = length(idx_f),
      size = 2L,
      prob = q
    )

    g[idx_m] <- stats::rbinom(
      n = length(idx_m),
      size = 1L,
      prob = q
    )

    pool[[snp]] <- g
  }
}

if (anyDuplicated(pool$ID) > 0L) {
  stop("El pool contiene identificadores duplicados.")
}

# ------------------------------------------------------------------------------
# 6.2. Control de codificación, AF y HWE
# ------------------------------------------------------------------------------

qc_af <- list()
qc_hwe <- list()

for (i in seq_len(nrow(variantes))) {

  snp <- variantes$SNP[[i]]
  chr <- variantes$CHR[[i]]
  q <- variantes$AF_ALT[[i]]
  p <- 1 - q

  g <- pool[[snp]]

  if (!es_chr_x(chr)) {

    if (!all(g %in% c(0L, 1L, 2L))) {
      stop(paste("Codificación inválida en", snp))
    }

    n_alelos <- 2 * nrow(pool)
    af_observada <- sum(g) / n_alelos

    observados <- c(
      sum(g == 0),
      sum(g == 1),
      sum(g == 2)
    )

    esperados <- nrow(pool) * c(
      p^2,
      2 * p * q,
      q^2
    )

    chi2 <- sum(
      (observados - esperados)^2 / esperados
    )

    p_hwe <- stats::pchisq(
      chi2,
      df = 1,
      lower.tail = FALSE
    )

    qc_hwe[[length(qc_hwe) + 1L]] <- data.frame(
      SNP = snp,
      CHR = chr,
      ESTRATO = "Todos",
      N = nrow(pool),
      CHI2_HWE = chi2,
      P_HWE = p_hwe
    )

  } else {

    idx_f <- pool$SEXO == "F"
    idx_m <- pool$SEXO == "M"

    if (!all(g[idx_f] %in% c(0L, 1L, 2L))) {
      stop(paste("Codificación femenina inválida en", snp))
    }

    if (!all(g[idx_m] %in% c(0L, 1L))) {
      stop(paste("Codificación masculina inválida en", snp))
    }

    n_alelos <- 2 * sum(idx_f) + sum(idx_m)

    af_observada <- (
      sum(g[idx_f]) +
        sum(g[idx_m])
    ) / n_alelos

    g_f <- g[idx_f]

    observados <- c(
      sum(g_f == 0),
      sum(g_f == 1),
      sum(g_f == 2)
    )

    esperados <- length(g_f) * c(
      p^2,
      2 * p * q,
      q^2
    )

    chi2 <- sum(
      (observados - esperados)^2 / esperados
    )

    p_hwe <- stats::pchisq(
      chi2,
      df = 1,
      lower.tail = FALSE
    )

    qc_hwe[[length(qc_hwe) + 1L]] <- data.frame(
      SNP = snp,
      CHR = chr,
      ESTRATO = "Mujeres",
      N = length(g_f),
      CHI2_HWE = chi2,
      P_HWE = p_hwe
    )
  }

  qc_af[[i]] <- data.frame(
    SNP = snp,
    CHR = chr,
    AF_ALT_OBJETIVO = q,
    AF_ALT_OBSERVADA = af_observada,
    DIFERENCIA = af_observada - q,
    DIFERENCIA_ABSOLUTA = abs(af_observada - q)
  )
}

qc_af <- dplyr::bind_rows(qc_af)
qc_hwe <- dplyr::bind_rows(qc_hwe)

guardar_csv(
  qc_af,
  "08_QC_pool_frecuencias_alelicas.csv",
  "resultados"
)

guardar_csv(
  qc_hwe,
  "09_QC_pool_HWE.csv",
  "resultados"
)

saveRDS(
  pool,
  file.path("datos", "pool_genotipico_10000.rds"),
  compress = TRUE
)


# ==============================================================================
# 7. CONSTRUCCIÓN DE LOS ESCENARIOS SIMULADOS
# ==============================================================================

# ------------------------------------------------------------------------------
# Función: proporción del alelo minoritario por locus
# Escala 0–1, respetando ploidía.
# ------------------------------------------------------------------------------

proporcion_minoritaria_snp <- function(datos, meta_fila) {

  snp <- meta_fila$SNP[[1]]
  chr <- meta_fila$CHR[[1]]
  alt_es_minor <- meta_fila$MINOR_IS_ALT[[1]]

  g_alt <- datos[[snp]]

  if (!es_chr_x(chr)) {

    dosis_minor <- if (alt_es_minor) {
      g_alt
    } else {
      2 - g_alt
    }

    return(dosis_minor / 2)
  }

  mujeres <- datos$SEXO == "F"
  varones <- datos$SEXO == "M"

  dosis_minor <- numeric(nrow(datos))

  if (alt_es_minor) {
    dosis_minor[mujeres] <- g_alt[mujeres]
    dosis_minor[varones] <- g_alt[varones]
  } else {
    dosis_minor[mujeres] <- 2 - g_alt[mujeres]
    dosis_minor[varones] <- 1 - g_alt[varones]
  }

  proporcion <- numeric(nrow(datos))
  proporcion[mujeres] <- dosis_minor[mujeres] / 2
  proporcion[varones] <- dosis_minor[varones]

  proporcion
}

pool_escenarios <- pool

for (i in seq_len(nrow(variantes))) {

  snp <- variantes$SNP[[i]]

  pool_escenarios[[paste0("M_", snp)]] <-
    proporcion_minoritaria_snp(
      pool_escenarios,
      variantes[i, , drop = FALSE]
    )
}

variables_minor <- paste0(
  "M_",
  variantes$SNP
)

# ------------------------------------------------------------------------------
# Función común: muestreo probabilístico por sexo, sin reemplazo
#
# score_z ya debe estar estandarizado.
# w_caso = exp(+lambda * score_z)
# w_control = exp(-lambda * score_z)
# ------------------------------------------------------------------------------

seleccionar_ponderado_por_sexo <- function(
    datos,
    score_z,
    semilla
) {

  stopifnot(length(score_z) == nrow(datos))

  set.seed(semilla)

  casos_idx <- integer()
  controles_idx <- integer()

  for (sx in c("F", "M")) {

    idx <- which(datos$SEXO == sx)

    pesos_caso <- exp(
      CONFIG$lambda * score_z[idx]
    )

    caso_local <- sample(
      seq_along(idx),
      size = CONFIG$n_casos_por_sexo,
      replace = FALSE,
      prob = pesos_caso
    )

    idx_caso <- idx[caso_local]

    restantes <- setdiff(
      idx,
      idx_caso
    )

    pesos_control <- exp(
      -CONFIG$lambda * score_z[restantes]
    )

    control_local <- sample(
      seq_along(restantes),
      size = CONFIG$n_controles_por_sexo,
      replace = FALSE,
      prob = pesos_control
    )

    idx_control <- restantes[control_local]

    casos_idx <- c(casos_idx, idx_caso)
    controles_idx <- c(controles_idx, idx_control)
  }

  if (length(intersect(casos_idx, controles_idx)) > 0L) {
    stop("Hay individuos compartidos entre casos y controles.")
  }

  casos <- datos[casos_idx, , drop = FALSE]
  controles <- datos[controles_idx, , drop = FALSE]

  casos$clase <- "caso"
  controles$clase <- "control"

  out <- dplyr::bind_rows(
    casos,
    controles
  )

  out$clase <- factor(
    out$clase,
    levels = c("control", "caso")
  )

  out
}

# ------------------------------------------------------------------------------
# 7.1. E1 — control negativo
# ------------------------------------------------------------------------------

generar_e1 <- function(datos) {

  set.seed(CONFIG$semilla_e1)

  elegidos <- list()

  for (sx in c("F", "M")) {

    idx <- which(datos$SEXO == sx)

    seleccion <- sample(
      idx,
      size = CONFIG$n_casos_por_sexo +
        CONFIG$n_controles_por_sexo,
      replace = FALSE
    )

    elegidos[[paste0("caso_", sx)]] <-
      seleccion[seq_len(CONFIG$n_casos_por_sexo)]

    elegidos[[paste0("control_", sx)]] <-
      seleccion[
        CONFIG$n_casos_por_sexo +
          seq_len(CONFIG$n_controles_por_sexo)
      ]
  }

  idx_casos <- c(
    elegidos$caso_F,
    elegidos$caso_M
  )

  idx_controles <- c(
    elegidos$control_F,
    elegidos$control_M
  )

  casos <- datos[idx_casos, , drop = FALSE]
  controles <- datos[idx_controles, , drop = FALSE]

  casos$clase <- "caso"
  controles$clase <- "control"

  out <- dplyr::bind_rows(
    casos,
    controles
  )

  out$clase <- factor(
    out$clase,
    levels = c("control", "caso")
  )

  out
}

E1_completo <- generar_e1(
  pool_escenarios
)

# ------------------------------------------------------------------------------
# 7.2. E2 — dependencia / combinación multilocus sintética
# ------------------------------------------------------------------------------

chr6 <- variantes |>
  dplyr::filter(.data$CHR == "6") |>
  dplyr::arrange(.data$POS)

pares_e2 <- data.frame(
  SNP_A = chr6$SNP[-nrow(chr6)],
  SNP_B = chr6$SNP[-1],
  POS_A = chr6$POS[-nrow(chr6)],
  POS_B = chr6$POS[-1],
  stringsAsFactors = FALSE
) |>
  dplyr::mutate(
    DISTANCIA_BP = .data$POS_B - .data$POS_A
  ) |>
  dplyr::filter(
    .data$DISTANCIA_BP <= 100000
  )

if (nrow(pares_e2) != 5L) {
  stop(
    paste0(
      "El criterio de proximidad <=100 kb produjo ",
      nrow(pares_e2),
      " pares; se esperaban 5."
    )
  )
}

guardar_csv(
  pares_e2,
  "10_E2_pares_chr6.csv",
  "resultados"
)

# Términos de coocurrencia I_jk = M_j * M_k.
nombres_interacciones <- character()

for (i in seq_len(nrow(pares_e2))) {

  a <- paste0("M_", pares_e2$SNP_A[[i]])
  b <- paste0("M_", pares_e2$SNP_B[[i]])
  nombre <- paste0("INT_", i)

  pool_escenarios[[nombre]] <-
    pool_escenarios[[a]] *
    pool_escenarios[[b]]

  nombres_interacciones <- c(
    nombres_interacciones,
    nombre
  )
}

# Estandarizar cada interacción para dar igual escala a los cinco pares.
for (nombre in nombres_interacciones) {

  media <- mean(
    pool_escenarios[[nombre]]
  )

  desviacion <- stats::sd(
    pool_escenarios[[nombre]]
  )

  pool_escenarios[[paste0(nombre, "_Z")]] <-
    (
      pool_escenarios[[nombre]] -
        media
    ) / desviacion
}

interacciones_z <- paste0(
  nombres_interacciones,
  "_Z"
)

pool_escenarios$SCORE_COMB_RAW <- rowMeans(
  pool_escenarios[, interacciones_z, drop = FALSE]
)

# Eliminar componente marginal:
# SCORE_COMB_RAW ~ M_1 + ... + M_11 + SEXO
formula_residual <- stats::as.formula(
  paste(
    "SCORE_COMB_RAW ~",
    paste(
      c(
        variables_minor,
        "SEXO"
      ),
      collapse = " + "
    )
  )
)

modelo_residual <- stats::lm(
  formula_residual,
  data = pool_escenarios
)

pool_escenarios$SCORE_COMB_RESIDUAL <-
  stats::residuals(modelo_residual)

pool_escenarios$Z_COMB <- as.numeric(
  scale(
    pool_escenarios$SCORE_COMB_RESIDUAL
  )
)

E2_completo <- seleccionar_ponderado_por_sexo(
  datos = pool_escenarios,
  score_z = pool_escenarios$Z_COMB,
  semilla = CONFIG$semilla_e2
)

# ------------------------------------------------------------------------------
# 7.3. E3 — señal aditiva / marginal de alelos minoritarios
# ------------------------------------------------------------------------------

pool_escenarios$SCORE_FREQ_RAW <- rowSums(
  pool_escenarios[
    ,
    variables_minor,
    drop = FALSE
  ]
)

pool_escenarios$Z_FREQ <- as.numeric(
  scale(
    pool_escenarios$SCORE_FREQ_RAW
  )
)

E3_completo <- seleccionar_ponderado_por_sexo(
  datos = pool_escenarios,
  score_z = pool_escenarios$Z_FREQ,
  semilla = CONFIG$semilla_e3
)

# ------------------------------------------------------------------------------
# 7.4. Dataset final de cada escenario
#
# Los modelos reciben solamente:
# ID, SEXO, 11 dosis ALT, clase.
#
# SEXO se conserva para estratificar, pero NO se usa como predictor.
# Scores e interacciones NO se proporcionan a los modelos.
# ------------------------------------------------------------------------------

columnas_escenario <- c(
  "ID",
  "SEXO",
  variantes$SNP,
  "clase"
)

E1 <- E1_completo[, columnas_escenario, drop = FALSE]
E2 <- E2_completo[, columnas_escenario, drop = FALSE]
E3 <- E3_completo[, columnas_escenario, drop = FALSE]

escenarios <- list(
  E1 = E1,
  E2 = E2,
  E3 = E3
)

for (nm in names(escenarios)) {

  esc <- escenarios[[nm]]

  if (nrow(esc) != 2000L) {
    stop(paste(nm, "no contiene 2.000 individuos."))
  }

  if (sum(esc$clase == "caso") != 1000L) {
    stop(paste(nm, "no contiene 1.000 casos."))
  }

  if (sum(esc$clase == "control") != 1000L) {
    stop(paste(nm, "no contiene 1.000 controles."))
  }

  if (anyDuplicated(esc$ID) > 0L) {
    stop(paste(nm, "contiene individuos duplicados."))
  }

  guardar_csv(
    esc,
    paste0("11_", nm, "_dataset_modelado.csv"),
    "datos"
  )
}


# ==============================================================================
# 8. PREPARACIÓN Y AJUSTE DE LOS DATOS PARA MACHINE LEARNING
# ==============================================================================

PREDICTORES <- variantes$SNP

# ------------------------------------------------------------------------------
# 8.1. Split 80/20 estratificado por clase × sexo
# ------------------------------------------------------------------------------

crear_split <- function(
    escenario,
    nombre
) {

  set.seed(
    CONFIG$semilla_split_base +
      match(nombre, c("E1", "E2", "E3"))
  )

  train_ids <- character()
  test_ids <- character()

  for (cl in c("control", "caso")) {
    for (sx in c("F", "M")) {

      ids <- escenario$ID[
        escenario$clase == cl &
          escenario$SEXO == sx
      ]

      ids <- sample(
        ids,
        length(ids),
        replace = FALSE
      )

      n_train <- floor(
        length(ids) *
          CONFIG$prop_train
      )

      train_ids <- c(
        train_ids,
        ids[seq_len(n_train)]
      )

      test_ids <- c(
        test_ids,
        ids[
          n_train +
            seq_len(length(ids) - n_train)
        ]
      )
    }
  }

  list(
    train_ids = train_ids,
    test_ids = test_ids
  )
}

# ------------------------------------------------------------------------------
# 8.2. CV estratificada repetida: 5 folds × 5 repeticiones
# ------------------------------------------------------------------------------

crear_resamples_cv <- function(
    train,
    escenario
) {

  res <- list()
  contador <- 1L

  e <- match(
    escenario,
    c("E1", "E2", "E3")
  )

  for (
    rep in seq_len(
      CONFIG$cv_repeticiones
    )
  ) {

    fold_id <- integer(
      nrow(train)
    )

    set.seed(
      CONFIG$semilla_cv_base +
        1000L * e +
        rep
    )

    estratos <- interaction(
      train$clase,
      train$SEXO,
      drop = TRUE,
      lex.order = TRUE
    )

    for (st in levels(estratos)) {

      idx <- which(
        estratos == st
      )

      idx <- sample(
        idx,
        length(idx),
        replace = FALSE
      )

      asignacion <- rep(
        seq_len(
          CONFIG$cv_folds
        ),
        length.out = length(idx)
      )

      fold_id[idx] <- asignacion
    }

    for (
      fold in seq_len(
        CONFIG$cv_folds
      )
    ) {

      res[[contador]] <- list(
        RESAMPLE_ID = sprintf(
          "R%02d_F%02d",
          rep,
          fold
        ),
        train_idx = which(
          fold_id != fold
        ),
        val_idx = which(
          fold_id == fold
        )
      )

      contador <- contador + 1L
    }
  }

  res
}

PREPARACION <- list()

for (nm in names(escenarios)) {

  esc <- escenarios[[nm]]
  split <- crear_split(
    esc,
    nm
  )

  train <- esc[
    match(
      split$train_ids,
      esc$ID
    ),
    ,
    drop = FALSE
  ]

  test <- esc[
    match(
      split$test_ids,
      esc$ID
    ),
    ,
    drop = FALSE
  ]

  train$clase <- factor(
    train$clase,
    levels = c("control", "caso")
  )

  test$clase <- factor(
    test$clase,
    levels = c("control", "caso")
  )

  resamples <- crear_resamples_cv(
    train,
    nm
  )

  PREPARACION[[nm]] <- list(
    train = train,
    test = test,
    resamples = resamples
  )
}

qc_split <- lapply(
  names(PREPARACION),
  function(nm) {

    x <- PREPARACION[[nm]]

    dplyr::bind_rows(
      x$train |>
        dplyr::count(
          .data$clase,
          .data$SEXO,
          name = "N"
        ) |>
        dplyr::mutate(
          ESCENARIO = nm,
          CONJUNTO = "train"
        ),

      x$test |>
        dplyr::count(
          .data$clase,
          .data$SEXO,
          name = "N"
        ) |>
        dplyr::mutate(
          ESCENARIO = nm,
          CONJUNTO = "test"
        )
    )
  }
) |>
  dplyr::bind_rows()

guardar_csv(
  qc_split,
  "12_QC_particiones_train_test.csv",
  "resultados"
)

saveRDS(
  lapply(
    PREPARACION,
    function(x) x$resamples
  ),
  file.path(
    "datos",
    "resamples_CV_5x5.rds"
  )
)


# ==============================================================================
# 9. MODELOS DE MACHINE LEARNING
# ==============================================================================

# Las funciones de esta sección DEFINEN los modelos.
# La ejecución se realiza después de definir las métricas comunes.

# ------------------------------------------------------------------------------
# 9.1. Regresión logística
# ------------------------------------------------------------------------------

ajustar_logistica <- function(
    train,
    validacion = NULL
) {

  formula <- stats::as.formula(
    paste(
      "clase ~",
      paste(
        PREDICTORES,
        collapse = " + "
      )
    )
  )

  fit <- stats::glm(
    formula,
    data = train,
    family = stats::binomial(
      link = "logit"
    )
  )

  if (is.null(validacion)) {
    return(fit)
  }

  list(
    fit = fit,
    prob = stats::predict(
      fit,
      newdata = validacion,
      type = "response"
    )
  )
}

# ------------------------------------------------------------------------------
# 9.2. Random Forest
# ------------------------------------------------------------------------------

GRID_RF <- expand.grid(
  num.trees = c(500L, 1000L),
  mtry = c(2L, 3L, 6L, 11L),
  min.node.size = c(1L, 10L, 20L),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

GRID_RF$CONFIG_ID <- seq_len(
  nrow(GRID_RF)
)

# ------------------------------------------------------------------------------
# 9.3. XGBoost
# ------------------------------------------------------------------------------

GRID_XGB <- expand.grid(
  nrounds = c(100L, 200L),
  max_depth = c(2L, 3L),
  eta = c(0.05, 0.10),
  min_child_weight = c(1, 5),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

GRID_XGB$subsample <- 0.8
GRID_XGB$colsample_bytree <- 0.8
GRID_XGB$CONFIG_ID <- seq_len(
  nrow(GRID_XGB)
)

# ------------------------------------------------------------------------------
# 9.4. LightGBM
# ------------------------------------------------------------------------------

GRID_LGB <- expand.grid(
  num_iterations = c(100L, 200L),
  learning_rate = c(0.03, 0.05),
  num_leaves = c(7L, 15L),
  min_data_in_leaf = c(20L, 40L),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

GRID_LGB$max_depth <- 3L
GRID_LGB$feature_fraction <- 0.8
GRID_LGB$bagging_fraction <- 0.8
GRID_LGB$bagging_freq <- 1L
GRID_LGB$CONFIG_ID <- seq_len(
  nrow(GRID_LGB)
)

# ------------------------------------------------------------------------------
# 9.5. Red neuronal feed-forward, una capa oculta
# ------------------------------------------------------------------------------

GRID_NN <- expand.grid(
  size = c(3L, 4L, 6L),
  decay = c(0, 0.01, 0.10),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

GRID_NN$maxit <- 500L
GRID_NN$CONFIG_ID <- seq_len(
  nrow(GRID_NN)
)

escalador_fit <- function(x) {

  media <- colMeans(x)
  sdv <- apply(
    x,
    2,
    stats::sd
  )

  sdv[
    !is.finite(sdv) |
      sdv == 0
  ] <- 1

  list(
    media = media,
    sd = sdv
  )
}

escalador_apply <- function(
    x,
    esc
) {

  x <- sweep(
    x,
    2,
    esc$media,
    "-"
  )

  sweep(
    x,
    2,
    esc$sd,
    "/"
  )
}


# ==============================================================================
# 10. MÉTRICAS DE EVALUACIÓN
# ==============================================================================

calcular_metricas <- function(
    truth,
    prob,
    umbral = 0.5
) {

  truth <- factor(
    truth,
    levels = c(
      CONFIG$clase_negativa,
      CONFIG$clase_positiva
    )
  )

  prob <- clip_prob(prob)

  pred <- factor(
    ifelse(
      prob >= umbral,
      CONFIG$clase_positiva,
      CONFIG$clase_negativa
    ),
    levels = c(
      CONFIG$clase_negativa,
      CONFIG$clase_positiva
    )
  )

  tp <- sum(
    pred == "caso" &
      truth == "caso"
  )

  tn <- sum(
    pred == "control" &
      truth == "control"
  )

  fp <- sum(
    pred == "caso" &
      truth == "control"
  )

  fn <- sum(
    pred == "control" &
      truth == "caso"
  )

  sensibilidad <- if (
    (tp + fn) > 0
  ) {
    tp / (tp + fn)
  } else {
    NA_real_
  }

  especificidad <- if (
    (tn + fp) > 0
  ) {
    tn / (tn + fp)
  } else {
    NA_real_
  }

  precision <- if (
    (tp + fp) > 0
  ) {
    tp / (tp + fp)
  } else {
    NA_real_
  }

  recall <- sensibilidad

  f1 <- if (
    is.finite(precision) &&
      is.finite(recall) &&
      (precision + recall) > 0
  ) {
    2 * precision * recall /
      (precision + recall)
  } else {
    NA_real_
  }

  roc_obj <- suppressMessages(
    pROC::roc(
      response = truth,
      predictor = prob,
      levels = c(
        "control",
        "caso"
      ),
      direction = "<",
      quiet = TRUE
    )
  )

  roc_auc <- as.numeric(
    pROC::auc(
      roc_obj
    )
  )

  pr_auc <- PRROC::pr.curve(
    scores.class0 = prob[
      truth == "caso"
    ],
    scores.class1 = prob[
      truth == "control"
    ],
    curve = FALSE
  )$auc.integral

  y <- as.integer(
    truth == "caso"
  )

  brier <- mean(
    (prob - y)^2
  )

  logloss <- -mean(
    y * log(prob) +
      (1 - y) * log(1 - prob)
  )

  c(
    ROC_AUC = roc_auc,
    PR_AUC = pr_auc,
    ACCURACY = mean(
      pred == truth
    ),
    BALANCED_ACCURACY = mean(
      c(
        sensibilidad,
        especificidad
      ),
      na.rm = TRUE
    ),
    SENSITIVITY = sensibilidad,
    SPECIFICITY = especificidad,
    PRECISION = precision,
    RECALL = recall,
    F1 = f1,
    BRIER = brier,
    LOGLOSS = logloss,
    TP = tp,
    TN = tn,
    FP = fp,
    FN = fn
  )
}

resumir_cv <- function(cv_raw) {

  metricas <- c(
    "ROC_AUC",
    "PR_AUC",
    "ACCURACY",
    "BALANCED_ACCURACY",
    "SENSITIVITY",
    "SPECIFICITY",
    "PRECISION",
    "RECALL",
    "F1",
    "BRIER",
    "LOGLOSS"
  )

  cv_raw |>
    dplyr::group_by(
      .data$ESCENARIO,
      .data$MODELO,
      .data$CONFIG_ID,
      .data$PARAMETROS
    ) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(metricas),
        list(
          MEDIA = ~mean(.x, na.rm = TRUE),
          SD = ~stats::sd(.x, na.rm = TRUE)
        )
      ),
      N_RESAMPLES = dplyr::n(),
      .groups = "drop"
    )
}

seleccionar_mejor_config <- function(
    cv_resumen
) {

  cv_resumen |>
    dplyr::arrange(
      dplyr::desc(
        .data$ROC_AUC_MEDIA
      ),
      dplyr::desc(
        .data$PR_AUC_MEDIA
      ),
      dplyr::desc(
        .data$BALANCED_ACCURACY_MEDIA
      ),
      .data$CONFIG_ID
    ) |>
    dplyr::slice(1)
}

param_string <- function(x) {
  paste(
    paste0(
      names(x),
      "=",
      unlist(
        x,
        use.names = FALSE
      )
    ),
    collapse = ";"
  )
}


# ==============================================================================
# 11. VALIDACIÓN CRUZADA, TUNING Y EVALUACIÓN FINAL
# ==============================================================================

# ------------------------------------------------------------------------------
# 11.1. Regresión logística
# ------------------------------------------------------------------------------

cv_logistica <- function(
    train,
    resamples,
    escenario
) {

  filas <- list()

  for (r in seq_along(resamples)) {

    rs <- resamples[[r]]

    tr <- train[
      rs$train_idx,
      ,
      drop = FALSE
    ]

    va <- train[
      rs$val_idx,
      ,
      drop = FALSE
    ]

    ajuste <- ajustar_logistica(
      tr,
      va
    )

    met <- calcular_metricas(
      va$clase,
      ajuste$prob
    )

    filas[[r]] <- data.frame(
      ESCENARIO = escenario,
      MODELO = "logistica",
      CONFIG_ID = 1L,
      PARAMETROS = "AJUSTE=MLE_LOGIT",
      RESAMPLE_ID = rs$RESAMPLE_ID,
      as.list(met),
      stringsAsFactors = FALSE
    )
  }

  dplyr::bind_rows(
    filas
  )
}

final_logistica <- function(
    train,
    test,
    escenario
) {

  fit <- ajustar_logistica(
    train
  )

  prob <- stats::predict(
    fit,
    newdata = test,
    type = "response"
  )

  met <- calcular_metricas(
    test$clase,
    prob
  )

  sm <- summary(fit)$coefficients

  coef_df <- data.frame(
    SNP = rownames(sm),
    BETA = sm[, "Estimate"],
    SE = sm[, "Std. Error"],
    Z = sm[, "z value"],
    P = sm[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(
      .data$SNP != "(Intercept)"
    ) |>
    dplyr::mutate(
      IMPORTANCIA = abs(.data$Z),
      ESCENARIO = escenario,
      MODELO = "logistica",
      METODO = "abs_z"
    )

  list(
    fit = fit,
    metricas = data.frame(
      ESCENARIO = escenario,
      MODELO = "logistica",
      CONFIG_ID = 1L,
      PARAMETROS = "AJUSTE=MLE_LOGIT",
      as.list(met),
      stringsAsFactors = FALSE
    ),
    predicciones = data.frame(
      ID = test$ID,
      clase = test$clase,
      PROB_CASO = prob,
      PRED_CLASE = ifelse(
        prob >= 0.5,
        "caso",
        "control"
      ),
      ESCENARIO = escenario,
      MODELO = "logistica",
      stringsAsFactors = FALSE
    ),
    importancia = coef_df
  )
}

# ------------------------------------------------------------------------------
# 11.2. Random Forest
# ------------------------------------------------------------------------------

cv_rf <- function(
    train,
    resamples,
    escenario
) {

  filas <- list()
  k <- 1L

  for (
    g in seq_len(
      nrow(GRID_RF)
    )
  ) {

    cfg <- GRID_RF[
      g,
      ,
      drop = FALSE
    ]

    params <- param_string(
      list(
        num.trees = cfg$num.trees[[1]],
        mtry = cfg$mtry[[1]],
        min.node.size = cfg$min.node.size[[1]]
      )
    )

    for (r in seq_along(resamples)) {

      rs <- resamples[[r]]

      tr <- train[
        rs$train_idx,
        c(PREDICTORES, "clase"),
        drop = FALSE
      ]

      va <- train[
        rs$val_idx,
        c(PREDICTORES, "clase"),
        drop = FALSE
      ]

      seed <- CONFIG$semilla_cv_base +
        g * 100L +
        r

      fit <- ranger::ranger(
        clase ~ .,
        data = tr,
        probability = TRUE,
        num.trees = cfg$num.trees[[1]],
        mtry = cfg$mtry[[1]],
        min.node.size = cfg$min.node.size[[1]],
        importance = "permutation",
        seed = seed
      )

      prob <- predict(
        fit,
        data = va
      )$predictions[, "caso"]

      met <- calcular_metricas(
        va$clase,
        prob
      )

      filas[[k]] <- data.frame(
        ESCENARIO = escenario,
        MODELO = "random_forest",
        CONFIG_ID = cfg$CONFIG_ID[[1]],
        PARAMETROS = params,
        RESAMPLE_ID = rs$RESAMPLE_ID,
        as.list(met),
        stringsAsFactors = FALSE
      )

      k <- k + 1L
    }
  }

  dplyr::bind_rows(
    filas
  )
}

final_rf <- function(
    train,
    test,
    escenario,
    mejor
) {

  id <- mejor$CONFIG_ID[[1]]

  cfg <- GRID_RF[
    GRID_RF$CONFIG_ID == id,
    ,
    drop = FALSE
  ]

  fit <- ranger::ranger(
    clase ~ .,
    data = train[
      ,
      c(PREDICTORES, "clase"),
      drop = FALSE
    ],
    probability = TRUE,
    num.trees = cfg$num.trees[[1]],
    mtry = cfg$mtry[[1]],
    min.node.size = cfg$min.node.size[[1]],
    importance = "permutation",
    seed = CONFIG$semilla_cv_base + id
  )

  prob <- predict(
    fit,
    data = test
  )$predictions[, "caso"]

  met <- calcular_metricas(
    test$clase,
    prob
  )

  importancia <- data.frame(
    SNP = names(
      fit$variable.importance
    ),
    IMPORTANCIA = as.numeric(
      fit$variable.importance
    ),
    ESCENARIO = escenario,
    MODELO = "random_forest",
    METODO = "permutacion_ranger",
    stringsAsFactors = FALSE
  ) |>
    dplyr::arrange(
      dplyr::desc(
        .data$IMPORTANCIA
      )
    )

  list(
    fit = fit,
    metricas = data.frame(
      ESCENARIO = escenario,
      MODELO = "random_forest",
      CONFIG_ID = id,
      PARAMETROS = param_string(
        as.list(cfg[
          ,
          c(
            "num.trees",
            "mtry",
            "min.node.size"
          )
        ])
      ),
      as.list(met),
      stringsAsFactors = FALSE
    ),
    predicciones = data.frame(
      ID = test$ID,
      clase = test$clase,
      PROB_CASO = prob,
      PRED_CLASE = ifelse(
        prob >= 0.5,
        "caso",
        "control"
      ),
      ESCENARIO = escenario,
      MODELO = "random_forest",
      stringsAsFactors = FALSE
    ),
    importancia = importancia
  )
}

# ------------------------------------------------------------------------------
# 11.3. XGBoost
# ------------------------------------------------------------------------------

xgb_matriz <- function(datos) {
  data.matrix(
    datos[
      ,
      PREDICTORES,
      drop = FALSE
    ]
  )
}

cv_xgb <- function(
    train,
    resamples,
    escenario
) {

  filas <- list()
  k <- 1L

  for (
    g in seq_len(
      nrow(GRID_XGB)
    )
  ) {

    cfg <- GRID_XGB[
      g,
      ,
      drop = FALSE
    ]

    params_texto <- param_string(
      as.list(
        cfg[
          ,
          c(
            "nrounds",
            "max_depth",
            "eta",
            "min_child_weight",
            "subsample",
            "colsample_bytree"
          )
        ]
      )
    )

    for (r in seq_along(resamples)) {

      rs <- resamples[[r]]

      tr <- train[
        rs$train_idx,
        ,
        drop = FALSE
      ]

      va <- train[
        rs$val_idx,
        ,
        drop = FALSE
      ]

      x_tr <- xgb_matriz(tr)
      y_tr <- as.integer(
        tr$clase == "caso"
      )

      x_va <- xgb_matriz(va)

      dtrain <- xgboost::xgb.DMatrix(
        x_tr,
        label = y_tr
      )

      seed <- CONFIG$semilla_cv_base +
        10000L +
        g * 100L +
        r

      params <- list(
        objective = "binary:logistic",
        eval_metric = "auc",
        max_depth = cfg$max_depth[[1]],
        eta = cfg$eta[[1]],
        min_child_weight = cfg$min_child_weight[[1]],
        subsample = cfg$subsample[[1]],
        colsample_bytree = cfg$colsample_bytree[[1]],
        seed = seed,
        nthread = 1L
      )

      fit <- xgboost::xgb.train(
        params = params,
        data = dtrain,
        nrounds = cfg$nrounds[[1]],
        verbose = 0
      )

      prob <- predict(
        fit,
        x_va
      )

      met <- calcular_metricas(
        va$clase,
        prob
      )

      filas[[k]] <- data.frame(
        ESCENARIO = escenario,
        MODELO = "xgboost",
        CONFIG_ID = cfg$CONFIG_ID[[1]],
        PARAMETROS = params_texto,
        RESAMPLE_ID = rs$RESAMPLE_ID,
        as.list(met),
        stringsAsFactors = FALSE
      )

      k <- k + 1L
    }
  }

  dplyr::bind_rows(
    filas
  )
}

final_xgb <- function(
    train,
    test,
    escenario,
    mejor
) {

  id <- mejor$CONFIG_ID[[1]]

  cfg <- GRID_XGB[
    GRID_XGB$CONFIG_ID == id,
    ,
    drop = FALSE
  ]

  x_tr <- xgb_matriz(train)
  y_tr <- as.integer(
    train$clase == "caso"
  )

  x_te <- xgb_matriz(test)

  dtrain <- xgboost::xgb.DMatrix(
    x_tr,
    label = y_tr
  )

  params <- list(
    objective = "binary:logistic",
    eval_metric = "auc",
    max_depth = cfg$max_depth[[1]],
    eta = cfg$eta[[1]],
    min_child_weight = cfg$min_child_weight[[1]],
    subsample = cfg$subsample[[1]],
    colsample_bytree = cfg$colsample_bytree[[1]],
    seed = CONFIG$semilla_cv_base + id,
    nthread = 1L
  )

  fit <- xgboost::xgb.train(
    params = params,
    data = dtrain,
    nrounds = cfg$nrounds[[1]],
    verbose = 0
  )

  prob <- predict(
    fit,
    x_te
  )

  met <- calcular_metricas(
    test$clase,
    prob
  )

  imp <- xgboost::xgb.importance(
    feature_names = PREDICTORES,
    model = fit
  )

  importancia <- data.frame(
    SNP = imp$Feature,
    IMPORTANCIA = imp$Gain,
    ESCENARIO = escenario,
    MODELO = "xgboost",
    METODO = "gain",
    stringsAsFactors = FALSE
  )

  list(
    fit = fit,
    metricas = data.frame(
      ESCENARIO = escenario,
      MODELO = "xgboost",
      CONFIG_ID = id,
      PARAMETROS = param_string(
        as.list(
          cfg[
            ,
            c(
              "nrounds",
              "max_depth",
              "eta",
              "min_child_weight",
              "subsample",
              "colsample_bytree"
            )
          ]
        )
      ),
      as.list(met),
      stringsAsFactors = FALSE
    ),
    predicciones = data.frame(
      ID = test$ID,
      clase = test$clase,
      PROB_CASO = prob,
      PRED_CLASE = ifelse(
        prob >= 0.5,
        "caso",
        "control"
      ),
      ESCENARIO = escenario,
      MODELO = "xgboost",
      stringsAsFactors = FALSE
    ),
    importancia = importancia
  )
}

# ------------------------------------------------------------------------------
# 11.4. LightGBM
# ------------------------------------------------------------------------------

cv_lgb <- function(
    train,
    resamples,
    escenario
) {

  filas <- list()
  k <- 1L

  for (
    g in seq_len(
      nrow(GRID_LGB)
    )
  ) {

    cfg <- GRID_LGB[
      g,
      ,
      drop = FALSE
    ]

    params_texto <- param_string(
      as.list(
        cfg[
          ,
          c(
            "num_iterations",
            "learning_rate",
            "num_leaves",
            "min_data_in_leaf",
            "max_depth",
            "feature_fraction",
            "bagging_fraction",
            "bagging_freq"
          )
        ]
      )
    )

    for (r in seq_along(resamples)) {

      rs <- resamples[[r]]

      tr <- train[
        rs$train_idx,
        ,
        drop = FALSE
      ]

      va <- train[
        rs$val_idx,
        ,
        drop = FALSE
      ]

      x_tr <- data.matrix(
        tr[
          ,
          PREDICTORES,
          drop = FALSE
        ]
      )

      y_tr <- as.integer(
        tr$clase == "caso"
      )

      x_va <- data.matrix(
        va[
          ,
          PREDICTORES,
          drop = FALSE
        ]
      )

      dtrain <- lightgbm::lgb.Dataset(
        data = x_tr,
        label = y_tr,
        free_raw_data = FALSE
      )

      seed <- CONFIG$semilla_cv_base +
        20000L +
        g * 100L +
        r

      params <- list(
        objective = "binary",
        metric = "auc",
        learning_rate = cfg$learning_rate[[1]],
        num_leaves = cfg$num_leaves[[1]],
        min_data_in_leaf = cfg$min_data_in_leaf[[1]],
        max_depth = cfg$max_depth[[1]],
        feature_fraction = cfg$feature_fraction[[1]],
        bagging_fraction = cfg$bagging_fraction[[1]],
        bagging_freq = cfg$bagging_freq[[1]],
        seed = seed,
        feature_fraction_seed = seed,
        bagging_seed = seed,
        num_threads = 1L,
        verbosity = -1L
      )

      fit <- lightgbm::lgb.train(
        params = params,
        data = dtrain,
        nrounds = cfg$num_iterations[[1]],
        verbose = -1
      )

      prob <- predict(
        fit,
        x_va
      )

      met <- calcular_metricas(
        va$clase,
        prob
      )

      filas[[k]] <- data.frame(
        ESCENARIO = escenario,
        MODELO = "lightgbm",
        CONFIG_ID = cfg$CONFIG_ID[[1]],
        PARAMETROS = params_texto,
        RESAMPLE_ID = rs$RESAMPLE_ID,
        as.list(met),
        stringsAsFactors = FALSE
      )

      k <- k + 1L
    }
  }

  dplyr::bind_rows(
    filas
  )
}

final_lgb <- function(
    train,
    test,
    escenario,
    mejor
) {

  id <- mejor$CONFIG_ID[[1]]

  cfg <- GRID_LGB[
    GRID_LGB$CONFIG_ID == id,
    ,
    drop = FALSE
  ]

  x_tr <- data.matrix(
    train[
      ,
      PREDICTORES,
      drop = FALSE
    ]
  )

  y_tr <- as.integer(
    train$clase == "caso"
  )

  x_te <- data.matrix(
    test[
      ,
      PREDICTORES,
      drop = FALSE
    ]
  )

  dtrain <- lightgbm::lgb.Dataset(
    data = x_tr,
    label = y_tr,
    free_raw_data = FALSE
  )

  seed <- CONFIG$semilla_cv_base + id

  params <- list(
    objective = "binary",
    metric = "auc",
    learning_rate = cfg$learning_rate[[1]],
    num_leaves = cfg$num_leaves[[1]],
    min_data_in_leaf = cfg$min_data_in_leaf[[1]],
    max_depth = cfg$max_depth[[1]],
    feature_fraction = cfg$feature_fraction[[1]],
    bagging_fraction = cfg$bagging_fraction[[1]],
    bagging_freq = cfg$bagging_freq[[1]],
    seed = seed,
    feature_fraction_seed = seed,
    bagging_seed = seed,
    num_threads = 1L,
    verbosity = -1L
  )

  fit <- lightgbm::lgb.train(
    params = params,
    data = dtrain,
    nrounds = cfg$num_iterations[[1]],
    verbose = -1
  )

  prob <- predict(
    fit,
    x_te
  )

  met <- calcular_metricas(
    test$clase,
    prob
  )

  imp <- lightgbm::lgb.importance(
    fit,
    percentage = TRUE
  )

  importancia <- data.frame(
    SNP = imp$Feature,
    IMPORTANCIA = imp$Gain,
    ESCENARIO = escenario,
    MODELO = "lightgbm",
    METODO = "gain",
    stringsAsFactors = FALSE
  )

  list(
    fit = fit,
    metricas = data.frame(
      ESCENARIO = escenario,
      MODELO = "lightgbm",
      CONFIG_ID = id,
      PARAMETROS = param_string(
        as.list(
          cfg[
            ,
            c(
              "num_iterations",
              "learning_rate",
              "num_leaves",
              "min_data_in_leaf",
              "max_depth",
              "feature_fraction",
              "bagging_fraction",
              "bagging_freq"
            )
          ]
        )
      ),
      as.list(met),
      stringsAsFactors = FALSE
    ),
    predicciones = data.frame(
      ID = test$ID,
      clase = test$clase,
      PROB_CASO = prob,
      PRED_CLASE = ifelse(
        prob >= 0.5,
        "caso",
        "control"
      ),
      ESCENARIO = escenario,
      MODELO = "lightgbm",
      stringsAsFactors = FALSE
    ),
    importancia = importancia
  )
}

# ------------------------------------------------------------------------------
# 11.5. Red neuronal
# ------------------------------------------------------------------------------

cv_nn <- function(
    train,
    resamples,
    escenario
) {

  filas <- list()
  k <- 1L

  for (
    g in seq_len(
      nrow(GRID_NN)
    )
  ) {

    cfg <- GRID_NN[
      g,
      ,
      drop = FALSE
    ]

    params_texto <- param_string(
      list(
        size = cfg$size[[1]],
        decay = cfg$decay[[1]],
        maxit = cfg$maxit[[1]]
      )
    )

    for (r in seq_along(resamples)) {

      rs <- resamples[[r]]

      tr <- train[
        rs$train_idx,
        ,
        drop = FALSE
      ]

      va <- train[
        rs$val_idx,
        ,
        drop = FALSE
      ]

      x_tr0 <- data.matrix(
        tr[
          ,
          PREDICTORES,
          drop = FALSE
        ]
      )

      x_va0 <- data.matrix(
        va[
          ,
          PREDICTORES,
          drop = FALSE
        ]
      )

      esc <- escalador_fit(
        x_tr0
      )

      x_tr <- escalador_apply(
        x_tr0,
        esc
      )

      x_va <- escalador_apply(
        x_va0,
        esc
      )

      y_tr <- as.integer(
        tr$clase == "caso"
      )

      seed <- CONFIG$semilla_cv_base +
        30000L +
        g * 100L +
        r

      set.seed(seed)

      fit <- nnet::nnet(
        x = x_tr,
        y = y_tr,
        size = cfg$size[[1]],
        decay = cfg$decay[[1]],
        maxit = cfg$maxit[[1]],
        trace = FALSE,
        linout = FALSE,
        rang = 0.1,
        MaxNWts = 10000
      )

      prob <- as.numeric(
        predict(
          fit,
          x_va,
          type = "raw"
        )
      )

      met <- calcular_metricas(
        va$clase,
        prob
      )

      filas[[k]] <- data.frame(
        ESCENARIO = escenario,
        MODELO = "red_neuronal",
        CONFIG_ID = cfg$CONFIG_ID[[1]],
        PARAMETROS = params_texto,
        RESAMPLE_ID = rs$RESAMPLE_ID,
        as.list(met),
        stringsAsFactors = FALSE
      )

      k <- k + 1L
    }
  }

  dplyr::bind_rows(
    filas
  )
}

final_nn <- function(
    train,
    test,
    escenario,
    mejor
) {

  id <- mejor$CONFIG_ID[[1]]

  cfg <- GRID_NN[
    GRID_NN$CONFIG_ID == id,
    ,
    drop = FALSE
  ]

  x_tr0 <- data.matrix(
    train[
      ,
      PREDICTORES,
      drop = FALSE
    ]
  )

  x_te0 <- data.matrix(
    test[
      ,
      PREDICTORES,
      drop = FALSE
    ]
  )

  esc <- escalador_fit(
    x_tr0
  )

  x_tr <- escalador_apply(
    x_tr0,
    esc
  )

  x_te <- escalador_apply(
    x_te0,
    esc
  )

  y_tr <- as.integer(
    train$clase == "caso"
  )

  set.seed(
    CONFIG$semilla_cv_base + id
  )

  fit <- nnet::nnet(
    x = x_tr,
    y = y_tr,
    size = cfg$size[[1]],
    decay = cfg$decay[[1]],
    maxit = cfg$maxit[[1]],
    trace = FALSE,
    linout = FALSE,
    rang = 0.1,
    MaxNWts = 10000
  )

  prob <- as.numeric(
    predict(
      fit,
      x_te,
      type = "raw"
    )
  )

  met <- calcular_metricas(
    test$clase,
    prob
  )

  list(
    fit = fit,
    metricas = data.frame(
      ESCENARIO = escenario,
      MODELO = "red_neuronal",
      CONFIG_ID = id,
      PARAMETROS = param_string(
        list(
          size = cfg$size[[1]],
          decay = cfg$decay[[1]],
          maxit = cfg$maxit[[1]]
        )
      ),
      as.list(met),
      stringsAsFactors = FALSE
    ),
    predicciones = data.frame(
      ID = test$ID,
      clase = test$clase,
      PROB_CASO = prob,
      PRED_CLASE = ifelse(
        prob >= 0.5,
        "caso",
        "control"
      ),
      ESCENARIO = escenario,
      MODELO = "red_neuronal",
      stringsAsFactors = FALSE
    )
  )
}


# ==============================================================================
# 12. EJECUCIÓN DE LOS CINCO MODELOS EN LOS TRES ESCENARIOS
# ==============================================================================

CV_RAW <- list()
CV_RESUMEN <- list()
METRICAS_TEST <- list()
PREDICCIONES_TEST <- list()
IMPORTANCIAS <- list()
MEJORES_CONFIG <- list()

for (escenario in names(PREPARACION)) {

  message(
    "Procesando ",
    escenario,
    "..."
  )

  train <- PREPARACION[[escenario]]$train
  test <- PREPARACION[[escenario]]$test
  resamples <- PREPARACION[[escenario]]$resamples

  # --------------------------------------------------------------------------
  # Logística
  # --------------------------------------------------------------------------

  cv_log <- cv_logistica(
    train,
    resamples,
    escenario
  )

  cv_log_res <- resumir_cv(
    cv_log
  )

  fin_log <- final_logistica(
    train,
    test,
    escenario
  )

  CV_RAW[[paste0(escenario, "_log")]] <- cv_log
  CV_RESUMEN[[paste0(escenario, "_log")]] <- cv_log_res
  METRICAS_TEST[[paste0(escenario, "_log")]] <- fin_log$metricas
  PREDICCIONES_TEST[[paste0(escenario, "_log")]] <- fin_log$predicciones
  IMPORTANCIAS[[paste0(escenario, "_log")]] <- fin_log$importancia
  MEJORES_CONFIG[[paste0(escenario, "_log")]] <- data.frame(
    ESCENARIO = escenario,
    MODELO = "logistica",
    CONFIG_ID = 1L,
    PARAMETROS = "AJUSTE=MLE_LOGIT"
  )

  saveRDS(
    fin_log$fit,
    file.path(
      "modelos",
      paste0(
        escenario,
        "_logistica.rds"
      )
    )
  )

  # --------------------------------------------------------------------------
  # Random Forest
  # --------------------------------------------------------------------------

  cv_rf_raw <- cv_rf(
    train,
    resamples,
    escenario
  )

  cv_rf_res <- resumir_cv(
    cv_rf_raw
  )

  best_rf <- seleccionar_mejor_config(
    cv_rf_res
  )

  fin_rf <- final_rf(
    train,
    test,
    escenario,
    best_rf
  )

  CV_RAW[[paste0(escenario, "_rf")]] <- cv_rf_raw
  CV_RESUMEN[[paste0(escenario, "_rf")]] <- cv_rf_res
  METRICAS_TEST[[paste0(escenario, "_rf")]] <- fin_rf$metricas
  PREDICCIONES_TEST[[paste0(escenario, "_rf")]] <- fin_rf$predicciones
  IMPORTANCIAS[[paste0(escenario, "_rf")]] <- fin_rf$importancia
  MEJORES_CONFIG[[paste0(escenario, "_rf")]] <- best_rf

  saveRDS(
    fin_rf$fit,
    file.path(
      "modelos",
      paste0(
        escenario,
        "_random_forest.rds"
      )
    )
  )

  # --------------------------------------------------------------------------
  # XGBoost
  # --------------------------------------------------------------------------

  cv_xgb_raw <- cv_xgb(
    train,
    resamples,
    escenario
  )

  cv_xgb_res <- resumir_cv(
    cv_xgb_raw
  )

  best_xgb <- seleccionar_mejor_config(
    cv_xgb_res
  )

  fin_xgb <- final_xgb(
    train,
    test,
    escenario,
    best_xgb
  )

  CV_RAW[[paste0(escenario, "_xgb")]] <- cv_xgb_raw
  CV_RESUMEN[[paste0(escenario, "_xgb")]] <- cv_xgb_res
  METRICAS_TEST[[paste0(escenario, "_xgb")]] <- fin_xgb$metricas
  PREDICCIONES_TEST[[paste0(escenario, "_xgb")]] <- fin_xgb$predicciones
  IMPORTANCIAS[[paste0(escenario, "_xgb")]] <- fin_xgb$importancia
  MEJORES_CONFIG[[paste0(escenario, "_xgb")]] <- best_xgb

  xgboost::xgb.save(
    fin_xgb$fit,
    file.path(
      "modelos",
      paste0(
        escenario,
        "_xgboost.ubj"
      )
    )
  )

  # --------------------------------------------------------------------------
  # LightGBM
  # --------------------------------------------------------------------------

  cv_lgb_raw <- cv_lgb(
    train,
    resamples,
    escenario
  )

  cv_lgb_res <- resumir_cv(
    cv_lgb_raw
  )

  best_lgb <- seleccionar_mejor_config(
    cv_lgb_res
  )

  fin_lgb <- final_lgb(
    train,
    test,
    escenario,
    best_lgb
  )

  CV_RAW[[paste0(escenario, "_lgb")]] <- cv_lgb_raw
  CV_RESUMEN[[paste0(escenario, "_lgb")]] <- cv_lgb_res
  METRICAS_TEST[[paste0(escenario, "_lgb")]] <- fin_lgb$metricas
  PREDICCIONES_TEST[[paste0(escenario, "_lgb")]] <- fin_lgb$predicciones
  IMPORTANCIAS[[paste0(escenario, "_lgb")]] <- fin_lgb$importancia
  MEJORES_CONFIG[[paste0(escenario, "_lgb")]] <- best_lgb

  lightgbm::lgb.save(
    fin_lgb$fit,
    file.path(
      "modelos",
      paste0(
        escenario,
        "_lightgbm.txt"
      )
    )
  )

  # --------------------------------------------------------------------------
  # Red neuronal
  # --------------------------------------------------------------------------

  cv_nn_raw <- cv_nn(
    train,
    resamples,
    escenario
  )

  cv_nn_res <- resumir_cv(
    cv_nn_raw
  )

  best_nn <- seleccionar_mejor_config(
    cv_nn_res
  )

  fin_nn <- final_nn(
    train,
    test,
    escenario,
    best_nn
  )

  CV_RAW[[paste0(escenario, "_nn")]] <- cv_nn_raw
  CV_RESUMEN[[paste0(escenario, "_nn")]] <- cv_nn_res
  METRICAS_TEST[[paste0(escenario, "_nn")]] <- fin_nn$metricas
  PREDICCIONES_TEST[[paste0(escenario, "_nn")]] <- fin_nn$predicciones
  MEJORES_CONFIG[[paste0(escenario, "_nn")]] <- best_nn

  saveRDS(
    fin_nn$fit,
    file.path(
      "modelos",
      paste0(
        escenario,
        "_red_neuronal.rds"
      )
    )
  )
}

cv_raw_todo <- dplyr::bind_rows(
  CV_RAW
)

cv_resumen_todo <- dplyr::bind_rows(
  CV_RESUMEN
)

metricas_test <- dplyr::bind_rows(
  METRICAS_TEST
)

predicciones_test <- dplyr::bind_rows(
  PREDICCIONES_TEST
)

importancias <- dplyr::bind_rows(
  IMPORTANCIAS
)

mejores_config <- dplyr::bind_rows(
  MEJORES_CONFIG
)

guardar_csv(
  cv_raw_todo,
  "13_CV_25_resamples_todos_modelos.csv",
  "resultados"
)

guardar_csv(
  cv_resumen_todo,
  "14_CV_resumen_todos_modelos.csv",
  "resultados"
)

guardar_csv(
  mejores_config,
  "15_hiperparametros_seleccionados.csv",
  "resultados"
)

guardar_csv(
  metricas_test,
  "16_metricas_test_final.csv",
  "resultados"
)

guardar_csv(
  predicciones_test,
  "17_predicciones_test.csv",
  "resultados"
)

guardar_csv(
  importancias,
  "18_importancia_variantes.csv",
  "resultados"
)


# ==============================================================================
# 13. COMPARACIONES PAREADAS COMPLEMENTARIAS DE VALIDACIÓN CRUZADA
# ==============================================================================

# Se realizan dentro de cada escenario sobre los mismos RESAMPLE_ID.
# Se consideran exploratorias/descriptivas porque los folds de una CV repetida
# no constituyen observaciones completamente independientes.

cv_mejor_config <- cv_raw_todo |>
  dplyr::inner_join(
    mejores_config |>
      dplyr::select(
        .data$ESCENARIO,
        .data$MODELO,
        .data$CONFIG_ID
      ) |>
      dplyr::distinct(),
    by = c(
      "ESCENARIO",
      "MODELO",
      "CONFIG_ID"
    )
  )

comparar_modelos_pareados <- function(
    cv_best,
    escenario,
    metrica
) {

  dat <- cv_best[
    cv_best$ESCENARIO == escenario,
    ,
    drop = FALSE
  ]

  modelos <- sort(
    unique(dat$MODELO)
  )

  pares <- utils::combn(
    modelos,
    2,
    simplify = FALSE
  )

  lapply(
    pares,
    function(par) {

      a <- dat[
        dat$MODELO == par[[1]],
        c(
          "RESAMPLE_ID",
          metrica
        )
      ]

      b <- dat[
        dat$MODELO == par[[2]],
        c(
          "RESAMPLE_ID",
          metrica
        )
      ]

      names(a)[[2]] <- "A"
      names(b)[[2]] <- "B"

      z <- merge(
        a,
        b,
        by = "RESAMPLE_ID"
      )

      dif <- z$A - z$B

      tt <- stats::t.test(
        dif,
        mu = 0
      )

      data.frame(
        ESCENARIO = escenario,
        METRICA = metrica,
        MODELO_A = par[[1]],
        MODELO_B = par[[2]],
        DIF_MEDIA_A_MENOS_B = mean(dif),
        IC95_INF = tt$conf.int[[1]],
        IC95_SUP = tt$conf.int[[2]],
        P_T_PAREADO = tt$p.value,
        N_RESAMPLES = nrow(z),
        stringsAsFactors = FALSE
      )
    }
  ) |>
    dplyr::bind_rows()
}

comparaciones_pareadas <- dplyr::bind_rows(
  lapply(
    c("E1", "E2", "E3"),
    function(e) {
      dplyr::bind_rows(
        comparar_modelos_pareados(
          cv_mejor_config,
          e,
          "ROC_AUC"
        ),
        comparar_modelos_pareados(
          cv_mejor_config,
          e,
          "BALANCED_ACCURACY"
        ),
        comparar_modelos_pareados(
          cv_mejor_config,
          e,
          "F1"
        )
      )
    }
  )
)

guardar_csv(
  comparaciones_pareadas,
  "19_comparaciones_pareadas_CV_exploratorias.csv",
  "resultados"
)


# ==============================================================================
# 14. REGISTRO DE REPRODUCIBILIDAD
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    "resultados",
    "sessionInfo.txt"
  )
)

cat(
  "\n============================================================\n",
  "PIPELINE COMPLETADO\n",
  "============================================================\n",
  "Resultados: ./resultados\n",
  "Figuras:    ./figuras\n",
  "Modelos:    ./modelos\n",
  "Datos derivados: ./datos\n",
  sep = ""
)
