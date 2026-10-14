# ============================================================================
# config.R — Configuration from environment variables
# ============================================================================

load_config <- function() {
  config <- list(
    input_dir      = Sys.getenv("INPUT_DIR", "/data/input"),
    output_dir     = Sys.getenv("OUTPUT_DIR", "/data/output"),
    output_mode    = tolower(Sys.getenv("OUTPUT_MODE", "expanded")),
    chunk_size     = as.integer(Sys.getenv("CHUNK_SIZE", "500")),
    environment    = Sys.getenv("ENVIRONMENT", "local"),
    # Peak-gene links (multiome only): compute with Signac::LinkPeaks when the
    # object has none. Genes in LINK_GENES go first, then the most variable
    # genes (so cell-type markers are included), in batches, until
    # LINK_MAX_GENES (0 = all) or the LINK_MAX_MINUTES budget (0 = no limit) is reached.
    link_peaks       = tolower(Sys.getenv("LINK_PEAKS", "true")) %in% c("true", "1", "yes"),
    link_max_genes   = as.integer(Sys.getenv("LINK_MAX_GENES", "5000")),
    link_genes       = parse_gene_list(Sys.getenv("LINK_GENES", "")),
    link_batch_size  = as.integer(Sys.getenv("LINK_BATCH_SIZE", "1000")),
    link_max_minutes = as.numeric(Sys.getenv("LINK_MAX_MINUTES", "0"))
  )

  # Validate output_mode
  if (!config$output_mode %in% c("expanded", "consolidated")) {
    stop(paste0("Invalid OUTPUT_MODE: '", config$output_mode,
                "'. Must be 'expanded' or 'consolidated'."))
  }

  config
}

print_config <- function(config) {
  log_info("Configuration:")
  log_info("  INPUT_DIR:       ", config$input_dir)
  log_info("  OUTPUT_DIR:      ", config$output_dir)
  log_info("  OUTPUT_MODE:     ", config$output_mode)
  log_info("  CHUNK_SIZE:      ", config$chunk_size)
  log_info("  ENVIRONMENT:     ", config$environment)
  log_info("  LINK_PEAKS:      ", config$link_peaks)
  log_info("  LINK_MAX_GENES:  ", config$link_max_genes)
  log_info("  LINK_GENES:      ", if (length(config$link_genes)) paste(config$link_genes, collapse = ", ") else "(none)")
  log_info("  LINK_BATCH_SIZE: ", config$link_batch_size)
  log_info("  LINK_MAX_MINUTES:", config$link_max_minutes)
}

#' "SCN10A, TAC1 PRPH" -> c("SCN10A", "TAC1", "PRPH")
parse_gene_list <- function(value) {
  genes <- trimws(unlist(strsplit(value, "[,;[:space:]]+")))
  unique(genes[genes != ""])
}
