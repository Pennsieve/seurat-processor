# ============================================================================
# config.R — Configuration from environment variables
# ============================================================================

load_config <- function() {
  config <- list(
    input_dir      = Sys.getenv("INPUT_DIR", "/data/input"),
    output_dir     = Sys.getenv("OUTPUT_DIR", "/data/output"),
    output_mode    = tolower(Sys.getenv("OUTPUT_MODE", "expanded")),
    chunk_size     = as.integer(Sys.getenv("CHUNK_SIZE", "500")),
    environment    = Sys.getenv("ENVIRONMENT", "local")
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
}
