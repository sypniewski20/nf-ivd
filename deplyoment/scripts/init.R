get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script_path <- sub("^--file=", "", file_arg)
  normalizePath(dirname(script_path))
}

PROJECT_ROOT <- normalizePath(file.path(get_script_dir(), ".."))
SCRIPTS_DIR  <- file.path(PROJECT_ROOT, "scripts")
DATA_DIR     <- file.path(PROJECT_ROOT, "data")
REF_DIR      <- file.path(PROJECT_ROOT, "reference")

LEDGER_PATH <- file.path(PROJECT_ROOT, "data", "_ledger.json")