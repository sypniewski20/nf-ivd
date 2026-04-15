suppressPackageStartupMessages({
  library(tidyverse)
  library(tools)
  library(optparse)
})

# ---------------------------
# CLI options
# ---------------------------
option_list <- list(
  make_option(c("-m", "--manifest"),
              type = "character",
              default = NULL,
              help = "Path to manifest CSV [default: %default]"
  ),
  make_option(c("-d", "--dir"),
              type = "character",
              default = NULL,
              help = "Download directory [default: %default]"
  ),
  make_option(c("-s", "--snapshot_dir"),
              type = "character",
              default = NULL,
              help = "Snapshot directory (optional)"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))


# ---------------------------
# Helper: validation
# ---------------------------
add_validation <- function(df) {
  df %>%
    mutate(
      exists = file.exists(path),
      local_md5 = map_chr(path, ~ {
        if (file.exists(.x)) {
          tools::md5sum(.x) %>% unname()
        } else {
          NA_character_
        }
      }),
      valid = exists & local_md5 == md5
    )
}

# ---------------------------
# Main function
# ---------------------------
download_and_validate <- function(manifest_path, dir, snapshot_dir = NULL) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  
  if (is.null(snapshot_dir)) {
    snapshot_dir <- file.path(dir, "_snapshot")
  }
  dir.create(snapshot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---------------------------
  # 1. Read manifest
  # ---------------------------
  manifest <- readr::read_csv(manifest_path, show_col_types = FALSE) %>%
    mutate(
      file = basename(url),
      path = file.path(dir, file)
    )
  
  # ---------------------------
  # 2. Initial validation
  # ---------------------------
  manifest <- manifest %>%
    add_validation()
  
  # ---------------------------
  # 3. Download missing/invalid
  # ---------------------------
  to_download <- manifest %>%
    filter(!valid)
  
  if (nrow(to_download) > 0) {
    message("Downloading ", nrow(to_download), " files...")
    
    pwalk(
      list(to_download$url, to_download$path),
      ~ {
        if (file.exists(.y)) file.remove(.y)
        system2("curl", c("-L", "--retry", "3", "--retry-delay", "5", "-o", .y, .x))
      }
    )
  } else {
    message("All files already present and valid.")
  }
  
  # ---------------------------
  # 4. Re-validation
  # ---------------------------
  manifest <- manifest %>%
    add_validation()
  
  # ---------------------------
  # 5. Summary
  # ---------------------------
  summary_tbl <- manifest %>%
    summarise(
      total = n(),
      passed = sum(valid),
      failed = sum(!valid)
    )
  
  passed <- manifest %>% filter(valid)
  failed <- manifest %>% filter(!valid)
  
  # ---------------------------
  # 6. Snapshot
  # ---------------------------
  snapshot <- list(
    timestamp = Sys.time(),
    directory = normalizePath(dir),
    summary = summary_tbl,
    manifest = manifest
  )
  
  snapshot_file <- file.path(
    snapshot_dir,
    paste0("snapshot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
  )
  
  saveRDS(snapshot, snapshot_file)
  
  # ---------------------------
  # 7. Logs
  # ---------------------------
  readr::write_csv(passed, file.path(snapshot_dir, "passed_files.csv"))
  readr::write_csv(failed, file.path(snapshot_dir, "failed_files.csv"))
  readr::write_csv(manifest, file.path(snapshot_dir, "full_manifest.csv"))
  
  # ---------------------------
  # 8. Reporting
  # ---------------------------
  message("Snapshot saved: ", snapshot_file)
  
  if (nrow(failed) > 0) {
    warning("Some files failed validation:")
    print(failed %>% select(dataset, file, md5, local_md5))
  } else {
    message("All files passed validation.")
  }
  
  invisible(snapshot)
}

# ---------------------------
# RUN FROM CLI ARGS
# ---------------------------
download_and_validate(
  manifest_path = opt$manifest,
  dir = opt$dir,
  snapshot_dir = opt$dir
)
