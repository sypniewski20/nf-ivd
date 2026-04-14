library(readr)
library(jsonlite)

# ----------------------------
# BOOTSTRAP
# ----------------------------

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script_path <- sub("^--file=", "", file_arg)
  normalizePath(dirname(script_path))
}

SCRIPT_DIR <- get_script_dir()

source(file.path(SCRIPT_DIR, "utils.R"))
source(file.path(SCRIPT_DIR, "init.R"))

# ----------------------------
# PATHS (from init.R)
# ----------------------------

manifest_path <- file.path(PROJECT_ROOT, "manifest.csv")
out_dir <- file.path(PROJECT_ROOT, "data")
ledger_path <- file.path(out_dir, "_ledger.json")

ensure_dir(out_dir)

manifest <- read_csv(manifest_path)

# ----------------------------
# LEDGER HELPER (local fallback if not in utils.R yet)
# ----------------------------

write_ledger <- function(record, path) {

  if (!file.exists(path)) {
    writeLines("[]", path)
  }

  ledger <- fromJSON(path)

  ledger <- append(ledger, list(record))

  write_json(ledger, path, pretty = TRUE, auto_unbox = TRUE)
}

# ----------------------------
# MAIN LOOP
# ----------------------------

for (i in seq_len(nrow(manifest))) {

  row <- manifest[i, ]
  id <- row$dataset

  dir <- file.path(out_dir, id)
  ensure_dir(dir)

  message("▶ READS: ", id)

  # ====================================================
  # SRA / ERR
  # ====================================================

  if (!is.na(row$accession) &&
      (startsWith(row$accession, "SRR") ||
       startsWith(row$accession, "ERR"))) {

    ok_flag <- file.path(dir, ".done_sra")

    if (file.exists(ok_flag)) {
      message("↻ Already processed: ", id)
      next
    }

    status <- "SUCCESS"

    tryCatch({

      system2(
        "fasterq-dump",
        c(row$accession, "--split-files", "--outdir", dir)
      )

      files <- list.files(dir, full.names = TRUE)

      if (length(files) == 0) {
        stop("SRA produced no files")
      }

      sha_list <- list()

      for (f in files) {
        sha <- sha256_file(f)
        write_sha(dir, sha, basename(f))
        sha_list[[basename(f)]] <- sha
      }

      writeLines("done", ok_flag)

      write_ledger(list(
        dataset = id,
        type = "SRA",
        accession = row$accession,
        files = basename(files),
        sha256 = sha_list,
        status = "SUCCESS",
        timestamp = as.character(Sys.time())
      ), ledger_path)

    }, error = function(e) {

      status <<- "FAILED"

      write_ledger(list(
        dataset = id,
        type = "SRA",
        accession = row$accession,
        status = "FAILED",
        error = conditionMessage(e),
        timestamp = as.character(Sys.time())
      ), ledger_path)

      stop(e)
    })

  } else {

    # ====================================================
    # HTTP / FTP
    # ====================================================

    filename <- basename(row$source_url)
    dest <- file.path(dir, filename)

    ok_flag <- file.path(dir, ".done_http")

    if (file.exists(ok_flag)) {
      message("↻ Already processed: ", id)
      next
    }

    ok <- safe_download(row$source_url, dest)

    if (!ok) {
      write_ledger(list(
        dataset = id,
        type = "HTTP",
        url = row$source_url,
        status = "FAILED",
        timestamp = as.character(Sys.time())
      ), ledger_path)

      stop("Download failed: ", row$source_url)
    }

    sha <- sha256_file(dest)
    write_sha(dir, sha, filename)

    writeLines("done", ok_flag)

    write_ledger(list(
      dataset = id,
      type = "HTTP",
      url = row$source_url,
      file = filename,
      sha256 = sha,
      status = "SUCCESS",
      timestamp = as.character(Sys.time())
    ), ledger_path)
  }
}

message("✔ READS download complete")