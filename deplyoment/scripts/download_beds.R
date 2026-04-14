library(readr)
library(jsonlite)

# ----------------------------
# BOOTSTRAP
# ----------------------------

source(file.path(get_script_dir(), "init.R"))
source(file.path(get_script_dir(), "utils.R"))

# ----------------------------
# INPUTS / OUTPUTS
# ----------------------------

bed_manifest <- read_csv(file.path(PROJECT_ROOT, "giab_bed_manifest.csv"))

out_dir <- DATA_DIR
ledger_path <- file.path(out_dir, "_ledger.json")

ensure_dir(out_dir)

# ----------------------------
# LEDGER FUNCTION (same contract as reads)
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

for (i in seq_len(nrow(bed_manifest))) {

  row <- bed_manifest[i, ]

  id <- paste0(row$dataset, "_", row$bed_type)

  dir <- file.path(out_dir, id)
  ensure_dir(dir)

  message("▶ BEDS: ", id)

  filename <- basename(row$source_url)
  dest <- file.path(dir, filename)

  status <- "SUCCESS"

  tryCatch({

    ok <- safe_download(row$source_url, dest)
    if (!ok) stop("BED download failed: ", row$source_url)

    sha <- sha256_file(dest)
    write_sha(dir, sha, filename)

    write_ledger(list(
      dataset = row$dataset,
      bed_type = row$bed_type,
      id = id,
      url = row$source_url,
      file = filename,
      sha256 = sha,
      status = "SUCCESS",
      timestamp = as.character(Sys.time())
    ), ledger_path)

  }, error = function(e) {

    status <<- "FAILED"

    write_ledger(list(
      dataset = row$dataset,
      bed_type = row$bed_type,
      id = id,
      url = row$source_url,
      status = "FAILED",
      error = conditionMessage(e),
      timestamp = as.character(Sys.time())
    ), ledger_path)

    stop(e)
  })
}

message("✔ BED download complete")