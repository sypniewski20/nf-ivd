library(digest)
library(jsonlite)

sha256_file <- function(path) {
  digest(path, algo = "sha256", file = TRUE)
}

write_sha <- function(dir, sha, file) {
  line <- paste0(sha, "  ", file)
  write(line, file = file.path(dir, "SHA256SUMS.txt"), append = TRUE)
}

safe_download <- function(url, dest) {

  if (file.exists(dest)) {
    message("↻ Exists, skipping: ", basename(dest))
    return(TRUE)
  }

  tmp <- paste0(dest, ".part")

  ok <- tryCatch({
    download.file(url, tmp, mode = "wb", quiet = TRUE)
    file.rename(tmp, dest)
    TRUE
  }, error = function(e) {
    file.remove(tmp)
    FALSE
  })

  return(ok)
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

write_ledger <- function(record, path) {

  if (!file.exists(path)) {
    writeLines("[]", path)
  }

  ledger <- jsonlite::fromJSON(path)

  ledger <- append(ledger, list(record))

  jsonlite::write_json(ledger, path, pretty = TRUE, auto_unbox = TRUE)
}