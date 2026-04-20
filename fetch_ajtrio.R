library(tidyverse)
library(data.table)

urls <- c("https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_indexes/AshkenazimTrio/sequence.index.AJtrio_Illumina300X_wgs_07292015_updated.HG002",
  "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_indexes/AshkenazimTrio/sequence.index.AJtrio_Illumina300X_wgs_07292015_updated.HG003",
  "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_indexes/AshkenazimTrio/sequence.index.AJtrio_Illumina300X_wgs_07292015_updated.HG004"
  )

index <- map(urls, ~ fread(.x)) %>%
  reduce(bind_rows)

parse_giab_index <- function(index) {
  index |>
    # work from R1 only (paired info already in PAIRED_FASTQ column)
    mutate(
      
      # ── from path: run folder ──────────────────────────────
      run_folder  = str_extract(FASTQ, "\\d{6}_[^/]+"),
      run_date    = str_extract(run_folder, "^\\d{6}"),
      instrument  = str_extract(run_folder, "(?<=\\d{6}_)[^_]+"),
      run_number  = str_extract(run_folder, "(?<=_)\\d{4}(?=_)"),
      flowcell    = str_extract(run_folder, "[^_]+$"),
      
      # ── from path: sample folder ───────────────────────────
      sample_dir  = str_extract(FASTQ, "(?<=Sample_)[^/]+"),
      
      # ── from filename ──────────────────────────────────────
      filename    = basename(FASTQ),
      sample_name = str_extract(filename, "^[^_]+"),
      barcode     = str_extract(filename, "(?<=_)[A-Z]+(?=_L)"),
      lane        = str_extract(filename, "(?<=L00)\\d"),
      chunk       = str_extract(filename, "(?<=R[12]_)\\d+"),
      
      # ── RG tags ────────────────────────────────────────────
      rg_sm       = NIST_SAMPLE_NAME,
      rg_lb       = str_glue("{sample_name}_{barcode}"),
      rg_pl       = "ILLUMINA",
      rg_id       = str_glue("{instrument}.{run_number}.{flowcell}.{lane}"),
      rg_pu       = str_glue("{flowcell}.{lane}.{barcode}"),
      
      # ── BWA -R string ──────────────────────────────────────
      rg_string   = str_glue("@RG\\tID:{rg_id}\\tSM:{rg_sm}\\tLB:{rg_lb}\\tPL:{rg_pl}\\tPU:{rg_pu}")
    )
}

result <- parse_giab_index(index)

manifest <- data.frame(
  SM = result$rg_sm,
  ID = result$rg_id,
  LB = result$rg_lb,
  PL = result$rg_pl,
  PU = result$rg_pu,
  R1 = result$FASTQ,
  R2 = result$PAIRED_FASTQ
)

manifest %>%
  write.csv("deployment/manifests/aj_trio_samplesheet.csv",
            row.names = FALSE,
            quote = FALSE)





