# 01_download_brfss2010.R
# ---------------------------------------------------------------------------
# Завантажує BRFSS 2010 (XPT format), розпаковує, читає в R, виконує базові
# перевірки і зберігає processed PHQ-8 subset як .rds для подальшої роботи.
#
# Запуск з кореня проєкту:
#   Rscript scripts/01_download_brfss2010.R
#
# Залежності: haven, here (встановляться автоматично якщо не знайдено)
# Час: ~3-5 хв на пристойному інтернеті (99 MB download)
# Дискове місце: ~250 MB (zip + extracted XPT + processed rds)
# ---------------------------------------------------------------------------

# 0. Setup --------------------------------------------------------------------
required_pkgs <- c("haven", "here")
new_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs) > 0) {
  message("Installing missing packages: ", paste(new_pkgs, collapse = ", "))
  install.packages(new_pkgs, repos = "https://cloud.r-project.org")
}
library(haven)
library(here)

# 1. Paths --------------------------------------------------------------------
data_dir <- here("source-documents", "data")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

zip_path <- file.path(data_dir, "CDBRFS10XPT.zip")
xpt_path <- file.path(data_dir, "CDBRFS10.XPT")  # canonical name after unzip
rds_path <- file.path(data_dir, "brfss2010_full.rds")
phq_rds_path <- file.path(data_dir, "brfss2010_phq8_highrisk.rds")

# 2. Download ----------------------------------------------------------------
url <- "https://www.cdc.gov/brfss/annual_data/2010/files/CDBRFS10XPT.zip"

if (!file.exists(zip_path)) {
  message("Downloading BRFSS 2010 XPT (~99 MB) from CDC...")
  options(timeout = 600)  # 10 min timeout
  download.file(url, zip_path, mode = "wb", quiet = FALSE)
  message("Download complete: ", round(file.info(zip_path)$size / 1024^2, 1), " MB")
} else {
  message("Zip already exists, skipping download: ", zip_path)
}

# 3. Unzip --------------------------------------------------------------------
if (!file.exists(xpt_path)) {
  message("Unzipping...")
  unzipped <- unzip(zip_path, exdir = data_dir)
  # CDC may use different extension casing — find actual file
  xpt_files <- list.files(data_dir, pattern = "\\.(XPT|xpt)$", full.names = TRUE)
  if (length(xpt_files) == 0) {
    stop("XPT file not found after unzip. Contents: ", paste(unzipped, collapse = ", "))
  }
  if (xpt_files[1] != xpt_path) {
    file.rename(xpt_files[1], xpt_path)
  }
  message("Unzipped: ", xpt_path, " (",
          round(file.info(xpt_path)$size / 1024^2, 1), " MB)")
} else {
  message("XPT already extracted, skipping unzip: ", xpt_path)
}

# 4. Load ---------------------------------------------------------------------
message("Reading XPT into R (this may take 30-60 s)...")
t0 <- Sys.time()
brfss <- haven::read_xpt(xpt_path)
t1 <- Sys.time()
message("Loaded in ", round(as.numeric(t1 - t0, units = "secs"), 1), " s")

# 5. Sanity checks ------------------------------------------------------------
cat("\n=== BRFSS 2010 sanity check ===\n")
cat("Rows:        ", nrow(brfss), "  (expected ~451075)\n")
cat("Cols:        ", ncol(brfss), "  (expected ~397)\n")
cat("Object size: ", format(object.size(brfss), units = "MB"), "\n")
cat("\nFirst 10 variable names:\n")
print(head(names(brfss), 10))

# 6. PHQ-8 змінні (звірити з codebook 2010!) ---------------------------------
# Очікувані імена за CDC BRFSS 2010 Anxiety and Depression module
# Якщо назви не співпадуть — відкрити codebook_10.pdf і виправити
phq_vars_candidates <- c(
  "ADPLEASR",  # (a) Little pleasure
  "ADDOWN",    # (b) Feeling down
  "ADSLEEP",   # (c) Sleep
  "ADENERGY",  # (d) Energy
  "ADEAT1",    # (e) Eating
  "ADFAIL",    # (f) Failure — "DAYS FELT LIKE FAILURE OR LET FAMILY DOWN"
  "ADTHINK",   # (g) Focus
  "ADMOVE"     # (h) Moving
)

phq_present <- phq_vars_candidates[phq_vars_candidates %in% names(brfss)]
phq_missing <- setdiff(phq_vars_candidates, phq_present)

cat("\n=== PHQ-8 variables check ===\n")
cat("Found    (", length(phq_present), "/8): ", paste(phq_present, collapse = ", "), "\n", sep = "")
if (length(phq_missing) > 0) {
  cat("MISSING  (", length(phq_missing), "/8): ", paste(phq_missing, collapse = ", "), "\n", sep = "")
  cat("\n>>> ACTION: відкрий", file.path(data_dir, "codebook_10.pdf"),
      "(окреме скачування) і знайди точні імена для відсутніх пунктів.\n",
      "Орієнтир: BRFSS 2010 Anxiety and Depression Module (стани AR, CO, DC, HI, LA, MD, MS, MT, NH, NJ, NV, NY, OK, PA, TN, VT, WV).\n")
}

if (length(phq_present) >= 6) {
  cat("\n=== Building high-risk subset ===\n")
  # У BRFSS PHQ-8 елементи кодуються:
  # 0-14 = number of days, 88 = none, 77 = don't know, 99 = refused
  # Перекодуємо 88→0, 77/99→NA
  phq_data <- brfss[, phq_present, drop = FALSE]
  phq_data <- as.data.frame(lapply(phq_data, function(x) {
    x <- as.numeric(x)
    x[x == 88] <- 0
    x[x == 77 | x == 99] <- NA
    x[x > 14] <- NA
    x
  }))

  # Sum score, тільки complete cases
  phq_data$sum_score <- rowSums(phq_data[, phq_present], na.rm = FALSE)
  complete_cases <- !is.na(phq_data$sum_score)
  cat("Complete PHQ-8 responses: ", sum(complete_cases), "\n")

  # High-risk per W&S p.33: "all eight symptoms occurred at least once
  # within the past 14 days (i.e., a sum score >= 8)"
  # KEY: парентеза оманлива — основна умова "all >= 1", не просто sum >= 8.
  # Сума ≥ 8 досягається тільки якщо всі 8 ≥ 1 (8 items × min 1 = 8).
  all_present <- complete_cases & apply(phq_data[, phq_present], 1,
                                         function(x) all(x >= 1))
  highrisk <- phq_data[all_present, ]
  cat("High-risk subset (all 8 symptoms >= 1): ", nrow(highrisk),
      "  (W&S reported N=2136)\n")

  # Зберегти full і high-risk
  saveRDS(brfss, rds_path)
  saveRDS(highrisk, phq_rds_path)
  cat("\nSaved:\n")
  cat("  Full BRFSS 2010:  ", rds_path, "\n")
  cat("  PHQ-8 high-risk:  ", phq_rds_path, "\n")

  # Розподільні характеристики (звірка з W&S Table 5)
  cat("\n=== Descriptive stats (compare with W&S Table 5) ===\n")
  desc <- sapply(highrisk[, phq_present], function(x) {
    c(M = mean(x, na.rm = TRUE), SD = sd(x, na.rm = TRUE),
      skew = mean((x - mean(x, na.rm = TRUE))^3, na.rm = TRUE) /
             sd(x, na.rm = TRUE)^3,
      exkurt = mean((x - mean(x, na.rm = TRUE))^4, na.rm = TRUE) /
               sd(x, na.rm = TRUE)^4 - 3)
  })
  print(round(desc, 2))
  cat("\nW&S reported: M ~7.3-10.2, SD ~4.3-4.9, skew ~[-0.66, 0.30], exkurt ~[-1.58, -1.06]\n")
} else {
  cat("\n>>> Зупиняюсь — PHQ-8 змінні не знайдено. Виправ імена після перегляду codebook.\n")
}

cat("\n=== DONE ===\n")
