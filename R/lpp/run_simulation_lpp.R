# =============================================================
# SKRIP UTAMA: EKSEKUSI SIMULASI LPP (LINEAR POINT PROCESS)
# =============================================================
Sys.setenv(RETICULATE_PYTHON = "C:/Users/Naufal/anaconda3/envs/lgbm-env/python.exe")

# --- 0.1 Setup Lingkungan ---
if (!require("reticulate")) install.packages("reticulate")
if (!require("dplyr")) install.packages("dplyr")
if (!require("spatstat")) install.packages("spatstat")
if (!require("spatstat.linnet")) install.packages("spatstat.linnet")
if (!require("openxlsx")) install.packages("openxlsx")
if (!require("viridis")) install.packages("viridis")
if (!require("sf")) install.packages("sf")
if (!require("sp")) install.packages("sp")

library(reticulate)
library(dplyr)
library(ggplot2)
library(spatstat)
library(spatstat.linnet)
library(spatstat.geom)
library(openxlsx)
library(viridis)
library(sf)
library(sp)

# --- 0.2 Koneksi Python & Impor Pustaka ---
tryCatch({
  use_condaenv(
    condaenv = "lgbm-env", 
    conda = "C:/Users/Naufal/anaconda3/Scripts/conda.exe",
    required = TRUE
  )
  cat("--- Berhasil terhubung ke Conda environment 'lgbm-env' ---\n")
}, error = function(e) {
  stop("Gagal terhubung ke Conda. Pastikan path ke conda.exe dan nama environment sudah benar.")
})

tryCatch({ 
  source_python("lgbpp.py"); 
  source_python("xgbpp.py") 
}, error = function(e) { stop("Pastikan file Python (lgbpp/xgbpp) ada.") })

cat("--- Mengimpor pustaka Python (XGBoost, LightGBM, Pandas) ---\n")
xgb <- reticulate::import("xgboost"); lgb <- reticulate::import("lightgbm"); pd <- reticulate::import("pandas")

# --- 1.1 Muat Fungsi-Fungsi R ---
source("sim-lpp.R") 
source("wrapper-sim-lpp.R") 

cat("--- Berhasil memuat fungsi R (simulasi LPP & analisis LPP) ---\n")

# --- 1.2 Fungsi Helper Training ---
train_lgbpp_fixed <- function(X, y, vol, loss, F_prime, base_params) {
  cat("--- Melatih LightGBM dengan parameter tetap...\n")
  # Pastikan parameter integer ditandai dengan L (misal 63L)
  final_params <- c(base_params, list(learning_rate = 0.001, lambda_l1 = 0, lambda_l2 = 50, num_leaves = 63L))
  
  train_set <- lgb$Dataset(data = as.matrix(X), label = pd$Series(y))
  callbacks <- list(lgb$early_stopping(stopping_rounds = 50L, verbose = FALSE))
  
  model <- lgbpp_py(data = train_set, vol = pd$Series(vol), params = final_params, loss = loss, F_prime = F_prime, num_boost_round = 2000L)
  return(model)
}

train_xgbpp_fixed <- function(X, y, vol, loss, F_prime, base_params) {
  cat("--- Melatih XGBoost dengan parameter tetap...\n")
  final_params <- c(base_params, list(eta = 0.001, alpha = 0, lambda = 50, max_depth = 6L))
  
  dtrain <- xgb$DMatrix(data = as.matrix(X), label = pd$Series(y))
  model <- xgbpp_py(dtrain = dtrain, vol = pd$Series(vol), params = final_params, loss = loss, F_prime = F_prime, num_boost_round = 2000L, verbose_eval = FALSE)
  return(model)
}

# =============================================================
# BAGIAN 2: EKSEKUSI UTAMA LPP
# =============================================================

# --- 2.1 Muat Data Jaringan & Definisikan Parameter Global ---
cat("--- Memuat data Jaringan (Nganjuk LPP)...\n")
nganjuk_ln <- readRDS('Data Analisis Tugas Akhir/nganjuk_linnet_rescaled.rds')
L <- nganjuk_ln 
datatl <- read.csv("Data Analisis Tugas Akhir/trafficlights.csv", header = TRUE, sep = ";")
coordinates(datatl) <- c("Long", "Lat")
proj4string(datatl) <- CRS("+proj=longlat +datum=WGS84")
datatl.UTM <- spTransform(datatl, CRS("+proj=utm +zone=49 ellps=WGS84"))
koor_tl <- datatl.UTM@coords
tl_rescaled <- 0.1 * koor_tl
tl_lpp <- lpp(tl_rescaled, L)
f_dist_tl <- distfun.lpp(tl_lpp) 

linmarks <- marks(nganjuk_ln$lines)
linmarks$Jenis_Jalan <- as.numeric(factor(linmarks$Jenis_Jalan, levels = c("Jalan Lokal", "Jalan Kolektor", "Jalan Arteri")))
linmarks$Kerb_L <- as.numeric(factor(ifelse(linmarks$Kerb_L == 0, "Tidak ada", "Ada"), levels = c("Tidak ada", "Ada")))
linmarks$Kerb_R <- as.numeric(factor(ifelse(linmarks$Kerb_R == 0, "Tidak ada", "Ada"), levels = c("Tidak ada", "Ada")))
linmarks$n_Lajur <- as.numeric(as.character(linmarks$n_Lajur)); linmarks$n_Lajur[is.na(linmarks$n_Lajur)] <- 0
linmarks$n_Jalur <- as.numeric(as.character(linmarks$n_Jalur)); linmarks$n_Jalur[is.na(linmarks$n_Jalur)] <- 0

covariate_names <- c("Jenis_Jalan", "n_Lajur", "n_Jalur", "Kerb_R", "Kerb_L")
linmarks_clean <- linmarks[, covariate_names]
funcListScaled <- lapply(linmarks_clean, function(z) { function(x,y,seg,tp) { z[seg] } })
linfunList <- lapply(funcListScaled, function(z, net) linfun(z, net), net = L)

# --- SETUP KOVARIAT DAN KOEFISIEN ---
ALL_LINFUNS_MASTER <- c(linfunList, list(Jarak_TL = f_dist_tl))
ALL_COEFFS_MASTER <- c(
  n_Lajur     = 0.9,
  n_Jalur     = 0.4,
  Kerb_R      = 0.2,
  Kerb_L      = 0.2,
  Jarak_TL    = -0.8
)

COVARIATES_TO_USE <- c("n_Lajur", "n_Jalur", "Jarak_TL")

all_covariates_LPP <- ALL_LINFUNS_MASTER[COVARIATES_TO_USE]
my_coeffs_LPP      <- ALL_COEFFS_MASTER[COVARIATES_TO_USE]

if (length(all_covariates_LPP) != length(my_coeffs_LPP)) stop("Error: Nama di COVARIATES_TO_USE tidak cocok.")
if (anyNA(all_covariates_LPP) || anyNA(my_coeffs_LPP)) stop("Error: Salah satu nama tidak ada di database master.")

cat(sprintf("--- Menggunakan %d kovariat untuk simulasi: %s ---\n", 
            length(COVARIATES_TO_USE), 
            paste(COVARIATES_TO_USE, collapse = ", ")))

N_SIMULATIONS <- 5 
num_sim_points <- 2000

# --- 2.2 Pengaturan Skenario, Model, dan Parameter ---
intensity_form_scenarios <- list(
  list(form = "linear"),
  list(form = "complex"),
  list(form = "complex_sparse")
)

models_to_run <- list(
  list(type = "xgb", loss = "poisson", analysis = "fixed", discretization = "pois"),
  list(type = "xgb", loss = "logistic", analysis = "fixed", discretization = "logi_nd2"),
  list(type = "lgb", loss = "poisson", analysis = "fixed", discretization = "pois"),
  list(type = "lgb", loss = "logistic", analysis = "fixed", discretization = "logi_nd2") 
)

# Parameter Dasar XGBoost
base_params_xgb <- list(booster = 'gbtree', subsample = 1.0, colsample_bytree = 1/3,
                        nthread = 1L, tree_method = 'exact', verbosity = 0L)

# Parameter Dasar LightGBM
base_params_lgb <-  list(boosting_type = 'goss',
                         top_rate = 0.2,
                         other_rate = 0.2,
                         is_enable_bundle = TRUE,
                         max_conflict_rate = 0.1,
                         max_bin = 255L,
                         bin_construct_sample_cnt = 200000L,
                         min_data_in_bin = 3L,
                         deterministic = TRUE,
                         verbose = -1L,
                         num_threads = 1L)

# --- 2.3 LOOP UTAMA PER BENTUK INTENSITAS (DIPERBAIKI) ---
for (scenario in intensity_form_scenarios) {
  
  intensity_form <- scenario$form
  
  cat(sprintf("\n\n========================================================\n"))
  cat(sprintf("===== MEMULAI SKENARIO LPP: INTENSITAS=%s =====\n", toupper(intensity_form)))
  cat(sprintf("========================================================\n"))
  
  scenario_runs_summary <- list() 
  
  for (i in 1:N_SIMULATIONS) {
    graphics.off()
    cat(sprintf("\n--- [Intensitas: %s] Memulai Run LPP #%d ---\n", toupper(intensity_form), i))
    set.seed(i)
    
    cat(sprintf("\n--- [Siklus %d] Menjalankan simulasi LPP...\n", i))
    true_simulation <- simulate_LPP_process(
      L = L,
      covariates_list = all_covariates_LPP,
      coefficients = my_coeffs_LPP,
      n_points = num_sim_points,
      intensity_form = intensity_form,
      n_pix = 128
    )
    
    for (model in models_to_run) {
      
      # 1. Pilih Data Discretization
      if (grepl("logi", model$discretization)) {
        current_sim_data <- true_simulation$sim_data_full_logi
      } else { # "pois"
        current_sim_data <- true_simulation$sim_data_full
      }
      
      base_folder_name <- file.path("simulation_results_LPP", intensity_form)
      
      # === PERBAIKAN DI SINI: Tambahkan model$type ke nama folder ===
      # Hasilnya akan menjadi: "run_1_xgb_pois" atau "run_1_lgb_pois"
      run_base_dir <- file.path(base_folder_name, paste0("run_", i, "_", model$type, "_", model$discretization))
      
      # 2. Pilih Parameter Berdasarkan Jenis Model
      if (model$type == "lgb") {
        current_base_params <- base_params_lgb
      } else {
        current_base_params <- base_params_xgb
      }
      
      # 3. Jalankan Analisis
      summary_df <- run_analysis_LPP(
        L = L, 
        model_type = model$type, 
        loss_type = model$loss, 
        sim_data = current_sim_data,
        sim_intensity_linim = true_simulation$intensity_linim_1D,
        sim_points_lpp = true_simulation$sim_points_lpp,
        base_output_dir = run_base_dir, 
        run_number = i,
        base_params = current_base_params 
      )
      
      summary_df$simulation_run <- i
      summary_df$model <- model$type
      summary_df$loss <- model$loss 
      summary_df$discretization <- model$discretization
      
      scenario_runs_summary[[length(scenario_runs_summary) + 1]] <- summary_df
    }
  } # Akhir loop N_SIMULATIONS
  
  # --- Simpan ringkasan ---
  if (length(scenario_runs_summary) > 0) {
    final_summary <- bind_rows(scenario_runs_summary) %>%
      tidyr::pivot_wider(names_from = Parameter, values_from = Value) %>%
      select(simulation_run, model, discretization, loss, 
             `left Log-Likelihood`,
             `Right Log-Likelihood`,
             `Log-Likelihood`,
             `Left Logistic Log-Lik`,
             `Right Logistic Log-Lik`,
             `Logistic Log-Lik`,
             `Predicted Events`, `True Events`, `MISE (Network)`,
             `Computation Time (mins)`)
    
    summary_folder <- file.path("simulation_results_LPP", intensity_form)
    if (!dir.exists(summary_folder)) dir.create(summary_folder, recursive = TRUE)
    summary_file_path <- file.path(summary_folder, "all_run_summary_LPP.xlsx")
    
    wb <- createWorkbook()
    addWorksheet(wb, "Summary_Runs_LPP")
    writeData(wb, "Summary_Runs_LPP", final_summary)
    saveWorkbook(wb, summary_file_path, overwrite = TRUE)
    cat(sprintf("\n--- Laporan ringkasan LPP untuk INTENSITAS=%s disimpan di: %s ---\n", toupper(intensity_form), summary_file_path))
  }
} # Akhir loop intensity_form

cat("\n--- SEMUA PROSES SIMULASI LPP SELESAI. ---\n")