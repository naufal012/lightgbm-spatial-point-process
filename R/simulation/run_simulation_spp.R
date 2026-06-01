# =============================================================
# SKRIP UTAMA: EKSEKUSI SIMULASI DAN ANALISIS (DENGAN 3 JENIS PROSES)
# =============================================================
Sys.setenv(RETICULATE_PYTHON = "C:/Users/Naufal/anaconda3/envs/lgbm-env/python.exe")
# --- 0.1 Setup Lingkungan ---
if (!require("reticulate")) install.packages("reticulate")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("spatstat")) install.packages("spatstat")
if (!require("openxlsx")) install.packages("openxlsx")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("viridis")) install.packages("viridis")

library(reticulate)
library(dplyr)
library(ggplot2)
library(spatstat)
library(openxlsx)
library(RColorBrewer)
library(viridis)

tryCatch({
  use_condaenv(
    condaenv = "lgbm-env", 
    conda = "C:/Users/Naufal/anaconda3/Scripts/conda.exe", # <-- TAMBAHKAN PATH INI
    required = TRUE
  )
  cat("--- Berhasil terhubung ke Conda environment 'lgbm-env' ---\n")
}, error = function(e) {
  stop("Gagal terhubung ke Conda. Pastikan path ke conda.exe dan nama environment sudah benar.")
})
tryCatch({ source_python("lgbpp.py"); source_python("xgbpp.py") }, error = function(e) { stop("Pastikan file Python ada.") })
cat("--- Mengimpor pustaka Python (XGBoost, LightGBM, Pandas) ---\n")
xgb <- reticulate::import("xgboost"); lgb <- reticulate::import("lightgbm"); pd <- reticulate::import("pandas")

# --- 1.1 Muat Fungsi-Fungsi R ---
source("strauss_pp.R") # Berisi simulate_poisson_process, simulate_thomas_process, simulate_lgcp_process
source("wrapper-sim-v2.R")  # Berisi run_analysis

cat("--- Berhasil memuat fungsi R dari file eksternal ---\n")

# --- 1.2 Fungsi Helper Training ---
train_lgbpp_fixed <- function(X, y, vol, loss, F_prime, base_params) {
  cat("--- Melatih LightGBM dengan parameter tetap...\n")
  final_params <- c(base_params, list(learning_rate = 0.001, lambda_l1 = 0, lambda_l2 = 0, num_leaves = 63L))
  train_set <- lgb$Dataset(data = as.matrix(X), label = pd$Series(y), feature_name = names(X))
  callbacks <- list(lgb$early_stopping(stopping_rounds = 50L, verbose = FALSE))
  model <- lgbpp_py(data = train_set, vol = pd$Series(vol), params = final_params, loss = loss, F_prime = F_prime, num_boost_round = 2000L)
  return(model)
}
train_xgbpp_fixed <- function(X, y, vol, loss, F_prime, base_params) {
  cat("--- Melatih XGBoost dengan parameter tetap...\n")
  final_params <- c(base_params, list(eta = 0.001, alpha = 0, lambda = 0, max_depth = 6L))
  dtrain <- xgb$DMatrix(data = as.matrix(X), label = pd$Series(y))
  dtrain$feature_names <- names(X)
  model <- xgbpp_py(dtrain = dtrain, vol = pd$Series(vol), params = final_params, loss = loss, F_prime = F_prime, num_boost_round = 2000L)
  return(model)
}

# =============================================================
# BAGIAN 2: EKSEKUSI UTAMA DENGAN LOOP
# =============================================================

# --- 2.1 Muat Data & Definisikan Parameter Global ---
cat("--- Memuat data BCI...\n")
load('dataset_bci/bci.covars.rda')
data(bei)
bci.covars$elev <- bei.extra$elev
win <- Window(bci.covars[[1]])
scale_factor <- 500
N_SIMULATIONS <- 1 
num_sim_points <- 2000

# --- 2.2 Pengaturan Skenario, Model, dan Parameter ---
# --- PERUBAHAN: Daftar Jenis Proses ---
process_scenarios <- list(
  list(name = "Poisson"),
  list(name = "Thomas"), # Contoh parameter Thomas
  list(name = "LGCP")      # Contoh parameter LGCP
)

intensity_form_scenarios <- list(
  list(form = "linear", n_cov = 8),
  list(form = "complex", n_cov = 8),
  list(form = "complex_sparse", n_cov = 2)
)

models_to_run <- list(
  list(type = "xgb", loss = "poisson", analysis = "fixed", discretization = "pois"),
  list(type = "xgb", loss = "logistic", analysis = "fixed", discretization = "logi_nd2"),
  list(type = "lgb", loss = "poisson", analysis = "fixed", discretization = "pois"),
  list(type = "lgb", loss = "logistic", analysis = "fixed", discretization = "logi_nd2")
)

base_params_xgb <- list(booster = 'gbtree', subsample = 1.0,
                        nthread = 1L, tree_method = 'exact', verbosity = 0L)

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

# --- 2.3 LOOP UTAMA PER JENIS PROSES (LOOP TERLUAR) ---
for (process_scenario in process_scenarios) {
  
  process_type <- process_scenario$name
  
  process_level_summary <- list() #NEW
  
  # --- LOOP KEDUA PER BENTUK INTENSITAS ---
  for (scenario in intensity_form_scenarios) {
    
    intensity_form <- scenario$form
    num_covariates <- scenario$n_cov
    
    cat(sprintf("\n\n========================================================\n"))
    cat(sprintf("===== MEMULAI SKENARIO: PROSES=%s, INTENSITAS=%s =====\n", toupper(process_type), toupper(intensity_form)))
    cat(sprintf("========================================================\n"))
    
    scenario_runs_summary <- list() 
    
    if (intensity_form == "linear") {
      sim_cov_names <- c("grad","elev","Cu","Nmin","P","pH","solar","twi")[1:num_covariates]
      sim_coefficients <- c(0.03, 0.03, -0.01, 0.02, 0.01, -0.01, 0.02, -0.01)[1:num_covariates]
    } else if (intensity_form == "complex") {
      sim_cov_names <- c("grad","elev","Cu","Nmin","P","pH","solar","twi")[1:num_covariates]
      sim_coefficients <- c(0.03, 0.03, -0.01, 0.02, 0.01, -0.01, 0.02, -0.01)[1:num_covariates]
    } else if (intensity_form == "complex_sparse") {
      sim_cov_names <- c("grad", "elev")[1:num_covariates]
      sim_coefficients <- c(0.03, 0.03)[1:num_covariates]
    }
    
    # --- Loop per run untuk skenario saat ini ---
    for (i in 1:N_SIMULATIONS) {
      graphics.off()
      cat(sprintf("\n--- [Proses: %s, Intensitas: %s] Memulai Run #%d ---\n", toupper(process_type), toupper(intensity_form), i))
      set.seed(i)
      
      # --- PERUBAHAN: PANGGIL FUNGSI SIMULASI SECARA DINAMIS ---
      cat(sprintf("\n--- [Siklus %d] Menjalankan simulasi point process jenis %s...\n", i, toupper(process_type)))
      if (process_type == "Poisson") {
        true_simulation <- simulate_poisson_process(
          covariate_names = sim_cov_names, coefficients = sim_coefficients, bci_covars = bci.covars, 
          bei_window = win, n_points = num_sim_points, scale_factor = scale_factor, intensity_form = intensity_form
        )
      } else if (process_type == "Thomas") {
        true_simulation <- simulate_thomas_process(
          covariate_names = sim_cov_names, coefficients = sim_coefficients, bci_covars = bci.covars, 
          bei_window = win, n_points = num_sim_points, scale_factor = scale_factor, intensity_form = intensity_form
        )
      } else if (process_type == "LGCP") {
        true_simulation <- simulate_lgc_process(
          covariate_names = sim_cov_names, coefficients = sim_coefficients, bci_covars = bci.covars, 
          bei_window = win, n_points = num_sim_points, scale_factor = scale_factor, intensity_form = intensity_form
        )
      }
      # --------------------------------------------------------
      
      for (model in models_to_run) {
        
        # 1. Determine Discretization Data
        if (model$discretization == "logi") {
          current_sim_data <- true_simulation$sim_data_full_logi
        } else if (model$discretization == "logi_nd2") {
          current_sim_data <- true_simulation$sim_data_full_logi_nd2
        } else { # "pois"
          current_sim_data <- true_simulation$sim_data_full_pois
        }
        
        features_cols <- setdiff(names(true_simulation$sim_data), c("x", "y", "label", "vol"))
        
        # 2. Define Output Paths
        base_folder_name <- file.path("simulation_results", process_type, intensity_form)
        run_base_dir <- file.path(base_folder_name, paste0("run_", i, "_", model$discretization))
        
        # 3. FIX: Select Parameters based on Model Type
        if (model$type == "lgb") {
          current_base_params <- base_params_lgb
        } else {
          current_base_params <- base_params_xgb
        }
        
        # 4. Run Analysis
        summary_df <- run_analysis(
          model_type = model$type, 
          loss_type = model$loss,
          sim_data = current_sim_data,
          sim_intensity = true_simulation$intensity,
          sim_points = true_simulation$sim_points,
          base_output_dir = run_base_dir, 
          run_number = i,
          base_params = current_base_params, # Now correctly points to lgb or xgb params
          analysis_type = model$analysis,
          scale_factor = scale_factor
        )
        
        summary_df$simulation_run <- i
        summary_df$model <- model$type
        summary_df$loss <- paste0(model$loss, "_", model$discretization)
        
        scenario_runs_summary[[length(scenario_runs_summary) + 1]] <- summary_df
      }
    }
    
    # --- Simpan ringkasan untuk kombinasi skenario saat ini ---
    if (length(scenario_runs_summary) > 0) {
      final_summary <- bind_rows(scenario_runs_summary) %>%
        tidyr::pivot_wider(names_from = Parameter, values_from = Value) %>%
        select(simulation_run, model, loss, `True Log-Likelihood`, `Scaled True Log-Likelihood`, `Left Scaled True Log-Likelihood`, `Right Scaled True Log-Likelihood`,
               `True Logistic Log-Likelihood`, `Left True Logistic Log-Likelihood`, `Right True Logistic Log-Likelihood`, `Log-Likelihood`,
               `Scaled Log-Likelihood`, `Left Scaled Log-Likelihood`, `Right Scaled Log-Likelihood`, `Left Logistic Log-Likelihood`,
               `Right Logistic Log-Likelihood`,`Logistic Log-Likelihood`, `Predicted Events`, `True Events`, `MISE`, `Scaled MISE`,
               `Computation Time (mins)`)
      
      # --- PERUBAHAN: Sesuaikan path penyimpanan ringkasan ---
      summary_folder <- file.path("simulation_results", process_type, intensity_form)
      if (!dir.exists(summary_folder)) dir.create(summary_folder, recursive = TRUE)
      summary_file_path <- file.path(summary_folder, "all_run_summary.xlsx")
      # ----------------------------------------------------
      
      wb <- createWorkbook()
      addWorksheet(wb, "Summary_Runs")
      writeData(wb, "Summary_Runs", final_summary)
      saveWorkbook(wb, summary_file_path, overwrite = TRUE)
      cat(sprintf("\n--- Laporan ringkasan untuk PROSES=%s, INTENSITAS=%s disimpan di: %s ---\n", toupper(process_type), toupper(intensity_form), summary_file_path))
      
      # 3. [BARU] Hitung Rata-Rata untuk Tabel Master
      # Kita menggunakan 'final_summary' yang baru saja dibuat di atas
      avg_summary <- final_summary %>%
        group_by(model, loss) %>%
        summarise(
          Avg_LogLik = mean(`Scaled Log-Likelihood`, na.rm = TRUE),
          Avg_MISE = mean(`Scaled MISE`, na.rm = TRUE),
          Avg_Time_Mins = mean(`Computation Time (mins)`, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(
          Intensity = intensity_form, 
          Model_Variant = paste0(toupper(model), "_", loss)
        ) %>%
        select(Intensity, Model_Variant, Avg_LogLik, Avg_MISE, Avg_Time_Mins)
      
      # Masukkan ke list akumulasi untuk digabungkan nanti
      process_level_summary[[length(process_level_summary) + 1]] <- avg_summary
    }
  } # Akhir loop intensity_form
    # --- [BARU] SIMPAN TABEL PERBANDINGAN FINAL PER PROSES ---
    # Setelah loop Intensity selesai (Linear, Complex, Sparse selesai semua),
    # Kita gabungkan hasilnya menjadi satu tabel besar mirip gambar referensi Anda.
    
    if (length(process_level_summary) > 0) {
      master_table <- bind_rows(process_level_summary)
      
      # Format agar rapi (opsional: scientific notation untuk MISE)
      master_table$Avg_MISE_Sci <- formatC(master_table$Avg_MISE, format = "e", digits = 2)
      
      output_master_file <- file.path("simulation_results", process_type, paste0("FINAL_SUMMARY_", process_type, ".xlsx"))
      
      wb <- createWorkbook()
      addWorksheet(wb, "Performance_Comparison")
      writeData(wb, "Performance_Comparison", master_table)
      saveWorkbook(wb, output_master_file, overwrite = TRUE)
      
      cat(sprintf("\n>>> TABEL PERBANDINGAN FINAL (%s) DISIMPAN DI: %s <<<\n", toupper(process_type), output_master_file))
      print(master_table)
    }
} # Akhir loop process_type

cat("\n--- SEMUA PROSES SIMULASI SELESAI. ---\n")