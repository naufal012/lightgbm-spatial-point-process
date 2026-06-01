# ================================
# 0. LIBRARIES
# ================================
library(dplyr)
library(readr)
library(openxlsx)
library(sp)
library(spatstat)
library(spatstat.linnet)
library(lubridate)
library(reticulate)
library(ggplot2)

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
# ================================
# 1. LOAD & CLEAN DATAALL (FINAL)
# ================================
dataall <- read_delim(
  "Data Analisis Tugas Akhir/dataall.csv",
  delim = ";",
  locale = locale(decimal_mark = ".")
)

# Fix column name
colnames(dataall)[1] <- "TanggalWaktu"

# ================================
# CLEAN DATETIME (FULL ROBUST)
# ================================
dataall$TanggalWaktu_clean <- dataall$TanggalWaktu

# Remove day names
dataall$TanggalWaktu_clean <- gsub(
  "Senin|Selasa|Rabu|Kamis|Jumat|Sabtu|Minggu",
  "",
  dataall$TanggalWaktu_clean
)

# Remove extra text
dataall$TanggalWaktu_clean <- gsub("Skj\\s*:\\s*", "", dataall$TanggalWaktu_clean)

# Remove "Jam" and "Wib"
dataall$TanggalWaktu_clean <- gsub("Jam\\.?", "", dataall$TanggalWaktu_clean)
dataall$TanggalWaktu_clean <- gsub("Wib|WIB", "", dataall$TanggalWaktu_clean)

# Fix typos
dataall$TanggalWaktu_clean <- gsub("Nopember", "November", dataall$TanggalWaktu_clean)
dataall$TanggalWaktu_clean <- gsub("Nopemeber", "November", dataall$TanggalWaktu_clean)

# Convert Indonesian full month → English
bulan_id <- c("Januari","Februari","Maret","April","Mei","Juni",
              "Juli","Agustus","September","Oktober","November","Desember")

bulan_en <- c("January","February","March","April","May","June",
              "July","August","September","October","November","December")

for(i in seq_along(bulan_id)){
  dataall$TanggalWaktu_clean <- gsub(
    paste0("(?i)", bulan_id[i]), 
    bulan_en[i], 
    dataall$TanggalWaktu_clean, 
    perl = TRUE
  )
}

# Convert Indonesian short month → English
bulan_short_id <- c("Jan","Feb","Mar","Apr","Mei","Jun","Jul","Agt","Sep","Okt","Nov","Des")
bulan_short_en <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

for(i in seq_along(bulan_short_id)){
  dataall$TanggalWaktu_clean <- gsub(
    paste0("(?i)", bulan_short_id[i]),
    bulan_short_en[i],
    dataall$TanggalWaktu_clean,
    perl = TRUE
  )
}

# Convert ALL dots to colon (time fix)
dataall$TanggalWaktu_clean <- gsub("\\.", ":", dataall$TanggalWaktu_clean)

# Remove commas and normalize spaces
dataall$TanggalWaktu_clean <- gsub(",", "", dataall$TanggalWaktu_clean)
dataall$TanggalWaktu_clean <- gsub("\\s+", " ", dataall$TanggalWaktu_clean)
dataall$TanggalWaktu_clean <- trimws(dataall$TanggalWaktu_clean)

# ================================
# PARSE DATETIME (FULL FLEXIBLE)
# ================================
dataall$datetime <- parse_date_time(
  dataall$TanggalWaktu_clean,
  orders = c("dmy HMS", "dmy HM", "dmy", "d-b-y")
)

# ================================
# CHECK PARSING
# ================================
cat("Unparsed rows:", sum(is.na(dataall$datetime)), "\n")

# Inspect remaining failures (if any)
failed <- dataall$TanggalWaktu_clean[is.na(dataall$datetime)]
head(failed, 10)

# ================================
# FIX NUMERIC
# ================================
dataall$JarakLampu <- as.numeric(dataall$JarakLampu)

# ================================
# FINAL STRUCTURE
# ================================
dataall <- dataall %>%
  mutate(
    x = Long * 0.1,
    y = Lat * 0.1
  ) %>%
  select(x, y, JarakLampu, Tahun, datetime)

# ================================
# 2. LOAD NETWORK & TRAFFIC LIGHT
# ================================
nganjuk_ln <- readRDS('Data Analisis Tugas Akhir/nganjuk_linnet_rescaled.rds')
L <- nganjuk_ln 

datatl <- read.csv("Data Analisis Tugas Akhir/trafficlights.csv", sep = ";")

coordinates(datatl) <- c("Long", "Lat")
proj4string(datatl) <- CRS("+proj=longlat +datum=WGS84")

datatl_utm <- spTransform(datatl, CRS("+proj=utm +zone=49 ellps=WGS84"))
tl_rescaled <- 0.1 * datatl_utm@coords

tl_lpp <- lpp(tl_rescaled, L)
f_dist_tl <- distfun.lpp(tl_lpp)

# ================================
# 3. FUNCTION FOR YEARLY DATA
# ================================
process_laka_year <- function(year, L, f_dist_tl) {
  
  file_path <- paste0("Data Analisis Tugas Akhir/LOKASI LAKA THN ", year, ".xlsx")
  
  cat("\nProcessing year:", year, "\n")
  
  laka_raw <- openxlsx::read.xlsx(file_path)
  
  print(colnames(laka_raw))
  
  # ================================
  # CLEAN & SELECT
  # ================================
  laka <- laka_raw %>%
    transmute(
      Lat  = as.numeric(`Koordinat.GPS.-.Lintang`),
      Long = as.numeric(`Koordinat.GPS.-.Bujur`),
      datetime = as.POSIXct(
        paste(`Tanggal.Kejadian`, `Waktu.Kejadian`),
        format = "%Y-%m-%d %H:%M:%S"
      ),
      Tahun = year
    ) %>%
    filter(!is.na(Lat), !is.na(Long))
  
  # ================================
  # SPATIAL
  # ================================
  coordinates(laka) <- c("Long", "Lat")
  proj4string(laka) <- CRS("+proj=longlat +datum=WGS84")
  
  laka_utm <- spTransform(laka, CRS("+proj=utm +zone=49 ellps=WGS84"))
  coords <- laka_utm@coords
  coords_rescaled <- 0.1 * coords
  
  # Build LPP
  laka_lpp <- lpp(coords_rescaled, L)
  
  # ================================
  # COMPUTE DISTANCE (SAFE)
  # ================================
  dist_tl <- f_dist_tl(laka_lpp)
  
  # ================================
  # ALIGN LENGTHS (KEY FIX)
  # ================================
  n <- min(nrow(coords_rescaled), length(dist_tl), length(laka$datetime))
  
  df <- data.frame(
    x = coords_rescaled[1:n,1],
    y = coords_rescaled[1:n,2],
    JarakLampu = dist_tl[1:n],
    Tahun = year,
    datetime = laka$datetime[1:n]
  )
  
  return(df)
}

# ================================
# 4. PROCESS ALL YEARS
# ================================
years <- 2022:2025

laka_all <- bind_rows(
  lapply(years, process_laka_year, L = L, f_dist_tl = f_dist_tl)
)

# ================================
# 5. COMBINE WITH DATAALL
# ================================
dataall_final <- bind_rows(dataall, laka_all)

# ================================
# 6. FINAL CHECK
# ================================
str(dataall_final)
summary(dataall_final)
head(dataall_final)


# ================================
# 1. LOAD NETWORK
# ================================
nganjuk_ln <- readRDS('Data Analisis Tugas Akhir/nganjuk_linnet_rescaled.rds')
ln_nganjuk <- nganjuk_ln[[1]]

# ================================
# 2. PREPARE EVENT DATA (FINAL)
# ================================
data_clean <- dataall_final %>%
  dplyr::select(x, y)

# Build LPP object
kcl_point <- lpp(data_clean, nganjuk_ln)

# ================================
# 3. TRAFFIC LIGHT DATA (COVARIATE)
# ================================
datatl <- read.csv("Data Analisis Tugas Akhir/trafficlights.csv", sep = ";")

coordinates(datatl) <- c("Long","Lat")
proj4string(datatl) <- CRS("+proj=longlat +datum=WGS84")

datatl_utm <- spTransform(datatl, CRS("+proj=utm +zone=49 ellps=WGS84"))
koor_tl <- datatl_utm@coords
tl_rescaled <- 0.1 * koor_tl

tl_lpp <- lpp(tl_rescaled, nganjuk_ln)

# Distance function
f <- distfun.lpp(tl_lpp)

# ================================
# 4. NETWORK COVARIATES
# ================================
linmarks <- marks(nganjuk_ln$lines)

linmarks <- linmarks %>% 
  mutate(
    Kerb_R = case_when(
      Kerb_R == 0 ~ "Tidak ada",
      Kerb_R == 1 ~ "Ada",
      TRUE ~ as.character(Kerb_R)
    ),
    Kerb_L = case_when(
      Kerb_L == 0 ~ "Tidak ada",
      Kerb_L == 1 ~ "Ada",
      TRUE ~ as.character(Kerb_L)
    )
  )

linmarks$Jenis_Jalan <- factor(linmarks$Jenis_Jalan,
                               levels = c("Jalan Lokal","Jalan Kolektor","Jalan Arteri"))

linmarks$Kerb_L <- factor(linmarks$Kerb_L, levels = c("Tidak ada","Ada"))
linmarks$Kerb_R <- factor(linmarks$Kerb_R, levels = c("Tidak ada","Ada"))

# Convert to linfun
funcListScaled <- lapply(linmarks, function(z){function(x,y,seg,tp){z[seg]}})
linfunList <- lapply(funcListScaled, function(z, net){linfun(z, net)}, net=nganjuk_ln)

Jenis_Jalan <- linfunList$Jenis_Jalan
n_Lajur     <- linfunList$n_Lajur
n_Jalur     <- linfunList$n_Jalur
Kerb_R      <- linfunList$Kerb_R
Kerb_L      <- linfunList$Kerb_L

# ================================
# 5. QUADRATURE SCHEME
# ================================
Q_lin <- spatstat.linnet::linequad(kcl_point)

# Data points
data_df <- as.data.frame(Q_lin$data)
data_df$type <- "data"
data_df$weight <- Q_lin$weights$data

# Dummy points
dummy_df <- as.data.frame(Q_lin$dummy)
dummy_df$type <- "dummy"
dummy_df$weight <- Q_lin$weights$dummy

# Combine
quad_df <- bind_rows(data_df, dummy_df)

# Convert to LPP
full_lpp <- lpp(quad_df[,c(1,2)], nganjuk_ln); head(full_lpp)

# ================================
# 6. EXTRACT COORDINATES + META
# ================================
network_coords <- as.data.frame(coords(full_lpp)); head(network_coords)

# network_coords$seg <- full_lpp$seg
# network_coords$tp  <- full_lpp$tp

network_coords$label <- spatstat.geom::is.data(Q_lin)
network_coords$vol   <- w.quad(Q_lin)

# ================================
# 7. BUILD FINAL DATASET
# ================================
final_data <- network_coords %>%
  mutate(
    # Distance covariate
    Jarak_Lampu = f(x, y, seg, tp),
    
    # Network covariates
    Jenis_Jalan = Jenis_Jalan(x, y, seg, tp),
    n_Lajur     = n_Lajur(x, y, seg, tp),
    n_Jalur     = n_Jalur(x, y, seg, tp),
    Kerb_R      = Kerb_R(x, y, seg, tp),
    Kerb_L      = Kerb_L(x, y, seg, tp)
  )

# Label: data = 1, dummy = -1
final_data$label <- ifelse(final_data$label, 1, -1)

# ================================
# 8. CHECK RESULT
# ================================
head(final_data)
nrow(final_data)
str(final_data)

# ================================
# 9. ANALYSIS
# ================================

train_lgbpp_fixed <- function(X, y, vol, loss = 'poisson', F_prime) {
  cat("--- Melatih LightGBM dengan parameter tetap...\n")
  
  base_params_lgb <- list(
    boosting_type = 'goss',
    top_rate = 0.2,
    other_rate = 0.2,
    is_enable_bundle = TRUE,
    max_conflict_rate = 0.1,
    max_bin = 255L,
    bin_construct_sample_cnt = 200000L,
    min_data_in_bin = 3L,
    deterministic = TRUE,
    verbose = -1L,
    num_threads = 1L
  )
  
  final_params <- c(base_params_lgb, list(
    learning_rate = 0.001,
    lambda_l1 = 0,
    lambda_l2 = 50,
    num_leaves = 63L
  ))
  
  # 🔥 KEY FIX
  train_set <- lgb$Dataset(
    data = as.matrix(X),
    label = pd$Series(y),
    feature_name = colnames(X)
  )
  
  model <- lgbpp_py(
    data = train_set,
    vol = pd$Series(vol),
    params = final_params,
    loss = loss,
    F_prime = F_prime,
    num_boost_round = 2000L,
    valid_sets = list(train_set),
    callbacks = list(lgb$log_evaluation(500L))
  )
  
  return(model)
}

train_xgbpp_fixed <- function(X, y, vol, loss = 'poisson', F_prime) {
  cat("--- Melatih XGBoost dengan parameter tetap...\n")
  
  base_params_xgb <- list(
    booster = 'gbtree',
    subsample = 1.0,
    colsample_bytree = 1/3,
    nthread = 1L,
    tree_method = 'exact',
    verbosity = 1L,
    min_child_weight = 1e-3,
    base_score = 0.001
  ) 
  
  final_params <- c(base_params_xgb, list(
    eta = 0.001,
    alpha = 0,
    lambda = 50,
    max_depth = 6L
  ))
  
  dtrain <- xgb$DMatrix(
    data = as.matrix(X),
    label = pd$Series(y)
  )
  
  # 🔥 KEY FIX
  dtrain$feature_names <- colnames(X)
  
  evals_list <- list(reticulate::tuple(dtrain, "train"))
  
  model <- xgbpp_py(
    dtrain = dtrain,
    vol = pd$Series(vol),
    params = final_params,
    loss = loss,
    F_prime = F_prime,
    num_boost_round = 2000L,
    evals = evals_list,
    verbose_eval = 500L
  )
  
  return(model)
}

run_models_LPP <- function(final_data, L, kcl_point, output_dir = "results_ml", type = 'squared', use_tuning = FALSE, n_trials = 20) {
  
  final_data_clean <- final_data %>%
    filter(
      is.finite(n_Lajur),
      is.finite(n_Jalur),
      is.finite(Jarak_Lampu),
      !is.na(n_Lajur),
      !is.na(n_Jalur),
      !is.na(Jarak_Lampu),
      !is.na(vol),
      vol > 0
    )
  
  cat("Before:", nrow(final_data), "\n")
  cat("After :", nrow(final_data_clean), "\n")
  
  final_data <- final_data_clean
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  model_list <- c("xgb", "lgb")
  
  png(file.path(output_dir, "original_points.png"),
      width = 800, height = 800, res = 120)
  
  plot(L, col = "grey70", lwd = 2,
       main = "Observed Accident Locations on Network")
  
  plot(kcl_point,
       add = TRUE,
       pch = 16,
       col = "red",
       cex = 0.6)
  
  dev.off()
  
  for (model_type in model_list) {
    
    cat("\n=============================\n")
    cat("Running model:", toupper(model_type), "\n")
    cat("=============================\n")
    
    start_time <- Sys.time()
    
    # ================================
    # PREPARE DATA
    # ================================
    sim_data <- final_data
    
    sf <- min(max(sim_data$y) - min(sim_data$y),
              max(sim_data$x) - min(sim_data$x))
    
    if(type == 'squared'){
      scale <- sf*sf
    }else{
      scale <- sf
    }
    
    sim_data$vol <- sim_data$vol / (scale)
    
    # ================================
    # FEATURE SELECTION (UPDATED)
    # ================================
    COVARIATES_TO_USE <- c("n_Lajur", "n_Jalur", "Jarak_Lampu")
    
    features_cols <- COVARIATES_TO_USE
    
    # Safety check
    missing_cols <- setdiff(features_cols, names(sim_data))
    if (length(missing_cols) > 0) {
      stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
    }
    
    X   <- sim_data[, features_cols]
    y   <- sim_data$label
    vol <- sim_data$vol
    
    # ================================
    # TRAIN MODEL
    # ================================
    F_prime <- 0
    
    cat("Feature columns used:\n")
    print(features_cols)
    
    cat("\nStructure of X:\n")
    str(X)
    
    if (model_type == "xgb") {
      
      if (use_tuning) {
        cat("--- Using Optuna tuning for XGBoost ---\n")
        
        tuning_result <- tune_xgbpp(
          X_df = pd$DataFrame(X),
          y_series = pd$Series(y),
          vol_series = pd$Series(vol),
          loss = "poisson",
          F_prime = F_prime,
          n_trials = as.integer(n_trials)
        )
        
        model <- tuning_result$final_model
        
        cat("Best params (XGB):\n")
        print(tuning_result$best_params)
      } else {
        cat("--- Using fixed XGBoost ---\n")
        model <- train_xgbpp_fixed(X, y, vol, "poisson", F_prime)
      }
      print(py_to_r(model$get_score()))
    } else {
      
      if (use_tuning) {
        cat("--- Using Optuna tuning for LightGBM ---\n")
        
        tuning_result <- tune_lgbpp(
          X_df = pd$DataFrame(X),
          y_series = pd$Series(y),
          vol_series = pd$Series(vol),
          loss = "poisson",
          F_prime = F_prime,
          n_trials = as.integer(n_trials)
        )
        
        model <- tuning_result$final_model
        
        cat("Best params (LGB):\n")
        print(tuning_result$best_params)
        
      }else {
        cat("--- Using fixed LightGBM ---\n")
        model <- train_lgbpp_fixed(X, y, vol, "poisson", F_prime)
      }
      print(py_to_r(model$feature_importance()))
    }
    
    
    # ================================
    # PREDICTION
    # ================================
    X_matrix <- as.matrix(X)
    
    dtest <- xgb$DMatrix(data = X_matrix)
    dtest$feature_names <- colnames(X)
    
    preds <- if (model_type == "xgb") {
      model$predict(dtest)
    } else {
      model$predict(X_matrix)
    }
    
    cat("======MODEL TRAINING SUCCESSFULLY EXECUTED======")
    
    # ================================
    # LOG-LIKELIHOOD
    # ================================
    left  <- sum(preds[y == 1] -  log(scale))
    right <- sum(exp(preds) * vol)
    
    loglik <- left - right
    
    comp_time <- Sys.time() - start_time
    
    # ================================
    # INTENSITY MAP
    # ================================
    sim_data$pred_intensity <- exp(preds)
    
    dummy <- sim_data[sim_data$label == -1, ]
    
    dummy_lpp <- lpp(dummy[, c("x","y")], L)
    marks(dummy_lpp) <- dummy$pred_intensity
    
    W <- as.owin(L)
    
    pp_quad <- ppp(
      x = dummy$x,
      y = dummy$y,
      window = W,
      marks = dummy$pred_intensity
    )
    
    # Smoothing (keep but slightly improve)
    im_pred <- Smooth(pp_quad, sigma = bw.diggle(pp_quad))
    
    intensity_map <- linim(L, im_pred)
    
    # ================================
    # SAVE INTENSITY PLOT
    # ================================
    
    png(file.path(output_dir, paste0("intensity_", model_type, ".png")),
        width = 800, height = 800, res = 120)
    
    plot(intensity_map,
         main = paste("Estimated Intensity -", toupper(model_type)),
         col = viridis::viridis(100),
         ribbon = TRUE)
    
    # Overlay network (thin)
    plot(L, add = TRUE, col = "grey40", lwd = 1)
    
    # Overlay observed points
    # plot(kcl_point,
    #      add = TRUE,
    #      pch = 16,
    #      col = "white",
    #      cex = 0.4)
    
    dev.off()
    
    # ================================
    # SAVE IMPORTANCE PLOT
    # ================================
    
    plot_feature_importance_xgb <- function(model, features_cols) {
      
      # SAFE extraction
      if (is.list(model) && "model" %in% names(model)) {
        model <- model$model
      }
      
      # ---- Get importance from model ----
      gain_dict  <- tryCatch(py_to_r(model$get_score(importance_type = "gain")), error = function(e) list())
      split_dict <- tryCatch(py_to_r(model$get_score(importance_type = "weight")), error = function(e) list())
      
      # ---- Initialize vectors ----
      gain_vec  <- setNames(rep(0, length(features_cols)), features_cols)
      split_vec <- setNames(rep(0, length(features_cols)), features_cols)
      
      if (length(gain_dict) > 0)
        gain_vec[names(gain_dict)] <- as.numeric(gain_dict)
      
      if (length(split_dict) > 0)
        split_vec[names(split_dict)] <- as.numeric(split_dict)
      
      # ---- DataFrames ----
      gain_df <- data.frame(feature = features_cols,
                            importance = gain_vec,
                            stringsAsFactors = FALSE)
      
      gain_df$importance_pct <- gain_df$importance / sum(gain_df$importance + 1e-12)
      gain_df <- gain_df[order(-gain_df$importance), ]
      
      split_df <- data.frame(feature = features_cols,
                             importance = split_vec,
                             stringsAsFactors = FALSE)
      
      split_df$importance_pct <- split_df$importance / sum(split_df$importance + 1e-12)
      split_df <- split_df[order(-split_df$importance), ]
      
      # ---- Plots ----
      p_gain <- ggplot(gain_df,
                       aes(x = reorder(feature, importance),
                           y = importance)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        coord_flip() +
        labs(title = "XGBoost Feature Importance (Gain)",
             x = "Feature", y = "Gain") +
        theme_minimal(base_size = 14)
      
      p_split <- ggplot(split_df,
                        aes(x = reorder(feature, importance),
                            y = importance)) +
        geom_bar(stat = "identity", fill = "darkorange") +
        coord_flip() +
        labs(title = "XGBoost Feature Importance (Split Count)",
             x = "Feature", y = "Split Count") +
        theme_minimal(base_size = 14)
      
      return(list(
        gain_plot  = p_gain,
        split_plot = p_split,
        gain_data  = gain_df,
        split_data = split_df
      ))
    }
    
    plot_feature_importance_lgb <- function(model, features_cols) {
      
      # SAFE extraction (for wrapper case)
      if (is.list(model) && "model" %in% names(model)) {
        model <- model$model
      }
      
      # ---- Get importance ----
      gain_vec  <- as.numeric(py_to_r(model$feature_importance(importance_type = "gain")))
      split_vec <- as.numeric(py_to_r(model$feature_importance(importance_type = "split")))
      
      # Ensure correct length
      if (length(gain_vec) != length(features_cols)) {
        stop("Mismatch between feature importance and feature columns")
      }
      
      # ---- DataFrames ----
      gain_df <- data.frame(
        feature = features_cols,
        importance = gain_vec
      )
      
      gain_df$importance_pct <- gain_df$importance / sum(gain_df$importance + 1e-12)
      gain_df <- gain_df[order(-gain_df$importance), ]
      
      split_df <- data.frame(
        feature = features_cols,
        importance = split_vec
      )
      
      split_df$importance_pct <- split_df$importance / sum(split_df$importance + 1e-12)
      split_df <- split_df[order(-split_df$importance), ]
      
      # ---- Plots ----
      p_gain <- ggplot(gain_df,
                       aes(x = reorder(feature, importance),
                           y = importance)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        coord_flip() +
        labs(title = "LightGBM Feature Importance (Gain)",
             x = "Feature", y = "Gain") +
        theme_minimal()
      
      p_split <- ggplot(split_df,
                        aes(x = reorder(feature, importance),
                            y = importance)) +
        geom_bar(stat = "identity", fill = "darkorange") +
        coord_flip() +
        labs(title = "LightGBM Feature Importance (Split Count)",
             x = "Feature", y = "Split Count") +
        theme_minimal()
      
      return(list(
        gain_plot  = p_gain,
        split_plot = p_split,
        gain_data  = gain_df,
        split_data = split_df
      ))
    }
    
    # ================================
    # FEATURE IMPORTANCE
    # ================================
    if (model_type == "xgb") {
      fi <- plot_feature_importance_xgb(model, features_cols)
      write.xlsx(sim_data, "results_dataframe_xgb.xlsx", rowNames = FALSE)
    } else {
      fi <- plot_feature_importance_lgb(model, features_cols)
      write.xlsx(sim_data, "results_dataframe_lgb.xlsx", rowNames = FALSE)
    }
    
    # Save plots
    ggsave(file.path(output_dir, paste0("fi_gain_", model_type, ".png")),
           fi$gain_plot, width = 8, height = 6)
    
    ggsave(file.path(output_dir, paste0("fi_split_", model_type, ".png")),
           fi$split_plot, width = 8, height = 6)
    
    # ================================
    # SAVE METRICS
    # ================================
    txt_path <- file.path(output_dir, paste0("summary_", model_type, ".txt"))
    
    cat(
      paste(
        "Model:", toupper(model_type),
        "\nLog-Likelihood:", round(loglik, 4),
        "\nPredicted Number of Events:", round(right, 4),
        "\nComputation Time (sec):", round(as.numeric(comp_time, units="secs"), 2),
        "\nScale Factor:",round(sf,4),
        "\nType:", type
      ),
      file = txt_path
    )
    
    cat("Finished:", model_type, "\n")
  }
}

run_models_LPP(final_data, L = nganjuk_ln, kcl_point = kcl_point, type = 'single', use_tuning = TRUE)