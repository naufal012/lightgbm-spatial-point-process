source('Code linquad_080524.R')

standardize_im_safe <- function(img) {
  # img: object of class 'im'
  if (!inherits(img, "im")) stop("standardize_im_safe: input must be an 'im' object")
  
  v <- img$v
  # replace non-finite values with NA
  v[!is.finite(v)] <- NA
  
  vals <- as.vector(v)
  mean_val <- mean(vals, na.rm = TRUE)
  sd_val <- sd(vals, na.rm = TRUE)
  
  # if sd is zero/undefined, force sd = 1 so we only center the image
  if (is.na(sd_val) || sd_val == 0) {
    warning("standardize_im_safe: zero or undefined sd -> centering only (sd forced = 1)")
    sd_val <- 1
  }
  
  standardized_values <- (vals - mean_val) / sd_val
  mat <- matrix(standardized_values, nrow = nrow(v), ncol = ncol(v))
  
  im(mat,
     xcol = img$xcol,
     yrow = img$yrow,
     xrange = img$xrange,
     yrange = img$yrange,
     unitname = img$unitname)
}

# Helper: safe exponentiation with clipping to avoid Inf
safe_exp_im <- function(logim, clip_lower = -50, clip_upper = 50) {
  if (!inherits(logim, "im")) stop("safe_exp_im: input must be 'im'")
  v <- logim$v
  v[!is.finite(v)] <- NA
  
  # clip using quantiles to avoid extreme values
  q_lo <- quantile(v, probs = 0.01, na.rm = TRUE)
  q_hi <- quantile(v, probs = 0.99, na.rm = TRUE)
  lower <- pmin(q_lo - 10, clip_lower)
  upper <- pmax(q_hi + 10, clip_upper)
  v_clipped <- v
  v_clipped[!is.na(v_clipped) & v_clipped < lower] <- lower
  v_clipped[!is.na(v_clipped) & v_clipped > upper] <- upper
  
  mat <- exp(v_clipped)
  mat[!is.finite(mat)] <- 0
  
  im(matrix(mat, nrow = nrow(logim$v), ncol = ncol(logim$v)),
     xcol = logim$xcol, yrow = logim$yrow,
     xrange = logim$xrange, yrange = logim$yrange,
     unitname = logim$unitname)
}

# safe as.im that forces window W and dimyx dims if needed
as.im_force <- function(x, W, dimyx = c(200,200)) {
  # x can be linfun or im; W is owin
  if (inherits(x, "im")) {
    # resample to window W if windows mismatch
    if (!identical(x$xrange, W$xrange) || !identical(x$yrange, W$yrange)) {
      return(as.im(x, W = W, dimyx = dimyx))
    } else return(x)
  }
  return(as.im(x, W = W, dimyx = dimyx))
}

simulate_LPP_process <- function(L,
                                      covariates_list,
                                      coefficients,
                                      n_points = 2000,
                                      n_pix = 200,
                                      intensity_form = 'linear',
                                      clip_quantiles = c(0.01, 0.99),
                                      smooth_sigma = NULL) {
  # L: linnet
  # covariates_list: named list of linfun objects (or functions with signature (x,y,seg,tp))
  # coefficients: numeric vector or named vector
  # n_pix: resolution for intermediate 'im'
  # intensity_form: 'linear' | 'complex' | 'complex_sparse'
  
  if (!inherits(L, "linnet")) stop("L must be a 'linnet' object")
  if (!is.list(covariates_list) || length(covariates_list) == 0) stop("covariates_list must be a non-empty list")
  
  # --- 1. Create 2D images for covariates covering the network window ---
  W <- as.owin(L)
  dims <- c(n_pix, n_pix)
  cov_images_2D <- lapply(covariates_list, function(f_loop) {
    # as.im will accept linfun or function; force window W and dims
    as.im_force(f_loop, W = W, dimyx = dims)
  })
  names(cov_images_2D) <- names(covariates_list)
  
  # --- 1b. Standardize images but keep NA where appropriate ---
  cov_images_2D <- lapply(cov_images_2D, standardize_im_safe)
  
  # --- 1c. Align coefficients ---
  cov_names <- names(cov_images_2D)
  if (is.null(names(coefficients))) {
    if (length(coefficients) != length(cov_images_2D)) stop(sprintf("coefficients must be length %d (or be a named vector).", length(cov_images_2D)))
    coefficients <- as.numeric(coefficients)
    names(coefficients) <- cov_names
  } else {
    if (!all(cov_names %in% names(coefficients))) stop("Named 'coefficients' must cover all covariate names: ", paste(cov_names, collapse = ", "))
    coefficients <- coefficients[cov_names]
  }
  
  # --- 2. Build log-intensity image robustly (do not propagate NA) ---
  template <- cov_images_2D[[1]]
  zero_mat <- matrix(0, nrow = nrow(template$v), ncol = ncol(template$v))
  logLambda_im <- im(zero_mat, xcol = template$xcol, yrow = template$yrow,
                     xrange = template$xrange, yrange = template$yrange,
                     unitname = template$unitname)
  
  safe_im_vals <- function(imobj) {
    v <- imobj$v
    v[!is.finite(v)] <- NA
    v[is.na(v)] <- 0   # treat missing covariate as 0-effect in linear combination
    im(matrix(v, nrow = nrow(imobj$v), ncol = ncol(imobj$v)),
       xcol = imobj$xcol, yrow = imobj$yrow,
       xrange = imobj$xrange, yrange = imobj$yrange,
       unitname = imobj$unitname)
  }
  
  if (intensity_form == 'linear') {
    for (i in seq_along(cov_images_2D)) {
      coef_i <- coefficients[i]
      im_i_clean <- safe_im_vals(cov_images_2D[[i]])
      logLambda_im <- logLambda_im + coef_i * im_i_clean
    }
  } else if (intensity_form == 'complex') {
    # implement cycle-based combination but use safe_im_vals for operands
    for (i in seq_along(cov_images_2D)) {
      coef_i <- coefficients[i]
      cycle_step <- (i - 1) %% 4
      im_i <- safe_im_vals(cov_images_2D[[i]])
      if (cycle_step == 0) {
        term_im <- coef_i * im_i
      } else if (cycle_step == 1) {
        im_prev <- safe_im_vals(cov_images_2D[[max(1, i-1)]])
        term_im <- coef_i * (im_prev * im_i)
      } else if (cycle_step == 2) {
        # apply exp on pixel values with clipping
        # compute exp on matrix level via safe_exp_im on im_i
        term_im <- coef_i * safe_exp_im(im_i)
      } else {
        # sin on pixel values (sin handles NA -> NA), replace NA with 0 prior
        mat <- im_i$v
        mat[is.na(mat)] <- 0
        term_im <- coef_i * im(matrix(sin(mat), nrow = nrow(mat), ncol = ncol(mat)),
                               xcol = im_i$xcol, yrow = im_i$yrow,
                               xrange = im_i$xrange, yrange = im_i$yrange,
                               unitname = im_i$unitname)
      }
      # replace NA by 0 in term_im before adding
      term_mat <- term_im$v
      term_mat[!is.finite(term_mat)] <- 0
      term_im$v <- term_mat
      logLambda_im <- logLambda_im + term_im
    }
  } else if (intensity_form == 'complex_sparse') {
    if (length(cov_images_2D) < 2 || length(coefficients) < 2) stop("For 'complex_sparse', need at least 2 covariates and 2 coefficients.")
    im1 <- safe_im_vals(cov_images_2D[[1]])
    im2 <- safe_im_vals(cov_images_2D[[2]])
    term1_mat <- (im1$v + im2$v) / 2
    term2_mat <- sqrt(pmax(0, exp(im1$v * im2$v)))
    term_mat <- coefficients[1] * term1_mat + coefficients[2] * term2_mat
    term_mat[!is.finite(term_mat)] <- 0
    logLambda_im <- im(matrix(term_mat, nrow = nrow(term_mat), ncol = ncol(term_mat)),
                       xcol = template$xcol, yrow = template$yrow,
                       xrange = template$xrange, yrange = template$yrange,
                       unitname = template$unitname)
  } else {
    stop("intensity_form must be 'linear', 'complex', or 'complex_sparse'.")
  }
  
  # optional smoothing / clipping before exponentiation
  if (!is.null(smooth_sigma)) logLambda_im <- Smooth(logLambda_im, sigma = smooth_sigma)
  
  # --- 3. Safe exponentiation to get Lambda_im_2D ---
  Lambda_im_2D <- safe_exp_im(logLambda_im)
  
  # ensure Lambda_im_2D covers window W
  if (!identical(Lambda_im_2D$xrange, W$xrange) || !identical(Lambda_im_2D$yrange, W$yrange)) {
    Lambda_im_2D <- as.im(Lambda_im_2D, W = W, dimyx = dims)
  }
  
  # fill NA with small positive value (median or eps)
  vals <- Lambda_im_2D$v
  if (all(is.na(vals))) {
    stop("Lambda_im_2D contains only NA values — check covariates and coefficients")
  }
  vals[is.na(vals)] <- median(vals, na.rm = TRUE)
  vals[vals < 0] <- 0
  Lambda_im_2D$v <- vals
  
  # --- 4. Convert to linim on network safely ---
  Lambda_linim_1D <- tryCatch(spatstat.linnet::linim(L, Lambda_im_2D), error = function(e) {
    warning("linim conversion failed: ", e$message, " -> attempting to build constant linim fallback")
    # fallback constant intensity: set to mean value on image
    const_val <- mean(Lambda_im_2D$v, na.rm = TRUE)
    im_const <- im(matrix(const_val, nrow = nrow(Lambda_im_2D$v), ncol = ncol(Lambda_im_2D$v)),
                   xcol = Lambda_im_2D$xcol, yrow = Lambda_im_2D$yrow,
                   xrange = Lambda_im_2D$xrange, yrange = Lambda_im_2D$yrange,
                   unitname = Lambda_im_2D$unitname)
    spatstat.linnet::linim(L, im_const)
  })
  
  # --- 5. Rescale to get approximately n_points expected ---
  total_integral <- tryCatch(integral(Lambda_linim_1D), error = function(e) NA_real_)
  beta0 <- NA_real_
  Lambda_final_linim <- NULL
  
  if (!is.finite(total_integral) || total_integral <= 0) {
    warning("Calculated intensity integral is zero or non-finite. Falling back to uniform intensity on network.")
    total_length <- sum(lengths(L))
    const_val <- n_points / total_length
    im_const <- im(matrix(const_val, nrow = nrow(Lambda_im_2D$v), ncol = ncol(Lambda_im_2D$v)),
                   xcol = Lambda_im_2D$xcol, yrow = Lambda_im_2D$yrow,
                   xrange = Lambda_im_2D$xrange, yrange = Lambda_im_2D$yrange,
                   unitname = Lambda_im_2D$unitname)
    Lambda_final_linim <- spatstat.linnet::linim(L, im_const)
    beta0 <- NA_real_
  } else {
    beta0 <- log(n_points / total_integral)
    # multiply linim by scalar exp(beta0)
    Lambda_final_linim <- Lambda_linim_1D * exp(beta0)
  }
  
  # --- 6. Simulate Poisson LPP from the final intensity and build quadrature ---
  sim_points_lpp <- rpoislpp(Lambda_final_linim)
  qd_scheme <- spatstat.linnet::linequad(sim_points_lpp)
  
  data_df <- as.data.frame(qd_scheme$data)
  data_df$weight <- qd_scheme$weights$data
  dummy_df <- as.data.frame(qd_scheme$dummy)
  dummy_df$weight <- qd_scheme$weights$dummy
  quad_df <- bind_rows(data_df, dummy_df)
  
  full_lpp <- lpp(quad_df[, c(1,2)], L)
  network_coords <- as.data.frame(coords(full_lpp))
  network_coords$label <- spatstat.geom::is.data(qd_scheme)
  network_coords$vol <- w.quad(qd_scheme)
  
  qd_scheme_logi <- quadscheme.logi.linnet(sim_points_lpp, dummytype = 'binomial')
  
  data_df_logi <- as.data.frame(qd_scheme_logi$data)
  dummy_df_logi <- as.data.frame(qd_scheme_logi$dummy)
  quad_df_logi <- bind_rows(data_df_logi, dummy_df_logi)
  
  full_lpp_logi <- lpp(quad_df_logi[, c(1,2)], L)
  network_coords_logi <- as.data.frame(coords(full_lpp_logi))
  network_coords_logi$weight <- qd_scheme_logi$w
  network_coords_logi$label <- spatstat.geom::is.data(qd_scheme_logi)
  network_coords_logi$vol <- w.quad(qd_scheme_logi)
  
  
  # --- 7. Evaluate covariates at quadrature locations and assemble dataframe ---
  cov_values_all <- lapply(covariates_list, function(f_loop) {
    vals <- f_loop(x = network_coords$x, y = network_coords$y, seg = network_coords$seg, tp = network_coords$tp)
    vals[!is.finite(vals)] <- NA
    vals
  })
  
  df_full <- as.data.frame(c(list(x = network_coords$x, y = network_coords$y), cov_values_all,
                             list(label = network_coords$label, vol = network_coords$vol)))
  
  cov_values_all_logi <- lapply(covariates_list, function(f_loop) {
    vals <- f_loop(x = network_coords_logi$x, y = network_coords_logi$y, seg = network_coords_logi$seg, tp = network_coords_logi$tp)
    vals[!is.finite(vals)] <- NA
    vals
  })
  
  df_full_logi <- as.data.frame(c(list(x = network_coords_logi$x, y = network_coords_logi$y), cov_values_all_logi,
                                  list(label = network_coords_logi$label, vol = network_coords_logi$vol)))
  
  # --- 8. Prepare outputs ---
  return(list(
    sim_points_lpp     = sim_points_lpp,
    intensity_im_2D    = Lambda_im_2D,
    intensity_linim_1D = Lambda_final_linim,
    sim_data_full      = df_full,
    sim_data_full_logi = df_full_logi,
    quad_scheme        = qd_scheme,
    quad_scheme_logi   = qd_scheme_logi,
    covs               = cov_images_2D,
    alpha              = beta0
  ))
}