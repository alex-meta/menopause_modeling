load("merged_with_fmp.rds")

## ===== helpers =====
.choose_first_present <- function(df, candidates) {
  for (nm in candidates) if (nm %in% names(df)) return(nm)
  stop("None of the columns were found: ", paste(candidates, collapse=", "))
}

## Euclidean distance on common (non-missing) features, adjusted for different vector length
.euclid_partial <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  p  <- sum(ok)
  if (p == 0) return(Inf)
  d  <- sqrt(sum((a[ok] - b[ok])^2))
  # normalization so distances are comparable with different numbers of shared features
  d * sqrt(length(a) / p)
}

## window-based standardization (z-score)
.zscore <- function(M) {
  mu <- colMeans(M, na.rm = TRUE)
  sd <- apply(M, 2, sd, na.rm = TRUE)
  sd[sd == 0 | !is.finite(sd)] <- 1
  Z  <- sweep(M, 2, mu, "-")
  Z  <- sweep(Z, 2, sd, "/")
  list(Z = Z, mu = mu, sd = sd)
}

## ===== main KNN predictor =====
# ref_df: long SWAN-by-visit table with known FMP age at participant level
# new_obs: data.frame with one row (age + same features) or a named vector
# feature_cols: features for proximity metric
# fmp_col: column in ref_df with estimated FMP age (e.g., est_age_fmp_amh_thr)
predict_fmp_knn <- function(
    ref_df,
    new_obs,
    feature_cols = c("FSH","AMH_pred_cal","GLUCRES","INSURES","TRIGRES","LDLRESU","HDLRESU","BMI","SHBG"),
    id_col   = "SWANID",
    age_col  = "AGE",
    visit_col= "VISIT",
    fmp_col  = "est_age_fmp_amh_thr",
    k        = 20,
    age_window = 1,
    floor_window_if_young = c(42, 44),   # if age < 43 → window [42,44]
    eps = 1e-6
) {
  stopifnot(id_col %in% names(ref_df), age_col %in% names(ref_df),
            visit_col %in% names(ref_df), fmp_col %in% names(ref_df))

  # keep only rows with known FMP in the reference
  ref_df <- ref_df[is.finite(ref_df[[fmp_col]]), , drop = FALSE]
  if (nrow(ref_df) == 0) stop("No rows in ref_df have a finite ", fmp_col)

  # coerce new_obs to a one-row data.frame
  if (is.null(dim(new_obs))) {
    new_obs <- as.data.frame(as.list(new_obs), stringsAsFactors = FALSE)
  }
  if (nrow(new_obs) != 1) stop("new_obs must contain exactly 1 row")

  # check age
  if (!age_col %in% colnames(new_obs)) stop("Missing age column in new_obs: ", age_col)
  new_age <- as.numeric(new_obs[[age_col]][1])
  if (!is.finite(new_age)) stop("Age in new_obs is non-numeric/NA")

  # age window
  if (new_age < 43) {
    age_min <- floor_window_if_young[1]
    age_max <- floor_window_if_young[2]
  } else {
    age_min <- new_age - age_window
    age_max <- new_age + age_window
  }

  # age-based candidates
  win_idx <- which(ref_df[[age_col]] >= age_min & ref_df[[age_col]] <= age_max)
  win <- ref_df[win_idx, , drop = FALSE]
  if (nrow(win) == 0) stop("No observations in age window [", age_min, ", ", age_max, "]")

  # align features list: keep only those present in both
  feats <- feature_cols[feature_cols %in% colnames(win) & feature_cols %in% colnames(new_obs)]
  if (length(feats) == 0) stop("No overlapping features between ref and new_obs")

  # drop observations with all-missing features
  keep <- rowSums(is.finite(as.matrix(win[, feats, drop=FALSE]))) > 0
  win  <- win[keep, , drop = FALSE]
  if (nrow(win) == 0) stop("No rows in the window have at least one non-missing predictor")

  # window standardization (robust to scale)
  zs <- .zscore(as.matrix(win[, feats, drop=FALSE]))
  Z_win <- zs$Z

  # prepare new vector (same columns), z-score using window stats
  x_new <- as.numeric(new_obs[1, feats, drop=TRUE])
  x_new <- (x_new - zs$mu) / zs$sd

  # compute distances for all rows in the window
  d <- vapply(seq_len(nrow(Z_win)), function(i) .euclid_partial(Z_win[i,], x_new), numeric(1))

  # de-dup by participant: keep minimum distance per SWANID
  tmp <- data.frame(
    SWANID = win[[id_col]],
    VISIT  = win[[visit_col]],
    AGE    = win[[age_col]],
    FMP    = win[[fmp_col]],
    dist   = d,
    stringsAsFactors = FALSE
  )

  # minimal distance per participant
  # (if you want to allow multiple visits per person — disable this block)
  ord <- order(tmp$SWANID, tmp$dist)
  tmp <- tmp[ord, ]
  min_by_id <- !duplicated(tmp$SWANID)
  pool <- tmp[min_by_id, , drop = FALSE]

  # pick k nearest participants
  pool <- pool[order(pool$dist), , drop = FALSE]
  k_eff <- min(k, nrow(pool))
  nn <- pool[seq_len(k_eff), , drop = FALSE]

  # weights = 1/(dist + eps), normalized
  w <- 1 / (nn$dist + eps)
  w <- w / sum(w)

  # weighted prediction of FMP age
  pred_fmp <- sum(w * nn$FMP)

  # return with diagnostics
  list(
    pred_age_fmp = pred_fmp,
    window = c(age_min = age_min, age_max = age_max),
    neighbors = transform(nn, weight = w),
    used_features = feats,
    new_age = new_age
  )
}



load("m_long_for_hackaton.bin")  # should provide m_long (long-format table)
## load("res2.bin")                 # optional: for variable descriptions
#
# ## Set of waves and variables to take from long:
# use_ds <- c("NHANES_0102","NHANES_0304","NHANES_0506","NHANES_0708","NHANES_0910",
#             "NHANES_1112","NHANES_1314","NHANES_1516","NHANES_1718","NHANES_1720",
#             "NHANES_2123","NHANES_9900","NHANES_III")
#
# use_variables <- c(
#   "RIAGENDR","RIDAGEYR","LBXAMH","LBXFSH","LBXEST","LBXLUH","LBXSHBG","LBXTST",
#   "LBDDHESI","LBXHSCRP","LBDLDLSI","LBDHDDSI","LBDTRSI","LBXTSH","BMXBMI","BMXHT",
#   "BPXPLS","RHQ060","RHQ010","LBDESOSI","LBDES1SI","LBDPG4SI","LBD17HSI","LBDANDSI",
#   "RHD143","RHD043","RHQ031","LBDINSI","LBDGLUSI","SSSHBG","BPXDI2"
# )

## Filter input long (if already filtered — this line is safe):
m_long <- subset(m_long, dataset %in% use_ds & variable %in% use_variables)

head(m_long)


## ---------- 1) Config: which waves/IDs to use ----------
use_ds <- c("NHANES_0102","NHANES_0304","NHANES_0506","NHANES_0708","NHANES_0910",
            "NHANES_1112","NHANES_1314","NHANES_1516","NHANES_1718","NHANES_1720",
            "NHANES_2123","NHANES_9900","NHANES_III")

## Map: canonical name -> NHANES ID
id_map <- c(
  AMH     = "LBXAMH",
  age     = "RIDAGEYR",
  fsh     = "LBXFSH",
  glucose = "LBDGLUSI",   # glucose (SI)
  insulin = "LBDINSI",    # insulin (SI)
  tg      = "LBDTRSI",    # triglycerides (SI)
  ldl     = "LBDLDLSI",   # LDL (SI)
  hdl     = "LBDHDDSI",   # HDL (SI)
  bmi     = "BMXBMI",
  shbg    = "LBXSHBG",
  sex         = "RIAGENDR",   # 1=male, 2=female
  preg_now    = "RHD143"      # 1=pregnant now
)

## Fallbacks (if primary ID is missing): feature name -> alternative IDs
fallbacks <- list(
  shbg = c("SSSHBG")  # if LBXSHBG is missing in a particular wave
)

## ---------- 2) Helpers ----------
winsor <- function(x, p=0.005){
  if (all(is.na(x))) return(x)
  ql <- suppressWarnings(quantile(x, p, na.rm=TRUE))
  qh <- suppressWarnings(quantile(x, 1-p, na.rm=TRUE))
  x[x < ql] <- ql; x[x > qh] <- qh; x
}
getcol <- function(df_wide, var_id){
  nm <- paste0("value.", var_id)
  if (nm %in% names(df_wide)) as.numeric(df_wide[[nm]]) else rep(NA_real_, nrow(df_wide))
}

## ---------- 3) Make wide from m_long and filter ----------
build_nhanes_wide <- function(m_long, id_map, fallbacks=list(), age_range=c(35,45), use_ds=NULL){
  # select required waves/variables
  if (!is.null(use_ds)) m_long <- subset(m_long, dataset %in% use_ds)
  keep_ids <- unique(c(unname(id_map), unlist(fallbacks, use.names = FALSE)))
  m_sub <- subset(m_long, variable %in% keep_ids)

  # females (RIAGENDR: 2)
  female_ids <- unique(subset(m_sub, variable==id_map["sex"] & value==2)$SEQN)
  m_sub <- m_sub[m_sub$SEQN %in% female_ids, ]

  # exclude currently pregnant (RHD143!=1)
  preg_ids <- unique(subset(m_sub, variable==id_map["preg_now"] & value==1)$SEQN)
  m_sub <- m_sub[!(m_sub$SEQN %in% preg_ids), ]

  # age filter
  age_rows <- subset(m_sub, variable==id_map["age"], select=c("SEQN","value"))
  ok_ids <- unique(age_rows$SEQN[age_rows$value >= age_range[1] & age_rows$value <= age_range[2]])
  m_sub <- m_sub[m_sub$SEQN %in% ok_ids, ]

  # aggregate duplicates (take first) and pivot to wide
  m_agg <- aggregate(value ~ SEQN + dataset + variable, data=m_sub, FUN=function(x) x[1])
  w <- reshape(m_agg, idvar=c("SEQN","dataset"), timevar="variable", direction="wide")
  w
}

## ---------- 4) Assemble canonical features from wide ----------
assemble_nhanes_features <- function(wide, id_map, fallbacks=list()){
  n <- nrow(wide)
  # take base IDs
  vals <- lapply(names(id_map), function(k) getcol(wide, id_map[[k]]))
  names(vals) <- names(id_map)
  # apply fallbacks if entire column is NA
  if (length(fallbacks)){
    for (nm in names(fallbacks)){
      if (all(is.na(vals[[nm]]))){
        for (fid in fallbacks[[nm]]){
          v_fb <- getcol(wide, fid)
          if (!all(is.na(v_fb))) { vals[[nm]] <- v_fb; break }
        }
      }
    }
  }
  df <- data.frame(
    SEQN    = wide$SEQN,
    dataset = wide$dataset,
    age     = vals$age,
    fsh     = vals$fsh,
    glucose = vals$glucose,
    insulin = vals$insulin,
    tg      = vals$tg,
    ldl     = vals$ldl,
    hdl     = vals$hdl,
    bmi     = vals$bmi,
    shbg    = vals$shbg,
    AMH     = vals$AMH,
    stringsAsFactors = FALSE
  )
  # light cleaning/winsorizing on non-missing
  for (v in c("age","fsh","glucose","insulin","tg","ldl","hdl","bmi","shbg")){
    if (v %in% names(df)) df[[v]] <- winsor(df[[v]])
  }
  df
}

## ---------- 5) Main workflow: m_long -> nhanes_df ----------
nhanes_from_mlong <- function(m_long, id_map, fallbacks=list(), use_ds=NULL,
                              age_range=c(35,45), require_complete=TRUE){
  w   <- build_nhanes_wide(m_long, id_map, fallbacks, age_range, use_ds)
  df  <- assemble_nhanes_features(w, id_map, fallbacks)

  # keep only rows with a complete set of keys
  req <- c("age","fsh","glucose","insulin","tg","ldl","hdl","bmi","shbg")
  if (require_complete){
    keep <- complete.cases(df[, req, drop=FALSE])
    df <- df[keep, , drop=FALSE]
  }
  rownames(df) <- NULL
  df
}

## ---------- 6) Call ----------
## Assumes m_long is already in memory (columns: SEQN, dataset, variable, value)
nhanes_df <- nhanes_from_mlong(
  m_long,
  id_map      = id_map,
  fallbacks   = fallbacks,
  use_ds      = use_ds,         # from section 1
  age_range   = c(35,45),
  require_complete = TRUE
)

## If AMH_pred_cal column is needed later:
## - if you have measured AMH -> simply copy it:
nhanes_df$AMH_pred_cal <- ifelse(!is.na(nhanes_df$AMH), nhanes_df$AMH, NA_real_)

## - OR (optionally) compute AMH from your two-stage GAM Gamma bundle:
## nhanes_df <- predict_amh_if_needed(nhanes_df, bundle)  # see function described earlier









## ---- Config ----
FEATURES <- c("AMH_pred_cal","FSH","GLUCRES","INSURES","TRIGRES",
              "LDLRESU","HDLRESU","BMI","SHBG")

K_NEIGHBORS <- 20
AGE_WIN <- 1   # window ±1 year by SWAN visit age
CLAMP_LOW <- 42  # if NHANES age < 43, use 42–44 (see below)

## ---- Helpers ----
zscale <- function(x, m, s) (x - m) / ifelse(s > 0, s, 1)

## Build SWAN reference: visits (visit features) + target fmp_age from personal field
## swan_all: one row per visit, must contain:
##   SWANID, VISIT, AGE, est_age_fmp_amh_thr (or your target), as well as FEATURES
build_swan_reference <- function(swan_all,
                                 features = FEATURES,
                                 fmp_col = "est_age_fmp_amh_thr") {
  stopifnot(all(c("SWANID","VISIT","AGE", fmp_col) %in% names(swan_all)))
  keep <- c("SWANID","VISIT","AGE", fmp_col, intersect(features, names(swan_all)))
  ref <- swan_all[ , keep]
  names(ref)[names(ref) == fmp_col] <- "fmp_age"
  # remove visits without target value
  ref <- ref[!is.na(ref$fmp_age), ]
  # light filtering of obviously absurd FMP (optional)
  ref <- ref[ref$fmp_age >= 40 & ref$fmp_age <= 60, ]
  ref
}

## If AMH is absent in NHANES — compute it via our two-stage GAM Gamma model
## bundle: list with model$early, model$late (+ iso, if present)
predict_amh_if_needed <- function(nh, bundle) {
  if (!"AMH_pred_cal" %in% names(nh)) {
    # check direct AMH
    if ("AMH" %in% names(nh)) {
      nh$AMH_pred_cal <- as.numeric(nh$AMH)
    } else {
      # need to derive from our features: age, fsh, glucose, insulin, tg, ldl, hdl, bmi, shbg
      req <- c("age","fsh","glucose","insulin","tg","ldl","hdl","bmi","shbg")
      miss <- setdiff(req, names(nh))
      if (length(miss)) stop("Missing features to compute AMH: ", paste(miss, collapse=", "))
      pred <- predict_two_stage_gam_gamma(bundle$model$early,
                                          bundle$model$late,
                                          nh, age_threshold = 37)
      if (!is.null(bundle$iso)) {
        pred <- iso_apply(bundle$iso, pred)  # calibration, if present
      }
      nh$AMH_pred_cal <- pmax(pred, 0)  # just in case
    }
  }
  nh
}

## kNN FMP forecast for a single NHANES woman
## nh_row: one row with age and FEATURES (in NHANES naming, see mapping below)
## ref   : swan_reference (from build_swan_reference)
## zstat : list with means/sds over SWAN for standardization
knn_fmp_one <- function(nh_row, ref, zstat, features = FEATURES,
                        k = K_NEIGHBORS, age_win = AGE_WIN, clamp_low = CLAMP_LOW) {
  a <- as.numeric(nh_row[["age"]])
  if (is.na(a)) return(NA_real_)

  # SWAN visit age window
  if (a < 43) {
    # rule from spec: if NHANES age < 43, use window 42–44
    ref_win <- ref[ref$AGE >= 42 & ref$AGE <= 44, , drop = FALSE]
  } else {
    ref_win <- ref[abs(ref$AGE - a) <= age_win, , drop = FALSE]
  }
  if (!nrow(ref_win)) return(NA_real_)

  # collect NHANES feature vector and reference features, z-normalize by SWAN
  fx <- features[features %in% names(ref_win) & features %in% names(nh_row)]
  if (!length(fx)) return(NA_real_)

  nh_vec <- as.numeric(nh_row[fx])
  ref_mat <- as.matrix(ref_win[ , fx, drop = FALSE])

  # if there are NAs — drop features with NA in the current row
  good <- !is.na(nh_vec)
  fx <- fx[good]
  nh_vec <- nh_vec[good]
  ref_mat <- ref_mat[ , good, drop = FALSE]
  if (!length(fx)) return(NA_real_)

  # z-normalization by SWAN stats
  m <- zstat$mean[fx]; s <- zstat$sd[fx]
  nh_z   <- zscale(nh_vec, m, s)
  ref_z  <- sweep(ref_mat, 2, m, FUN = "-")
  ref_z  <- sweep(ref_z,  2, s, FUN = "/")

  # distances
  d <- sqrt(rowSums((ref_z - matrix(nh_z, nrow(ref_z), ncol(ref_z), byrow=TRUE))^2))
  idx <- order(d)[seq_len(min(k, length(d)))]
  d_k <- d[idx]; w <- 1 / (d_k + 1e-6)

  # target: donor’s FMP age (personal!), NOT visit age
  fmp <- ref_win$fmp_age[idx]
  # weighted estimate
  sum(w * fmp) / sum(w)
}

## Wrapper: predict FMP age for NHANES table
## nhanes_df: NHANES table of women 35–45 with required features (see mapping below)
## swan_all : SWAN visits table containing fmp_age and visit features
predict_fmp_nhanes_from_swan <- function(nhanes_df, swan_all,
                                         features = FEATURES,
                                         fmp_col = "est_age_fmp_amh_thr",
                                         k = K_NEIGHBORS, age_win = AGE_WIN) {

  # 1) prepare SWAN reference
  ref <- build_swan_reference(swan_all, features, fmp_col)

  # 2) z-normalization stats over SWAN (all visits)
  zmean <- sapply(features, function(v) mean(ref[[v]], na.rm=TRUE))
  zsd   <- sapply(features, function(v)   sd(ref[[v]], na.rm=TRUE))
  zstat <- list(mean = zmean, sd = zsd)

  # 3) predict AMH in NHANES if missing
  # nhanes_df <- predict_amh_if_needed(nhanes_df, bundle)  # if AMH must be computed by model

  # 4) iterate over rows
  res <- apply(nhanes_df, 1, function(r) {
    knn_fmp_one(as.list(r), ref, zstat, features = features, k = k, age_win = age_win)
  })
  as.numeric(res)
}


predict_fmp_nhanes<-predict_fmp_nhanes_from_swan(nhanes_df = nhanes_df,swan_all = swan_all,k = 10,age_win = 2)

plot(density(predict_fmp_nhanes,na.rm = T))
abline(v=52,lty=2)

plot(density(swan_all$est_age_fmp_amh_thr,na.rm = T))
abline(v=52,lty=2)


# problem: we predict ~48 years on average, but it should be ~52
summary(nh_pred$pred_age_fmp)










# check distribution of values between SWAN and NHANES

## safe numeric
safe_num <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    x <- trimws(x)
    x <- gsub(",", ".", x, fixed = TRUE)
    x <- gsub("[^0-9eE+\\-\\.]", "", x)
  }
  suppressWarnings(as.numeric(x))
}

## SWAN → SI, canonical names as in nhanes_df
swan_to_canon_si <- function(swan_df) {
  stopifnot(is.data.frame(swan_df))
  out <- data.frame(
    SWANID = swan_df$SWANID,
    VISIT  = swan_df$VISIT,
    age    = safe_num(swan_df$AGE),                       # years (matches)
    fsh    = safe_num(swan_df$FSH),                       # mIU/mL — matches NHANES
    ## SWAN units → SI:
    ## glucose mg/dL → mmol/L: *0.0555
    glucose = safe_num(swan_df$GLUCRES) * 0.0555,
    ## insulin uIU/mL → pmol/L: *6.945   (NHANES LBDINSI — pmol/L)
    insulin = safe_num(swan_df$INSURES) * 6.945,
    ## TG mg/dL → mmol/L: *0.01129
    tg  = safe_num(swan_df$TRIGRES) * 0.01129,
    ## LDL/HDL mg/dL → mmol/L: *0.02586
    ldl = safe_num(swan_df$LDLRESU) * 0.02586,
    hdl = safe_num(swan_df$HDLRESU) * 0.02586,
    ## BMI (as is), SHBG — nM == nmol/L (matches)
    bmi = safe_num(swan_df$BMI),
    shbg = safe_num(swan_df$SHBG)
  )
  out
}

safe_num <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    x <- trimws(x)
    x <- gsub(",", ".", x, fixed = TRUE)
    x <- gsub("[^0-9eE+\\-\\.]", "", x)
  }
  suppressWarnings(as.numeric(x))
}

qsum <- function(x) {
  x <- safe_num(x)
  nn  <- sum(!is.na(x))
  med <- if (nn > 0) median(x, na.rm = TRUE) else NA_real_
  p10 <- if (nn > 0) quantile(x, 0.10, na.rm = TRUE) else NA_real_
  p90 <- if (nn > 0) quantile(x, 0.90, na.rm = TRUE) else NA_real_
  c(n = nn, med = med, p10 = p10, p90 = p90)
}

print_block <- function(df, cols, title = "Summary") {
  cat("\n---", title, "---\n")
  present <- intersect(cols, names(df))
  if (!length(present)) {
    cat("None of these columns are present: ", paste(cols, collapse = ", "), "\n")
    return(invisible(NULL))
  }
  stats_list <- lapply(present, function(v) qsum(df[[v]]))
  mat <- do.call(rbind, stats_list)
  mat <- as.data.frame(mat, stringsAsFactors = FALSE)
  # ensure column names
  colnames(mat) <- c("n","med","p10","p90")
  mat$col <- present
  mat <- mat[, c("col","n","med","p10","p90")]
  mat$med <- round(as.numeric(mat$med), 3)
  mat$p10 <- round(as.numeric(mat$p10), 3)
  mat$p90 <- round(as.numeric(mat$p90), 3)
  print(mat, row.names = FALSE)
  invisible(mat)
}

common_cols <- c("age","fsh","glucose","insulin","tg","ldl","hdl","bmi","shbg")
print_block(swan_canon, common_cols, "SWAN (converted to SI, canonical)")
print_block(nhanes_df, common_cols, "NHANES (SI, canonical)")




## =========================
## SWAN → canonical SI (without overwriting originals)
## expected input columns in merged_by_fpm:
## SWANID, VISIT, AGE, FSH (mIU/mL), GLUCRES (mg/dL), INSURES (uIU/mL),
## TRIGRES (mg/dL), LDLRESU (mg/dL), HDLRESU (mg/dL), BMI, SHBG (nmol/L),
## (opt.) AMH_pred (ng/mL) or AMH_pred_cal (ng/mL)
## =========================

make_swan_canon <- function(merged_by_fpm,
                            keep_original = TRUE,
                            add_amh_si = TRUE,
                            amh_ngml_to_pmol = 7.14) {
  stopifnot(is.data.frame(merged_by_fpm))
  df <- merged_by_fpm

  # convenient getter (if column is missing — returns rep(NA, n))
  take <- function(nm) if (nm %in% names(df)) df[[nm]] else rep(NA_real_, nrow(df))

  # conversions to SI
  # glucose: mg/dL → mmol/L
  mgdl_to_mmol_glu <- function(x) x * 0.05551
  # triglycerides: mg/dL → mmol/L
  mgdl_to_mmol_tg  <- function(x) x * 0.01129
  # cholesterol: mg/dL → mmol/L
  mgdl_to_mmol_ch  <- function(x) x * 0.02586
  # insulin: uIU/mL → pmol/L
  uiul_to_pmol_ins <- function(x) x * 6.945

  # assemble canonical (SI) with unified names
  swan_canon <- data.frame(
    SWANID  = take("SWANID"),
    VISIT   = take("VISIT"),

    # age — as is
    age     = suppressWarnings(as.numeric(take("AGE"))),

    # FSH in SWAN — mIU/mL; 1 mIU/mL = 1 IU/L → numerically identical to SI.
    fsh     = suppressWarnings(as.numeric(take("FSH"))),

    # metabolic markers in SI
    glucose = mgdl_to_mmol_glu(suppressWarnings(as.numeric(take("GLUCRES")))),
    insulin = uiul_to_pmol_ins(suppressWarnings(as.numeric(take("INSURES")))),
    tg      = mgdl_to_mmol_tg (suppressWarnings(as.numeric(take("TRIGRES")))),
    ldl     = mgdl_to_mmol_ch (suppressWarnings(as.numeric(take("LDLRESU")))),
    hdl     = mgdl_to_mmol_ch (suppressWarnings(as.numeric(take("HDLRESU")))),

    # anthropometrics / SHBG (in SWAN it’s already nmol/L)
    bmi     = suppressWarnings(as.numeric(take("BMI"))),
    shbg    = suppressWarnings(as.numeric(take("SHBG")))
  )

  # AMH (if model predictions in ng/mL are present)
  amh_ng <- NULL
  if ("AMH_pred_cal" %in% names(df)) amh_ng <- suppressWarnings(as.numeric(df$AMH_pred_cal))
  if (is.null(amh_ng) && "AMH_pred" %in% names(df)) amh_ng <- suppressWarnings(as.numeric(df$AMH_pred))

  if (!is.null(amh_ng)) {
    swan_canon$amh_pred_ng_ml <- amh_ng
    if (add_amh_si) {
      swan_canon$amh_pred_pmol_l <- amh_ng * amh_ngml_to_pmol
    }
  }

  # optionally — append back all original columns (for debugging/comparison convenience)
  if (keep_original) {
    # do not duplicate already created names
    extra <- df[, setdiff(names(df), names(swan_canon)), drop = FALSE]
    swan_canon <- cbind(swan_canon, extra)
  }

  # light sanity-check of ranges (can be commented out)
  rng <- function(x) if (all(is.na(x))) c(NA, NA) else range(x, na.rm = TRUE)
  msg <- sprintf(
    paste(
      "SWAN canon — ranges:",
      "age:   [%.1f, %.1f]",
      "fsh:   [%.1f, %.1f] (IU/L)",
      "glu:   [%.2f, %.2f] (mmol/L)",
      "ins:   [%.1f, %.1f] (pmol/L)",
      "tg:    [%.2f, %.2f] (mmol/L)",
      "ldl:   [%.2f, %.2f] (mmol/L)",
      "hdl:   [%.2f, %.2f] (mmol/L)",
      "bmi:   [%.1f, %.1f]",
      "shbg:  [%.1f, %.1f] (nmol/L)",
      sep = "\n  "
    ),
    rng(swan_canon$age)[1],     rng(swan_canon$age)[2],
    rng(swan_canon$fsh)[1],     rng(swan_canon$fsh)[2],
    rng(swan_canon$glucose)[1], rng(swan_canon$glucose)[2],
    rng(swan_canon$insulin)[1], rng(swan_canon$insulin)[2],
    rng(swan_canon$tg)[1],      rng(swan_canon$tg)[2],
    rng(swan_canon$ldl)[1],     rng(swan_canon$ldl)[2],
    rng(swan_canon$hdl)[1],     rng(swan_canon$hdl)[2],
    rng(swan_canon$bmi)[1],     rng(swan_canon$bmi)[2],
    rng(swan_canon$shbg)[1],    rng(swan_canon$shbg)[2]
  )
  cat(msg, "\n")

  return(swan_canon)
}

## usage example:
swan_canon <- make_swan_canon(merged_with_fmp, keep_original = TRUE)
head(swan_canon)
















## =========================
## 1) Matching settings
## =========================

# Which features to use for similarity (without FSH/TSH):
FEATURES_MATCH <- c("age","AMH","glucose","insulin","tg","ldl","hdl","bmi","shbg")

# Feature weights in the metric (emphasize age and AMH)
FEATURE_WEIGHTS <- c(
  age    = 2.0,
  AMH    = 2.0,
  glucose= 1.0,
  insulin= 1.0,
  tg     = 1.0,
  ldl    = 1.0,
  hdl    = 1.0,
  bmi    = 1.0,
  shbg   = 1.0
)

AGE_WINDOW <- 1      # ±1 year around the request age
K_NEIGHBORS <- 20    # number of nearest trajectories

## =========================
## 2) SWAN reference preparation
## =========================

# Expected swan_canon: SWANID, VISIT, age, fsh, glucose, insulin, tg, ldl, hdl, bmi, shbg, AMH(if available) and fmp_age (menopause age)
# If measured AMH is absent, allow NA — it will just drop from used features; if AMH_pred_cal exists — rename to AMH
prepare_swan_reference <- function(swan_df, amh_col = c("AMH","AMH_pred_cal")) {
  df <- swan_df
  # choose AMH column: measured preferred, else model
  amh_name <- amh_col[amh_col %in% names(df)][1]
  if (!is.na(amh_name)) {
    if (!"AMH" %in% names(df)) df$AMH <- df[[amh_name]]
  }

  df$fmp_age <- df$est_age_fmp_amh_thr

  # remove rows without fmp_age (we need the “truth” on menopause)
  if (!"fmp_age" %in% names(df)) {
    stop("Column fmp_age is missing in swan_df — compute it first (FMP age per your logic).")
  }
  swan_ref <- df[, c("SWANID","VISIT","age","glucose","insulin","tg","ldl","hdl","bmi","shbg","AMH","fmp_age")]
  swan_ref
}

## =========================
## 3) Robust normalization within age window
## =========================

robust_z <- function(x) {
  m <- median(x, na.rm=TRUE)
  s <- stats::mad(x, constant = 1.4826, na.rm=TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

# Convert to robust z-scores per column within window
zscore_window <- function(ref_win, cols) {
  out <- ref_win
  for (c in cols) {
    if (c %in% names(out)) out[[c]] <- robust_z(out[[c]])
  }
  out
}

## =========================
## 4) Distance and kNN within age window
## =========================

# Euclidean distance with weights on robust z-scores
dist_weighted <- function(x, M, weights) {
  # x — vector of length p, M — matrix NxP, weights — named vector of length p
  w <- weights[colnames(M)]
  if (any(is.na(w))) {
    # if some columns have no weight (unexpected), set to 1
    w[is.na(w)] <- 1.0
  }
  # scale by weights
  dif <- sweep(M, 2, x, FUN = "-")
  difw <- sweep(dif, 2, sqrt(w), FUN = "*")
  sqrt(rowSums(difw^2))
}

# Select age window (with adjustment for “young”)
age_window_bounds <- function(age, window = AGE_WINDOW) {
  if (age < 43) return(c(42, 44))
  c(age - window, age + window)
}

## =========================
## 5) FMP prediction via nearest SWAN trajectories
## =========================

predict_fmp_knn <- function(nhanes_row, swan_ref, features = FEATURES_MATCH, weights = FEATURE_WEIGHTS,
                            k = K_NEIGHBORS, age_window = AGE_WINDOW, min_neighbors = 5) {
  a <- as.numeric(nhanes_row[["age"]])
  if (is.na(a)) return(NA_real_)

  # age window in SWAN
  bnds <- age_window_bounds(a, window = age_window)
  ref_win <- subset(swan_ref, age >= bnds[1] & age <= bnds[2])
  if (nrow(ref_win) < min_neighbors) return(NA_real_)

  # which features are available for this NHANES row and in SWAN window
  avail <- features[features %in% names(ref_win) & features %in% names(nhanes_row)]
  # drop features where the query has NA
  avail <- avail[!is.na(nhanes_row[avail])]
  if (length(avail) == 0) return(NA_real_)

  # robust standardization within window for avail
  Z <- zscore_window(ref_win, avail)
  # query vector to same scales: compute its z using window (median and MAD)
  z_of <- function(x, ref_vec) {
    m <- median(ref_vec, na.rm=TRUE)
    s <- stats::mad(ref_vec, constant = 1.4826, na.rm=TRUE)
    if (is.na(s) || s==0) return(0)
    (x - m)/s
  }
  xz <- sapply(avail, function(cn) z_of(as.numeric(nhanes_row[[cn]]), ref_win[[cn]]))

  # feature matrix of window (z-scores)
  M <- as.matrix(Z[, avail, drop=FALSE])
  colnames(M) <- avail

  # weights for available features
  w <- weights[avail]; names(w) <- avail

  # distances
  d <- dist_weighted(xz, M, w)
  ord <- order(d)
  k_use <- min(k, sum(is.finite(d)))
  if (k_use < min_neighbors) return(NA_real_)

  nn <- ref_win[ord[seq_len(k_use)], , drop=FALSE]
  dd <- d[ord[seq_len(k_use)]]

  # neighbor weights ~ 1 / (d + eps) so nearer matter more
  eps <- 1e-6
  ww <- 1 / (dd + eps)
  ww <- ww / sum(ww)

  # FMP age prediction
  pred_fmp <- sum(nn$fmp_age * ww, na.rm=TRUE)
  as.numeric(pred_fmp)
}

## =========================
## 6) Apply to NHANES
## =========================

# nhanes_df: data.frame with columns age, AMH (or AMH_pred_cal→AMH), glucose, insulin, tg, ldl, hdl, bmi, shbg
# swan_canon: as in your pipeline, + column fmp_age
# returns nhanes_df with column fmp_pred_knn (without using FSH/TSH)
run_fmp_matching <- function(nhanes_df, swan_canon,
                             features = FEATURES_MATCH,
                             weights  = FEATURE_WEIGHTS,
                             k = K_NEIGHBORS,
                             age_window = AGE_WINDOW) {
  # prepare SWAN reference
  swan_ref <- prepare_swan_reference(swan_canon, amh_col = c("AMH","AMH_pred_cal"))
  # bring AMH in NHANES: if no measured AMH, take modelled
  if (!"AMH" %in% names(nhanes_df)) {
    if ("AMH_pred_cal" %in% names(nhanes_df)) nhanes_df$AMH <- nhanes_df$AMH_pred_cal
  }
  # only women 35–45 (if needed)
  # nhanes_df <- subset(nhanes_df, age >= 35 & age <= 45)

  preds <- mapply(function(i) {
    predict_fmp_knn(nhanes_df[i, , drop=FALSE], swan_ref, features, weights, k, age_window)
  }, seq_len(nrow(nhanes_df)))

  nhanes_df$fmp_pred_knn <- as.numeric(preds)
  nhanes_df
}


nhanes_with_fmp <- run_fmp_matching(nhanes_df, swan_canon)
summary(nhanes_with_fmp$fmp_pred_knn)



names(swan_canon)



## ===== 1) Prepare SWAN reference with target fmp_age (YEARS) =====
## swan_canon: your canonical SWAN table (SI), must contain age and required biomarkers
## merged_by_fpm (or analogue): table where menopause age is already computed by your rule (in YEARS)
## fmp_col: column name with menopause age (YEARS!), e.g. "fmp_age_amh" or "fmp_age_status"

build_swan_reference <- function(swan_canon, merged_by_fpm, fmp_col,
                                 feature_keys = c("age","fsh","glucose","insulin","tg","ldl","hdl","bmi","shbg")) {
  stopifnot(is.data.frame(swan_canon), is.data.frame(merged_by_fpm))
  if (!fmp_col %in% names(merged_by_fpm))
    stop(sprintf("Column %s (FMP age in years) is missing in merged_by_fpm.", fmp_col))

  # take per-person target fmp_age and a “portrait” around the age neighborhood
  # simplest way: join by identifier, then keep only visits with non-NA age and features
  ref <- merge(swan_canon, merged_by_fpm[, c("SWANID", fmp_col)], by = "SWANID", all.x = TRUE)
  names(ref)[names(ref) == fmp_col] <- "fmp_age"

  # filter rows where target age (in years) is present
  ref <- ref[!is.na(ref$fmp_age), , drop = FALSE]

  # keep only columns needed in KNN
  keep <- unique(c("SWANID","VISIT","fmp_age", feature_keys))
  keep <- keep[keep %in% names(ref)]
  ref  <- ref[, keep, drop = FALSE]

  # basic cleaning
  ok_feats <- complete.cases(ref[, feature_keys, drop = FALSE]) & !is.na(ref$fmp_age)
  ref <- ref[ok_feats, , drop = FALSE]

  if (nrow(ref) == 0) stop("Empty SWAN reference after filtering: no rows with a full feature set and fmp_age (years).")
  ref
}

## ===== 2) Predict FMP age (YEARS) via KNN on SWAN =====
## nhanes_df: NHANES table (SI, canonical) with the same feature_keys
## window = ±1 year; for age < 42 use window 42–44

predict_fmp_knn <- function(nhanes_df, swan_ref,
                            feature_keys = c("age","fsh","glucose","insulin","tg","ldl","hdl","bmi","shbg"),
                            k = 20, age_window = 1, force_min_age = 42, force_max_age = 44) {

  stopifnot(is.data.frame(nhanes_df), is.data.frame(swan_ref))

  # keep only rows with complete feature set
  need <- feature_keys
  nh <- nhanes_df[complete.cases(nhanes_df[, need, drop = FALSE]), , drop = FALSE]
  if (nrow(nh) == 0) stop("No NHANES rows have a complete set of features.")

  # standardization: separately for each “local” subsample (by age window)
  # (to avoid dragging global normalization across very different ages)
  knn_pred <- rep(NA_real_, nrow(nh))

  for (i in seq_len(nrow(nh))) {
    a <- nh$age[i]

    # age window
    if (a < force_min_age) {
      min_age <- force_min_age
      max_age <- force_max_age
    } else {
      min_age <- a - age_window
      max_age <- a + age_window
    }

    pool <- swan_ref[swan_ref$age >= min_age & swan_ref$age <= max_age, , drop = FALSE]
    if (nrow(pool) < 5L) {
      # expand window if too few neighbors
      pool <- swan_ref[swan_ref$age >= (min_age - 1) & swan_ref$age <= (max_age + 1), , drop = FALSE]
    }
    if (nrow(pool) < 5L) next  # leave as NA

    # standardize by local pool
    X_pool <- as.matrix(pool[, feature_keys, drop = FALSE])
    X_q    <- as.numeric(nh[i, feature_keys, drop = FALSE])

    mu <- colMeans(X_pool, na.rm = TRUE)
    sdv <- apply(X_pool, 2, sd, na.rm = TRUE); sdv[sdv == 0 | is.na(sdv)] <- 1

    Z_pool <- sweep(X_pool, 2, mu, "-"); Z_pool <- sweep(Z_pool, 2, sdv, "/")
    Z_q    <- (X_q - mu) / sdv

    # distances and top-K
    d <- sqrt(rowSums((Z_pool - matrix(Z_q, nrow=nrow(Z_pool), ncol=length(Z_q), byrow=TRUE))^2))
    ord <- order(d)
    k_use <- min(k, sum(is.finite(d)))
    if (k_use == 0) next

    idx  <- ord[seq_len(k_use)]
    dd   <- d[idx]
    ww   <- 1 / (dd + 1e-6)          # inverse-distance weights
    ww   <- ww / sum(ww)

    fmp_hat <- sum(pool$fmp_age[idx] * ww)
    knn_pred[i] <- fmp_hat
  }

  # assemble result into a data.frame of the same length as input nhanes_df
  out <- nhanes_df
  out$fmp_pred_knn <- NA_real_
  # fill only where features were complete
  out$fmp_pred_knn[match(rownames(nh), rownames(out))] <- knn_pred

  out
}

## ===== 3) Usage example =====
swan_ref <- build_swan_reference(swan_canon, merged_with_fmp, fmp_col = "est_age_fmp_amh_thr")
# nhanes_with_fmp <- predict_fmp_knn(nhanes_df, swan_ref, k=20, age_window=1)

## sanity-check: these should be YEARS, not months.
# summary(nhanes_with_fmp$fmp_pred_knn)
# hist(nhanes_with_fmp$fmp_pred_knn, breaks=30)


names(merged_with_fmp)
