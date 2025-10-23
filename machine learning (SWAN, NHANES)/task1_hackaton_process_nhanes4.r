## =========================================================
## 0) Imports and global settings
## =========================================================
suppressPackageStartupMessages({
          library(glmnet)
          library(mgcv)
          library(ggplot2)
          # we'll call splines::ns with a prefix (no attach)
})

## --- (uncomment if needed) ---
load("nhanes_data_long.rds")  # should create m_long (long-format table)

m_long<-nhanes_data_long

## NHANES waves and variables to keep from long:
use_ds <- c("NHANES_0102","NHANES_0304","NHANES_0506","NHANES_0708","NHANES_0910",
            "NHANES_1112","NHANES_1314","NHANES_1516","NHANES_1718","NHANES_1720",
            "NHANES_2123","NHANES_9900","NHANES_III")

use_variables <- c(
          "RIAGENDR","RIDAGEYR","LBXAMH","LBXFSH","LBXEST","LBXLUH","LBXSHBG","LBXTST",
          "LBDDHESI","LBXHSCRP","LBDLDLSI","LBDHDDSI","LBDTRSI","LBXTSH","BMXBMI","BMXHT",
          "BPXPLS","RHQ060","RHQ010","LBDESOSI","LBDES1SI","LBDPG4SI","LBD17HSI","LBDANDSI",
          "RHD143","RHD043","RHQ031","LBDINSI","LBDGLUSI","SSSHBG","BPXDI2"
)

## Filter incoming long (safe to call even if already filtered):
m_long <- subset(m_long, dataset %in% use_ds & variable %in% use_variables)

## =========================================================
## 1) Unified configuration
## =========================================================
## Map: canonical name → NHANES id
id_map <- c(
          AMH     = "LBXAMH",
          age     = "RIDAGEYR",
          fsh     = "LBXFSH",
          glucose = "LBDGLUSI",
          insulin = "LBDINSI",
          tg      = "LBDTRSI",
          ldl     = "LBDLDLSI",
          hdl     = "LBDHDDSI",
          shbg    = "LBXSHBG",
          bmi     = "BMXBMI",
          # service fields:
          sex              = "RIAGENDR",
          preg_now         = "RHD143",
          last_menses_age  = "RHQ060",
          had_period_12m   = "RHQ031",
          reason_no_period = "RHD043"
)

## Fallbacks (if the primary id is empty):
fallbacks <- list(
          shbg = c("SSSHBG")
)

## Unified list of predictors (used everywhere)
PREDICTOR_KEYS <- c("age","fsh","glucose","insulin","tg","ldl","hdl","bmi","shbg")

## Threshold for two-stage models
AGE_THRESHOLD <- 37

## =========================================================
## 2) Utilities
## =========================================================
getcol <- function(df, var_id) {
          nm <- paste0("value.", var_id)
          if (nm %in% names(df)) as.numeric(df[[nm]]) else rep(NA_real_, nrow(df))
}
winsor <- function(x, p=0.005) {
          if (all(is.na(x))) return(x)
          ql <- suppressWarnings(stats::quantile(x, p, na.rm=TRUE))
          qh <- suppressWarnings(stats::quantile(x, 1-p, na.rm=TRUE))
          x[x < ql] <- ql; x[x > qh] <- qh; x
}
not_all_na <- function(x) !all(is.na(x))

## =========================================================
## 3) Long → Wide + filters (sex/pregnancy/age)
## =========================================================
long_to_wide <- function(m_long, keep_ids, age_range=c(30,55)) {
          keep_ids <- unique(keep_ids[!is.na(keep_ids) & nzchar(keep_ids)])
          m_sub <- subset(m_long, variable %in% keep_ids)

          female_ids <- unique(subset(m_sub, variable==id_map["sex"] & value==2)$SEQN)
          m_sub <- m_sub[m_sub$SEQN %in% female_ids, ]

          preg_ids <- unique(subset(m_sub, variable==id_map["preg_now"] & value==1)$SEQN)
          m_sub <- m_sub[!(m_sub$SEQN %in% preg_ids), ]

          age_rows <- subset(m_sub, variable==id_map["age"], select=c("SEQN","value"))
          ok_ids <- unique(age_rows$SEQN[age_rows$value >= age_range[1] & age_rows$value <= age_range[2]])
          m_sub <- m_sub[m_sub$SEQN %in% ok_ids, ]

          m_agg <- stats::aggregate(value ~ SEQN + dataset + variable, data=m_sub, FUN=function(x) x[1])
          stats::reshape(m_agg, idvar=c("SEQN","dataset"), timevar="variable", direction="wide")
}

## =========================================================
## 4) Status labeling (optional; for debugging/stratification)
## =========================================================
label_status <- function(w) {
          RIDAGEYR <- getcol(w, id_map["age"])
          RHQ031   <- getcol(w, id_map["had_period_12m"])
          RHD043   <- getcol(w, id_map["reason_no_period"])
          RHQ060   <- getcol(w, id_map["last_menses_age"])
          RHD143   <- getcol(w, id_map["preg_now"])

          RHQ060_clean <- ifelse(is.na(RHQ060), NA,
                                 ifelse(RHQ060 %in% c(777,999), NA,
                                        ifelse(RHQ060 >= 0 & RHQ060 <= 150, RHQ060, NA)))
          years_since_FMP <- ifelse(!is.na(RHQ060_clean) & !is.na(RIDAGEYR),
                                    RIDAGEYR - RHQ060_clean, NA)
          years_since_FMP[years_since_FMP < 0 | years_since_FMP > RIDAGEYR] <- NA

          is_menopause_reason <- RHD043 == 7
          is_hysterectomy     <- RHD043 == 3
          is_preg_bf_reason   <- RHD043 %in% c(1,2)

          status <- rep(NA_character_, nrow(w))
          exclude_preg_now <- RHD143 == 1
          status[exclude_preg_now] <- NA_character_

          status[is.na(status) & RHQ031 == 1] <- "pre"
          status[is.na(status) & RHQ031 == 2 & is_menopause_reason] <- "post"

          surgical_flag <- RHQ031 == 2 & is_hysterectomy
          status[is.na(status) & surgical_flag] <- NA_character_

          exclude_pregbf <- RHQ031 == 2 & is_preg_bf_reason
          status[is.na(status) & exclude_pregbf] <- NA_character_

          status[is.na(status) & !is.na(years_since_FMP) & years_since_FMP >= 1] <- "post"
          status[is.na(status) & !is.na(years_since_FMP) & years_since_FMP >= 0 & years_since_FMP < 1] <- "peri"

          w$years_since_FMP  <- years_since_FMP
          w$surgical_oophyst <- as.integer(surgical_flag)
          w$exclude_preg_now <- as.integer(exclude_preg_now | exclude_pregbf)
          w$menopausal_status<- status
          w
}

## =========================================================
## 5) Feature assembly from wide (uses PREDICTOR_KEYS)
## =========================================================
assemble_features <- function(w, id_map, fallbacks=list(), predictor_keys=PREDICTOR_KEYS, drop_incomplete=FALSE) {
          n <- nrow(w)
          take_by_id <- function(id) {
                    if (is.null(id) || is.na(id) || length(id)==0) return(rep(NA_real_, n))
                    getcol(w, id)
          }
          ALL_KEYS <- unique(c("AMH", predictor_keys))
          vals <- setNames(vector("list", length(ALL_KEYS)), ALL_KEYS)

          for (k in ALL_KEYS) {
                    base_id <- unname(id_map[k])
                    v <- take_by_id(base_id)
                    if (k %in% names(fallbacks)) {
                              for (fid in as.character(fallbacks[[k]])) {
                                        if (all(is.na(v))) {
                                                  v_fb <- take_by_id(fid)
                                                  if (!all(is.na(v_fb))) { v <- v_fb; break }
                                        }
                              }
                    }
                    vals[[k]] <- v
          }

          df <- as.data.frame(vals, stringsAsFactors = FALSE)
          df <- df[!is.na(df$AMH) & df$AMH > 0 & df$AMH < 10, ]
          for (v in predictor_keys) if (v %in% names(df)) df[[v]] <- winsor(df[[v]])

          if (drop_incomplete) {
                    needed <- unique(c("AMH", predictor_keys))
                    df <- df[stats::complete.cases(df[, needed, drop=FALSE]), , drop=FALSE]
          }
          df
}

## =========================================================
## 6) Design matrix and train/test split
## =========================================================
mk_design <- function(df, predictor_keys=PREDICTOR_KEYS) {
          for (v in predictor_keys) if (!v %in% names(df)) df[[v]] <- NA_real_
          age_s <- if ("age" %in% predictor_keys) splines::ns(df$age, df=4) else NULL
          fsh_s <- if ("fsh" %in% predictor_keys) splines::ns(df$fsh, df=4) else NULL
          if (!is.null(age_s)) colnames(age_s) <- paste0("s_age_", 1:ncol(age_s))
          if (!is.null(fsh_s)) colnames(fsh_s) <- paste0("s_fsh_", 1:ncol(fsh_s))
          inter1 <- if (!is.null(age_s) && !is.null(fsh_s)) age_s[,1] * fsh_s[,1] else NULL
          lin_keys <- setdiff(predictor_keys, c("age","fsh"))
          X_lin <- if (length(lin_keys)) as.matrix(df[, lin_keys, drop=FALSE]) else NULL
          parts <- list(age_s, fsh_s, inter1=inter1, X_lin)
          parts <- parts[!sapply(parts, is.null)]
          as.matrix(do.call(cbind, parts))
}

train_test_split <- function(df, seed=42, frac=0.8, predictor_keys=PREDICTOR_KEYS) {
          needed <- unique(c("AMH", predictor_keys))
          keep <- stats::complete.cases(df[, needed, drop=FALSE])
          df2 <- df[keep, , drop=FALSE]
          set.seed(seed)
          n <- nrow(df2); if (n < 10) stop(sprintf("Too few observations after filtering: n=%d", n))
          idx_tr <- sample.int(n, floor(frac*n))
          list(train=df2[idx_tr, , drop=FALSE], test=df2[-idx_tr, , drop=FALSE])
}

## =========================================================
## 7) Baseline models (linear scale)
## =========================================================
fit_lm <- function(train, predictor_keys=PREDICTOR_KEYS) {
          have <- predictor_keys[predictor_keys %in% names(train) & sapply(train[predictor_keys], not_all_na)]
          if (length(have)==0) have <- intersect(c("age","fsh"), names(train))
          form <- stats::as.formula(paste0("AMH ~ ", paste(have, collapse=" + ")))
          list(model = stats::lm(form, data=train), predictors=have)
}
fit_enet <- function(train, predictor_keys=PREDICTOR_KEYS) {
          X <- mk_design(train, predictor_keys); y <- train$AMH
          alphas <- c(0.0, 0.25, 0.5, 0.75, 1.0)
          cvs <- lapply(alphas, function(a) glmnet::cv.glmnet(X, y, family="gaussian", alpha=a, standardize=TRUE, nfolds=5))
          best_i <- which.min(sapply(cvs, function(cv) min(cv$cvm)))
          list(model = glmnet::glmnet(X, y, family="gaussian", alpha=alphas[best_i],
                                      lambda=cvs[[best_i]]$lambda.min, standardize=TRUE),
               alpha=alphas[best_i], lambda=cvs[[best_i]]$lambda.min, cols=colnames(X))
}
fit_gam <- function(train, predictor_keys=PREDICTOR_KEYS) {
          have <- intersect(predictor_keys, names(train))
          rhs <- c()
          if ("age" %in% have) rhs <- c(rhs, "s(age, k=6)")
          if ("fsh" %in% have) rhs <- c(rhs, "s(fsh, k=6)")
          if (all(c("age","fsh") %in% have)) rhs <- c(rhs, "ti(age, fsh, k=c(5,5))")
          rhs <- c(rhs, setdiff(have, c("age","fsh")))
          form <- stats::as.formula(paste("AMH ~", paste(rhs, collapse=" + ")))
          mgcv::gam(form, data=train, method="REML")
}

## =========================================================
## 8) Positive-scale models (AMH > 0)
## =========================================================
## EN on log(AMH) + Duan smearing
duan_smearing <- function(residuals_log) mean(exp(residuals_log), na.rm = TRUE)

## === helpers for spline df selection ===
.unique_n <- function(x) length(unique(stats::na.omit(x)))

.pick_k <- function(x, k_max = 10) {
          u <- .unique_n(x)
          if (is.na(u) || u <= 2) return(0)      # too few for a spline → use linear term
          # conservative cap: no more than u-1 and k_max
          k <- min(k_max, max(3, u - 1))
          k
}

## === adaptive Gamma-GAM with log link ===
fit_gam_gamma <- function(train, predictor_keys = PREDICTOR_KEYS) {
          have <- intersect(predictor_keys, names(train))

          # adaptive k based on data
          k_age <- if ("age" %in% have) .pick_k(train$age, k_max = 10) else 0
          k_fsh <- if ("fsh" %in% have) .pick_k(train$fsh, k_max = 10) else 0

          rhs <- c()

          # age: spline s(age,k) or linear age or nothing
          if ("age" %in% have) {
                    if (k_age >= 3) rhs <- c(rhs, sprintf("s(age, k=%d)", k_age)) else rhs <- c(rhs, "age")
          }
          # fsh: spline or linear
          if ("fsh" %in% have) {
                    if (k_fsh >= 3) rhs <- c(rhs, sprintf("s(fsh, k=%d)", k_fsh)) else rhs <- c(rhs, "fsh")
          }

          # tensor interaction only if both splines are valid
          if (k_age >= 3 && k_fsh >= 3) {
                    rhs <- c(rhs, sprintf("ti(age, fsh, k=c(%d,%d))", k_age, k_fsh))
          }

          # remaining covariates — linear
          rhs <- c(rhs, setdiff(have, c("age","fsh")))

          # fallback if nothing left
          if (length(rhs) == 0) rhs <- c("1")

          form <- stats::as.formula(paste("AMH ~", paste(rhs, collapse = " + ")))

          mgcv::gam(
                    form,
                    data    = train,
                    family  = stats::Gamma(link = "log"),
                    method  = "REML",
                    select  = TRUE   # penalize redundant smooths
          )
}

fit_enet_log <- function(train, predictor_keys=PREDICTOR_KEYS) {
          X <- mk_design(train, predictor_keys)
          y_log <- log(train$AMH)
          alphas <- c(0.0, 0.25, 0.5, 0.75, 1.0)
          cvs <- lapply(alphas, function(a) glmnet::cv.glmnet(X, y_log, family="gaussian",
                                                              alpha=a, standardize=TRUE, nfolds=5))
          best_i <- which.min(sapply(cvs, function(cv) min(cv$cvm)))
          fit <- glmnet::glmnet(X, y_log, family="gaussian", alpha=alphas[best_i],
                                lambda=cvs[[best_i]]$lambda.min, standardize=TRUE)
          yhat_tr <- as.numeric(stats::predict(fit, newx = X))
          smear <- duan_smearing(y_log - yhat_tr)
          list(model=fit, alpha=alphas[best_i], lambda=cvs[[best_i]]$lambda.min,
               cols=colnames(X), smear=smear)
}
predict_enet_linear <- function(fit_enet_log, newdf, predictor_keys=PREDICTOR_KEYS) {
          Xte <- mk_design(newdf, predictor_keys)
          pred_log <- as.numeric(stats::predict(fit_enet_log$model, newx = Xte))
          p <- exp(pred_log) * fit_enet_log$smear
          p[p < 0] <- 0
          p
}

predict_gam_gamma <- function(fit_gam, newdf) {
          p <- as.numeric(stats::predict(fit_gam, newdata=newdf, type="response"))
          p[p < 0] <- 0
          p
}

## (optional) Quantile regression τ=0.9 + blend with EN — keep as helper
fit_rq90 <- function(train, predictor_keys=PREDICTOR_KEYS) {
          if (!requireNamespace("quantreg", quietly = TRUE))
                    stop("Package 'quantreg' is not installed")
          have <- intersect(predictor_keys, names(train))
          form <- stats::as.formula(paste("AMH ~", paste(have, collapse=" + ")))
          quantreg::rq(form, data=train, tau=0.9)
}
predict_rq <- function(fit_rq, newdf) as.numeric(stats::predict(fit_rq, newdata=newdf))
blend_high_tail <- function(median_like, upper_like, w = 0.3) (1 - w) * median_like + w * upper_like

## =========================================================
## 9) Two-stage models (threshold 37 years)
## =========================================================
split_by_age <- function(df, thr = AGE_THRESHOLD) {
          list(early = df[df$age <  thr, , drop=FALSE],
               late  = df[df$age >= thr, , drop=FALSE])
}

## === two-stage Gamma-GAM with adaptive k ===
fit_two_stage_gam_gamma <- function(train, predictor_keys = PREDICTOR_KEYS, thr = AGE_THRESHOLD) {
          parts <- list(
                    early = train[train$age <  thr, , drop = FALSE],
                    late  = train[train$age >= thr, , drop = FALSE]
          )

          # minimal checks to avoid failures on tiny subsets
          if (nrow(parts$early) < 20 || nrow(parts$late) < 20) {
                    warning(sprintf(
                              "Few observations in a subgroup (<20): early=%d, late=%d. Training anyway; quality may degrade.",
                              nrow(parts$early), nrow(parts$late)
                    ))
          }

          m_early <- fit_gam_gamma(parts$early, predictor_keys)
          m_late  <- fit_gam_gamma(parts$late,  predictor_keys)

          list(type="two_stage_gam_gamma", thr=thr, early=m_early, late=m_late, predictors=predictor_keys)
}
predict_two_stage_gam_gamma <- function(fit2, newdf) {
          stopifnot(fit2$type == "two_stage_gam_gamma")
          thr <- fit2$thr
          idx_e <- which(newdf$age <  thr)
          idx_l <- which(newdf$age >= thr)
          out <- rep(NA_real_, nrow(newdf))
          if (length(idx_e)) out[idx_e] <- predict_gam_gamma(fit2$early, newdf[idx_e, , drop=FALSE])
          if (length(idx_l)) out[idx_l] <- predict_gam_gamma(fit2$late,  newdf[idx_l, , drop=FALSE])
          out[out < 0] <- 0
          out
}

## Two-stage EN (log) + Duan
fit_two_stage_enet_log <- function(train, predictor_keys = PREDICTOR_KEYS, thr = AGE_THRESHOLD) {
          parts <- split_by_age(train, thr)
          if (nrow(parts$early) < 20 || nrow(parts$late) < 20)
                    stop("Too few observations in one age subgroup for a two-stage model.")
          m_early <- fit_enet_log(parts$early, predictor_keys)
          m_late  <- fit_enet_log(parts$late,  predictor_keys)
          list(type="two_stage_enet_log", thr=thr, early=m_early, late=m_late, predictors=predictor_keys)
}
predict_two_stage_enet_linear <- function(fit2, newdf, predictor_keys = PREDICTOR_KEYS) {
          stopifnot(fit2$type == "two_stage_enet_log")
          thr <- fit2$thr
          idx_e <- which(newdf$age <  thr)
          idx_l <- which(newdf$age >= thr)
          out <- rep(NA_real_, nrow(newdf))
          if (length(idx_e)) out[idx_e] <- predict_enet_linear(fit2$early, newdf[idx_e, , drop=FALSE], predictor_keys)
          if (length(idx_l)) out[idx_l] <- predict_enet_linear(fit2$late,  newdf[idx_l, , drop=FALSE], predictor_keys)
          out[out < 0] <- 0
          out
}

## Evaluation of two-stage models
evaluate_two_stage <- function(models2, test, predictor_keys = PREDICTOR_KEYS) {
          mae  <- function(a,b) mean(abs(a-b))
          rmse <- function(a,b) sqrt(mean((a-b)^2))
          p_gam2 <- predict_two_stage_gam_gamma(models2$gam2,  test)
          p_en2  <- predict_two_stage_enet_linear(models2$enet2, test, predictor_keys)
          cat(sprintf("[2-stage GAM Γ(log)] MAE=%.4f  RMSE=%.4f\n", mae(p_gam2, test$AMH), rmse(p_gam2, test$AMH)))
          cat(sprintf("[2-stage EN log→lin] MAE=%.4f  RMSE=%.4f\n", mae(p_en2,  test$AMH), rmse(p_en2, test$AMH)))
          invisible(list(gam2 = p_gam2, en2 = p_en2))
}

## =========================================================
## 10) Metrics and plots
## =========================================================
mae  <- function(a,b) mean(abs(a-b))
rmse <- function(a,b) sqrt(mean((a-b)^2))

evaluate_models_linear <- function(models, test, predictor_keys=PREDICTOR_KEYS) {
          p_lm <- as.numeric(stats::predict(models$lm$model,   newdata=test))
          p_en <- as.numeric(stats::predict(models$enet$model, newx = mk_design(test, predictor_keys)))
          p_gm <- as.numeric(stats::predict(models$gam,        newdata=test))
          cat(sprintf("[LM ] Test: MAE=%.4f  RMSE=%.4f\n", mae(p_lm, test$AMH), rmse(p_lm, test$AMH)))
          cat(sprintf("[EN ] Test: MAE=%.4f  RMSE=%.4f\n", mae(p_en, test$AMH), rmse(p_en, test$AMH)))
          cat(sprintf("[GAM] Test: MAE=%.4f  RMSE=%.4f\n", mae(p_gm, test$AMH), rmse(p_gm, test$AMH)))
          invisible(list(lm=p_lm, en=p_en, gam=p_gm))
}

evaluate_models_positive <- function(models, test, predictor_keys=PREDICTOR_KEYS) {
          p_en  <- predict_enet_linear(models$enet_log, test, predictor_keys)
          p_gam <- predict_gam_gamma(models$gam_gamma, test)
          p_out <- list(en=p_en, gam=p_gam)
          if (!is.null(models$rq90)) {
                    p_q90 <- predict_rq(models$rq90, test)
                    p_blend <- blend_high_tail(p_en, p_q90, w = 0.3)
                    p_out$rq90  <- p_q90
                    p_out$blend <- p_blend
          }
          cat(sprintf("[EN  log→lin+smear]  MAE=%.4f  RMSE=%.4f\n", mae(p_en,  test$AMH), rmse(p_en,  test$AMH)))
          cat(sprintf("[GAM Gamma(log)]     MAE=%.4f  RMSE=%.4f\n", mae(p_gam, test$AMH), rmse(p_gam, test$AMH)))
          if (!is.null(p_out$blend))
                    cat(sprintf("[Blend EN + RQ90]    MAE=%.4f  RMSE=%.4f\n", mae(p_out$blend, test$AMH), rmse(p_out$blend, test$AMH)))
          invisible(p_out)
}

plot_obs_vs_pred <- function(test, preds, title, n_sample=100) {
          dfp <- data.frame(AMH_obs=test$AMH, AMH_pred=preds, age=test$age)
          set.seed(123)
          idx <- sample(seq_len(nrow(dfp)), size=min(n_sample, nrow(dfp)), replace=FALSE)
          dfp <- dfp[idx, ]
          ggplot2::ggplot(dfp, ggplot2::aes(x=AMH_obs, y=AMH_pred, color=age)) +
                    ggplot2::geom_abline(slope=1, intercept=0, linetype="dashed") +
                    ggplot2::geom_point(size=2, alpha=0.85) +
                    ggplot2::labs(title=title, x="AMH (observed), ng/mL", y="AMH (predicted), ng/mL", color="Age, years") +
                    ggplot2::theme_minimal()
}

## =========================================================
## 11) PIPELINE
## =========================================================
run_pipeline <- function(m_long, id_map, fallbacks=list(),
                         predictor_keys=PREDICTOR_KEYS, age_range=c(30,55),
                         use_quantile_blend=TRUE) {

          keep_ids <- unique(c(unname(id_map), unlist(fallbacks, use.names = FALSE)))
          w  <- long_to_wide(m_long, keep_ids, age_range)
          w  <- label_status(w)

          df <- assemble_features(w, id_map, fallbacks, predictor_keys, drop_incomplete = TRUE)

          split <- train_test_split(df, seed=42, frac=0.8, predictor_keys=predictor_keys)
          train <- split$train;  test <- split$test

          ## --- Linear-scale models (for comparison)
          mdl_lm   <- fit_lm(train, predictor_keys)
          mdl_enet <- fit_enet(train, predictor_keys)
          mdl_gam  <- fit_gam(train, predictor_keys)

          cat("\n=== HOLD-OUT (linear scale) ===\n")
          preds_lin <- evaluate_models_linear(list(lm=mdl_lm, enet=mdl_enet, gam=mdl_gam), test, predictor_keys)
          print(plot_obs_vs_pred(test, preds_lin$lm,  "[LM]   AMH: obs vs pred"))
          print(plot_obs_vs_pred(test, preds_lin$en,  "[EN]   AMH: obs vs pred"))
          print(plot_obs_vs_pred(test, preds_lin$gam, "[GAM]  AMH: obs vs pred"))

          ## --- Positive-scale models
          mdl_enet_log  <- fit_enet_log(train, predictor_keys)
          mdl_gam_gamma <- fit_gam_gamma(train, predictor_keys)
          mdl_rq90      <- NULL
          if (use_quantile_blend && requireNamespace("quantreg", quietly = TRUE)) {
                    mdl_rq90 <- fit_rq90(train, predictor_keys)
          }

          cat("\n=== HOLD-OUT (positive-scale) ===\n")
          preds_pos <- evaluate_models_positive(
                    list(enet_log=mdl_enet_log, gam_gamma=mdl_gam_gamma, rq90=mdl_rq90),
                    test, predictor_keys
          )
          print(plot_obs_vs_pred(test, preds_pos$en,    "[EN log→lin+smear] AMH: obs vs pred"))
          print(plot_obs_vs_pred(test, preds_pos$gam,   "[GAM Gamma(log)]  AMH: obs vs pred"))
          if (!is.null(preds_pos$blend))
                    print(plot_obs_vs_pred(test, preds_pos$blend, "[Blend EN+RQ90]   AMH: obs vs pred"))

          ## --- Two-stage (threshold 37 years)
          mdl_gam2  <- fit_two_stage_gam_gamma(train, predictor_keys, thr = AGE_THRESHOLD)
          mdl_enet2 <- fit_two_stage_enet_log(train,  predictor_keys, thr = AGE_THRESHOLD)

          cat("\n=== HOLD-OUT (two-stage @ 37y) ===\n")
          preds_2stage <- evaluate_two_stage(list(gam2 = mdl_gam2, enet2 = mdl_enet2), test, predictor_keys)
          print(plot_obs_vs_pred(test, preds_2stage$gam2, "[2-stage GAM Γ(log)] AMH: obs vs pred"))
          print(plot_obs_vs_pred(test, preds_2stage$en2,  "[2-stage EN log→lin] AMH: obs vs pred"))

          invisible(list(
                    wide=w, features=df, train=train, test=test,
                    lm=mdl_lm, enet=mdl_enet, gam=mdl_gam,
                    enet_log=mdl_enet_log, gam_gamma=mdl_gam_gamma, rq90=mdl_rq90,
                    gam2=mdl_gam2, enet2=mdl_enet2,
                    preds_linear=preds_lin, preds_positive=preds_pos, preds_2stage=preds_2stage,
                    predictor_keys=predictor_keys
          ))
}

## =========================================================
## 12) Run
## =========================================================
res <- run_pipeline(
          m_long,
          id_map      = id_map,
          fallbacks   = fallbacks,
          predictor_keys = PREDICTOR_KEYS,
          age_range   = c(30,55),
          use_quantile_blend = TRUE
)

##########################

# After training:
two_stage <- fit_two_stage_gam_gamma(train, predictor_keys = PREDICTOR_KEYS, thr = AGE_THRESHOLD)

# Make a validation split from train (e.g., 80/20 inside train)
set.seed(7)
idx_val <- sample.int(nrow(train), floor(0.2*nrow(train)))
train_sub <- train[-idx_val, , drop=FALSE]
val       <- train[ idx_val, , drop=FALSE]

# Refit on train_sub, calibrate on val
two_stage_cal <- fit_two_stage_gam_gamma(train_sub, predictor_keys = PREDICTOR_KEYS, thr = AGE_THRESHOLD)
pred_val_raw  <- predict_two_stage_gam_gamma(two_stage_cal, val)

# Isotonic calibration
if (!requireNamespace("isotone", quietly=TRUE)) install.packages("isotone")
iso_fit <- function(pred, obs) {
          o <- order(pred); pred_s <- pred[o]; obs_s <- obs[o]
          cal <- isotone::gpava(pred_s, obs_s)$x
          list(pred_sorted = pred_s, cal_sorted = cal)
}
iso <- iso_fit(pred_val_raw, val$AMH)

# Final refit on FULL train (as requested)
two_stage_final <- fit_two_stage_gam_gamma(train, predictor_keys = PREDICTOR_KEYS, thr = AGE_THRESHOLD)

# Save everything as a single bundle
bundle <- list(
          model        = two_stage_final,          # two-stage GAM Γ(log)
          threshold    = AGE_THRESHOLD,            # age threshold
          predictor_keys = PREDICTOR_KEYS,         # predictor list
          iso          = iso,                      # isotonic calibration
          version      = "amh_two_stage_gam_gamma_v1"
)
dir.create("models", showWarnings = FALSE)
saveRDS(bundle, file = "models/amh_two_stage_gam_gamma_bundle.rds")
