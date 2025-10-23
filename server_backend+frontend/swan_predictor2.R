#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(mgcv)   # predict.gam
})

# --- только JSON наружу ---
json_fail <- function(msg, status = 1) {
  cat(toJSON(list(error = msg), auto_unbox = TRUE))
  quit(status = status)
}

# ВОНИ отключаем; будем избегать их источников, но на всякий
options(warn = -1, show.error.locations = FALSE)


# --- безопасное логирование ---
log_dir <- Sys.getenv("SWAN_LOG_DIR", unset = "logs")
writable <- function(d) dir.exists(d) && suppressWarnings(file.access(d, 2) == 0)
if (!dir.exists(log_dir)) suppressWarnings(dir.create(log_dir, recursive = TRUE, showWarnings = FALSE))
if (!writable(log_dir)) log_dir <- tempdir()  # silent fallback

ts_now <- function() format(Sys.time(), "%Y%m%d_%H%M%S")   # <--- ДОБАВЬ ЭТО

logwrite <- function(name, txt) {
  fn <- file.path(log_dir, sprintf("%s_%s.log", name, ts_now()))
  invisible(try(writeLines(txt, fn), silent = TRUE))
}






# --- глобальный обработчик ошибок: лог внутрь, наружу — JSON ---
options(error = function(e) {
  logwrite("r_error", paste0(conditionMessage(e), "\n"))
  json_fail("internal_error")  # чистый JSON для PHP
})


suppressPackageStartupMessages({
  library(jsonlite)
  library(mgcv)  # нужен для predict.gam
})

# 1) Считывание JSON безопасно (NULL -> NA)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript swan_predictor2.R '<json>'")
raw_json <- args[1]
payload <- tryCatch(fromJSON(raw_json, simplifyVector = FALSE), error = function(e) NULL)
if (is.null(payload)) stop("Invalid JSON")

keys <- c("age","fsh","glucose","insulin","tg","ldl","hdl","bmi","shbg","amh")
to_row <- function(lst, keys) {
  out <- setNames(vector("list", length(keys)), keys)
  for (k in keys) {
    v <- lst[[k]]
    if (is.null(v) || length(v) == 0) {
      out[[k]] <- NA_real_
    } else {
      nv <- suppressWarnings(as.numeric(v))
      out[[k]] <- if (is.finite(nv)) nv else NA_real_
    }
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}
person <- to_row(payload, keys)

# Лог входа
logwrite("input_raw", raw_json)

# 2) Приведение к NHANES SI (без «умных» эвристик, только мягко)
safe_fin <- function(x) isTRUE(is.finite(x))
normalize_to_nhanes_si <- function(df) {
  x <- as.list(df)
  if (safe_fin(x$glucose) && x$glucose > 20) x$glucose <- x$glucose / 18.0
  if (safe_fin(x$tg)      && x$tg      > 10) x$tg      <- x$tg      / 88.57
  if (safe_fin(x$ldl)     && x$ldl     > 15) x$ldl     <- x$ldl     / 38.67
  if (safe_fin(x$hdl)     && x$hdl     > 15) x$hdl     <- x$hdl     / 38.67
  if (safe_fin(x$insulin) && x$insulin < 60) x$insulin <- x$insulin * 6.0
  as.data.frame(x, stringsAsFactors = FALSE)
}
person <- normalize_to_nhanes_si(person)

# 3) Если AMH = NA → предсказать из бандла
predict_gam_gamma <- function(fit, newdf) {
  as.numeric(mgcv::predict.gam(fit, newdata=newdf, type="response", newdata.guaranteed=TRUE))
}
predict_two_stage_gam_gamma <- function(fit2, newdf) {
  stopifnot(fit2$type == "two_stage_gam_gamma")
  thr <- fit2$thr
  idx_e <- which(newdf$age <  thr)
  idx_l <- which(newdf$age >= thr)
  out <- rep(NA_real_, nrow(newdf))
  if (length(idx_e)) out[idx_e] <- predict_gam_gamma(fit2$early, newdf[idx_e,,drop=FALSE])
  if (length(idx_l)) out[idx_l] <- predict_gam_gamma(fit2$late,  newdf[idx_l,,drop=FALSE])
  out[out < 0] <- 0
  out
}
predict_amh <- function(bundle, newdf) {
  pred_raw <- predict_two_stage_gam_gamma(bundle$model, newdf)
  # изотонная калибровка
  approx(x=bundle$iso$pred_sorted, y=bundle$iso$cal_sorted, xout=pred_raw,
         ties="ordered", rule=2)$y |> pmax(0)
}

bundle_path <- file.path("models","amh_two_stage_gam_gamma_bundle.rds")
bundle <- if (file.exists(bundle_path)) readRDS(bundle_path) else NULL

amh_imputed <- NA_real_
if (is.na(person$amh)) {
  if (is.null(bundle)) stop("No AMH in input and bundle missing")
  amh_imputed <- as.numeric(predict_amh(bundle, person))
  person$amh <- amh_imputed
}

# Лог нормализованного входа
logwrite("input_norm", toJSON(as.list(person), pretty=TRUE, auto_unbox=TRUE))

# 4) Загрузить swan_ref и взять VISIT==5
if (!file.exists("swan_ref.bin")) stop("swan_ref.bin not found")
load("swan_ref.bin")               # объект swan_ref
if (!exists("swan_ref")) stop("swan_ref not found in bin")
ref <- subset(swan_ref, VISIT == 5)
if (!all(c("age","AMH","fmp_age") %in% names(ref))) {
  stop("swan_ref must have age, AMH, fmp_age")
}
ref <- ref[is.finite(ref$age)&is.finite(ref$AMH)&is.finite(ref$fmp_age), c("age","AMH","fmp_age")]
if (nrow(ref) < 50) stop("swan_ref too small after filtering")

# 5) Простой двухстадийный kNN (жёстко и коротко)
predict_fmp_two_stage <- function(q, ref, top_k=300, k_final=5, sigma=1.5) {
  za_m <- mean(ref$age); za_s <- sd(ref$age); if (!is.finite(za_s) || za_s == 0) za_s <- 1
  zm_m <- mean(ref$AMH); zm_s <- sd(ref$AMH); if (!is.finite(zm_s) || zm_s == 0) zm_s <- 1

  z_age <- (ref$age - za_m)/za_s
  z_amh <- (ref$AMH - zm_m)/zm_s
  qz_age <- (q$age - za_m)/za_s
  qz_amh <- (q$amh - zm_m)/zm_s

  d <- sqrt((z_age - qz_age)^2 + (z_amh - qz_amh)^2)
  ord <- order(d)
  idx <- ord[seq_len(min(length(ord), top_k))]

  d2 <- d[idx]
  y  <- ref$fmp_age[idx]
  ord2 <- order(d2)
  idx2 <- idx[ord2[seq_len(min(length(idx), k_final))]]
  d3 <- d[idx2]; y2 <- ref$fmp_age[idx2]

  w <- exp(-0.5*(d3/sigma)^2); w <- ifelse(is.finite(w), w, 0)
  if (sum(w)==0) w <- rep(1, length(w))
  w <- w/sum(w)

  yhat <- sum(w*y2)
  varw <- sum(w*(y2 - yhat)^2)
  neff <- (sum(w))^2/sum(w^2)
  se <- sqrt(varw/max(neff,1))
  z975 <- qnorm(0.975)

  list(est=yhat, lo=yhat - z975*se, hi=yhat + z975*se,
       neigh=data.frame(age=ref$age[idx2], AMH=ref$AMH[idx2], fmp_age=y2, dist=d3, weight=w))
}

fmp <- predict_fmp_two_stage(person, ref, top_k=300, k_final=5, sigma=1.5)
est_age <- as.numeric(fmp$est)
ci_low  <- as.numeric(fmp$lo)
ci_high <- as.numeric(fmp$hi)

# Лог соседей
invisible(try(
  write.csv(fmp$neigh, file = file.path(log_dir, paste0("neighbors_", ts_now(), ".csv")), row.names = FALSE),
  silent = TRUE
))


# 6) Симптомы (минимум кода)
predict_with_ci <- function(model, newdata) {
  pr <- predict(model, newdata=newdata, type="link", se.fit=TRUE)
  z <- qnorm(0.975); lo <- plogis(pr$fit - z*pr$se.fit); hi <- plogis(pr$fit + z*pr$se.fit)
  data.frame(p=plogis(pr$fit), low=lo, high=hi)
}
risk_level <- function(p) ifelse(p<0.33,"Low", ifelse(p<0.66,"Medium","High"))

person$est_age_fmp_amh_thr <- est_age
person$age_minus_thr <- person$age - est_age
person$is_post_thr   <- as.integer(person$age >= est_age)

symptom_files <- list.files("models", pattern="^model_.*\\.rds$", full.names=TRUE)
models <- setNames(lapply(symptom_files, readRDS),
                   gsub("^model_|\\.rds$","", basename(symptom_files)))

symptom_results <- list()
for (nm in names(models)) {
  ci <- predict_with_ci(models[[nm]], person)
  symptom_results[[nm]] <- list(
    probability = round(ci$p,3),
    ci_low      = round(ci$low,3),
    ci_high     = round(ci$high,3),
    risk_level  = unname(risk_level(ci$p))
  )
}

# 7) Ответ
out <- list(
  input = c(as.list(person[, keys]), amh_imputed = if (is.finite(amh_imputed)) amh_imputed else NULL),
  predicted_menopause_age = round(est_age,2),
  ci_low  = round(ci_low,2),
  ci_high = round(ci_high,2),
  symptoms = symptom_results
)

# Лог и вывод
j <- toJSON(out, pretty=TRUE, auto_unbox=TRUE)
logwrite("output", j)
cat(j)
