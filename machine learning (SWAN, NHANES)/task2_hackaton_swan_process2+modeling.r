## Universal variable label collector
make_meta <- function(df) {
          # 1) per-column labels (haven / labelled)
          per_col <- vapply(df, function(x) {
                    lab <- attr(x, "label", exact = TRUE)
                    if (is.null(lab)) NA_character_ else as.character(lab)
          }, character(1))

          # 2) data frame–level labels (foreign::read.spss, etc.)
          df_lab <- attr(df, "variable.labels", exact = TRUE)
          if (!is.null(df_lab)) {
                    # align to column names order
                    df_lab <- as.character(df_lab[ names(df) ])
          } else {
                    df_lab <- rep(NA_character_, length(df))
          }

          # priority: per_col, then df_lab
          desc <- per_col
          need  <- is.na(desc) | desc == ""
          desc[need] <- df_lab[need]

          data.frame(
                    column_name = names(df),
                    description = desc,
                    stringsAsFactors = FALSE
          )
}

m1<-make_meta(da28762.0001)

write.csv(m1,file = "m1.csv",row.names = F)


#    AMH age    fsh glucose  insulin    tg     ldl    hdl    bmi    shbg

id_id<-c("SWANID","VISIT","INTDAY0","AGE0")
id_blood_markers<-c("FSH0","GLUCRES0","INSURES0","TRIGRES0","LDLRESU0","HDLRESU0","BMI0","SHBG0")
id_mentsr_status<-c("FLOWDAY0","TENDAFL0","FLOWAMT0","FLOODIN0","CLOTS0","STARTDA0","USUALDA0","FLOAGE20","LMPDAY0")
id_symptoms<-c("HOTFLAS0","RESTLES0","TRBLSLE0","TYPNIGH0","MOODCHN0","MOODCHG0","LEAKURI0","DAYSLEA0","DEPRESS0","BREASTP0","PELVIC0","AVCIGDA0","ALCHWK0","ALCHMON0")

id_all<-c(id_id,id_blood_markers,id_symptoms)

di<-da28762.0001[,id_all]

mi<-m1[id_all,]

row.names(mi)<-1:nrow(mi)


dim(di)


## =========================
## SWAN: merge all .rda visits -> one long table
## =========================

# ---- config ----
dir_in  <- "Rdata"   # CHANGE THIS FOLDER
save_rds <- "swan_merged.rds"
save_csv <- "swan_merged.csv"   # set to NULL if not needed

## SWAN: build long table
##  - keep SWANID and VISIT as-is from source df
##  - strip numeric suffixes for other fields (FSH0→FSH, TRIGRES2→TRIGRES)
##  - bind together only the “whitelist” of variables

KEEP_KEYS <- c(
          # id/visit/time
          "SWANID","VISIT","INTDAY","AGE",
          # hormones/metabolism
          "FSH","GLUCRES","INSURES","TRIGRES","LDLRESU","HDLRESU","BMI","SHBG",
          # symptoms
          "HOTFLAS","RESTLES","TRBLSLE","TYPNIGH",
          "MOODCHN","MOODCHG","DEPRESS",
          "LEAKURI","DAYSLEA","BREASTP","PELVIC",
          # behavioral
          "AVCIGDA","ALCHWK","ALCHMON",
          # menstrual characteristics (if needed)
          "FLOWDAY","TENDAFL","FLOWAMT","FLOODIN","CLOTS","STARTDA","USUALDA","FLOAGE2","LMPDAY"
)

## Which fields to strip numeric suffixes for
STRIP_SUFFIX_FOR <- setdiff(KEEP_KEYS, c("SWANID","VISIT"))

load_first_df_from_rda <- function(f) {
          env <- new.env(parent = emptyenv())
          objs <- load(f, envir = env)
          for (nm in objs) {
                    x <- get(nm, envir = env)
                    if (is.data.frame(x)) return(x)
          }
          stop(sprintf("No data.frame found in %s", basename(f)))
}

strip_suffix <- function(nms, keys = STRIP_SUFFIX_FOR) {
          sapply(nms, function(x) {
                    # only touch names that start with any of the keys
                    if (!any(startsWith(x, keys))) return(x)
                    # two-digit suffix at the end?
                    if (grepl("\\d{2}$", x)) {
                              substr(x, 1L, nchar(x) - 2L)
                    } else if (grepl("\\d$", x)) {
                              substr(x, 1L, nchar(x) - 1L)
                    } else {
                              x
                    }
          }, USE.NAMES = FALSE)
}


rda_files <- list.files(dir_in, pattern="\\.rda$", full.names = TRUE)
stopifnot(length(rda_files) > 0)

pieces <- vector("list", length(rda_files))

for (i in seq_along(rda_files)) {
          f <- rda_files[i]
          cat(sprintf(">>> %s\n", basename(f)))
          df0 <- load_first_df_from_rda(f)

          # SWANID and VISIT must exist
          if (!all(c("SWANID","VISIT") %in% names(df0))) {
                    stop(sprintf("No SWANID or VISIT in %s", basename(f)))
          }

          # keep only needed columns + possible suffixed variants
          # 1) take columns that match KEEP_KEYS OR their “stripped” versions
          nms0 <- names(df0)
          nms_stripped <- strip_suffix(nms0, STRIP_SUFFIX_FOR)

          # mask: name in KEEP_KEYS OR stripped(name) in KEEP_KEYS
          keep_mask <- (nms0 %in% KEEP_KEYS) | (nms_stripped %in% KEEP_KEYS)
          df <- df0[, keep_mask, drop = FALSE]

          # rename ONLY non-ID fields to stripped names;
          # keep SWANID and VISIT unchanged
          nms <- names(df)
          nms_new <- strip_suffix(nms, STRIP_SUFFIX_FOR)
          nms_new[nms %in% c("SWANID","VISIT")] <- nms[nms %in% c("SWANID","VISIT")]
          names(df) <- nms_new

          # if stripping created duplicate columns (e.g., FSH from different sources),
          # keep the one with fewer NAs
          if (any(duplicated(names(df)))) {
                    dups <- unique(names(df)[duplicated(names(df))])
                    for (nm in dups) {
                              cols <- which(names(df) == nm)
                              if (length(cols) > 1) {
                                        na_counts <- sapply(cols, function(j) sum(is.na(df[[j]])))
                                        keep_col <- cols[which.min(na_counts)]
                                        drop_cols <- setdiff(cols, keep_col)
                                        df <- df[, -drop_cols, drop = FALSE]
                              }
                    }
          }

          # fix types for ID/visit
          df$SWANID <- as.character(df$SWANID)
          df$VISIT  <- as.integer(as.character(df0$VISIT))

          # keep only whitelist (+ create missing NA-columns if needed)
          miss <- setdiff(KEEP_KEYS, names(df))
          if (length(miss)) df[miss] <- NA
          df <- df[, KEEP_KEYS, drop = FALSE]

          # source file
          df$SOURCE_FILE <- basename(f)

          pieces[[i]] <- df
}

# bind all
swan_core <- do.call(rbind, pieces)

# report
cat("\nDone.\nRows: ", nrow(swan_core), "  Cols: ", ncol(swan_core), "\n", sep = "")
print(table(swan_core$VISIT, useNA = "ifany"))

# save
saveRDS(swan_core, file = save_rds)
utils::write.csv(swan_core, file = save_csv, row.names = FALSE)




## --- small robust helpers ---
.as_num <- function(x) suppressWarnings(as.numeric(x))

.extract_code <- function(x) {
  # Input: vector (factor/char/num/NA). Output: numeric codes (first integer found in string).
  if (is.factor(x)) x <- as.character(x)
  if (is.numeric(x)) return(as.numeric(x))
  if (!length(x)) return(numeric(0))

  out <- rep(NA_real_, length(x))
  ok  <- !is.na(x)  # work only over non-NA
  if (!any(ok)) return(out)

  xx <- x[ok]
  m  <- regexpr("\\d+", xx, perl = TRUE)        # position of the first number
  has <- m > 0L                                 # where a number is found

  if (any(has)) {
    starts <- m[has]
    lens   <- attr(m, "match.length")[has]
    substrs <- substring(xx[has], starts, starts + lens - 1L)
    out[ok][has] <- suppressWarnings(as.numeric(substrs))
  }
  out
}


.make_ordered_factor <- function(raw, code) {
  if (all(is.na(code))) return(NULL)
  # unique (code, raw) ordered by code
  o <- order(code, na.last = NA)
  tab <- unique(data.frame(code = code[o], raw = raw[o], stringsAsFactors = FALSE))
  # clean labels: drop leading "(k)" and "k:" prefixes if present
  clean_lab <- function(s) {
    s <- sub("^\\s*\\(\\d+\\)\\s*", "", s)
    s <- sub("^\\s*\\d+\\s*:\\s*", "", s)
    trimws(s)
  }
  levs <- sort(unique(tab$code))
  labs <- vapply(levs, function(k) {
    i <- which(tab$code == k)[1]
    if (length(i) == 1L) clean_lab(tab$raw[i]) else as.character(k)
  }, character(1))
  factor(code, levels = levs, labels = labs, ordered = TRUE)
}

# estimate days for FLOWDAY (from textual range or code)
.parse_flowday_days <- function(raw, code) {
  if (is.factor(raw)) raw <- as.character(raw)
  out <- rep(NA_real_, length(raw))
  # 1) range like "1-5" or "9–13"
  rng <- regexec("(\\d+)\\s*[-–]\\s*(\\d+)", raw %||% "", perl = TRUE)
  mt  <- regmatches(raw, rng)
  has <- lengths(mt) > 0
  if (any(has)) {
    lo <- as.numeric(vapply(mt[has], function(z) z[2], ""))
    hi <- as.numeric(vapply(mt[has], function(z) z[3], ""))
    out[has] <- (lo + hi)/2
  }
  # 2) fall back to code
  # empirical midpoints for SWAN:
  # 1: 1–2 days → 1.5; 2: 3–5 → 4; 3: 6–8 → 7; 4: 9–13 → 11; 5: every day → 14
  mid <- c(`1`=1.5, `2`=4, `3`=7, `4`=11, `5`=14)
  need <- is.na(out) & !is.na(code)
  if (any(need)) out[need] <- unname(mid[as.character(code[need])])
  out
}
`%||%` <- function(x, y) if (is.null(x)) y else x


## ---------- Helpers (place near .extract_code and .make_ordered_factor) ----------
.map_yesno <- function(code_vec) {
  # expect codes 1/2 -> 0/1
  out <- ifelse(is.na(code_vec), NA_real_,
                ifelse(code_vec %in% c(1, 2), code_vec - 1, NA_real_))
  out
}

.make_yesno_factor <- function(code_vec) {
  labs <- c("No","Yes")
  out <- factor(ifelse(is.na(code_vec), NA, code_vec), levels = c(1,2), labels = labs, ordered = TRUE)
  out
}

# FLOWAMT: (1) Light, (2) Moderate, (3) Heavy
.make_flowamt <- function(raw, code) {
  lv <- c(1,2,3)
  labs <- c("Light","Moderate","Heavy")
  fct <- factor(ifelse(is.na(code), NA, code), levels = lv, labels = labs, ordered = TRUE)
  list(code = ifelse(code %in% lv, code, NA_real_), fct = fct)
}

# FLOWDAY: text categories -> estimated days (interval midpoint), else use numeric as-is
.flowday_to_days <- function(raw, code) {
  # SWAN mapping:
  # 1: 1–2 days -> ~1.5
  # 2: 3–5 days -> ~4
  # 3: 6–8 days -> ~7
  # 4: 9–13 days -> ~11
  # 5: Every day -> NA (better leave NA than force a value)
  mp <- c(`1`=1.5, `2`=4.0, `3`=7.0, `4`=11.0, `5`=NA_real_)
  out <- rep(NA_real_, length(code))

  # if there is a code from parentheses
  has_code <- !is.na(code)
  out[has_code] <- unname(mp[as.character(code[has_code])])

  # if no code but looks like a bare number (e.g. "24"), try to parse
  raw_num <- suppressWarnings(as.numeric(raw))
  use_raw_num <- is.na(out) & !is.na(raw_num)
  out[use_raw_num] <- raw_num[use_raw_num]

  out
}

## ---------- Updated recode_swan (minimal, without extra “magic”) ----------
recode_swan <- function(df,
                        numeric_vars = c("AGE","FSH","GLUCRES","INSURES","TRIGRES","LDLRESU","HDLRESU","BMI","SHBG"),
                        survey_vars  = c("HOTFLAS","RESTLES","TRBLSLE","TYPNIGH",
                                         "MOODCHN","MOODCHG","DEPRESS",
                                         "LEAKURI","DAYSLEA","BREASTP","PELVIC",
                                         "AVCIGDA","ALCHWK","ALCHMON",
                                         # new:
                                         "TENDAFL","FLOWAMT","FLOODIN","CLOTS","STARTDA","FLOWDAY"),
                        make_factors = TRUE) {

  stopifnot(is.data.frame(df))

  # 1) numeric markers
  for (v in intersect(numeric_vars, names(df))) {
    df[[paste0(v, "_raw")]] <- df[[v]]
    df[[v]] <- .as_num(df[[v]])
  }

  # 2) survey fields (unified code parser)
  for (v in intersect(survey_vars, names(df))) {
    raw  <- df[[v]]
    code <- .extract_code(raw)
    # if absolutely nothing was parsed — leave column unchanged
    if (all(is.na(code)) && all(is.na(raw))) next

    df[[paste0(v, "_raw")]] <- raw

    # special handling per field:
    if (v %in% c("TENDAFL","FLOODIN","CLOTS","STARTDA")) {
      # Yes/No -> 0/1 and factor
      df[[v]] <- .map_yesno(code)
      if (make_factors) df[[paste0(v, "_fct")]] <- .make_yesno_factor(code)

    } else if (v == "FLOWAMT") {
      am <- .make_flowamt(raw, code)
      df[[v]] <- am$code
      if (make_factors) df[[paste0(v, "_fct")]] <- am$fct

    } else if (v == "FLOWDAY") {
      # both code and “days estimate”
      df[[v]] <- code
      df[[paste0(v, "_days")]] <- .flowday_to_days(raw, code)
      if (make_factors) {
        # neat ordered factor exactly per SWAN coding
        lv <- c(1,2,3,4,5)
        labs <- c("1–2 days","3–5 days","6–8 days","9–13 days","Every day")
        df[[paste0(v, "_fct")]] <- factor(ifelse(is.na(code), NA, code), levels = lv, labels = labs, ordered = TRUE)
      }

    } else {
      # default: just codes + ordered factor over observed codes
      df[[v]] <- code
      if (make_factors) {
        lev <- sort(unique(na.omit(code)))
        if (length(lev)) {
          df[[paste0(v, "_fct")]] <- factor(ifelse(is.na(code), NA, code),
                                            levels = lev, ordered = TRUE)
        }
      }
    }
  }

  df
}


swan_core2<-recode_swan(swan_core)



drop_raw_cols <- function(df) {
  stopifnot(is.data.frame(df))
  df2<- df[, !grepl("_raw$", names(df)), drop = FALSE]
  df2[, !grepl("_fct$", names(df2)), drop = FALSE]
}

## ---- unit conversions helpers (SWAN -> SI, aligned with NHANES) ----
swan_to_si <- function(df) {
  # Expected raw SWAN names:
  # AGE (years), FSH (mIU/mL), GLUCRES (mg/dL), INSURES (uIU/mL),
  # TRIGRES (mg/dL), LDLRESU (mg/dL), HDLRESU (mg/dL),
  # BMI (kg/m^2), SHBG (nM == nmol/L)

  out <- df

  # Glucose: mg/dL -> mmol/L (÷ ~18.0)
  if ("GLUCRES" %in% names(out)) out$glucose <- as.numeric(out$GLUCRES) / 18.0

  # Insulin: uIU/mL -> pmol/L (× ~6.0)
  if ("INSURES" %in% names(out)) out$insulin <- as.numeric(out$INSURES) * 6.0

  # Triglycerides: mg/dL -> mmol/L (÷ 88.57)
  if ("TRIGRES" %in% names(out)) out$tg <- as.numeric(out$TRIGRES) / 88.57

  # Lipids (cholesterol): mg/dL -> mmol/L (÷ 38.67)
  if ("LDLRESU" %in% names(out)) out$ldl <- as.numeric(out$LDLRESU) / 38.67
  if ("HDLRESU" %in% names(out)) out$hdl <- as.numeric(out$HDLRESU) / 38.67

  # Direct carries (units already match)
  if ("AGE"  %in% names(out)) out$age  <- as.numeric(out$AGE)
  if ("FSH"  %in% names(out)) out$fsh  <- as.numeric(out$FSH)
  if ("BMI"  %in% names(out)) out$bmi  <- as.numeric(out$BMI)
  if ("SHBG" %in% names(out)) out$shbg <- as.numeric(out$SHBG)  # nM = nmol/L — OK

  # Keep only what we need
  keep <- c("SWANID","VISIT","age","fsh","glucose","insulin","tg","ldl","hdl","bmi","shbg")
  keep <- keep[keep %in% names(out)]
  out[, keep, drop = FALSE]
}

## ---- two-stage GAM Gamma(log) prediction on new data ----
predict_two_stage_gam_gamma <- function(model_early, model_late, newdata, thr = 37) {
  stopifnot(is.data.frame(newdata),
            "age" %in% names(newdata))
  # split by age
  idx_early <- which(newdata$age <  thr)
  idx_late  <- which(newdata$age >= thr)

  pred <- rep(NA_real_, nrow(newdata))

  if (length(idx_early)) {
    pred[idx_early] <- as.numeric(predict(model_early,
                                          newdata = newdata[idx_early, , drop = FALSE],
                                          type = "response"))   # Gamma(link=log) => direct ng/mL
  }
  if (length(idx_late)) {
    pred[idx_late] <- as.numeric(predict(model_late,
                                         newdata = newdata[idx_late, , drop = FALSE],
                                         type = "response"))
  }
  # guard against negatives (Gamma(log) should not, but keep just in case)
  pred[pred < 0] <- 0
  pred
}

## ---- apply isotonic calibrator if available (e.g., list with $apply) ----
iso_apply_safe <- function(iso, yhat) {
  if (is.null(iso)) return(yhat)
  if (is.function(iso)) return(iso(yhat))
  if (is.list(iso) && is.function(iso$apply)) return(iso$apply(yhat))
  # unknown iso format -> no calibration
  yhat
}

## ---- main function to apply AMH model to SWAN ----
apply_amh_model_swan <- function(swan_df, bundle,
                                 required_keys = c("age","fsh","glucose","insulin","tg","ldl","hdl","bmi","shbg"),
                                 age_threshold = 37) {
  stopifnot(is.data.frame(swan_df),
            is.list(bundle),
            !is.null(bundle$model$early),
            !is.null(bundle$model$late))

  # 1) convert SWAN to SI + align feature names
  feats <- swan_to_si(swan_df)

  # 2) keep only complete cases for predictors
  need <- intersect(required_keys, names(feats))
  if (!all(required_keys %in% names(feats))) {
    miss <- setdiff(required_keys, names(feats))
    warning("Missing predictors in SWAN: ", paste(miss, collapse = ", "),
            ". Predictions will be returned only for complete rows.")
  }
  keep <- complete.cases(feats[, need, drop = FALSE])
  feats_ok <- feats[keep, , drop = FALSE]
  if (!nrow(feats_ok)) {
    warning("No rows with a complete set of predictors — returning empty.")
    return(feats[0, c("SWANID","VISIT"), drop=FALSE])
  }

  # 3) predict with two-stage GAM Γ(log)
  yhat <- predict_two_stage_gam_gamma(
    model_early = bundle$model$early,
    model_late  = bundle$model$late,
    newdata     = feats_ok,
    thr         = age_threshold
  )

  # 4) isotonic calibration (if any)
  yhat_cal <- iso_apply_safe(bundle$iso, yhat)

  # 5) assemble result
  res <- data.frame(
    SWANID        = if ("SWANID" %in% names(feats_ok)) feats_ok$SWANID else NA,
    VISIT         = if ("VISIT"  %in% names(feats_ok)) feats_ok$VISIT  else NA,
    AMH_pred      = yhat,
    AMH_pred_cal  = yhat_cal,
    feats_ok[, intersect(c("age","fsh","glucose","insulin","tg","ldl","hdl","bmi","shbg"),
                         names(feats_ok)), drop = FALSE],
    row.names = NULL
  )

  # 6) merge back to the original table by index (only where complete cases existed)
  out <- merge(
    swan_df[, intersect(c("SWANID","VISIT"), names(swan_df)), drop = FALSE],
    res,
    by = c("SWANID","VISIT"),
    all.x = TRUE,
    sort = FALSE
  )

  out
}


# 1) drop auxiliary *_raw columns if needed
swan_clean <- drop_raw_cols(swan_core2)

# 2) load trained bundle (what you saved earlier)
#load("models/amh_two_stage_gam_gamma_bundle.rds")  # bundle$model$early, bundle$model$late, bundle$iso

load_bundle <- function(path) {
  stopifnot(file.exists(path))
  ext <- tolower(tools::file_ext(path))
  obj <- NULL

  if (ext == "rds") {
    # standard path for a single object
    obj <- readRDS(path)
  } else if (ext %in% c("rda","rdata")) {
    # load into isolated environment
    e <- new.env(parent = emptyenv())
    nm <- load(path, envir = e)
    if (length(nm) == 1) {
      obj <- e[[nm]]
    } else {
      # multiple objects — collect into a list
      obj <- mget(nm, envir = e, inherits = FALSE)
    }
  } else {
    stop("Unsupported extension: ", ext)
  }

  # verify it looks like the expected bundle.
  # If not list(model=..., iso=...), but early/late/iso exist — normalize to common format.
  if (is.list(obj) && all(c("model","iso") %in% names(obj))) {
    return(obj)
  }

  # “raw” variants:
  if (is.list(obj) && all(c("early","late") %in% names(obj))) {
    return(list(model = list(early = obj$early, late = obj$late),
                iso   = obj$iso %||% NULL,
                meta  = obj$meta %||% list()))
  }

  # list loaded from .rda:
  if (is.list(obj)) {
    # try to find components
    early <- obj$early %||% obj$model$early %||% NULL
    late  <- obj$late  %||% obj$model$late  %||% NULL
    iso   <- obj$iso   %||% NULL
    if (!is.null(early) && !is.null(late)) {
      return(list(model = list(early = early, late = late),
                  iso = iso, meta = obj$meta %||% list()))
    }
  }

  stop("Does not look like a two-stage GAM bundle. Contents: ", paste(names(obj), collapse=", "))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# Example usage:
bundle <- load_bundle("models/amh_two_stage_gam_gamma_bundle.rds")


library(mgcv)




# 3) predict AMH on SWAN
pred_swan <- apply_amh_model_swan(
  swan_df = swan_clean,
  bundle  = bundle,             # your saved bundle
  required_keys = c("age","fsh","glucose","insulin","tg","ldl","hdl","bmi","shbg"),
  age_threshold = 37
)

pred_swan_sample<-pred_swan[sample(1:nrow(pred_swan),1000),]

plot(pred_swan_sample$age,pred_swan_sample$AMH_pred)

plot(pred_swan_sample$age,pred_swan_sample$fsh)

id<-10464
pred_swan_tst<-pred_swan[pred_swan$SWANID==id,]



## ========= Helpers =========
nz_num <- function(x) {    # convert to numeric; invalid/negative to NA
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x)] <- NA_real_
  x
}
coalesce_num <- function(...) {
  args <- list(...)
  for (a in args) {
    if (!all(is.na(a))) return(a)
  }
  rep(NA_real_, length(args[[1]]))
}

## ========= Staging function =========
## Input: swan_clean with already recoded numeric columns:
##   AGE, FSH, LMPDAY  (see your recode_swan)
## Output: two new columns:
##   is_postmeno (TRUE/FALSE/NA)
##   meno_stage ∈ {pre, early_peri, late_peri, post, unknown}
stage_swan <- function(df,
                       lmp_thresh_post  = 365,  # days without menses for post
                       lmp_thresh_late  = 90,   # days for late perimenopause
                       age_fallback_post = 55,  # if no LMPDAY, age fallback
                       fsh_early = 25,          # early peri support
                       fsh_post  = 40) {        # high FSH → late/post

  stopifnot(is.data.frame(df))
  if (!all(c("AGE","FSH","LMPDAY") %in% names(df))) {
    stop("Expect columns AGE, FSH, LMPDAY (numeric) in df.")
  }

  age    <- nz_num(df$AGE)
  fsh    <- nz_num(df$FSH)
  lmp    <- nz_num(df$LMPDAY)           # days since last menses
  lmp[lmp < 0] <- NA_real_              # drop negatives

  is_post <- rep(NA, nrow(df))

  ## 1) hard rule by LMPDAY
  is_post[!is.na(lmp) & lmp >= lmp_thresh_post] <- TRUE
  is_post[!is.na(lmp) & lmp <  lmp_thresh_post] <- FALSE

  ## 2) fallback: if LMPDAY is unknown
  idx_na_lmp <- which(is.na(lmp))
  if (length(idx_na_lmp)) {
    # age fallback
    is_post[idx_na_lmp[age[idx_na_lmp] >= age_fallback_post]] <- TRUE
    # hormonal fallback
    is_post[idx_na_lmp[is.na(is_post[idx_na_lmp]) & !is.na(fsh[idx_na_lmp]) & fsh[idx_na_lmp] >= fsh_post]] <- TRUE
    is_post[idx_na_lmp[is.na(is_post[idx_na_lmp]) & !is.na(fsh[idx_na_lmp]) & fsh[idx_na_lmp] <  fsh_post]] <- FALSE
  }

  ## 3) stage label
  stage <- rep(NA_character_, nrow(df))

  # post via LMPDAY/fallbacks
  stage[is_post %in% TRUE]  <- "post"

  # late peri: 90–364 days without bleeding
  stage[is.na(stage) & !is.na(lmp) & lmp >= lmp_thresh_late & lmp < lmp_thresh_post] <- "late_peri"

  # early peri: menses (<90 days) + elevated FSH (or just use FSH as early transition marker)
  stage[is.na(stage) & !is.na(lmp) & lmp < lmp_thresh_late & !is.na(fsh) & fsh >= fsh_early] <- "early_peri"

  # pre: recent menses (<90 days) and FSH not elevated
  stage[is.na(stage) & !is.na(lmp) & lmp < lmp_thresh_late & !is.na(fsh) & fsh < fsh_early] <- "pre"

  # if no LMPDAY: FSH hints stage
  stage[is.na(stage) & is.na(lmp) & !is.na(fsh) & fsh >= fsh_post]  <- "late_peri"
  stage[is.na(stage) & is.na(lmp) & !is.na(fsh) & fsh >= fsh_early & fsh < fsh_post] <- "early_peri"
  stage[is.na(stage) & is.na(lmp) & !is.na(fsh) & fsh < fsh_early] <- "pre"

  # unclear -> unknown
  stage[is.na(stage)] <- "unknown"

  df$is_postmeno <- is_post
  df$meno_stage  <- factor(stage, levels = c("pre","early_peri","late_peri","post","unknown"))
  df
}



swan_flagged <- stage_swan(swan_clean)
table(swan_flagged$meno_stage, useNA = "ifany")


## Forward-fixing postmenopause:
## if a participant first shows ≥ k_consec consecutive visits with status "post",
## then ALL subsequent (including that run) visits are labeled "post" and is_postmeno=TRUE.
## Previous visits are NOT modified.
enforce_post_forward <- function(df,
                                 id_col      = "SWANID",
                                 visit_col   = "VISIT",
                                 stage_col   = "meno_stage",
                                 is_post_col = "is_postmeno",
                                 k_consec    = 2) {
  stopifnot(is.data.frame(df),
            id_col %in% names(df),
            visit_col %in% names(df),
            stage_col %in% names(df),
            is_post_col %in% names(df))

  # preserve factor levels if present
  stage_levels <- NULL
  if (is.factor(df[[stage_col]])) stage_levels <- levels(df[[stage_col]])

  # helper: return START index of the first run with ≥ k consecutive TRUE; NA if none
  first_run_start <- function(x, k) {
    if (!length(x)) return(NA_integer_)
    r <- rle(x %in% TRUE)
    ends   <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1
    i <- which(r$values & r$lengths >= k)[1]
    if (is.na(i)) return(NA_integer_)
    starts[i]
  }

  # mark that from this row (in visit order) everything is post
  to_force <- logical(nrow(df))

  split_idx <- split(seq_len(nrow(df)), df[[id_col]])
  for (rows in split_idx) {
    ord <- rows[order(df[[visit_col]][rows], na.last = TRUE)]
    st  <- as.character(df[[stage_col]][ord])
    is_post_here <- !is.na(st) & st == "post"

    start_run <- first_run_start(is_post_here, k_consec)
    if (!is.na(start_run)) {
      # from start_run forward — enforce post
      to_force[ord[start_run:length(ord)]] <- TRUE
    }
  }

  # apply changes
  if (!is.null(stage_levels)) {
    if (!("post" %in% stage_levels)) stage_levels <- c(stage_levels, "post")
    df[[stage_col]] <- as.character(df[[stage_col]])
    df[[stage_col]][to_force] <- "post"
    df[[stage_col]] <- factor(df[[stage_col]], levels = stage_levels)
  } else {
    df[[stage_col]][to_force] <- "post"
  }

  # is_postmeno flag: forward only
  if (!is.logical(df[[is_post_col]])) {
    df[[is_post_col]] <- as.logical(df[[is_post_col]])
  }
  df[[is_post_col]][to_force] <- TRUE
  df[[is_post_col]][is.na(df[[is_post_col]])] <- FALSE

  # optional: service marker for rows modified by forward pass
  df$._post_forced_forward <- FALSE
  df$._post_forced_forward[to_force] <- TRUE

  df
}





id<-10245
swan_flagged_tst<-swan_flagged[swan_flagged$SWANID==id,]


patients_stat<-as.data.frame(table(swan_flagged$SWANID))


swan_flagged_corrected_status<-enforce_post_forward(swan_flagged)


## join datasets

## 1) build keys
prep_keys <- function(df, id_col="SWANID", visit_col="VISIT") {
  df[[id_col]]    <- as.character(df[[id_col]])
  df[[visit_col]] <- as.integer(df[[visit_col]])
  df$key <- paste(df[[id_col]], df[[visit_col]], sep = "-")
  df
}

pred_swan_k <- prep_keys(pred_swan)
flags_k     <- prep_keys(swan_flagged_corrected_status)

## 2) Left join by key (simple)
# preserve pred_swan order:
idx <- match(pred_swan_k$key, flags_k$key)
merged <- cbind(
  pred_swan_k,
  flags_k[ idx, setdiff(names(flags_k), c("SWANID","VISIT","key")), drop = FALSE ]
)

names(merged)

## 3) (optional) drop auxiliary key
# merged$key <- NULL

## Quick summary
cat("Rows in pred_swan:", nrow(pred_swan_k), "\n")
cat("Matched:", sum(!is.na(idx)), "\n")
cat("Unmatched:", sum(is.na(idx)), "\n")



## Determine age at menopause by visits
## df         — merged table (after join), one row = visit
## id_col     — participant id column
## visit_col  — visit number (integer, monotonic over time)
## age_col    — age at visit (years, numeric)
## post_col   — boolean TRUE/FALSE whether visit is post;
##              if not available, you may use status_col == "post"
## status_col — (opt.) categorical status ("pre"/"peri"/"post"); used if post_col is absent
menopause_onset_by_person <- function(
    df,
    id_col    = "SWANID",
    visit_col = "VISIT",
    age_col   = "AGE",
    post_col  = "is_postmeno",     # set your boolean flag column here
    status_col = NULL              # e.g., "status" or "menopausal_status"
) {
  stopifnot(all(c(id_col, visit_col, age_col) %in% names(df)))

  # build logical post_flag
  if (!is.null(post_col) && post_col %in% names(df)) {
    post_flag <- as.logical(df[[post_col]])
  } else if (!is.null(status_col) && status_col %in% names(df)) {
    post_flag <- !is.na(df[[status_col]]) & tolower(as.character(df[[status_col]])) == "post"
  } else {
    stop("Need either post_col (boolean) or status_col ('post' as a value).")
  }

  # local copies
  id    <- df[[id_col]]
  visit <- df[[visit_col]]
  age   <- suppressWarnings(as.numeric(df[[age_col]]))

  tmp <- data.frame(
    id    = as.character(id),
    visit = as.integer(visit),
    age   = age,
    post  = post_flag,
    stringsAsFactors = FALSE
  )

  # drop rows missing id/visit/age
  tmp <- tmp[!is.na(tmp$id) & !is.na(tmp$visit) & !is.na(tmp$age), , drop = FALSE]

  # sort by id, then visit
  o <- order(tmp$id, tmp$visit, tmp$age)
  tmp <- tmp[o, ]

  # split by participant
  sp <- split(tmp, tmp$id)

  res_list <- lapply(sp, function(dd) {
    # dd: one participant
    # first visit with post = TRUE
    idx_post <- which(dd$post %in% TRUE)
    ever_post <- length(idx_post) > 0

    if (!ever_post) {
      return(data.frame(
        SWANID            = dd$id[1],
        ever_post         = FALSE,
        visit_first_post  = NA_integer_,
        age_first_post    = NA_real_,
        age_last_nonpost  = if (any(!dd$post)) max(dd$age[!dd$post], na.rm = TRUE) else NA_real_,
        age_mid_transition= NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    i_first <- idx_post[1]
    age_first <- dd$age[i_first]
    vis_first <- dd$visit[i_first]

    # age at last “non-post” before first “post”
    if (i_first > 1) {
      prev_nonpost_idx <- max(which(!dd$post[1:(i_first-1)]), na.rm = TRUE)
      age_prev <- if (is.finite(prev_nonpost_idx)) dd$age[prev_nonpost_idx] else NA_real_
    } else {
      age_prev <- NA_real_
    }

    # midpoint of transition (if pre-post ages available)
    age_mid <- if (!is.na(age_prev)) mean(c(age_prev, age_first)) else age_first

    data.frame(
      SWANID            = dd$id[1],
      ever_post         = TRUE,
      visit_first_post  = as.integer(vis_first),
      age_first_post    = as.numeric(age_first),
      age_last_nonpost  = as.numeric(age_prev),
      age_mid_transition= as.numeric(age_mid),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, res_list)
  rownames(out) <- NULL
  out
}

# merged — your data frame after join (contains SWANID, VISIT, AGE, is_postmeno OR status)
onset_df <- menopause_onset_by_person(
  df         = merged,
  id_col     = "SWANID",
  visit_col  = "VISIT",
  age_col    = "AGE",
  post_col   = "is_postmeno",         # set your boolean flag column here
  status_col = "status"               # or text status if no boolean flag
)

plot(density(onset_df$age_mid_transition,na.rm = T))

# quick check
summary(onset_df$age_first_post)
head(onset_df, 10)

# somehow ~48 — seems low

## Normalize LMPDAY: extract numbers from strings; negative and unrealistic -> NA
extract_num <- function(x) {
  if (is.numeric(x)) return(x)
  m <- regexpr("[-]?[0-9]+", as.character(x))
  out <- rep(NA_real_, length(x))
  has <- m > 0
  out[has] <- suppressWarnings(as.numeric(substr(x[has], m[has], m[has] + attr(m, "match.length")[has] - 1L)))
  out
}

clean_lmp <- function(x) {
  v <- extract_num(x)
  v[!is.finite(v)] <- NA_real_
  v[v < 0] <- NA_real_           # drop negatives
  v[v > 365*20] <- NA_real_      # trim obvious artifacts
  v
}

## Compute FMP age by LMPDAY (12-month rule) and monotonicity check
## df: by visits (rows) with columns SWANID, VISIT, AGE, LMPDAY (as is)
## Returns one row per woman
compute_fmp_age <- function(df, id_col="SWANID", visit_col="VISIT",
                            age_col="AGE", lmp_col="LMPDAY") {
  stopifnot(all(c(id_col, visit_col, age_col, lmp_col) %in% names(df)))
  tmp <- df[, c(id_col, visit_col, age_col, lmp_col)]
  names(tmp) <- c("id","visit","age","lmp")
  tmp$age <- suppressWarnings(as.numeric(tmp$age))
  tmp$lmp <- clean_lmp(tmp$lmp)

  tmp <- tmp[is.finite(tmp$age) & !is.na(tmp$visit), ]
  tmp <- tmp[order(tmp$id, tmp$visit), ]
  sp  <- split(tmp, tmp$id)

  res <- lapply(sp, function(dd){
    # candidates: visits with LMPDAY >= 365
    cand <- which(is.finite(dd$lmp) & dd$lmp >= 365)
    if (length(cand) == 0) {
      return(data.frame(SWANID = dd$id[1],
                        ever_fmp = FALSE,
                        age_fmp  = NA_real_,
                        visit_fmp= NA_integer_,
                        stringsAsFactors = FALSE))
    }

    # “reset” check after a candidate: if later lmp < 365, candidate is doubtful
    good <- NA
    for (i in cand) {
      after <- dd$lmp[(i+1):nrow(dd)]
      if (length(after) == 0 || all(is.na(after) | after >= 365)) {
        good <- i; break
      }
    }
    if (is.na(good)) {
      # no stable candidates — take the first one (you may mark as doubtful if needed)
      good <- cand[1]
    }

    # FMP age = age at visit minus (days without menses)/365.25
    age_fmp <- dd$age[good] - dd$lmp[good]/365.25

    data.frame(SWANID   = dd$id[1],
               ever_fmp = TRUE,
               age_fmp  = as.numeric(age_fmp),
               visit_fmp= as.integer(dd$visit[good]),
               stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, res); rownames(out) <- NULL
  out
}

## Example:
## merged_df — your merged SWAN table (after join), has SWANID, VISIT, AGE, LMPDAY
## Get FMP age estimate:
fmp_df <- compute_fmp_age(merged, id_col="SWANID", visit_col="VISIT",
                          age_col="AGE", lmp_col="LMPDAY")


# again ~46 — odd

## Distribution diagnostics:
hist(fmp_df$age_fmp, breaks=40, main="FMP age (LMPDAY ≥ 365)", xlab="years")
summary(fmp_df$age_fmp, na.rm=TRUE)


all_ids<-patients_stat[patients_stat$Freq>8,]$Var1





# ---- 1) column-chooser helpers ----
.choose_first_present <- function(df, candidates) {
  for (nm in candidates) if (nm %in% names(df)) return(nm)
  stop("None of the columns found: ", paste(candidates, collapse=", "))
}

# ---- 2) main function ----
estimate_fmp_by_amh_threshold <- function(
    df,
    id_col              = "SWANID",
    age_cols_candidates = c("AGE","AGE0","age"),
    visit_cols_cand     = c("VISIT","visit"),
    day_cols_cand       = c("INTDAY","INTDAY0","intday"),
    amh_cols_cand       = c("AMH_pred_cal","AMH_pred","AMH","amh"),
    threshold           = 0.15,
    return_if_all_below = c("NA","first_visit_age")  # behavior if never > threshold
) {
  return_if_all_below <- match.arg(return_if_all_below)

  stopifnot(id_col %in% names(df))
  age_col   <- .choose_first_present(df, age_cols_candidates)
  amh_col   <- .choose_first_present(df, amh_cols_cand)
  visit_col <- if (any(visit_cols_cand %in% names(df))) {
    .choose_first_present(df, visit_cols_cand)
  } else NULL
  day_col   <- if (any(day_cols_cand %in% names(df))) {
    .choose_first_present(df, day_cols_cand)
  } else NULL

  # sort keys
  order_keys <- c(age_col, visit_col, day_col)
  order_keys <- order_keys[!is.null(order_keys) & nzchar(order_keys)]

  ids <- unique(df[[id_col]])
  out_list <- vector("list", length(ids))

  for (i in seq_along(ids)) {
    sid <- ids[i]
    sub <- df[df[[id_col]] == sid, , drop = FALSE]

    # keep rows with both age and AMH present
    keep <- !is.na(sub[[age_col]]) & !is.na(sub[[amh_col]])
    sub  <- sub[keep, , drop = FALSE]
    if (nrow(sub) == 0) {
      out_list[[i]] <- data.frame(
        SWANID = sid,
        est_age_fmp_amh_thr = NA_real_,
        last_above_visit_idx = NA_integer_,
        n_visits_used = 0L,
        stringsAsFactors = FALSE
      )
      next
    }

    # sort: age -> visit -> intday
    ord <- do.call(order, c(unname(sub[order_keys]), list(na.last = TRUE)))
    sub <- sub[ord, , drop = FALSE]

    # indices where AMH > threshold
    idx_above <- which(sub[[amh_col]] > threshold)

    if (length(idx_above) == 0) {
      # never above threshold
      est_age <- if (return_if_all_below == "first_visit_age") {
        as.numeric(sub[[age_col]][1])
      } else {
        NA_real_
      }
      last_above <- NA_integer_
    } else {
      last_above <- max(idx_above)
      # if there is a next visit — take its age; else NA (right-censored)
      if (last_above < nrow(sub)) {
        est_age <- as.numeric(sub[[age_col]][last_above + 1L])
      } else {
        est_age <- NA_real_
      }
    }

    out_list[[i]] <- data.frame(
      SWANID = sid,
      est_age_fmp_amh_thr = est_age,
      last_above_visit_idx = last_above,
      n_visits_used = nrow(sub),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out_list)
}

# ---- 3) example: apply to merged and join back ----
result_fmp <- estimate_fmp_by_amh_threshold(merged, threshold = 0.15) #0.15

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#43.00   47.00   48.00   48.42   50.00   58.00   14576

merged_with_fmp <- merge(merged, result_fmp, by = "SWANID", all.x = TRUE)

table(merged_with_fmp$VISIT)

merged_with_fmp_visit0<-merged_with_fmp[merged_with_fmp$VISIT==0,]
merged_with_fmp_visit10<-merged_with_fmp[merged_with_fmp$VISIT==10,]

summary(merged_with_fmp$est_age_fmp_amh_thr)
summary(merged_with_fmp_visit0$est_age_fmp_amh_thr)
summary(merged_with_fmp_visit10$est_age_fmp_amh_thr)

plot(density(merged_with_fmp_visit10$est_age_fmp_amh_thr,na.rm = T))
abline(v=52,lty=2)

id<-sample(all_ids,1)

# 70179 43861

merged_tst<-merged[merged$SWANID==id,]
merged_tst<-merged_tst[order(merged_tst$VISIT),]


age_fmp<-merged_with_fmp[merged_with_fmp$SWANID==id,]$est_age_fmp_amh_thr[1]

colori<-ifelse(merged_tst$AGE==age_fmp,"red","black")


plot(merged_tst$AMH_pred,merged_tst$fsh,type = "n",xlim = c(0,1),xlab = "AHM, predicted",ylab = "TSH, SWAN data",main = paste0("SWANID = ",id))
text(merged_tst$AMH_pred,merged_tst$fsh,labels = as.character(merged_tst$AGE),col=colori)
abline(v=0.15,lty=2)


summary(merged_with_fmp$est_age_fmp_amh_thr)

save(merged_with_fmp,file = "merged_with_fmp.rds")
