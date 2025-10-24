# symptom prediction using SWAN

load("merged_with_fmp.bin")


# === Simple SWAN symptom predictor ===
# Only requires the pROC package for AUC calculation
if(!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
library(pROC)

# --- Symptoms and features ---
symptoms <- c(
  "HOTFLAS","RESTLES","TRBLSLE","TYPNIGH",
  "MOODCHN","MOODCHG","DEPRESS",
  "LEAKURI","BREASTP","PELVIC","TENDAFL"
)

features <- c("age","fsh","glucose","insulin","tg","ldl","hdl","bmi","shbg")#,"est_age_fmp_amh_thr")

# --- Engineered features ---
merged_with_fmp$age_minus_thr <- merged_with_fmp$age - merged_with_fmp$est_age_fmp_amh_thr
merged_with_fmp$is_post_thr   <- as.integer(merged_with_fmp$age >= merged_with_fmp$est_age_fmp_amh_thr)
#features <- c(features, "age_minus_thr", "is_post_thr")

# --- Symptom binarization ---
binarize <- function(x) {
  # if a 1–5 numeric scale, 1 = no, >=2 = yes
  if (is.numeric(x)) {
    return(ifelse(is.na(x), NA, ifelse(x >= 2, 1L, 0L)))
  } else {
    x <- tolower(trimws(as.character(x)))
    yes  <- c("yes","y","true","present","often","sometimes","moderate","severe","mild")
    no   <- c("no","n","false","absent","never","none","not at all")
    return(ifelse(x %in% yes, 1L,
                  ifelse(x %in% no, 0L, NA_integer_)))
  }
}

# --- Main training loop ---
results <- data.frame()
models <- list()

for (sym in symptoms) {
  if (!sym %in% names(merged_with_fmp)) next
  cat("=== Symptom:", sym, "===\n")

  df <- merged_with_fmp[, c(sym, features)]
  df[[sym]] <- binarize(df[[sym]])
  df <- df[complete.cases(df), ]

  # skip if there is too little data
  if (nrow(df) < 50 || length(unique(df[[sym]])) < 2) {
    cat("Skip:", sym, "- not enough data or single class\n")
    next
  }

  form <- as.formula(paste(sym, "~", paste(features, collapse = " + ")))
  m <- glm(form, data = df, family = binomial())
  models[[sym]] <- m

  pr <- predict(m, type = "response")
  auc <- tryCatch({
    as.numeric(pROC::roc(df[[sym]], pr, quiet = TRUE)$auc)
  }, error = function(e) NA)

  results <- rbind(
    results,
    data.frame(symptom = sym, n = nrow(df),
               pos_rate = mean(df[[sym]]),
               AUC = auc, stringsAsFactors = FALSE)
  )
}

cat("\n--- Summary ---\n")
print(results)

# --- Prediction for a specific person ---
predict_symptoms_for_person <- function(person, models) {
  base <- as.data.frame(person)
  thr <- base$est_age_fmp_amh_thr
  base$age_minus_thr <- base$age - thr
  base$is_post_thr   <- as.integer(base$age >= thr)

  pre  <- base; pre$age <- thr - 0.1
  pre$age_minus_thr <- pre$age - thr
  pre$is_post_thr   <- 0L

  post <- base; post$age <- thr + 0.1
  post$age_minus_thr <- post$age - thr
  post$is_post_thr   <- 1L

  out <- data.frame(symptom=character(), p_pre_thr=numeric(), p_post_thr=numeric(), delta=numeric())

  for (sym in names(models)) {
    m <- models[[sym]]
    p_pre  <- predict(m, newdata = pre,  type = "response")
    p_post <- predict(m, newdata = post, type = "response")
    out <- rbind(out, data.frame(symptom=sym,
                                 p_pre_thr=p_pre,
                                 p_post_thr=p_post,
                                 delta=p_post - p_pre))
  }
  out
}

# --- Example prediction ---
# person <- list(
#   age=49.2, fsh=34.0, glucose=92.0, insulin=8.1,
#   tg=120.0, ldl=130.0, hdl=55.0, bmi=27.4,
#   shbg=45.0, est_age_fmp_amh_thr=50.1
# )
# predict_symptoms_for_person(person, models)



summary(models[["HOTFLAS"]])
summary(models[["PELVIC"]])


summary(models[["HOTFLAS"]])


#
# # === Volcano plots by predictors ===
# if(!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
# library(ggplot2)
#
# make_volcano <- function(model, title) {
#   s <- summary(model)$coefficients
#   df <- data.frame(
#     term = rownames(s),
#     estimate = s[, "Estimate"],
#     pvalue = s[, "Pr(>|z|)"]
#   )
#   df <- df[df$term != "(Intercept)", ]
#   df$logp <- -log10(df$pvalue)
#
#   ggplot(df, aes(x = estimate, y = logp, label = term)) +
#     geom_point(color = "firebrick", size = 3, alpha = 0.7) +
#     geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
#     geom_text(vjust = -0.7, size = 6, color = "black") +
#     labs(
#       title = paste("Volcano plot:", title),
#       x = "Model coefficient (logit)",
#       y = "-log10(p-value)"
#     ) +
#     theme_minimal(base_size = 14) +
#     theme(
#       plot.title = element_text(face = "bold", hjust = 0.5),
#       panel.grid.minor = element_blank()
#     )
# }
#
# # --- Generate for all models ---
# for (sym in names(models)) {
#   cat("Drawing volcano for:", sym, "\n")
#   print(make_volcano(models[[sym]], sym))
# }

# === Volcano plots for each symptom model (white background) ===
if(!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

# --- readable English names for symptoms ---
symptom_labels <- c(
  HOTFLAS = "Hot flashes / Night sweats",
  RESTLES = "Restless sleep",
  TRBLSLE = "Trouble sleeping",
  TYPNIGH = "Night awakenings",
  MOODCHN = "Mood changes (general)",
  MOODCHG = "Mood swings",
  DEPRESS = "Depressive symptoms",
  LEAKURI = "Urinary leakage",
  BREASTP = "Breast pain / tenderness",
  PELVIC  = "Pelvic pain / discomfort",
  TENDAFL = "Tendon / muscle ache"
)

# --- create a folder for images ---
plot_dir <- "plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir)

# --- function to build a volcano plot ---
make_volcano <- function(model, title) {
  s <- summary(model)$coefficients
  df <- data.frame(
    term = rownames(s),
    estimate = s[, "Estimate"],
    pvalue = s[, "Pr(>|z|)"]
  )
  df <- df[df$term != "(Intercept)", ]
  df$logp <- -log10(df$pvalue)

  ggplot(df, aes(x = estimate, y = logp, label = term)) +
    geom_point(color = "firebrick", size = 3, alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_text(vjust = -0.8, size = 3, color = "black") +
    labs(
      title = paste("Volcano plot:", title),
      x = "Coefficient (log-odds)",
      y = "-log10(p-value)"
    ) +
    theme_bw(base_size = 14) +                # <-- white background
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

# --- loop over models ---
for (sym in names(models)) {
  if (is.null(models[[sym]])) next

  # Title for the plot and file
  title <- ifelse(sym %in% names(symptom_labels),
                  symptom_labels[[sym]],
                  sym)
  file_name <- paste0(plot_dir, "/", sym, "_volcano.png")

  cat("Saving volcano plot for:", sym, "→", file_name, "\n")

  # Build the plot
  p <- make_volcano(models[[sym]], title)

  # Save PNG
  ggsave(filename = file_name, plot = p, width = 7, height = 5, dpi = 300, bg = "white")
}

cat("\n✅ Volcano plots saved to folder:", plot_dir, "with white background.\n")


# if (!dir.exists("models")) dir.create("models")
# for (sym in names(models)) {
#   saveRDS(models[[sym]], file = file.path("models", paste0("model_", sym, ".rds")))
# }
