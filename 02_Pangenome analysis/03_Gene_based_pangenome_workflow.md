
# Gene based pangenome using Orthofinder 

## Running orthofinder with all longest isoforms proteins sequneces

```bash

#!/bin/bash

#SBATCH -D /path/to/directory/orthofinder_pb
#SBATCH -J pb_orth
#SBATCH -o Pb-%j.out
#SBATCH -c 10
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=150G

# Increase open files limit for this job
ulimit -n 65536

#Load the software with module if applicable:

module load orthofinder/2.5.5 diamond/2.1.11 ncbiblast/2.16.0 mcl/14-137 fastme/2.1.5 python/2.7 mafft/7.525 

# command 

orthofinder -f /path/to/directory/orthofinder_pb/pb_final_longest_isoform_proteins -t 10

```

## Pan-core genome accumulation and new gene family discovery rate with Heaps Law classification

```R

# ─────────────────────────────────────────────────────────────────────────────
# Two separate figures:
#   Figure 1: Pan/Core accumulation curves (no inline labels, legend instead)
#   Figure 2: New gene family discovery curve (standalone, full panel border)
# ─────────────────────────────────────────────────────────────────────────────

library(ggplot2)
library(scales)


# ── 0. Parameters ─────────────────────────────────────────────────────────────
INPUT_FILE <- "/path/to/directory/pangene_PAV_binary.tsv"   # <-- change to your file path
N_ITER     <- 1000
SEED          <- 42
OUT_MAIN_PDF  <- "./path/to/directory/pancore_curves.pdf"
OUT_MAIN_PNG  <- "./path/to/directory/pancore_curves.png"
OUT_INSET_PDF <- "./path/to/directory/new_gene_discovery.pdf"
OUT_INSET_PNG <- "./path/to/directory/new_gene_discovery.png"


# ── 1. Load matrix ─────────────────────────────────────────────────────────────
cat("Loading matrix...\n")
mat_raw  <- read.delim(INPUT_FILE, row.names = 1, check.names = FALSE)
mat      <- as.matrix((mat_raw > 0) * 1L)
outliers <- c("cl02_proteins", "on04_proteins", "nl03_proteins")
mat      <- mat[, !colnames(mat) %in% outliers]
genomes   <- colnames(mat)
n_genomes <- length(genomes)
cat(sprintf("  %d orthogroups x %d genomes\n", nrow(mat), n_genomes))

# ── 2. Simulate ────────────────────────────────────────────────────────────────
cat(sprintf("Running %d simulations per n...\n", N_ITER))
set.seed(SEED)

pan_mat  <- matrix(NA_integer_, nrow = n_genomes, ncol = N_ITER)
core_mat <- matrix(NA_integer_, nrow = n_genomes, ncol = N_ITER)

for (n in seq_len(n_genomes)) {
  for (iter in seq_len(N_ITER)) {
    sampled <- sample(genomes, size = n, replace = FALSE)
    rs      <- rowSums(mat[, sampled, drop = FALSE])
    pan_mat[n,  iter] <- sum(rs >= 1)
    core_mat[n, iter] <- sum(rs == n)
  }
  if (n %% 10 == 0) cat(sprintf("  n = %d / %d\n", n, n_genomes))
}

# ── 3. Summarise ──────────────────────────────────────────────────────────────
calc_stats <- function(m) {
  data.frame(
    med = apply(m, 1, median),
    q1  = apply(m, 1, quantile, 0.25),
    q3  = apply(m, 1, quantile, 0.75)
  )
}

pan_s  <- calc_stats(pan_mat)
core_s <- calc_stats(core_mat)

main_df <- data.frame(
  n        = seq_len(n_genomes),
  pan_med  = pan_s$med,  pan_q1  = pan_s$q1,  pan_q3  = pan_s$q3,
  core_med = core_s$med, core_q1 = core_s$q1, core_q3 = core_s$q3
)

# LOESS smooth
for (col in names(main_df)[-1]) {
  fit <- loess(main_df[[col]] ~ main_df$n, span = 0.25)
  main_df[[paste0(col,"_s")]] <- predict(fit, main_df$n)
}

# ── 4. New gene discovery ─────────────────────────────────────────────────────
pan_sd_v     <- apply(pan_mat, 1, sd)
new_genes    <- diff(pan_s$med)
new_genes_sd <- sqrt(pan_sd_v[-1]^2 + pan_sd_v[-n_genomes]^2)
n_new        <- seq(2, n_genomes)

fit_mask <- new_genes > 0
n_fit    <- n_new[fit_mask]
ng_fit   <- new_genes[fit_mask]

new_df <- data.frame(
  n       = n_new,
  new_med = pmax(new_genes, 0),
  new_lo  = pmax(new_genes - new_genes_sd, 0),
  new_hi  = new_genes + new_genes_sd
)

# ── 5. Fit new gene model ─────────────────────────────────────────────────────
log_fit <- lm(log(ng_fit) ~ log(n_fit))
a_start <- -coef(log_fit)[2]
k_start <- exp(coef(log_fit)[1])
k_fit <- k_start; a_fit <- a_start; k_se <- NA; a_se <- NA; r2 <- NA

tryCatch({
  nls_new <- nls(ng_fit ~ k * n_fit^(-a),
                 start   = list(k = k_start, a = a_start),
                 control = nls.control(maxiter = 500))
  k_fit <- coef(nls_new)["k"]
  a_fit <- coef(nls_new)["a"]
  k_se  <- summary(nls_new)$coefficients["k", "Std. Error"]
  a_se  <- summary(nls_new)$coefficients["a", "Std. Error"]
  r2    <- 1 - sum(residuals(nls_new)^2) /
    sum((ng_fit - mean(ng_fit))^2)
  cat(sprintf("  k=%.1f+/-%.1f, a=%.3f+/-%.3f, R2=%.3f\n",
              k_fit, k_se, a_fit, a_se, r2))
}, error = function(e) cat("NLS failed:", conditionMessage(e), "\n"))

inset_max   <- 10
inset_df    <- new_df[new_df$n <= inset_max, ]
n_curve     <- seq(2, inset_max, by = 0.2)
inset_curve <- data.frame(n = n_curve, new = k_fit * n_curve^(-a_fit))

# ── 6. Colours ────────────────────────────────────────────────────────────────
col_pan  <- "#88694c"
col_core <- "#5f9fc6"

# ── 7. Y-axis limits ──────────────────────────────────────────────────────────
y_min  <- round(min(main_df$core_q1_s, na.rm = TRUE) - 300, -2)
y_max  <- round(max(main_df$pan_q3_s,  na.rm = TRUE) + 300, -2)
y_step <- max(round((y_max - y_min) / 6, -2), 500)

# ── 8. FIGURE 1: Pan/Core curves — legend instead of inline labels ────────────
# Reshape to long format for legend
long_df <- rbind(
  data.frame(n        = main_df$n,
             med_s    = main_df$pan_med_s,
             q1_s     = main_df$pan_q1_s,
             q3_s     = main_df$pan_q3_s,
             curve    = "Pan-genome"),
  data.frame(n        = main_df$n,
             med_s    = main_df$core_med_s,
             q1_s     = main_df$core_q1_s,
             q3_s     = main_df$core_q3_s,
             curve    = "Core-genome")
)
long_df$curve <- factor(long_df$curve,
                        levels = c("Pan-genome", "Core-genome"))

p_main <- ggplot(long_df, aes(x = n, colour = curve, fill = curve)) +
  
  # ── IQR whiskers (bold, both curves) ──
  geom_linerange(aes(ymin = q1_s, ymax = q3_s),
                 linewidth = 0.55) +
  
  # ── Median lines ──
  geom_line(aes(y = med_s), linewidth = 0.7) +
  
  # ── Median squares ──
  geom_point(aes(y = med_s), shape = 22, size = 2.0) +
  
  # ── Manual colour/fill scales ──
  scale_colour_manual(
    name   = NULL,
    values = c("Pan-genome" = col_pan, "Core-genome" = col_core)
  ) +
  scale_fill_manual(
    name   = NULL,
    values = c("Pan-genome" = col_pan, "Core-genome" = col_core)
  ) +
  
  scale_x_continuous(
    breaks = c(1, seq(5, n_genomes, by = 5)),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    labels = comma,
    limits = c(y_min, y_max),
    breaks = seq(y_min, y_max, by = y_step),
    expand = c(0, 0)
  ) +
  
  labs(
    x       = "Number of genomes",
    y       = "Number of genes",
    title   = "Simulations of pangenome expansion\nand core-genome contraction",
    caption = paste0("For each given number of accessions, ", N_ITER,
                     " random combinations of accessions were sampled.")
  ) +
  
  theme_classic(base_size = 11) +
  theme(
    plot.title       = element_text(size = 10, face = "bold",
                                    hjust = 0.5, lineheight = 1.3),
    plot.caption     = element_text(size = 7.5, colour = "grey50", hjust = 0),
    axis.title       = element_text(size = 10),
    axis.text        = element_text(size = 9, colour = "black"),
    axis.line        = element_line(colour = "black", linewidth = 0.5),
    axis.ticks       = element_line(colour = "black"),
    panel.grid       = element_blank(),
    # Legend: top-right inside panel
    legend.position  = c(0.82, 0.55),
    legend.background = element_rect(fill = "white", colour = "grey80",
                                     linewidth = 0.3),
    legend.key.size  = unit(0.5, "cm"),
    legend.text      = element_text(size = 9),
    plot.margin      = margin(8, 8, 8, 8)
  )

# ── 9. Save Figure 1 ──────────────────────────────────────────────────────────
ggsave(OUT_MAIN_PDF, plot = p_main, width = 7, height = 5)
ggsave(OUT_MAIN_PNG, plot = p_main, width = 7, height = 5, dpi = 600)
cat(sprintf("\nFigure 1 saved:\n  %s\n  %s\n", OUT_MAIN_PDF, OUT_MAIN_PNG))

# ── 10. FIGURE 2: New gene discovery — fixed whiskers + no top/right border ────
# FIX 1: Use per-iteration differences for IQR whiskers (tighter, more accurate)
new_iter_mat <- pan_mat[-1, ] - pan_mat[-n_genomes, ]  # differences per iteration
new_q1_iter  <- apply(new_iter_mat, 1, quantile, 0.25)
new_q3_iter  <- apply(new_iter_mat, 1, quantile, 0.75)
new_med_iter <- apply(new_iter_mat, 1, median)

inset_df2 <- data.frame(
  n       = n_new,
  new_med = pmax(new_med_iter, 0),
  new_lo  = pmax(new_q1_iter,  1),
  new_hi  = pmax(new_q3_iter,  2)
)
inset_df2 <- inset_df2[inset_df2$n <= inset_max, ]

y_inset_max <- max(inset_df2$new_hi, na.rm = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
Y_FLOOR      <- 5       # y-axis starts this far above zero → whiskers float above x-axis
Y_TOP_MULT   <- 2.2     # upper y limit = max whisker × this
LABEL_Y_MULT <- 0.62    # annotation box y = y_top × this
LABEL_X      <- 5.2     # annotation box x position
PT_SIZE      <- 3.0     # triangle point size
LINE_WIDTH   <- 1.1     # fitted curve thickness
WHISKER_W    <- 0.7     # IQR whisker line thickness
BASE_SIZE    <- 11      # base font size
FIG_W        <- 5       # figure width (inches)
FIG_H        <- 4       # figure height (inches)

y_top   <- y_inset_max * Y_TOP_MULT
label_y <- y_top * LABEL_Y_MULT

inset_label <- sprintf(
  "n = k x N^(-a)\nk = %.0f +/- %.0f\na = %.2f +/- %.2f\nR2 = %.3f",
  k_fit, k_se, a_fit, a_se, r2)

p_discovery <- ggplot(inset_df2, aes(x = n)) +
  geom_linerange(aes(ymin = new_lo, ymax = new_hi),
                 colour = col_pan, linewidth = WHISKER_W, alpha = 0.9) +
  geom_point(aes(y = new_med),
             colour = col_pan, shape = 17, size = PT_SIZE) +
  geom_line(data = inset_curve, aes(x = n, y = new),
            colour = col_pan, linewidth = LINE_WIDTH) +
  annotate("label",
           x          = LABEL_X,
           y          = label_y,
           label      = inset_label,
           colour     = "black",
           fill       = "white",
           size       = 3.2,
           hjust      = 0,
           label.size = 0.3,
           family     = "mono") +
  scale_x_continuous(
    breaks = c(2, 4, 6, 8, 10),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  scale_y_continuous(
    labels = comma,
    limits = c(Y_FLOOR, y_top),   # Y_FLOOR lifts axis above zero
    expand = c(0, 0)
  ) +
  labs(
    x       = "Sample (N)",
    y       = "New gene families (n)",
    title   = "New gene family discovery rate",
    caption = paste0("Whiskers = IQR across ", N_ITER,
                     " samplings per N  |  Model: n = k x N^(-a)")
  ) +
  theme_classic(base_size = BASE_SIZE) +
  theme(
    plot.title         = element_text(size = 10, face = "bold", hjust = 0.5),
    plot.caption       = element_text(size = 7.5, colour = "grey50", hjust = 0),
    axis.title         = element_text(size = 10),
    axis.text          = element_text(size = 9, colour = "black"),
    axis.line.x.bottom = element_line(colour = "black", linewidth = 0.5),
    axis.line.y.left   = element_line(colour = "black", linewidth = 0.5),
    axis.line.x.top    = element_blank(),
    axis.line.y.right  = element_blank(),
    panel.border       = element_blank(),
    axis.ticks         = element_line(colour = "black"),
    panel.grid         = element_blank(),
    plot.margin        = margin(8, 8, 8, 8)
  )

# ── 11. Save Figure 2 ─────────────────────────────────────────────────────────
ggsave(OUT_INSET_PDF, plot = p_discovery, width = FIG_W, height = FIG_H)
ggsave(OUT_INSET_PNG, plot = p_discovery, width = FIG_W, height = FIG_H, dpi = 600)
cat(sprintf("\nFigure 2 saved:\n  %s\n  %s\n", OUT_INSET_PDF, OUT_INSET_PNG))

# ── 12. Summary ───────────────────────────────────────────────────────────────
cat("\n===================================================\n")
cat("  FINAL RESULTS SUMMARY\n")
cat("===================================================\n")
cat(sprintf("  Pan-genome size  (n=%d): %.0f\n", n_genomes, pan_s$med[n_genomes]))
cat(sprintf("  Core-genome size (n=%d): %.0f\n", n_genomes, core_s$med[n_genomes]))
cat(sprintf("  New gene model k:  %.1f +/- %.1f\n", k_fit, k_se))
cat(sprintf("  New gene model a:  %.4f +/- %.4f\n", a_fit, a_se))
cat(sprintf("  New gene model R2: %.4f\n", r2))
cat(sprintf("  Pangenome type:    %s\n",
            ifelse(a_fit > 1, "CLOSED (a > 1)", "OPEN (a < 1)")))
cat("===================================================\n")

```
