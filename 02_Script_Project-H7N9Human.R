# =============================================================================
# PROJECT B - H7N9 HUMAN (bulk transcriptomics)
# Design: 5 conditions (mock, H3N2, H7N7, H7N9, H5N1/VN1203)
#         x 4 time points (3h, 7h, 12h, 24h) x 4 replicates
# =============================================================================


# 0. SETUP ---------------------------------------------------------------------

## 0.1 Install packages (run once, then keep commented) ----
# install.packages(c("ggplot2", "ggrepel", "enrichR", "MASS", "BiocManager"))
# BiocManager::install("limma")

## 0.2 Load packages ----
library(limma)
library(ggplot2)
library(ggrepel)
library(enrichR)
library(MASS)

## 0.3 Working directory & output folder ----
project_dir <- "~/Desktop/M2/AIDA/Data analysis in immunology workshop/Project Virus"
setwd(project_dir)
out_dir <- file.path(project_dir, "Results")          # all figures/tables go here
dir.create(out_dir, showWarnings = FALSE)

# helper: save a ggplot into Results/ (PNG + PDF)
save_plot <- function(p, name, width = 7, height = 5) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = width, height = height)
  invisible(p)
}

## 0.4 Global settings (colours) ----
virus_cols <- c("mock" = "grey50", "H3N2" = "#1f77b4", "H7N7" = "#2ca02c",
                "H7N9" = "#d62728", "H5N1 (VN1203)" = "#ff7f0e")


# 1. DATA IMPORT ---------------------------------------------------------------

data <- read.delim("01_Data.txt")
data <- as.matrix(data)
dim(data)   # genes x samples


# 2. SAMPLE METADATA -----------------------------------------------------------

## 2.1 Build metadata from column names (virus_time_replicate) ----
parts <- do.call(rbind, strsplit(colnames(data), "_"))
meta  <- data.frame(sample = colnames(data),
                    virus  = parts[, 1],
                    hour   = as.integer(sub("h$", "", parts[, 2])),
                    rep    = parts[, 3],
                    stringsAsFactors = FALSE)
meta$virus[meta$virus == "VN1203"] <- "H5N1 (VN1203)"
meta$virus  <- factor(meta$virus, levels = c("mock", "H3N2", "H7N7", "H7N9", "H5N1 (VN1203)"))
meta$hour_f <- factor(meta$hour, levels = c(3, 7, 12, 24), labels = paste0(c(3, 7, 12, 24), " h"))
rownames(meta) <- meta$sample

## 2.2 Check missing samples ----
# 4 replicates expected per group; H3N2 and VN1203 have a missing sample at 3h
table(meta$virus, meta$hour_f)


# 3. FILTERING & TRANSFORMATION ------------------------------------------------

data <- data[apply(data, 1, max) > 0, ]   # remove genes never expressed
logx <- log2(data + 1)                    # log2 transform
dim(logx)


# 4. PCA -----------------------------------------------------------------------

## 4.1 Functions ----
run_pca <- function(mat, meta, n_top = 2000, scale. = TRUE) {
  top <- names(sort(apply(mat, 1, var), decreasing = TRUE))[seq_len(n_top)]
  p   <- prcomp(t(mat[top, ]), center = TRUE, scale. = scale.)
  list(scores = cbind(as.data.frame(p$x[, 1:10]), meta[rownames(p$x), ]),
       var    = 100 * p$sdev^2 / sum(p$sdev^2),
       load   = p$rotation[, 1:10])
}

plot_pca <- function(res, x = "PC1", y = "PC2", title = "") {
  xi <- as.integer(sub("PC", "", x)); yi <- as.integer(sub("PC", "", y))
  ggplot(res$scores, aes(.data[[x]], .data[[y]], colour = virus, shape = hour_f)) +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_vline(xintercept = 0, colour = "grey85") +
    geom_point(size = 3, alpha = .85) +
    scale_colour_manual(values = virus_cols) +
    labs(title = title, colour = "Virus", shape = "Time",
         x = sprintf("%s (%.1f%%)", x, res$var[xi]),
         y = sprintf("%s (%.1f%%)", y, res$var[yi])) +
    theme_bw()
}

## 4.2 PCA - all samples ----
res <- run_pca(logx, meta)
p12 <- plot_pca(res, "PC1", "PC2", "PCA - all samples"); print(p12)
p13 <- plot_pca(res, "PC1", "PC3", "PCA - all samples"); print(p13)
p23 <- plot_pca(res, "PC2", "PC3", "PCA - all samples"); print(p23)

save_plot(p12, "PCA_all_PC1_PC2")
save_plot(p13, "PCA_all_PC1_PC3")
save_plot(p23, "PCA_all_PC2_PC3")
write.csv(res$scores, file.path(out_dir, "PCA_all_scores.csv"))


# 5. DIFFERENTIAL EXPRESSION (virus vs mock, per time point) -------------------
# Data are already normalised (non-integer values) -> limma on log2 data
# (eBayes with trend = TRUE, i.e. "limma-trend"), not DESeq2/edgeR.

de_dir <- file.path(out_dir, "DE")
dir.create(de_dir, showWarnings = FALSE)

## 5.1 Parameters ----
lfc_thr  <- 1      # |log2FC| > 1  (fold change > 2)
padj_thr <- 0.05   # BH-adjusted p-value

## 5.2 Remove lowly expressed genes ----
expr <- logx[rowMeans(logx) > 1, rownames(meta)]
dim(expr)

## 5.3 Design: one group per virus x time ----
meta$group <- factor(paste(sub("_.*", "", meta$sample), meta$hour, sep = "_"))  # e.g. H7N9_24
design <- model.matrix(~ 0 + group, data = meta)
colnames(design) <- levels(meta$group)

fit <- lmFit(expr, design)

## 5.4 Contrasts: each virus vs mock at the same time point ----
viruses <- c("H3N2", "H7N7", "H7N9", "VN1203")
hours   <- c(3, 7, 12, 24)
grid    <- expand.grid(virus = viruses, hour = hours, stringsAsFactors = FALSE)
contr_names <- paste0(grid$virus, "_", grid$hour, "h_vs_mock")
contr_def   <- paste0(grid$virus, "_", grid$hour, " - mock_", grid$hour)

cm <- makeContrasts(contrasts = contr_def, levels = design)
colnames(cm) <- contr_names
fit2 <- eBayes(contrasts.fit(fit, cm), trend = TRUE)

## 5.5 Volcano plot function ----
plot_volcano <- function(tt, title, n_lab = 15) {
  lab <- head(tt[tt$status != "NS", ], n_lab)   # tt is sorted by p-value
  ggplot(tt, aes(logFC, -log10(P.Value), colour = status)) +
    geom_point(size = .8, alpha = .6) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = 2, colour = "grey50") +
    geom_text_repel(data = lab, aes(label = gene), size = 3, colour = "black",
                    max.overlaps = Inf) +
    scale_colour_manual(values = c(Up = "#d62728", Down = "#1f77b4", NS = "grey75")) +
    labs(title = title, x = "log2 fold change (virus / mock)",
         y = "-log10 p-value", colour = NULL) +
    theme_bw()
}

## 5.6 Tables + volcanos for every contrast (saved in Results/DE) ----
de_list <- list()
for (cn in contr_names) {
  tt <- topTable(fit2, coef = cn, number = Inf, sort.by = "P")
  tt$gene   <- rownames(tt)
  tt$status <- ifelse(tt$adj.P.Val < padj_thr & tt$logFC >  lfc_thr, "Up",
               ifelse(tt$adj.P.Val < padj_thr & tt$logFC < -lfc_thr, "Down", "NS"))
  de_list[[cn]] <- tt

  write.csv(tt, file.path(de_dir, paste0("DE_", cn, ".csv")), row.names = FALSE)
  save_plot(plot_volcano(tt, gsub("_", " ", cn)), file.path("DE", paste0("Volcano_", cn)))
}

## 5.7 Number of DEGs per virus and time point ----
deg_summary <- data.frame(
  contrast = contr_names,
  virus    = grid$virus,
  hour     = grid$hour,
  Up       = sapply(de_list, function(t) sum(t$status == "Up")),
  Down     = sapply(de_list, function(t) sum(t$status == "Down"))
)
deg_summary
write.csv(deg_summary, file.path(de_dir, "DEG_summary.csv"), row.names = FALSE)

deg_long <- rbind(data.frame(deg_summary[, c("virus", "hour")], direction = "Up",   n = deg_summary$Up),
                  data.frame(deg_summary[, c("virus", "hour")], direction = "Down", n = deg_summary$Down))
deg_long$virus  <- factor(sub("VN1203", "H5N1 (VN1203)", deg_long$virus), levels = names(virus_cols)[-1])
deg_long$hour_f <- factor(paste0(deg_long$hour, " h"), levels = paste0(hours, " h"))

p_deg <- ggplot(deg_long, aes(hour_f, n, fill = virus)) +
  geom_col(position = "dodge") +
  facet_wrap(~ direction) +
  scale_fill_manual(values = virus_cols) +
  labs(title = "Number of DEGs vs mock",
       subtitle = sprintf("adj. p < %.2f, |log2FC| > %g", padj_thr, lfc_thr),
       x = "Time post-infection", y = "Number of genes", fill = "Virus") +
  theme_bw()
print(p_deg)
save_plot(p_deg, "DEG_number_per_virus_time", width = 9, height = 4.5)


# 6. FUNCTIONAL ENRICHMENT (enrichR) -------------------------------------------
# TODO


# 7. FIGURES / EXPORT ----------------------------------------------------------
# TODO
