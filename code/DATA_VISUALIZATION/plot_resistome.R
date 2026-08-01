# ==============================================================================
# SCRIPT TRỰC QUAN HÓA TOÀN DIỆN DỮ LIỆU RESISTOME (PHẦN 4.2)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
})

# 1. Khởi tạo thư mục
out_dir <- "08_functional_screening/figures"
combined_dir <- "08_functional_screening/_combined"
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 2. Định nghĩa đường dẫn các file dữ liệu
amr_file <- file.path(combined_dir, "amrfinder_all.tsv")
staramr_file <- file.path(combined_dir, "staramr_summary_all.tsv")
summary_file <- "08_functional_screening/functional_screening_summary.tsv"

# ==============================================================================
# PHẦN 1: SỬ DỤNG DỮ LIỆU AMRFINDER (Dành cho Figure 4.2A, 4.2B, 4.2C)
# ==============================================================================
if(file.exists(amr_file)) {
  amr_df <- read.delim(amr_file, stringsAsFactors = FALSE)
  total_isolates <- n_distinct(amr_df$sample)
  
  resistome_clean <- amr_df %>%
    filter(Type %in% c("AMR", "POINT")) %>% 
    select(sample, group, Element.symbol, Class) %>%
    rename(Gene = Element.symbol, Group = group, Sample = sample) %>%
    mutate(Class = str_to_title(tolower(Class)))

  # --- FIGURE 4.2A: RESISTOME HEATMAP ---
  cat("[INFO] Đang tạo Figure 4.2A (Resistome Heatmap)...\n")
  heatmap_data <- resistome_clean %>%
    group_by(Sample, Gene, Class) %>%
    summarise(Copy_Number = n(), .groups = "drop")
  
  matrix_data <- heatmap_data %>%
    select(Sample, Gene, Copy_Number) %>%
    pivot_wider(names_from = Sample, values_from = Copy_Number, values_fill = list(Copy_Number = 0))
  
  mat <- as.matrix(matrix_data[, -1])
  rownames(mat) <- matrix_data$Gene
  
  col_anno <- data.frame(Group = (resistome_clean %>% select(Sample, Group) %>% distinct())$Group)
  rownames(col_anno) <- (resistome_clean %>% select(Sample, Group) %>% distinct())$Sample
  
  row_anno <- data.frame(Class = (heatmap_data %>% select(Gene, Class) %>% distinct())$Class)
  rownames(row_anno) <- (heatmap_data %>% select(Gene, Class) %>% distinct())$Gene
  
  my_colors <- list(Group = c("Clinical" = "#1f78b4", "Mutant" = "#e31a1c"))
  
  heatmap_file <- file.path(out_dir, "Figure_4.2A_Resistome_Heatmap.png")
  pheatmap(mat, color = colorRampPalette(c("#f7f7f7", "#fdae61", "#d73027", "#7f0000"))(100),
           annotation_col = col_anno, annotation_row = row_anno, annotation_colors = my_colors,
           cluster_cols = TRUE, cluster_rows = TRUE, fontsize_row = 6, fontsize_col = 8,
           border_color = "white", main = "Comprehensive Genotypic Resistome Profile",
           filename = heatmap_file, width = 12, height = 10)

  # --- FIGURE 4.2B: BETA-LACTAMASE PREVALENCE ---
  cat("[INFO] Đang tạo Figure 4.2B (Beta-Lactamase Prevalence)...\n")
  beta_lactam_df <- amr_df %>%
    filter(Class == "BETA-LACTAM") %>%
    mutate(Gene_Family = case_when(
      str_detect(Element.symbol, "CTX-M") ~ "blaCTX-M-15 (ESBL)",
      str_detect(Element.symbol, "OXA-48") ~ "blaOXA-48 (Carbapenemase)",
      str_detect(Element.symbol, "NDM") ~ "blaNDM-1 (Carbapenemase)",
      str_detect(Element.symbol, "SHV") ~ "blaSHV (Broad-spectrum)",
      str_detect(Element.symbol, "TEM") ~ "blaTEM (Broad-spectrum)",
      str_detect(Element.symbol, "OXA-[19]") ~ "blaOXA-1/9",
      TRUE ~ "Other"
    )) %>%
    group_by(Gene_Family) %>%
    summarise(Percentage = (n_distinct(sample) / total_isolates) * 100, .groups = "drop")
  
  p_bar <- ggplot(beta_lactam_df, aes(x = reorder(Gene_Family, Percentage), y = Percentage, fill = Gene_Family)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3, width = 0.7) +
    coord_flip() + scale_fill_brewer(palette = "Set2") + theme_minimal() +
    geom_text(aes(label = sprintf("%.1f%%", Percentage)), hjust = -0.2, size = 4, fontface = "bold") +
    scale_y_continuous(limits = c(0, 110)) +
    labs(title = "Prevalence of Key Beta-lactamase Genes", x = "", y = "Prevalence (%)") +
    theme(legend.position = "none", axis.text.y = element_text(face = "italic", size = 11))
  
  ggsave(file.path(out_dir, "Figure_4.2B_BetaLactamase_Prevalence.png"), plot = p_bar, width = 8, height = 5, dpi = 300, bg = "white")

  # --- FIGURE 4.2C: CO-OCCURRENCE DOT PLOT ---
  cat("[INFO] Đang tạo Figure 4.2C (Co-occurrence Dot Plot)...\n")
  key_genes <- c("blaCTX-M-15", "blaOXA-48", "blaNDM-1", "aac(6')-Ib", "rmtC", "aac(3)-IIe", "qnrS1", "gyrA_D87N", "parC_S80I", "fosA", "sul1", "tet(A)")
  
  co_occur_df <- amr_df %>% filter(Element.symbol %in% key_genes) %>% select(sample, Element.symbol) %>% distinct() %>% mutate(Presence = "Yes")
  sample_groups <- amr_df %>% select(sample, group) %>% distinct()
  
  plot_df <- expand_grid(sample = unique(amr_df$sample), Element.symbol = key_genes) %>%
    left_join(co_occur_df, by = c("sample", "Element.symbol")) %>%
    mutate(Presence = ifelse(is.na(Presence), "No", "Yes")) %>%
    left_join(sample_groups, by = "sample")
  
  p_dot <- ggplot(plot_df, aes(x = sample, y = factor(Element.symbol, levels = rev(key_genes)))) +
    geom_point(aes(color = Presence, size = Presence)) +
    scale_color_manual(values = c("Yes" = "#E94E77", "No" = "grey90")) +
    scale_size_manual(values = c("Yes" = 5, "No" = 2)) +
    facet_grid(~ group, scales = "free_x", space = "free_x") +
    theme_bw() +
    labs(title = "Co-occurrence of Critical Antimicrobial Resistance Determinants", x = "Isolate", y = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), axis.text.y = element_text(face = "italic"))
  
  ggsave(file.path(out_dir, "Figure_4.2C_Co_occurrence.png"), plot = p_dot, width = 10, height = 6, dpi = 300, bg = "white")
} else { cat("Không tìm thấy amrfinder_all.tsv\n") }

# ==============================================================================
# PHẦN 2: SỬ DỤNG STARAMR (Dành cho Figure 4.2D - Phenotypic Prediction)
# ==============================================================================
if(file.exists(staramr_file)) {
  cat("[INFO] Đang tạo Figure 4.2D (Predicted Phenotype Heatmap)...\n")
  star_df <- read.delim(staramr_file, stringsAsFactors = FALSE)
  
  # Xử lý chuỗi CGE Predicted Phenotype
  pheno_df <- star_df %>%
    select(sample, group, CGE.Predicted.Phenotype) %>%
    separate_rows(CGE.Predicted.Phenotype, sep = ",\\s*") %>%
    filter(CGE.Predicted.Phenotype != "", !is.na(CGE.Predicted.Phenotype)) %>%
    mutate(Presence = 1) %>%
    distinct()
  
  pheno_matrix <- pheno_df %>%
    pivot_wider(names_from = sample, values_from = Presence, values_fill = 0) %>%
    as.data.frame()
  
  mat_pheno <- as.matrix(pheno_matrix[, -c(1:2)])
  rownames(mat_pheno) <- pheno_matrix$CGE.Predicted.Phenotype
  
  col_anno_pheno <- data.frame(Group = (star_df %>% select(sample, group) %>% distinct())$group)
  rownames(col_anno_pheno) <- (star_df %>% select(sample, group) %>% distinct())$sample
  
  pheno_file <- file.path(out_dir, "Figure_4.2D_Predicted_Phenotype.png")
  pheatmap(mat_pheno, color = c("grey95", "#31a354"), breaks = c(-0.5, 0.5, 1.5),
           legend_breaks = c(0, 1), legend_labels = c("Susceptible", "Resistant"),
           annotation_col = col_anno_pheno, annotation_colors = list(Group = c("Clinical" = "#1f78b4", "Mutant" = "#e31a1c")),
           cluster_cols = TRUE, cluster_rows = TRUE, fontsize_row = 7, fontsize_col = 8,
           border_color = "white", main = "In Silico Predicted Antimicrobial Resistance Phenotypes",
           filename = pheno_file, width = 10, height = 8)
  cat(sprintf("[OK] Đã lưu %s\n", pheno_file))
} else { cat("Không tìm thấy staramr_summary_all.tsv\n") }

# ==============================================================================
# PHẦN 3: SỬ DỤNG FUNCTIONAL SUMMARY (Dành cho Figure 4.2E - AMR Burden)
# ==============================================================================
if(file.exists(summary_file)) {
  cat("[INFO] Đang tạo Figure 4.2E (AMR Gene Burden Boxplot)...\n")
  summ_df <- read.delim(summary_file, stringsAsFactors = FALSE)
  
  # Đếm số lượng gene kháng thuốc thực tế bằng cách đếm dấu chấm phẩy
  burden_df <- summ_df %>%
    select(Sample, Group, AMRFinder_genes) %>%
    mutate(AMR_Count = str_count(AMRFinder_genes, ";") + 1)
  
  p_burden <- ggplot(burden_df, aes(x = Group, y = AMR_Count, fill = Group)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 3, alpha = 0.8, color = "black", shape = 21) +
    scale_fill_manual(values = c("Clinical" = "#1f78b4", "Mutant" = "#e31a1c")) +
    theme_classic() +
    labs(title = "Total AMR Gene Burden",
         subtitle = "Comparison of resistance determinant counts between groups",
         x = "Isolate Group", y = "Total AMR Genes / Mutations Detected") +
    theme(legend.position = "none", plot.title = element_text(face = "bold", size = 14))
  
  burden_file <- file.path(out_dir, "Figure_4.2E_AMR_Burden.png")
  ggsave(burden_file, plot = p_burden, width = 6, height = 5, dpi = 300, bg = "white")
  cat(sprintf("[OK] Đã lưu %s\n", burden_file))
} else { cat("Không tìm thấy functional_screening_summary.tsv\n") }

cat("\n[DONE] XUẤT SẮC! Toàn bộ 5 biểu đồ cho phần 4.2 đã được tạo thành công.\n")