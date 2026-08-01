# Nạp thư viện
library(ggplot2)
library(dplyr)
library(tidyr) # Thêm tidyr để xử lý dữ liệu cho biểu đồ Boxplot

# ==========================================
# 1. NẠP DỮ LIỆU THỰC TẾ TỪ PIPELINE
# ==========================================

# Dữ liệu Assembly QC (Bao gồm cả CheckM2)
qc_data_text <- "Sample\tGroup\tStatus\tQUAST_N50\tQUAST_Contigs\tBUSCO_Completeness\tCheckM2_Completeness\tCheckM2_Contamination\tMessage
Kpn-17\tClinical\tSUCCESS\t5269722\t8\t98.9%\t100.0\t0.89\tOK
Kpn-12\tClinical\tSUCCESS\t5290859\t8\t98.9%\t100.0\t0.08\tOK
Kpn-15\tClinical\tSUCCESS\t3562395\t20\t98.7%\t100.0\t0.22\tOK
Kpn-26\tClinical\tSUCCESS\t5259925\t4\t98.9%\t100.0\t0.5\tOK
Kpn-24\tClinical\tSUCCESS\t5130949\t11\t98.9%\t100.0\t0.52\tOK
Kpn-28\tClinical\tSUCCESS\t5511943\t9\t98.7%\t100.0\t0.15\tOK
Kpn-29\tClinical\tSUCCESS\t5312522\t12\t98.7%\t100.0\t0.14\tOK
Kpn-3\tClinical\tSUCCESS\t5266793\t9\t98.9%\t100.0\t0.13\tOK
Kpn-35\tClinical\tSUCCESS\t5263749\t9\t98.9%\t100.0\t0.12\tOK
Kpn-59\tClinical\tSUCCESS\t5267070\t9\t98.9%\t100.0\t0.11\tOK
Kpn-62\tClinical\tSUCCESS\t5256629\t10\t98.9%\t100.0\t0.19\tOK
Kpn-7\tClinical\tSUCCESS\t5312301\t10\t98.9%\t100.0\t0.11\tOK
Kpn-72\tClinical\tSUCCESS\t4019176\t19\t98.9%\t100.0\t0.25\tOK
Kpn-9\tClinical\tSUCCESS\t5406725\t8\t98.9%\t100.0\t0.11\tOK
Kpn-73\tClinical\tSUCCESS\t5174875\t10\t98.9%\t100.0\t0.9\tOK
Kpn-17.11\tMutant\tSUCCESS\t5269147\t20\t98.9%\t100.0\t0.85\tOK
Kpn-24.4\tMutant\tSUCCESS\t5266403\t10\t98.9%\t100.0\t0.52\tOK
Kpn-72.11\tMutant\tSUCCESS\t5402003\t6\t98.7%\t100.0\t0.16\tOK
Kpn-73.9\tMutant\tSUCCESS\t5271629\t8\t98.9%\t100.0\t0.89\tOK"

st_data_text <- "Sample\tKleborate_ST
Kpn-12\tST101
Kpn-15\tST15
Kpn-17\tST15
Kpn-24\tST147
Kpn-26\tST525
Kpn-28\tST395
Kpn-29\tST395
Kpn-3\tST101
Kpn-35\tST101
Kpn-59\tST101
Kpn-62\tST101
Kpn-7\tST101
Kpn-72\tST101
Kpn-73\tST15
Kpn-9\tST101
Kpn-17.11\tST15
Kpn-24.4\tST147
Kpn-72.11\tST101
Kpn-73.9\tST15"

# Đọc dữ liệu
df_qc <- read.table(text = qc_data_text, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
df_st <- read.table(text = st_data_text, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Tiền xử lý dữ liệu QC
df_qc$QUAST_N50_Mb <- df_qc$QUAST_N50 / 1000000
# Xóa dấu '%' trong cột BUSCO và ép kiểu sang dạng số (numeric)
df_qc$BUSCO_Completeness <- as.numeric(gsub("%", "", df_qc$BUSCO_Completeness))

# ==========================================
# 2. FIGURE 1A & 1B (Giữ nguyên như cũ)
# ==========================================
plot_1a <- ggplot(df_qc, aes(x = QUAST_N50_Mb, y = QUAST_Contigs, color = Group, shape = Group)) +
  geom_point(size = 5, alpha = 0.85) +
  scale_color_manual(values = c("Clinical" = "#2c7bb6", "Mutant" = "#d7191c")) +
  theme_bw(base_size = 14) +
  labs(title = "Assembly Quality: N50 vs Contig Count", x = "QUAST N50 (Mb)", y = "Number of Contigs") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16), axis.title = element_text(face = "bold"), legend.position = "top", legend.title = element_blank(), panel.grid.minor = element_blank())

ggsave("Figure_1A_Assembly_Quality.png", plot = plot_1a, width = 8, height = 6, dpi = 300)

st_counts <- df_st %>% count(Kleborate_ST, name = "Count") %>% arrange(desc(Count))
plot_1b <- ggplot(st_counts, aes(x = reorder(Kleborate_ST, -Count), y = Count, fill = Kleborate_ST)) +
  geom_col(color = "black", width = 0.6) +
  geom_text(aes(label = Count), vjust = -0.5, fontface = "bold", size = 5) +
  scale_fill_viridis_d(option = "D") + 
  theme_classic(base_size = 14) +
  labs(title = "Distribution of Sequence Types (ST)", x = "Sequence Type (ST)", y = "Number of Isolates") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16), axis.title = element_text(face = "bold"), axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"), legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

ggsave("Figure_1B_ST_Distribution.png", plot = plot_1b, width = 8, height = 6, dpi = 300)

# ==========================================
# 3. VẼ FIGURE 1C: GENOME COMPLETENESS
# ==========================================
# Gộp 2 cột BUSCO và CheckM2 Completeness lại để vẽ Boxplot chung
df_comp <- df_qc %>%
  select(Sample, BUSCO_Completeness, CheckM2_Completeness) %>%
  pivot_longer(cols = c("BUSCO_Completeness", "CheckM2_Completeness"), 
               names_to = "Tool", values_to = "Completeness") %>%
  mutate(Tool = ifelse(Tool == "BUSCO_Completeness", "BUSCO", "CheckM2"))

plot_1c <- ggplot(df_comp, aes(x = Tool, y = Completeness, fill = Tool)) +
  geom_boxplot(width = 0.4, alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 3, alpha = 0.6, color = "black") +
  scale_fill_manual(values = c("BUSCO" = "#1b9e77", "CheckM2" = "#7570b3")) +
  theme_bw(base_size = 14) +
  labs(
    title = "Genome Completeness Evaluation",
    x = "Evaluation Tool",
    y = "Completeness (%)"
  ) +
  scale_y_continuous(limits = c(98, 100.5), breaks = seq(98, 100, by = 0.5)) + # Focus vào dải 98-100%
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave("Figure_1C_Genome_Completeness.png", plot = plot_1c, width = 6, height = 6, dpi = 300)
cat("[OK] Da luu Figure_1C_Genome_Completeness.png\n")

# ==========================================
# 4. VẼ FIGURE 1D: GENOME CONTAMINATION
# ==========================================
plot_1d <- ggplot(df_qc, aes(x = Group, y = CheckM2_Contamination, fill = Group)) +
  geom_boxplot(width = 0.4, alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 3, alpha = 0.6, color = "black") +
  scale_fill_manual(values = c("Clinical" = "#2c7bb6", "Mutant" = "#d7191c")) +
  geom_hline(yintercept = 5.0, linetype = "dashed", color = "red", size = 1) + # Đường giới hạn chuẩn 5%
  annotate("text", x = 1.5, y = 4.8, label = "MIMAG High-Quality Cutoff (< 5%)", color = "red", fontface = "italic") +
  theme_bw(base_size = 14) +
  labs(
    title = "Genome Contamination (CheckM2)",
    x = "Sample Group",
    y = "Contamination (%)"
  ) +
  scale_y_continuous(limits = c(0, 5.5)) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave("Figure_1D_Genome_Contamination.png", plot = plot_1d, width = 6, height = 6, dpi = 300)
cat("[OK] Da luu Figure_1D_Genome_Contamination.png\n")