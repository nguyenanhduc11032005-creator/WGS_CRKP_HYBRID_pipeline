# ==============================================================================
# SCRIPT HOÀN CHỈNH CHO PHẦN 4.3: PLASMID COPY NUMBER (PCN)
# ==============================================================================

# 1. THIẾT LẬP THƯ MỤC LƯU TRỮ (WORKING DIRECTORY)
setwd("/mnt/d18t/ANHDUC/WGS_CRKP_HYBRID")

# 2. LOAD LIBRARIES
library(ggplot2)
library(dplyr)
library(gt)
library(webshot2) # Đã cài qua conda

# ==============================================================================
# 3. KHAI BÁO DỮ LIỆU
# ==============================================================================

# Dữ liệu cho Figure 4A & Table 1 (Clinical Baseline)
clinical_data <- data.frame(
  Isolate = factor(c("Kpn-15", "Kpn-35", "Kpn-59", "Kpn-24", "Kpn-9", "Kpn-72", 
                     "Kpn-3", "Kpn-29", "Kpn-73", "Kpn-17", "Kpn-28", "Kpn-7", "Kpn-62"),
                   levels = c("Kpn-15", "Kpn-35", "Kpn-59", "Kpn-24", "Kpn-9", "Kpn-72", 
                              "Kpn-3", "Kpn-29", "Kpn-73", "Kpn-17", "Kpn-28", "Kpn-7", "Kpn-62")),
  PCN = c(0.7, 0.8, 0.8, 1.0, 1.0, 1.2, 1.6, 1.8, 2.5, 2.6, 5.3, 6.0, 6.6),
  Gene = c("blaOXA-48", "blaOXA-48", "blaOXA-48", "blaOXA-48", "blaOXA-48", "blaCTX-M-15", 
           "blaOXA-48", "blaOXA-48", "blaCTX-M-15", "blaCTX-M-15", "blaOXA-48", "blaOXA-48", "blaOXA-48")
)

table1_data <- data.frame(
  Isolate = c("Kpn-15", "Kpn-35", "Kpn-59", "Kpn-24", "Kpn-9", "Kpn-72", 
              "Kpn-3", "Kpn-29", "Kpn-73", "Kpn-17", "Kpn-28", "Kpn-7", 
              "Kpn-62", "Kpn-12", "Kpn-26"),
  BetaLactamase = c("*bla*~OXA-48~ / *bla*~NDM-1~", "*bla*~OXA-48~", "*bla*~OXA-48~", 
                    "*bla*~OXA-48~", "*bla*~OXA-48~", "*bla*~CTX-M-15~", 
                    "*bla*~OXA-48~", "*bla*~OXA-48~", "*bla*~CTX-M-15~", 
                    "*bla*~CTX-M-15~", "*bla*~OXA-48~", "*bla*~OXA-48~", 
                    "*bla*~OXA-48~", "*bla*~CTX-M-15~", "*bla*~CTX-M-15~"),
  Replicon = c("IncL", "IncL", "IncL", "IncL", "IncL", "IncR", "IncL", 
               "IncL", "IncFII", "IncFII", "IncL", "IncL", "IncL", 
               "IncFII / IncR", "IncFII"),
  PCN = c("0.7X", "0.8X", "0.8X", "1.0X", "1.0X", "1.2X", "1.6X", 
          "1.8X", "2.5X", "2.6X", "5.3X", "6.0X", "6.6X", 
          "Chromosomally Integrated", "Fragmented Assembly")
)

# Dữ liệu cho Figure 4B & Table 2 (Mutant Amplification)
mutant_data <- data.frame(
  Pair = factor(c("Kpn-24", "Kpn-24", "Kpn-17", "Kpn-17", "Kpn-72", "Kpn-72", "Kpn-73", "Kpn-73"),
                levels = c("Kpn-24", "Kpn-17", "Kpn-72", "Kpn-73")),
  Status = factor(c("Native", "Mutant", "Native", "Mutant", 
                    "Native", "Mutant", "Native", "Mutant"),
                  levels = c("Native", "Mutant")),
  PCN = c(1.0, 3.2, 2.6, 9.4, 1.2, 6.1, 2.5, 16.5)
)

table2_data <- data.frame(
  Pair = c("Kpn-24 vs 24.4", "Kpn-17 vs 17.11", "Kpn-72 vs 72.11", "Kpn-73 vs 73.9"),
  Replicon = c("IncL", "IncFII", "IncR", "IncFII"),
  BetaLactamase = c("*bla*~OXA-48~", "*bla*~CTX-M-15~", "*bla*~CTX-M-15~", "*bla*~CTX-M-15~"),
  Native_PCN = c("1.0X", "2.6X", "1.2X", "2.5X"),
  Mutant_PCN = c("3.2X", "9.4X", "6.1X", "16.5X"),
  Fold_Increase = c("**3.2x**", "**3.6x**", "**5.1x**", "**6.6x**")
)

# ==============================================================================
# 4. VẼ VÀ LƯU FILE
# ==============================================================================

# --- Figure 4A ---
plot_4A <- ggplot(clinical_data, aes(x = Isolate, y = PCN, fill = Gene)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  scale_fill_manual(values = c("blaCTX-M-15" = "#F4B183", "blaOXA-48" = "#8FAADC")) +
  geom_hline(yintercept = 2.0, linetype = "dashed", color = "red", alpha = 0.6) +
  labs(title = "Figure 4A. Baseline Plasmid Copy Numbers in Clinical Isolates",
       x = "Clinical Isolate", y = "Plasmid Copy Number (Relative to Chromosome)", fill = "Resistance Gene") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
        axis.text.y = element_text(face = "bold", color = "black"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 15)),
        legend.position = c(0.2, 0.8),
        legend.background = element_rect(color = "black", linewidth = 0.5))

ggsave("Figure_4A_Baseline_PCN.png", plot = plot_4A, width = 10, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 4. VẼ VÀ LƯU FIGURE 4B (ĐÃ FIX LỖI ĐÈ CHỮ)
# ------------------------------------------------------------------------------
plot_4B <- ggplot(mutant_data, aes(x = Pair, y = PCN, fill = Status)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +
  geom_text(aes(label = sprintf("%.1fX", PCN)), 
            position = position_dodge(width = 0.8), 
            vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Native" = "#9BC2E6", "Mutant" = "#E06666")) +
  
  # CÁCH FIX: Định dạng lại nhãn trục X bằng cách thêm \n (xuống dòng) cho Replicon
  scale_x_discrete(labels = c(
    "Kpn-24" = "Kpn-24\n(IncL)", 
    "Kpn-17" = "Kpn-17\n(IncFII)", 
    "Kpn-72" = "Kpn-72\n(IncR)", 
    "Kpn-73" = "Kpn-73\n(IncFII)"
  )) +
  
  coord_cartesian(ylim = c(0, 18)) +
  labs(
    title = "Figure 4B. PCN Amplification Following Meropenem Exposure",
    x = "Isolate Pair (Replicon Type)",
    y = "Plasmid Copy Number",
    fill = "Isolate Status"
  ) +
  theme_classic(base_size = 14) +
  theme(
    # Tinh chỉnh lại lề của trục X cho cân đối
    axis.text.x = element_text(face = "bold", color = "black", margin = margin(t = 5, b = 5)),
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.title = element_text(face = "bold", margin = margin(t = 10)),
    plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 15)),
    legend.position = c(0.15, 0.85),
    legend.background = element_rect(color = "black", linewidth = 0.5)
  )

print(plot_4B)
ggsave("Figure_4B_PCN_Amplification.png", plot = plot_4B, width = 8, height = 6, dpi = 300)

# --- Table 1 ---
table1_gt <- table1_data %>%
  gt() %>%
  tab_header(title = md("**Table 1. Baseline Copy Numbers of Primary Resistance Plasmids in 15 Clinical Isolates**")) %>%
  cols_label(Isolate = md("**Clinical Isolate**"), BetaLactamase = md("**Primary Resistance Driver**"),
             Replicon = md("**Target Plasmid Replicon**"), PCN = md("**Computed PCN (Relative to Chromosome)**")) %>%
  fmt_markdown(columns = BetaLactamase) %>%
  tab_options(table.border.top.color = "black", table.border.top.width = px(2),
              table.border.bottom.color = "black", table.border.bottom.width = px(2),
              heading.align = "left",
              column_labels.border.bottom.color = "black", column_labels.border.bottom.width = px(2),
              table_body.border.bottom.color = "black", table_body.border.bottom.width = px(2)) %>%
  cols_align(align = "center", columns = c(Replicon, PCN)) %>%
  cols_align(align = "left", columns = c(Isolate, BetaLactamase))

gtsave(table1_gt, "Table_1_Baseline_PCN.png", vwidth = 800)

# --- Table 2 ---
table2_gt <- table2_data %>%
  gt() %>%
  tab_header(title = md("**Table 2. Comparative Amplification of Plasmid Copy Numbers (PCN) Following Meropenem Exposure**")) %>%
  cols_label(Pair = md("**Pairwise Comparison**"), Replicon = md("**Target Plasmid Replicon**"),
             BetaLactamase = md("**Primary Beta-lactamase**"), Native_PCN = md("**Native Strain PCN**"),
             Mutant_PCN = md("**Mutant Strain PCN**"), Fold_Increase = md("**PCN Fold Increase**")) %>%
  fmt_markdown(columns = c(BetaLactamase, Fold_Increase)) %>%
  tab_options(table.border.top.color = "black", table.border.top.width = px(2),
              table.border.bottom.color = "black", table.border.bottom.width = px(2),
              heading.align = "left",
              column_labels.border.bottom.color = "black", column_labels.border.bottom.width = px(2),
              table_body.border.bottom.color = "black", table_body.border.bottom.width = px(2)) %>%
  cols_align(align = "center", columns = c(Replicon, Native_PCN, Mutant_PCN, Fold_Increase)) %>%
  cols_align(align = "left", columns = c(Pair, BetaLactamase)) %>%
  tab_style(style = cell_fill(color = "#f2f2f2"), locations = cells_body(columns = Fold_Increase))

gtsave(table2_gt, "Table_2_PCN_Amplification.png", vwidth = 800)

print("Đã lưu thành công 4 file (2 biểu đồ và 2 bảng) vào thư mục /mnt/d18t/ANHDUC/WGS_CRKP_HYBRID")