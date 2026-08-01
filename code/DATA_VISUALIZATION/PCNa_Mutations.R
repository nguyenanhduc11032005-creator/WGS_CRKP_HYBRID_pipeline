library(ggplot2)
library(dplyr)

# 1. TẠO DỮ LIỆU TỪ DEEP_DIVE_SUMMARY CHO CẢ 4 CẶP ISOLATE
depth_data <- data.frame(
  Pair = c(rep("Kpn-24 (IncL)", 4), rep("Kpn-17 (IncFII)", 4), 
           rep("Kpn-72 (IncR)", 4), rep("Kpn-73 (IncFII)", 4)),
  Status = rep(c("Native", "Native", "Mutant", "Mutant"), 4),
  Replicon = rep(c("Chromosome", "Target Plasmid", "Chromosome", "Target Plasmid"), 4),
  Absolute_Depth = c(
    # Kpn-24 (Native Chr: 132.3 | Mutant Chr: 134.9 | PCN: 1.0X -> 3.2X)
    132.3, 132.3 * 1.0,   
    134.9, 134.9 * 3.2,   
    
    # Kpn-17 (Native Chr: 114.1 | Mutant Chr: 152.9 | PCN: 2.6X -> 9.4X)
    114.1, 114.1 * 2.6,   
    152.9, 152.9 * 9.4,   
    
    # Kpn-72 (Native Chr: 303.5 | Mutant Chr: 135.7 | PCN: 1.2X -> 6.1X)
    303.5, 303.5 * 1.2,   
    135.7, 135.7 * 6.1,   
    
    # Kpn-73 (Native Chr: 158.2 | Mutant Chr: 117.0 | PCN: 2.5X -> 16.5X)
    158.2, 158.2 * 2.5,   
    117.0, 117.0 * 16.5   
  )
)

# 2. XỬ LÝ ĐỊNH DẠNG VÀ THỨ TỰ
depth_data$Status <- factor(depth_data$Status, levels = c("Native", "Mutant"))
depth_data$Replicon <- factor(depth_data$Replicon, levels = c("Chromosome", "Target Plasmid"))
depth_data$Pair <- factor(depth_data$Pair, levels = c("Kpn-24 (IncL)", "Kpn-17 (IncFII)", "Kpn-72 (IncR)", "Kpn-73 (IncFII)"))

# 3. VẼ BIỂU ĐỒ 6B VỚI 4 PANEL
fig6b <- ggplot(depth_data, aes(x = Status, y = Absolute_Depth, fill = Replicon)) +
  geom_col(position = position_dodge(width = 0.8), color = "black", width = 0.7, alpha = 0.9) +
  
  # Chia thành 4 panel (sử dụng nrow = 1 để dàn ngang cho đẹp)
  facet_wrap(~Pair, scales = "free_x", nrow = 1) +
  
  geom_text(aes(label = paste0(round(Absolute_Depth, 0), "X")),
            position = position_dodge(width = 0.8), 
            vjust = -0.6, size = 4, fontface = "bold") +
            
  scale_fill_manual(values = c("Chromosome" = "#ABB2B9", "Target Plasmid" = "#F39C12")) +
  
  # Mở rộng trục Y để text không bị cắt
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  
  labs(
    title = "Figure 6B: Absolute Sequencing Depth Verification",
    subtitle = "Confirming true plasmid amplification independent of chromosomal coverage",
    x = "Isolate Status",
    y = "Mean Sequencing Depth (X)",
    fill = "Genomic Element"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30", margin = margin(b = 15)),
    axis.text.x = element_text(face = "bold", color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 11),
    axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
    axis.title.y = element_text(margin = margin(r = 10), face = "bold"),
    strip.text = element_text(face = "bold", size = 13),
    strip.background = element_rect(fill = "#EBF5FB", color = "black"),
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )

print(fig6b)

# 4. XUẤT ẢNH CHẤT LƯỢNG CAO (Nên tăng width lên 12 để ảnh chữ nhật dài vừa 4 cột)
ggsave("Figure_6B_Absolute_Depth.png", plot = fig6b, width = 12, height = 6, dpi = 300, bg = "white")