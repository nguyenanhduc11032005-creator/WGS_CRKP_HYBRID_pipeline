# Tải các thư viện cần thiết
library(ggplot2)
library(dplyr)

# 1. TẠO DATAFRAME
pcn_data <- data.frame(
  Pair = factor(c("Kpn-24\n(IncL)", "Kpn-17\n(IncFII)", "Kpn-72\n(IncR)", "Kpn-73\n(IncFII)"),
                levels = c("Kpn-24\n(IncL)", "Kpn-17\n(IncFII)", "Kpn-72\n(IncR)", "Kpn-73\n(IncFII)")),
  Native = c(1.0, 2.6, 1.2, 2.5),
  Mutant = c(3.2, 9.4, 6.1, 16.5),
  Mutation = c("No local mutation", "copA: C59A", "481-bp insertion", "copA: C58A")
)

# 2. VẼ BIỂU ĐỒ ARROW/DUMBBELL PLOT 
fig6 <- ggplot(pcn_data) +
  
  # Vẽ mũi tên kết nối từ Native lên Mutant (Thể hiện sức bật/Amplification)
  geom_segment(aes(x = Pair, xend = Pair, y = Native, yend = Mutant),
               arrow = arrow(length = unit(0.4, "cm"), type = "closed"),
               color = "grey60", linewidth = 1.2) +
  
  # Điểm cho Native (Hình tròn nhỏ hơn, màu xanh)
  geom_point(aes(x = Pair, y = Native, fill = "Native"),
             shape = 21, size = 5, color = "black", stroke = 1) +
  
  # Điểm cho Mutant (Hình tròn to hơn, màu đỏ gạch)
  geom_point(aes(x = Pair, y = Mutant, fill = "Mutant"),
             shape = 21, size = 6.5, color = "black", stroke = 1) +
  
  # Hiển thị số PCN của Native (nằm dưới điểm xanh)
  geom_text(aes(x = Pair, y = Native - 0.8, label = sprintf("%.1fX", Native)), 
            size = 4, fontface = "bold", color = "#2874A6") +
  
  # Hiển thị số PCN của Mutant (nằm trên điểm đỏ)
  geom_text(aes(x = Pair, y = Mutant + 1.2, label = sprintf("%.1fX", Mutant)), 
            size = 4, fontface = "bold", color = "#B03A2E") +
  
  # Gắn nhãn Đột biến (Nằm ngay giữa thân mũi tên)
  geom_label(aes(x = Pair, y = Native + (Mutant - Native)/2, label = Mutation),
             fill = "#E8F8F5", color = "#0E6251", fontface = "bold.italic", size = 4,
             label.padding = unit(0.3, "lines"), label.r = unit(0.15, "lines")) +
  
  # Tùy chỉnh màu
  scale_fill_manual(values = c("Native" = "#85C1E9", "Mutant" = "#E74C3C"), 
                    name = "Isolate Status") +
  
  # Tùy chỉnh trục Y
  scale_y_continuous(limits = c(0, 19), breaks = seq(0, 18, by = 3)) +
  
  # Tên biểu đồ và các trục
  labs(
    title = "Figure 6A: Mutation-Driven PCN Amplification",
    subtitle = "Fold-change shift from native clinical isolates to meropenem-exposed mutants",
    x = "Isolate Pair (Target Plasmid Replicon)",
    y = "Plasmid Copy Number (Relative to Chromosome)"
  ) +
  
  # Tùy chỉnh Theme gọn gàng, tinh tế cho bài báo khoa học
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30", margin = margin(b = 20)),
    axis.text.x = element_text(face = "bold", color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title.x = element_text(margin = margin(t = 15), face = "bold"),
    axis.title.y = element_text(margin = margin(r = 15), face = "bold"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(), # Bỏ đường kẻ dọc để biểu đồ thoáng hơn
    panel.grid.minor.y = element_blank()
  )

# Hiển thị biểu đồ trong console (không bắt buộc nhưng tốt cho log)
print(fig6)

# 3. LƯU THẲNG RA FILE PNG (Bật ggsave)
# bg = "white" đảm bảo nền ảnh màu trắng thay vì trong suốt (khi chèn vào Word/PDF sẽ không bị lỗi)
ggsave("Figure_6_PCN_Mutations.png", plot = fig6, width = 9, height = 7, dpi = 300, bg = "white")