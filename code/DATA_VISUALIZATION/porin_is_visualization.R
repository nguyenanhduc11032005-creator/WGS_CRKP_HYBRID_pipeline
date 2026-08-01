library(ggplot2)
library(dplyr)
library(tidyr)

# 1. NHẬP DỮ LIỆU TRỰC TIẾP TỪ VĂN BẢN (Copy-paste toàn bộ bảng của bạn)
raw_data <- "Sample  Group   Status  Chromosome_Depth    Plasmids_PCN    ompK35_Status   ompK36_Status   IS_Interrupting_Porin   Time
Kpn-12  Clinical    SUCCESS 139.9   2:0.9X;3:1.1X;4:1.5X;5:2.3X;6:0.9X;7:2.7X;8:3.3X    Missing_or_Truncated    Intact  None    800s
Kpn-24  Clinical    SUCCESS 132.3   2:1.1X;3:1.1X;4:0.9X;5:1.6X;6:1.0X;7:0.9X;8:1.1X;9:3.5X;10:3.7X;11:2.7X Intact  Missing_or_Truncated    None    803s
Kpn-15  Clinical    SUCCESS 124.5   2:1.0X;3:1.2X;4:1.6X;5:1.0X;6:0.7X;7:4.0X;8:0.9X;9:0.9X;10:1.2X;11:2.0X;12:4.8X;13:0.9X;14:1.0X;15:13.8X;16:32.2X;17:0.9X;18:0.8X   Intact  Missing_or_Truncated    None    816s
Kpn-26  Clinical    SUCCESS 131.4   2:1.1X;3:0.7X;4:1.6X    Missing_or_Truncated    Intact  None    817s
Kpn-17  Clinical    SUCCESS 114.1   2:0.8X;3:1.2X;4:2.6X;5:2.3X;6:9.5X;7:17.3X;8:21.1X  Missing_or_Truncated    Missing_or_Truncated    None    828s
Kpn-29  Clinical    SUCCESS 131.2   2:0.9X;3:1.8X;4:1.3X;5:3.1X;6:0.9X;7:3.1X;8:3.0X;9:4.4X;10:4.6X;11:5.5X;12:1.1X Intact  Missing_or_Truncated    None    795s
Kpn-28  Clinical    SUCCESS 157.6   2:5.3X;3:5.0X;4:15.2X;5:18.5X;6:9.5X;7:25.6X;8:25.2X;9:22.6X    Intact  Missing_or_Truncated    None    809s
Kpn-3   Clinical    SUCCESS 136.1   2:0.9X;3:1.6X;4:1.4X;5:1.6X;6:2.0X;7:0.9X;8:6.1X;9:6.2X Missing_or_Truncated    Intact  None    803s
Kpn-35  Clinical    SUCCESS 278.9   2:0.9X;3:0.8X;4:1.0X;5:2.5X;6:0.9X;7:3.6X;8:4.1X;9:0.9X Missing_or_Truncated    Intact  None    811s
Kpn-59  Clinical    SUCCESS 238.5   2:0.9X;3:1.1X;4:0.8X;5:2.5X;6:1.0X;7:0.9X;8:3.7X;9:3.9X Missing_or_Truncated    Intact  None    817s
Kpn-62  Clinical    SUCCESS 193.8   2:1.1X;3:4.2X;4:3.9X;5:6.6X;6:0.6X;7:1.1X;8:7.2X;9:19.8X;10:21.3X   Missing_or_Truncated    Missing_or_Truncated    None    817s
Kpn-7   Clinical    SUCCESS 120.2   2:1.1X;3:4.4X;4:5.4X;5:6.0X;6:0.7X;7:1.1X;8:19.0X;9:20.0X;10:32.4X  Missing_or_Truncated    Missing_or_Truncated    None    812s
Kpn-9   Clinical    SUCCESS 93.9    2:1.3X;3:1.0X;4:2.4X;5:1.2X;6:4.5X;7:4.4X;8:0.8X    Missing_or_Truncated    Intact  None    809s
Kpn-73  Clinical    SUCCESS 158.2   2:1.0X;3:1.2X;4:2.5X;5:2.0X;6:0.9X;7:0.9X;8:9.7X;9:16.5X;10:24.3X   Intact  Missing_or_Truncated    None    828s
Kpn-72  Clinical    SUCCESS 303.5   2:0.7X;3:1.3X;4:0.8X;5:0.8X;6:0.8X;7:1.9X;8:1.2X;9:2.4X;10:0.7X;11:3.6X;12:0.8X;13:0.7X;14:0.5X Missing_or_Truncated    Intact  None    841s
Kpn-24.4    Mutant  SUCCESS 134.9   2:1.4X;3:2.6X;4:3.3X;5:3.2X;6:2.5X;7:11.4X;8:14.2X;9:10.5X  Intact  Missing_or_Truncated    None    815s
Kpn-17.11   Mutant  SUCCESS 152.9   2:0.4X;3:1.2X;4:9.4X;5:0.4X;6:0.4X;7:2.3X;8:0.6X;9:0.4X;10:0.4X;11:0.4X;12:8.9X;13:0.4X;14:17.1X;15:0.7X;16:0.4X;17:0.4X;18:19.6X   Missing_or_Truncated    Missing_or_Truncated    None    837s
Kpn-72.11   Mutant  SUCCESS 135.7   2:1.9X;3:1.9X;4:6.1X;5:7.6X;6:5.9X  Missing_or_Truncated    Missing_or_Truncated    None    808s
Kpn-73.9    Mutant  SUCCESS 117.0   2:0.9X;3:1.7X;4:16.5X;5:4.0X;6:19.1X;7:41.7X;8:54.8X    Missing_or_Truncated    Missing_or_Truncated    None    838s"

df <- read.table(text = raw_data, header = TRUE, stringsAsFactors = FALSE)

# 2. XỬ LÝ DỮ LIỆU
porin_data <- df %>%
  select(Sample, Group, ompK35_Status, ompK36_Status) %>%
  # Xoay dữ liệu từ ngang sang dọc để vẽ biểu đồ ngói (Tile Plot)
  pivot_longer(cols = c(ompK35_Status, ompK36_Status), 
               names_to = "Porin_Gene", 
               values_to = "Status") %>%
  # Sửa tên hiển thị cho đẹp
  mutate(Porin_Gene = ifelse(Porin_Gene == "ompK35_Status", "ompK35", "ompK36"),
         Status = gsub("_", " ", Status))

# Sắp xếp thứ tự Sample để các Mutant nằm cạnh nhau ở cuối
porin_data$Sample <- factor(porin_data$Sample, levels = df$Sample)

# 3. VẼ BIỂU ĐỒ HEATMAP PORIN
fig5 <- ggplot(porin_data, aes(x = Sample, y = Porin_Gene, fill = Status)) +
  geom_tile(color = "white", linewidth = 1.5) +
  
  # Tạo lưới chia theo Clinical và Mutant để dễ quan sát
  facet_grid(~Group, scales = "free_x", space = "free_x") +
  
  # Tùy chỉnh màu sắc (Xanh = Bình thường, Cam = Đột biến/Mất)
  scale_fill_manual(values = c("Intact" = "#7DCEA0", "Missing or Truncated" = "#F1948A")) +
  
  labs(
    title = "Figure 5. Genomic Status of Major Porins (ompK35/ompK36)",
    subtitle = "In silico prediction of structural integrity across clinical and mutant isolates",
    x = "Klebsiella pneumoniae Isolate",
    y = "Porin Gene",
    fill = "Genomic Status"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30", margin = margin(b=15)),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "italic", size = 14, color = "black"),
    strip.text = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "grey90", color = NA),
    legend.position = "bottom"
  )

print(fig5)
ggsave("Figure_5_Porin_Status.png", plot = fig5, width = 11, height = 5, dpi = 300, bg="white")