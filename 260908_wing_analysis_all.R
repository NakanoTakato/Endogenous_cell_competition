# 全体を上から実行してください。計算 → グラフ → 統計解析の順に進みます。
# ============================================================
# 1. 設定（通常はここだけ変更）
# ============================================================
Control <- "nubGFP"
Sample1 <- "nubp35"
Sample2 <- "nubmiRHG"
CommonName <- "260908 adult nub p35 test"

DataPath <- "/Volumes/IODATASSD/2024 Igaki lab/Data/2609/260902 nub p35 miRHG egg30 wing size"
InputFile <- "results.csv"
SavePath <- file.path(DataPath, "R")
Width <- 6

# 描画する群。Sample2も表示する場合は c(Control, Sample1, Sample2)
plot_limits <- c(Control, Sample1)
genotype_colors <- setNames(c("black", "magenta", "magenta"),
                            c(Control, Sample1, Sample2))
# 表示名は必要に応じて変更
ControlLabel <- expression(italic("nub > GFP"))
Sample1Label <- expression(italic("nub > GFP, p35"))
Sample2Label <- expression(italic("nub > GFP, miRHG"))
genotype_labels <- setNames(c(ControlLabel, Sample1Label, Sample2Label),
                            c(Control, Sample1, Sample2))

# LRペアのキーは元コードと同じ。日ごとにsample番号が重複するなら
# c("date", "Genotype", "sample") に変更してください。
PairKeys <- c("Genotype", "sample")

# ============================================================
# 2. パッケージ・データ読み込み
# ============================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbreak)
library(ggbeeswarm)
library(gghalves)
library(exactRankTests)


raw_dataset <- read.csv(file.path(DataPath, InputFile))

required_columns <- unique(c("Genotype", "date", PairKeys, "LR", "whole",
  "anterior", "posterior", "CellSizeW", "CellSizeA", "CellSizeP",
  "TotalCellW", "TotalCellA", "TotalCellP"))

missing_columns <- setdiff(required_columns, names(raw_dataset))

if (length(missing_columns)) stop("不足している列: ", paste(missing_columns, collapse = ", "))

print(table(raw_dataset$Genotype))

# ============================================================
# 3. wing size：日ごとのControl平均を求め、各測定値を割る
# ============================================================
day_means_wing_size <- raw_dataset %>%
  filter(Genotype == Control) %>%
  group_by(date) %>%
  summarise(
    mean_W = mean(whole, na.rm = TRUE),
    mean_A = mean(anterior, na.rm = TRUE),
    mean_P = mean(posterior, na.rm = TRUE),
    .groups = "drop"
  )

wing_size_df <- raw_dataset %>%
  left_join(day_means_wing_size, by = "date") %>%
  mutate(
    RelativeW = whole / mean_W,
    RelativeA = anterior / mean_A,
    RelativeP = posterior / mean_P
  )
print(day_means_wing_size)

# ============================================================
# 4. cell size：日ごとのControl平均を求め、各測定値を割る
# ============================================================
day_means_cell_size <- raw_dataset %>%
  filter(Genotype == Control) %>%
  group_by(date) %>%
  summarise(
    mean_W = mean(CellSizeW, na.rm = TRUE),
    mean_A = mean(CellSizeA, na.rm = TRUE),
    mean_P = mean(CellSizeP, na.rm = TRUE),
    .groups = "drop"
  )

cell_size_df <- raw_dataset %>%
  left_join(day_means_cell_size, by = "date") %>%
  mutate(
    RelativeW = CellSizeW / mean_W,
    RelativeA = CellSizeA / mean_A,
    RelativeP = CellSizeP / mean_P
  )
print(day_means_cell_size)

# ============================================================
# 5. cell number：日ごとのControl平均を求め、各測定値を割る
# ============================================================
day_means_cell_number <- raw_dataset %>%
  filter(Genotype == Control) %>%
  group_by(date) %>%
  summarise(
    mean_W = mean(TotalCellW, na.rm = TRUE),
    mean_A = mean(TotalCellA, na.rm = TRUE),
    mean_P = mean(TotalCellP, na.rm = TRUE),
    .groups = "drop"
  )

cell_number_df <- raw_dataset %>%
  left_join(day_means_cell_number, by = "date") %>%
  mutate(
    RelativeW = TotalCellW / mean_W,
    RelativeA = TotalCellA / mean_A,
    RelativeP = TotalCellP / mean_P
  )
print(day_means_cell_number)

# ============================================================
# 6. AP ratio：posterior / anterior を日ごとのControl平均で正規化
# ============================================================
AP_df <- raw_dataset %>% mutate(ratio = posterior / anterior)
day_means_AP <- AP_df %>%
  filter(Genotype == Control) %>%
  group_by(date) %>%
  summarise(mean_ratio = mean(ratio, na.rm = TRUE), .groups = "drop")
AP_df <- AP_df %>%
  left_join(day_means_AP, by = "date") %>%
  mutate(Relativeratio = ratio / mean_ratio)
print(day_means_AP)
# 元コードと同じく、描画にはRelativeratio、検定にはratioを使用します。

# ============================================================
# 7. LR difference：左右両方がある個体だけを残し、左右の絶対差を計算
# ============================================================
paired_data <- raw_dataset %>%
  group_by(across(all_of(PairKeys))) %>%
  filter(all(c("Left", "Right") %in% LR)) %>%
  ungroup() %>%
  filter(LR %in% c("Left", "Right"))
if (anyDuplicated(paired_data[c(PairKeys, "LR")])) {
  stop("LRペアが重複しています。PairKeysにdateなど個体を区別する列を追加してください。")
}
if (nrow(paired_data) == 0) stop("左右がそろった個体がありません。")
LR_df <- paired_data %>%
  pivot_wider(
    id_cols = all_of(PairKeys),
    names_from = LR,
    values_from = c(whole, anterior, posterior)
  ) %>%
  mutate(
    LR_w = abs(whole_Left - whole_Right),
    LR_a = abs(anterior_Left - anterior_Right),
    LR_p = abs(posterior_Left - posterior_Right)
  )

# ============================================================
# 8. 共通のggplot設定（すべてのグラフがこの設定を使用）
# ============================================================
# data: 描画するデータ、y_column: 縦軸の列名、y_label: 縦軸の表示名
# file_name: 保存名、relative: 相対値ならTRUE、jitter: 点の追加表示
save_plot <- function(data, y_column, y_label, file_name,
                      relative = TRUE, height = 10, jitter = TRUE) {
  p <- ggplot(data, aes(x = Genotype, y = .data[[y_column]], colour = Genotype)) +
    stat_summary(fun = "mean", geom = "crossbar", width = 0.3) +
    geom_half_boxplot(errorbar.length = 0.125, lwd = 1, width = 0.85,
      nudge = 0.2, outlier.shape = NA, outlier.color = NA, alpha = 0.6) +
    geom_beeswarm(cex = 1, orientation = "x", size = 1, alpha = 0.3)
  if (jitter) p <- p + geom_jitter(size = 1, alpha = 0.3, width = 0.1)
  p <- p + scale_colour_manual(values = genotype_colors) +
    labs(y = y_label, x = "") +
    scale_x_discrete(labels = genotype_labels, limits = plot_limits) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 15, angle = 60, hjust = 1.1),
      axis.title.x = element_text(size = 0),
      axis.title.y = element_text(size = 15, hjust = 0.7),
      axis.text.y.right = element_blank(), axis.ticks.y.right = element_blank(),
      legend.position = "none", axis.text = element_text(size = 15))
  if (relative) {
    p <- p + scale_y_continuous(limits = c(0, 1.2),
      breaks = c(0, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6),
      labels = c("0", "0.4", "0.6", "0.8", "1.0", "1.2", "1.4", "1.6")) +
      scale_y_break(c(0.1, 0.6), scales = 10)
  } else {
    p <- p + scale_y_continuous(limits = c(0, 0.1))
  }
  print(p)
  ggsave(filename = paste0(file_name, ".png"), plot = p, path = SavePath,
         width = Width, height = height, units = "cm", dpi = 300)
  invisible(p)
}

# 各グラフは「データ・列・縦軸名・保存名」だけを指定します。
# wing size
file_name <- paste(CommonName, "wing size")
plot_wing_size_W <- save_plot(wing_size_df, "RelativeW",
  "Relative wing size", paste0(file_name, "_W"))
plot_wing_size_A <- save_plot(wing_size_df, "RelativeA",
  "Relative anterior wing size", paste0(file_name, "_A"))
plot_wing_size_P <- save_plot(wing_size_df, "RelativeP",
  "Relative posterior wing size", paste0(file_name, "_P"))

# cell size
file_name <- paste(CommonName, "cell size")
plot_cell_size_W <- save_plot(cell_size_df, "RelativeW",
  "Relative wing cell size", paste0(file_name, "_W"), jitter = FALSE)
plot_cell_size_A <- save_plot(cell_size_df, "RelativeA",
  "Relative anterior wing cell size", paste0(file_name, "_A"))
plot_cell_size_P <- save_plot(cell_size_df, "RelativeP",
  "Relative posterior wing cell size", paste0(file_name, "_P"))

# cell number
file_name <- paste(CommonName, "cell number")
plot_cell_number_W <- save_plot(cell_number_df, "RelativeW",
  "Relative wing cell number", paste0(file_name, "_W"))
plot_cell_number_A <- save_plot(cell_number_df, "RelativeA",
  "Relative anterior wing cell number", paste0(file_name, "_A"))
plot_cell_number_P <- save_plot(cell_number_df, "RelativeP",
  "Relative posterior wing cell number", paste0(file_name, "_P"))

# AP ratio
file_name <- paste(CommonName, "AP ratio")
plot_AP <- save_plot(AP_df, "Relativeratio", "P/A wing size ratio", file_name)

# LR difference
file_name <- paste(CommonName, "LR difference")
plot_LR_W <- save_plot(LR_df, "LR_w",
  expression("abs. L-R difference (mm"^"2"~")"),
  paste0(file_name, "_W"), relative = FALSE, height = 12)
plot_LR_A <- save_plot(LR_df, "LR_a",
  expression("abs. anterior L-R difference (mm"^"2"~")"),
  paste0(file_name, "_A"), relative = FALSE, height = 12)
plot_LR_P <- save_plot(LR_df, "LR_p",
  expression("abs. posterior L-R difference (mm"^"2"~")"),
  paste0(file_name, "_P"), relative = FALSE, height = 12)

# ============================================================
# 9. 統計解析：各解析の比較をここにまとめて指定
# ============================================================
# 同じ検定処理だけを関数にまとめます。
# Controlと指定した群を、指定した列で比較（対応なしWilcoxon）。
# 欠損・非有限値を除外し、有効nとp値を返します。
# 詳細な検定オブジェクトはstat_testsに保存します。
stat_tests <- list()
compare_group <- function(data, column, sample_group, analysis, region) {
  x <- data[[column]][which(data$Genotype == Control)]
  y <- data[[column]][which(data$Genotype == sample_group)]
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  result <- if (!length(x) || !length(y)) {
    simpleError("比較する群の有効データがありません")
  } else {
    tryCatch(wilcox.exact(x = x, y = y, paired = FALSE), error = identity)
  }
  ok <- inherits(result, "htest")
  key <- paste(analysis, region, sample_group, sep = "__")
  if (ok) stat_tests[[key]] <<- result
  data.frame(analysis = analysis, region = region, value_column = column,
    control = Control, sample = sample_group,
    n_control = length(x), n_sample = length(y),
    statistic = if (ok) unname(result$statistic) else NA_real_,
    p_value = if (ok) result$p.value else NA_real_,
    status = if (ok) "OK" else conditionMessage(result))
}

# 元コードの比較群・検定対象を維持しています。
stat_results <- bind_rows(
  # wing size：Control vs Sample1 / Sample2
  compare_group(wing_size_df, "RelativeW", Sample1, "wing_size", "W"),
  compare_group(wing_size_df, "RelativeA", Sample1, "wing_size", "A"),
  compare_group(wing_size_df, "RelativeP", Sample1, "wing_size", "P"),
  # cell size：Control vs Sample1
  compare_group(cell_size_df, "RelativeW", Sample1, "cell_size", "W"),
  compare_group(cell_size_df, "RelativeA", Sample1, "cell_size", "A"),
  compare_group(cell_size_df, "RelativeP", Sample1, "cell_size", "P"),
  # cell number：Control vs Sample1
  compare_group(cell_number_df, "RelativeW", Sample1, "cell_number", "W"),
  compare_group(cell_number_df, "RelativeA", Sample1, "cell_number", "A"),
  compare_group(cell_number_df, "RelativeP", Sample1, "cell_number", "P"),
  # AP ratio：未正規化ratioを検定
  compare_group(AP_df, "ratio", Sample1, "AP", "P/A"),
  # LR difference：Control vs Sample1 / Sample2
  compare_group(LR_df, "LR_w", Sample1, "LR", "W"),
  compare_group(LR_df, "LR_a", Sample1, "LR", "A"),
  compare_group(LR_df, "LR_p", Sample1, "LR", "P")
)
print(stat_results, row.names = FALSE)
write.csv(stat_results, file.path(SavePath, paste(CommonName, "statistics.csv")),
          row.names = FALSE)
# 一覧を再表示：stat_results
# 検定の詳細：stat_tests[[paste("wing_size", "W", Sample1, sep = "__")]]
