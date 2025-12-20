# Filter Laplacian with N=128, aux=0.5, varying tol
filter_laplacian_n128_aux05_varying_tol <- function(df) {
  df %>% filter(matrix_type == "laplacian", N == 128, aux == 0.5)
}
# Plot functions for p-varying (aux) Laplacian/ER plots
plot_time_comparison_p <- function(df) {
  long <- df %>%
    dplyr::select(aux, matrix_type, any_of("iteration"), time.rbd, time.brute) %>%
    pivot_longer(c(time.rbd, time.brute), names_to = "method", values_to = "time_seconds") %>%
    mutate(method = recode(method, "time.rbd" = "RBD", "time.brute" = "Standard"))

  ggplot(long, aes(x = factor(aux), y = time_seconds, fill = method)) +
    geom_boxplot(outlier.alpha = 0.2, width = 0.6, linewidth = 0.25,
                 position = position_dodge2(width = 0.7, preserve = "single")) +
    stat_summary(aes(group = method), fun = median, geom = "point", shape = 21, size = 1.8,
                 position = position_dodge(width = 0.7), color = "black", fill = "white") +
    scale_fill_manual(values = c("RBD" = "#2E86AB", "Standard" = "#A23B72")) +
    scale_y_log10(labels = label_number()) +
    labs(x = "Edge Probability (p)", y = "Time (seconds)", fill = "Method") +
    theme_minimal(base_size = 12) +
    theme(axis.title = element_text(size = 14.4), legend.position = "top")
}

plot_speedup_p <- function(df) {
  ggplot(df, aes(x = factor(aux), y = speedup)) +
    geom_boxplot(fill = "#06A77D", alpha = 0.7, outlier.alpha = 0.3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    labs(x = "Edge Probability (p)", y = "Standard Time Taken / RBD Time Taken") +
    scale_y_log10() +
    theme_minimal(base_size = 12) +
    theme(axis.title = element_text(size = 14.4))
}

plot_kl_p <- function(df) {
  ggplot(df, aes(x = factor(aux), y = kl_abs)) +
    geom_boxplot(fill = "#F18F01", alpha = 0.7, outlier.alpha = 0.3) +
    labs(x = "Edge Probability (p)", y = "|KL|") +
    scale_y_log10() +
    theme_minimal(base_size = 12) +
    theme(axis.title = element_text(size = 14.4))
}
# Utility plotting helpers for simulation study results
# Load results from data/combined_simulation_study_results.csv and provide
# functions to visualize timing and KL divergence comparisons between RBD and brute.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

results_path_default <- "data/combined_simulation_study_results.csv"

load_results <- function(path = results_path_default) {
  if (!file.exists(path)) {
    stop(sprintf("Results file not found at %s", path))
  }
  df <- read.csv(path)
  df %>%
    mutate(
      speedup = time.brute / time.rbd,
      time_diff_sec = time.brute - time.rbd,
      kl_abs = abs(kl_divergence)
    )
}

# Generic helper: filter by matrix_type and index IDs
filter_by_type_ids <- function(df, matrix_type_value, ids) {
  df %>% filter(matrix_type == matrix_type_value, index %in% ids)
}

# Convenience wrappers for specific studies
filter_toeplitz <- function(df, ids = 511:515) {
  filter_by_type_ids(df, "toeplitz", ids)
}

filter_inverse_wishart <- function(df, ids = 481:485) {
  filter_by_type_ids(df, "inverse_wishart", ids)
}

filter_er_gnp <- function(df, ids = 121:125) {
  filter_by_type_ids(df, "gnp", ids)
}

filter_laplacian <- function(df, ids = 391:395) {
  filter_by_type_ids(df, "laplacian", ids)
}

# Filter Laplacian with tol=1e-06, aux=0.5 across all N
filter_laplacian_tol1e6_aux05 <- function(df) {
  df %>% filter(matrix_type == "laplacian", rbd_tol == 1e-06, aux == 0.5)
}

# Filter Toeplitz with tol=1e-06, aux=0.5 across all N
filter_toeplitz_tol1e6_aux05 <- function(df) {
  df %>% filter(matrix_type == "toeplitz", rbd_tol == 1e-06, is.na(aux) | aux == 0.5)
}

# Filter ER (gnp) with N=64, tol=1e-06, varying aux (edge probability p)
filter_er_n64_varying_p <- function(df) {
  df %>% filter(matrix_type == "gnp", N == 64, rbd_tol == 1e-06)
}

# Filter ER (gnp) with N specified, tol=1e-06, varying aux (edge probability p)
filter_er_n_varying_p <- function(df, n_value) {
  df %>% filter(matrix_type == "gnp", N == n_value, rbd_tol == 1e-06)
}

# Filter ER (gnp) with N=64, aux=0.5, varying tol
filter_er_n64_aux05_varying_tol <- function(df) {
  df %>% filter(matrix_type == "gnp", N == 64, aux == 0.5)
}

# Filter Laplacian with N=64, tol=1e-06, varying aux (edge probability p)
filter_laplacian_n64_varying_p <- function(df) {
  df %>% filter(matrix_type == "laplacian", N == 64, rbd_tol == 1e-06)
}

# Filter Laplacian with N specified, tol=1e-06, varying aux (edge probability p)
filter_laplacian_n_varying_p <- function(df, n_value) {
  df %>% filter(matrix_type == "laplacian", N == n_value, rbd_tol == 1e-06)
}

plot_time_comparison <- function(df) {
  long <- df %>%
    dplyr::select(N, matrix_type, any_of("iteration"), time.rbd, time.brute) %>%
    pivot_longer(c(time.rbd, time.brute), names_to = "method", values_to = "time_seconds") %>%
    mutate(method = recode(method, "time.rbd" = "RBD", "time.brute" = "Standard"))

  p <- ggplot(long, aes(x = factor(N), y = time_seconds, fill = method)) +
    geom_boxplot(outlier.alpha = 0.2, width = 0.6, linewidth = 0.25,
                 position = position_dodge2(width = 0.7, preserve = "single")) +
    stat_summary(aes(group = method), fun = median, geom = "point", shape = 21, size = 1.8,
                 position = position_dodge(width = 0.7), color = "black", fill = "white") +
    scale_fill_manual(values = c("RBD" = "#2E86AB", "Standard" = "#A23B72")) +
    scale_y_log10(labels = label_number()) +
    labs(x = "Number of Study Objects (N)", y = "Time (seconds)", fill = "Method") +
    theme_minimal(base_size = 12) +
    theme(axis.title = element_text(size = 14.4), legend.position = "top")

  if (length(unique(df$matrix_type)) > 1) p <- p + facet_wrap(~matrix_type)
  p
}

# Fit linear model to log(time) vs N for both RBD and Standard methods
fit_log_time_vs_log_N <- function(df) {
  long <- df %>%
    dplyr::select(N, time.rbd, time.brute) %>%
    pivot_longer(c(time.rbd, time.brute), names_to = "method", values_to = "time_seconds") %>%
    mutate(method = recode(method, "time.rbd" = "RBD", "time.brute" = "Standard"))

  # Remove non-positive times to avoid log(0) or log(negative)
  long <- long %>% filter(time_seconds > 0)

  models <- long %>%
    group_by(method) %>%
    do(model = lm(log(time_seconds) ~ log(N), data = .))

  # Print summary for each model
  for (i in seq_len(nrow(models))) {
    cat("\nLinear model for", models$method[i], ":\n")
    print(summary(models$model[[i]]))
  }

  invisible(models)
}

plot_speedup <- function(df) {
  p <- ggplot(df, aes(x = factor(N), y = speedup)) +
    geom_boxplot(fill = "#06A77D", alpha = 0.7, outlier.alpha = 0.3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    labs(x = "Number of Study Objects (N)", y = "Standard Time Taken / RBD Time Taken") +
    scale_y_log10() +
    theme_minimal(base_size = 12) +
    theme(axis.title = element_text(size = 14.4))

  if (length(unique(df$matrix_type)) > 1) p <- p + facet_wrap(~matrix_type)
  p
}

plot_kl <- function(df) {
  p <- ggplot(df, aes(x = factor(N), y = kl_abs)) +
    geom_boxplot(fill = "#F18F01", alpha = 0.7, outlier.alpha = 0.3) +
    labs(x = "Number of Study Objects (N)", y = "|KL|") +
    scale_y_log10() +
    theme_minimal(base_size = 12) +
    theme(axis.title = element_text(size = 14.4))

  if (length(unique(df$matrix_type)) > 1) p <- p + facet_wrap(~matrix_type)
  p
}

# Plot functions for tolerance sensitivity (x-axis = tolerance)
plot_time_comparison_tol <- function(df) {
  long <- df %>%
    select(rbd_tol, matrix_type, any_of("iteration"), time.rbd, time.brute) %>%
    pivot_longer(c(time.rbd, time.brute), names_to = "method", values_to = "time_seconds") %>%
    mutate(method = recode(method, "time.rbd" = "RBD", "time.brute" = "Standard"))

  p <- ggplot(long, aes(x = factor(rbd_tol), y = time_seconds, fill = method)) +
    geom_boxplot(outlier.alpha = 0.2, width = 0.6, linewidth = 0.25,
                 position = position_dodge2(width = 0.7, preserve = "single")) +
    stat_summary(aes(group = method), fun = median, geom = "point", shape = 21, size = 1.8,
                 position = position_dodge(width = 0.7), color = "black", fill = "white") +
    scale_fill_manual(values = c("RBD" = "#2E86AB", "Standard" = "#A23B72")) +
    scale_y_log10(labels = label_number()) +
    labs(x = "Tolerance", y = "Time (seconds)", fill = "Method") +
    theme_minimal(base_size = 12) +
    theme(axis.title = element_text(size = 14.4), legend.position = "top")

  if (length(unique(df$matrix_type)) > 1) p <- p + facet_wrap(~matrix_type)
  p
}

plot_speedup_tol <- function(df) {
  p <- ggplot(df, aes(x = factor(rbd_tol), y = speedup)) +
    geom_boxplot(fill = "#06A77D", alpha = 0.7, outlier.alpha = 0.3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    labs(x = "Tolerance", y = "Standard Time Taken / RBD Time Taken") +
    scale_y_log10() +
    theme_minimal(base_size = 12) +
    theme(axis.title = element_text(size = 14.4))

  if (length(unique(df$matrix_type)) > 1) p <- p + facet_wrap(~matrix_type)
  p
}

plot_kl_tol <- function(df) {
  p <- ggplot(df, aes(x = factor(rbd_tol), y = kl_abs)) +
    geom_boxplot(fill = "#F18F01", alpha = 0.7, outlier.alpha = 0.3) +
    labs(x = "Tolerance", y = "|KL|") +
    scale_y_log10() +
    theme_minimal(base_size = 12) +
    theme(axis.title = element_text(size = 14.4))

  if (length(unique(df$matrix_type)) > 1) p <- p + facet_wrap(~matrix_type)
  p
}

# Helper to save a plot
save_plot <- function(plot_obj, filename, width = 10, height = 6, dpi = 300) {
  dir.create("plots", showWarnings = FALSE)
  ggsave(file.path("plots", filename), plot_obj, width = width, height = height, dpi = dpi)
}

# Create summary table for ER with tol=1e-6, varying p across all N values
summarize_er_n_varying_p <- function(df) {
  summary_data <- df %>%
    filter(matrix_type == 'gnp', rbd_tol == 1e-06) %>%
    group_by(N, aux) %>%
    summarise(
      RBD_time_mean = round(mean(time.rbd, na.rm = TRUE), 4),
      Standard_time_mean = round(mean(time.brute, na.rm = TRUE), 4),
      Speedup_mean = round(mean(speedup, na.rm = TRUE), 2),
      KL_mean = formatC(mean(kl_abs, na.rm = TRUE), format = "e", digits = 2),
      .groups = 'drop'
    ) %>%
    rename('N' = N, 'p' = aux, 'RBD (s)' = RBD_time_mean, 'Standard (s)' = Standard_time_mean, 'Speedup' = Speedup_mean, '|KL|' = KL_mean) %>%
    arrange(N, p)
  
  summary_data
}

# Create summary table for Laplacian with tol=1e-6, varying p across all N values
summarize_laplacian_n_varying_p <- function(df) {
  summary_data <- df %>%
    filter(matrix_type == 'laplacian', rbd_tol == 1e-06) %>%
    group_by(N, aux) %>%
    summarise(
      RBD_time_mean = round(mean(time.rbd, na.rm = TRUE), 4),
      Standard_time_mean = round(mean(time.brute, na.rm = TRUE), 4),
      Speedup_mean = round(mean(speedup, na.rm = TRUE), 2),
      KL_mean = formatC(mean(kl_abs, na.rm = TRUE), format = "e", digits = 2),
      .groups = 'drop'
    ) %>%
    rename('N' = N, 'p' = aux, 'RBD (s)' = RBD_time_mean, 'Standard (s)' = Standard_time_mean, 'Speedup' = Speedup_mean, '|KL|' = KL_mean) %>%
    arrange(N, p)
  
  summary_data
}

# Create summary table for Toeplitz with tol=1e-6, aux=0.5 across all N values
summarize_toeplitz_tol1e6_aux05 <- function(df) {
  summary_data <- df %>%
    filter_toeplitz_tol1e6_aux05() %>%
    group_by(N) %>%
    summarise(
      RBD_min = round(min(time.rbd, na.rm = TRUE), 4),
      RBD_mean = round(mean(time.rbd, na.rm = TRUE), 4),
      RBD_median = round(median(time.rbd, na.rm = TRUE), 4),
      RBD_max = round(max(time.rbd, na.rm = TRUE), 4),
      Standard_min = round(min(time.brute, na.rm = TRUE), 4),
      Standard_mean = round(mean(time.brute, na.rm = TRUE), 4),
      Standard_median = round(median(time.brute, na.rm = TRUE), 4),
      Standard_max = round(max(time.brute, na.rm = TRUE), 4),
      .groups = 'drop'
    ) %>%
    rename(
      'N' = N,
      'RBD min (s)' = RBD_min,
      'RBD mean (s)' = RBD_mean,
      'RBD median (s)' = RBD_median,
      'RBD max (s)' = RBD_max,
      'Standard min (s)' = Standard_min,
      'Standard mean (s)' = Standard_mean,
      'Standard median (s)' = Standard_median,
      'Standard max (s)' = Standard_max
    ) %>%
    arrange(N)

  summary_data
}

# Create summary table for Inverse Wishart with tol=1e-6, aux=0.5 across all N values
summarize_inverse_wishart_tol1e6_aux05 <- function(df) {
  summary_data <- df %>%
    filter_inverse_wishart_tol1e6_aux05() %>%
    group_by(N) %>%
    summarise(
      RBD_min = round(min(time.rbd, na.rm = TRUE), 4),
      RBD_mean = round(mean(time.rbd, na.rm = TRUE), 4),
      RBD_median = round(median(time.rbd, na.rm = TRUE), 4),
      RBD_max = round(max(time.rbd, na.rm = TRUE), 4),
      Standard_min = round(min(time.brute, na.rm = TRUE), 4),
      Standard_mean = round(mean(time.brute, na.rm = TRUE), 4),
      Standard_median = round(median(time.brute, na.rm = TRUE), 4),
      Standard_max = round(max(time.brute, na.rm = TRUE), 4),
      .groups = 'drop'
    ) %>%
    rename(
      'N' = N,
      'RBD min (s)' = RBD_min,
      'RBD mean (s)' = RBD_mean,
      'RBD median (s)' = RBD_median,
      'RBD max (s)' = RBD_max,
      'Standard min (s)' = Standard_min,
      'Standard mean (s)' = Standard_mean,
      'Standard median (s)' = Standard_median,
      'Standard max (s)' = Standard_max
    ) %>%
    arrange(N)

  summary_data
}

# Create summary table for Laplacian with tol=1e-6, aux=0.5 across all N values
summarize_laplacian_tol1e6_aux05 <- function(df) {
  summary_data <- df %>%
    filter_laplacian_tol1e6_aux05() %>%
    group_by(N) %>%
    summarise(
      RBD_min = round(min(time.rbd, na.rm = TRUE), 4),
      RBD_mean = round(mean(time.rbd, na.rm = TRUE), 4),
      RBD_median = round(median(time.rbd, na.rm = TRUE), 4),
      RBD_max = round(max(time.rbd, na.rm = TRUE), 4),
      Standard_min = round(min(time.brute, na.rm = TRUE), 4),
      Standard_mean = round(mean(time.brute, na.rm = TRUE), 4),
      Standard_median = round(median(time.brute, na.rm = TRUE), 4),
      Standard_max = round(max(time.brute, na.rm = TRUE), 4),
      .groups = 'drop'
    ) %>%
    rename(
      'N' = N,
      'RBD min (s)' = RBD_min,
      'RBD mean (s)' = RBD_mean,
      'RBD median (s)' = RBD_median,
      'RBD max (s)' = RBD_max,
      'Standard min (s)' = Standard_min,
      'Standard mean (s)' = Standard_mean,
      'Standard median (s)' = Standard_median,
      'Standard max (s)' = Standard_max
    ) %>%
    arrange(N)

  summary_data
}

# Export table to LaTeX format
table_to_latex <- function(summary_table, caption = "", label = "") {
  # Check if xtable is available
  if (!requireNamespace("xtable", quietly = TRUE)) {
    message("Installing xtable package for LaTeX export...")
    install.packages("xtable", repos = "https://cloud.r-project.org/")
  }
  
  library(xtable)
  
  xtab <- xtable(summary_table, caption = caption, label = label)
  print(xtab, 
        include.rownames = FALSE,
        booktabs = TRUE,
        sanitize.text.function = identity)
}

# Example usage (uncomment to run interactively):
df <- load_results()
df <- filter_laplacian(df)
p1 <- plot_time_comparison(df); print(p1); save_plot(p1, "time_comparison.png")
p2 <- plot_speedup(df); print(p2); save_plot(p2, "speedup.png")
p3 <- plot_kl(df); print(p3); save_plot(p3, "kl_divergence.png")


# # Generate Laplacian plots for varying N and p:
df <- load_results()
summary_table_lap <- summarize_laplacian_n_varying_p(df)
print(summary_table_lap, n = Inf)
table_to_latex(summary_table_lap, caption = "Laplacian Matrix Performance", label = "tab:laplacian")
#
# # Generate Toeplitz table for tol=1e-6, aux=0.5 across N values:
# df <- load_results()
# summary_table_toep <- summarize_toeplitz_tol1e6_aux05(df)
# print(summary_table_toep, n = Inf)
# table_to_latex(summary_table_toep, caption = "Toeplitz Performance (tol=1e-6, aux=0.5)", label = "tab:toeplitz_tol1e6_aux05")
# 
# # Generate ER plots for varying tolerance (N=64, aux=0.5):
# df <- load_results()
# df_tol <- filter_er_n64_aux05_varying_tol(df)
# df_tol$rbd_tol <- factor(df_tol$rbd_tol, labels = paste0('tol=', sort(unique(df_tol$rbd_tol))))
# save_plot(plot_time_comparison(df_tol), 'er_n64_aux05_tol_time_comparison.png')
# save_plot(plot_speedup(df_tol), 'er_n64_aux05_tol_speedup.png')
# save_plot(plot_kl(df_tol), 'er_n64_aux05_tol_kl.png')



df<-load_results();
df_tol <- filter_laplacian_n128_aux05_varying_tol(df);

df<-load_results();
df <- filter_laplacian(df)
fit_log_time_vs_log_N(df)
