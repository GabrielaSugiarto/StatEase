# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   https://r-pkgs.org
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#=============== MULTIVARIATE
#----------- Generate Data-----------
generate_dataset <- function(n_samples, n_features, multivariate_normal = TRUE, data_range = c(0, 1)) {
  if (multivariate_normal) {
    # Generate a random covariance matrix
    random_matrix <- matrix(runif(n_features^2, -1, 1), nrow = n_features)
    cov_matrix <- t(random_matrix) %*% random_matrix  # Ensure symmetric positive-definite matrix
    mean_vector <- runif(n_features, data_range[1], data_range[2])

    # Generate multivariate normal data
    library(MASS)
    data <- mvrnorm(n = n_samples, mu = mean_vector, Sigma = cov_matrix)
  } else {
    # Generate uniformly distributed data
    data <- matrix(runif(n_samples * n_features, data_range[1], data_range[2]),
                   nrow = n_samples, ncol = n_features)
  }

  # Convert to data frame
  colnames(data) <- paste0("Variable_", seq_len(n_features))
  return(as.data.frame(data))
}

#----------- Col Median ----------
colMedians <- function(data) {
  # Initialize a vector to store column medians
  manual_col_medians <- numeric(ncol(data))

  # Loop through each column
  for (i in seq_along(data)) {
    if (is.numeric(data[[i]])) {
      # Calculate the median for numeric columns (remove NA values)
      manual_col_medians[i] <- median(data[[i]], na.rm = TRUE)
    } else {
      # Assign NA for non-numeric columns
      manual_col_medians[i] <- NA
    }
  }

  # Assign column names to the result
  names(manual_col_medians) <- colnames(data)

  # Return the result
  return(manual_col_medians)
}

#----------- Col Modes --------------
colModes <- function(data) {
  # Helper function to calculate mode
  calculate_mode <- function(x) {
    # Remove NA values
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA)  # Handle empty case
    # Tabulate frequencies and return the most frequent value
    freq_table <- table(x)
    mode_value <- names(freq_table)[which.max(freq_table)]
    return(mode_value)
  }

  # Initialize a vector to store column modes
  manual_col_modes <- vector("list", ncol(data))  # Use a list for mixed types

  # Loop through each column
  for (i in seq_along(data)) {
    # Calculate mode for all column types
    manual_col_modes[[i]] <- calculate_mode(data[[i]])
  }

  # Assign column names to the result
  names(manual_col_modes) <- colnames(data)

  # Return the result
  return(unlist(manual_col_modes))  # Convert list to a named vector
}

#--------- critical value
# Fungsi untuk menghitung nilai kritis
hitung_nilai_kritis <- function(n, alpha) {
  # Tabel nilai kritis
  tabel_kritis <- data.frame(
    n = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 75, 100, 150, 200, 300),
    a0.01 = c(0.8299, 0.8801, 0.9126, 0.9269, 0.9410, 0.9479, 0.9538, 0.9599, 0.9632, 0.9671, 0.9695, 0.9720,
              0.9771, 0.9822, 0.9879, 0.9905, 0.9935),
    a0.05 = c(0.8788, 0.9198, 0.9389, 0.9508, 0.9591, 0.9652, 0.9682, 0.9726, 0.9749, 0.9768, 0.9787, 0.9801,
              0.9838, 0.9873, 0.9913, 0.9931, 0.9953),
    a0.1 = c(0.9032, 0.9351, 0.9503, 0.9604, 0.9665, 0.9715, 0.9740, 0.9771, 0.9792, 0.9809, 0.9822, 0.9836,
             0.9866, 0.9895, 0.9925, 0.9942, 0.9960)
  )

  # Buat nama kolom alpha
  kolom_alpha <- paste0("a", alpha)

  # Validasi alpha
  if (!kolom_alpha %in% colnames(tabel_kritis)) {
    stop(paste("Nilai alpha tidak valid. Gunakan 0.01, 0.05, atau 0.1. Nilai alpha saat ini:", alpha))
  }

  # Validasi n
  if (n < min(tabel_kritis$n) || n > max(tabel_kritis$n)) {
    stop(paste("Nilai n di luar jangkauan tabel kritis. Nilai n saat ini:", n))
  }

  # Cari batas bawah dan atas untuk interpolasi
  lower_row <- max(which(tabel_kritis$n <= n))
  upper_row <- min(which(tabel_kritis$n > n))

  # Ambil nilai untuk interpolasi
  x1 <- tabel_kritis$n[lower_row]
  x2 <- tabel_kritis$n[upper_row]
  y1 <- tabel_kritis[[kolom_alpha]][lower_row]
  y2 <- tabel_kritis[[kolom_alpha]][upper_row]

  # Lakukan interpolasi linier
  nilai_kritis <- y1 + ((n - x1) * (y2 - y1) / (x2 - x1))

  return(nilai_kritis)
}

#----------- Normality Test with Mahalanobis ------
mulmah <- function(data) {
  # Hanya gunakan kolom numerik
  dependen_vars <- data[sapply(data, is.numeric)]

  # Hapus baris dengan nilai NA
  dependen_vars <- na.omit(dependen_vars)

  # Hitung rata-rata kolom dan matriks kovarian
  means <- colMeans(dependen_vars)
  cov_matrix <- cov(dependen_vars)

  # Tambahkan kolom Mahalanobis Distance ke data asli
  data_clean <- na.omit(data)
  data_clean$MahalanobisDistance <- mahalanobis(dependen_vars, center = means, cov = cov_matrix)

  # Urutkan berdasarkan Mahalanobis Distance
  data_clean <- data_clean[order(data_clean$MahalanobisDistance), ]

  # Hitung nilai probabilitas dan kuantil chi-square
  n <- nrow(data_clean)
  df <- ncol(dependen_vars)  # Degrees of freedom
  data_clean$obs <- seq_len(n)
  data_clean$prob_value <- (data_clean$obs - 0.5) / n
  data_clean$chi_square <- qchisq(data_clean$prob_value, df = df)

  # --- Q-Q Plot ---
  qqplot(data_clean$MahalanobisDistance, data_clean$chi_square,
         main = "Q-Q Plot of Mahalanobis Distance vs Chi-square Distribution",
         xlab = "Sample Quantiles (Mahalanobis Distance)",
         ylab = "Theoretical Quantiles (Chi-square)")
  abline(0, 1, col = 'red', lwd = 2)  # Garis y = x
  grid()

  # --- Correlation Test ---
  correlation_test <- cor.test(data_clean$MahalanobisDistance, data_clean$chi_square, method = "pearson")
  r_squared <- correlation_test$estimate^2

  # --- Evaluasi Normalitas Berdasarkan Nilai Kritis ---
  hasil <- list()
  for (alpha in c(0.01, 0.05, 0.1)) {
    # Hitung nilai kritis berdasarkan fungsi hitung_nilai_kritis
    nilai_kritis <- hitung_nilai_kritis(n, alpha)

    kesimpulan <- if (r_squared < nilai_kritis) {
      "Tidak terdistribusi normal multivariat"
    } else {
      "Terdistribusi normal multivariat"
    }

    hasil[[paste0("α = ", alpha)]] <- list(
      nilai_kritis = round(nilai_kritis, 5),
      r_squared = round(r_squared, 3),
      kesimpulan = kesimpulan
    )
  }

  # Cetak hasil
  for (alpha in names(hasil)) {
    res <- hasil[[alpha]]
    cat("Signifikansi", alpha, "\n")
    cat("  Nilai kritis =", res$nilai_kritis, "\n")
    cat("  R² =", res$r_squared, "\n")
    cat("  Kesimpulan:", res$kesimpulan, "\n\n")
  }

  return(data_clean)
}


#----------- Detect Outlier ----------
detect_outliers <- function(data, alpha = 0.05) {
  # Cek apakah input adalah data frame
  if (!is.data.frame(data)) stop("Input harus berupa data frame.")

  # Hitung Mahalanobis Distance
  mean_vec <- colMeans(data, na.rm = TRUE)
  cov_matrix <- cov(data, use = "complete.obs")
  mahal_dist <- mahalanobis(data, mean_vec, cov_matrix)

  # Hitung threshold berdasarkan distribusi Chi-square
  threshold <- qchisq(1 - alpha, df = ncol(data))

  # Tandai apakah outlier
  is_outlier <- mahal_dist > threshold

  # Kembalikan hasil dalam data frame
  return(data.frame(Distance = mahal_dist, Outlier = is_outlier))
}
#----------- PCA Plot -----------
plot_pca_with_groups <- function(data, group_var, scale_data = TRUE) {
  # Validasi input
  if (!is.data.frame(data)) stop("Data harus berupa data frame.")
  if (is.null(group_var) || !group_var %in% colnames(data)) {
    stop("group_var harus diisi dan ada dalam data.")
  }

  # Ekstrak data numerik (tanpa kolom grup)
  numeric_data <- data[sapply(data, is.numeric)]

  # Standarisasi data jika diperlukan
  if (scale_data) numeric_data <- scale(numeric_data)

  # Lakukan PCA
  pca_result <- prcomp(numeric_data, center = TRUE, scale. = scale_data)

  # Ambil komponen utama
  pca_data <- as.data.frame(pca_result$x[, 1:2]) # PC1 dan PC2
  colnames(pca_data) <- c("PC1", "PC2")

  # Tambahkan kolom grup
  pca_data$Group <- as.factor(data[[group_var]])

  # Plot hasil PCA
  library(ggplot2)
  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(title = "PCA Scatterplot with Groups", x = "PC1", y = "PC2") +
    theme_minimal() +
    theme(legend.position = "right")

  # Tampilkan plot
  print(p)
}

#============= UNIVARIATE
#----------- confidence Interval--------------
confidence_interval_plot <- function(data, conf_level = 0.95) {
  # Memastikan data adalah vektor numerik
  if (!is.numeric(data)) {
    stop("Data harus berupa vektor numerik.")
  }

  # Menghitung rata-rata dan standar deviasi
  n <- length(data)
  mean_value <- mean(data)
  sd_value <- sd(data)

  # Menghitung error standar
  se <- sd_value / sqrt(n)

  # Menghitung batas bawah dan atas interval kepercayaan
  alpha <- 1 - conf_level
  critical_value <- qt(1 - alpha / 2, df = n - 1)
  margin_of_error <- critical_value * se

  lower_bound <- mean_value - margin_of_error
  upper_bound <- mean_value + margin_of_error

  # Menampilkan hasil
  cat("Rata-rata:", mean_value, "\n")
  cat("Interval kepercayaan (", conf_level * 100, "%): [", lower_bound, ", ", upper_bound, "]\n", sep = "")

  # Membuat plot
  plot(data, main = "Data Sampel dan Interval Kepercayaan",
       xlab = "Indeks", ylab = "Nilai", pch = 19, col = "blue")
  abline(h = mean_value, col = "red", lwd = 2) # Rata-rata
  segments(1, lower_bound, length(data), lower_bound, col = "green", lwd = 2) # Batas bawah
  segments(1, upper_bound, length(data), upper_bound, col = "green", lwd = 2) # Batas atas
  legend("topright", legend = c("Rata-rata", "Interval Kepercayaan"),
         col = c("red", "green"), lwd = 2)
}

#----------- Anova-----------
anova_sampling <- function(data, group_col, value_col, n_samples = 1000) {
  # Data split berdasarkan grup
  groups <- split(data[[value_col]], data[[group_col]])

  # Hitung F-statistik asli
  overall_mean <- mean(data[[value_col]])
  group_means <- sapply(groups, mean)
  n_group <- sapply(groups, length)

  ss_between <- sum(n_group * (group_means - overall_mean)^2)
  ss_within <- sum(sapply(groups, function(g) sum((g - mean(g))^2)))

  df_between <- length(groups) - 1
  df_within <- nrow(data) - length(groups)

  f_obs <- (ss_between / df_between) / (ss_within / df_within)

  # Sampling distribusi F-statistik
  f_samples <- numeric(n_samples)
  for (i in 1:n_samples) {
    permuted_data <- sample(data[[value_col]])
    permuted_groups <- split(permuted_data, data[[group_col]])
    permuted_means <- sapply(permuted_groups, mean)
    ss_between_perm <- sum(n_group * (permuted_means - overall_mean)^2)
    ss_within_perm <- sum(sapply(permuted_groups, function(g) sum((g - mean(g))^2)))
    f_samples[i] <- (ss_between_perm / df_between) / (ss_within_perm / df_within)
  }

  # Nilai p-value
  p_value <- mean(f_samples >= f_obs)

  # Visualisasi
  hist(
    f_samples, breaks = 30, probability = TRUE, col = "skyblue",
    main = "Distribusi Sampling F-statistik",
    xlab = "F-statistik"
  )
  abline(v = f_obs, col = "red", lwd = 2, lty = 2)
  legend("topright", legend = paste("F-observasi =", round(f_obs, 2)), col = "red", lwd = 2, lty = 2)

  # Interpretasi hasil
  if (p_value < 0.05) {
    interpretation <- paste(
      "F-observasi =", round(f_obs, 2),
      "terletak di area ekstrem distribusi sampling.\n",
      "Ini menunjukkan bahwa terdapat perbedaan signifikan antara kelompok (p-value =", round(p_value, 4), ")."
    )
  } else {
    interpretation <- paste(
      "F-observasi =", round(f_obs, 2),
      "tidak terletak di area ekstrem distribusi sampling.\n",
      "Ini menunjukkan bahwa tidak terdapat perbedaan signifikan antara kelompok (p-value =", round(p_value, 4), ")."
    )
  }

  cat("\nInterpretasi Hasil:\n", interpretation, "\n")

  # Output hasil
  list(
    f_obs = f_obs,
    p_value = p_value,
    f_distribution = f_samples,
    interpretation = interpretation
  )
}

#----------- Handle Missing Value -----------
handle_data_issues <- function(data, numeric_vars) {
  library(dplyr)
  library(ggplot2)

  # === Bagian 1: Menghitung dan Mengisi Missing Values ===
  cat("===== Menghitung Missing Values =====\n")

  # Hitung jumlah dan persentase missing value
  missing_count <- sapply(data, function(x) sum(is.na(x)))
  total_rows <- nrow(data)
  missing_percentage <- (missing_count / total_rows) * 100
  missing_summary <- data.frame(
    Var = names(missing_count),
    MissingCount = missing_count,
    MissingPercentage = missing_percentage
  )

  # Tampilkan hasil
  print(missing_summary)

  # Mengisi missing value dengan rata-rata (untuk variabel numerik)
  data <- data %>%
    mutate_if(is.numeric, ~ifelse(is.na(.), mean(., na.rm = TRUE), .))

  # Verifikasi setelah pengisian
  remaining_missing <- sapply(data, function(x) sum(is.na(x)))
  cat("\nJumlah missing value setelah diisi rata-rata:\n")
  print(remaining_missing)

  # === Bagian 2: Deteksi dan Visualisasi Outliers ===
  cat("\n===== Deteksi dan Visualisasi Outliers =====\n")

  # Standardisasi variabel numerik
  standardized_data <- data %>%
    mutate(across(all_of(numeric_vars), scale, .names = "z{col}"))

  # Scatterplot (menggunakan variabel pertama dan kedua, jika ada)
  if (length(numeric_vars) >= 2) {
    scatterplot <- ggplot(data, aes_string(x = numeric_vars[2], y = numeric_vars[1])) +
      geom_point() +
      labs(x = numeric_vars[2], y = numeric_vars[1]) +
      ggtitle(paste("Scatterplot", numeric_vars[1], "vs.", numeric_vars[2])) +
      theme_minimal()
    print(scatterplot)
  }

  # Boxplot untuk semua variabel numerik
  for (var in numeric_vars) {
    boxplot <- ggplot(data, aes_string(y = var)) +
      geom_boxplot(fill = "skyblue", color = "blue") +
      labs(y = var) +
      ggtitle(paste("Boxplot", var)) +
      theme_minimal()
    print(boxplot)
  }

  return(list(
    cleaned_data = data,
    standardized_data = standardized_data,
    missing_summary = missing_summary
  ))
}


