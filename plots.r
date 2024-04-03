dev.off()
rm(list = ls())

serinc5 <- "ENSG00000164300"

if (rstudioapi::isAvailable()) {
  if (require("rstudioapi") != TRUE) {
    install.packages("rstudioapi")
  }else {
    library(rstudioapi)
  }
  wdir <- dirname(getActiveDocumentContext()$path)
}
library("reporter")

# Setting Working Directory
setwd(paste0(wdir, "/Data"))

# Importing nef data
nef_ratio <- read.csv("infectivity_with_rpm.csv", header = TRUE)

# Plotting bar plot (Fig 1a)
color_vector <- c("red", "red", "red", "red", "red",
                  "red", "red", "red", "green", "green",
                  "green", "green", "green", "blue", "blue")


barplot(nef_ratio$infectivity_ratio,
        horiz = TRUE,
        col = color_vector,
        names.arg = nef_ratio$cell_type,
        xlab = "Nef+/Nef - infectivity ratio",
        ylab = "Cell Type",
        xlim = c(0, 40), las = 1, cex.axis = 0.75, cex.names = 0.6,
        main = "Bar Plot of cell lines according to their infectivity ratio")

# Plotting scatter plot of SERINC5 (Fig 1d)

rpm_serinc5 <- nef_ratio[, paste0("RPM_", serinc5)]

corr_spearman <- cor.test(nef_ratio$infectivity_ratio, rpm_serinc5 , method = "spearman")

corr <- cor.test(nef_ratio$infectivity_ratio, rpm_serinc5 , method = "pearson")

linear_model <- lm(rpm_serinc5 ~ nef_ratio$infectivity_ratio)
summary_linear_model <- summary(linear_model)

r <- corr$estimate
p_value <- corr$p.value
r_squared <- summary_linear_model$r.squared


plot(nef_ratio$infectivity_ratio, rpm_serinc5, pch = 19, col = color_vector, xlab = "infectivity ratio", ylab = "RPM")
abline(linear_model)

text(4, max(rpm_serinc5),
     paste("r", ":", round(r, 3)), pos = 2)
text(5, max(rpm_serinc5 - 20),
     paste("R" %p% supsc("2"), ":", round(r_squared, 4)), pos = 2)
text(4.5, max(rpm_serinc5 - 40),
     "P < 0.001", pos = 2)

# Hypothesis testing

shapiro.test(rpm_serinc5)
shapiro.test(nef_ratio$infectivity_ratio)
# Create QQ plot
qqnorm(rpm_serinc5, main = "QQ Plot")  # QQ plot with data points
qqline(rpm_serinc5, col = 2)  # Add a reference line for a theoretical normal distribution

# Add labels and legend
legend("topleft", legend = "Data", pch = 1, col = 1)
legend("topright", legend = "Theoretical Normal", lty = 1, col = 2)


proteins_below_threshold <- read.csv(file = "proteins_below_threshold.csv",
                                     header = TRUE)


plot_correlation <- function(number, df, method) {
  iterations <- 0
  i <- 1
  while (iterations < number) {
    protein_id <- df$id[i]
    protein_name <- df$hgnc_symbol[i]
    rpm <- nef_ratio[, paste0("RPM_", protein_id)]

    shapiro <- shapiro.test(rpm)

    method_here <- ifelse(shapiro$p.value < 0.05, "spearman", "pearson")

    if (method_here != method) {
      i <- i + 1
      next
    }

    corr <- cor.test(nef_ratio$infectivity_ratio, rpm, method = method)
    linear_model <- lm(rpm ~ nef_ratio$infectivity_ratio)
    # summary_linear_model <- summary(linear_model)

    r <- corr$estimate
    p_value <- corr$p.value


    plot(nef_ratio$infectivity_ratio, rpm, pch = 19, col = color_vector,
         xlab = "infectivity ratio", ylab = "RPM",
         main = paste("Protein:", protein_name, "Method: ", method,
                      "P:", round(p_value, 4),
                      "r:", r))

    abline(linear_model)
    i <- i + 1
    iterations <- iterations + 1
  }
}

# Plot Positive and Negative Correlation top 5
proteins_num <- 5

num_rows <- ceiling(sqrt(proteins_num))  # Round up to the nearest integer
num_cols <- ceiling(proteins_num / num_rows)

neg_corr_df <- proteins_below_threshold[proteins_below_threshold$corr < 0, ]
neg_corr_df_spearman <- proteins_below_threshold[proteins_below_threshold$corr_spearman < 0, ]
neg_corr_df_spearman <- neg_corr_df_spearman[order(neg_corr_df_spearman$corr_spearman), ]

proteins_below_threshold <- proteins_below_threshold[order(proteins_below_threshold$corr, decreasing = TRUE), ]

par(mfrow = c(num_rows, num_cols))
plot_correlation(proteins_num, proteins_below_threshold, "pearson")
mtext("Top 5 Positive Correlation with normally distributed data",
      side = 3, line = -1.5, outer = TRUE, cex = 1.5, font = 2)

par(mfrow = c(num_rows, num_cols))
plot_correlation(proteins_num, neg_corr_df, "pearson")
mtext("Top 5 Negative Correlation with normally distributed data",
      side = 3, line = -1.5, outer = TRUE, cex = 1.5, font = 2)

proteins_below_threshold <- proteins_below_threshold[order(proteins_below_threshold$p_value_spearman), ]
par(mfrow = c(num_rows, num_cols))
plot_correlation(proteins_num, proteins_below_threshold, "spearman")
mtext("Top 5 Positive Correlation with non-normally distributed data",
      side = 3, line = -1.5, outer = TRUE, cex = 1.5, font = 2)


par(mfrow = c(num_rows, num_cols))
plot_correlation(proteins_num, neg_corr_df_spearman, "spearman")
mtext("Top 5 - Negative Correlation with non-normally distributed data",
      side = 3, line = -1.5, outer = TRUE, cex = 1.5, font = 2)

# Plot PCA
columns_to_remove <- c("filename", "cell_type")
PCA_filter_df <- nef_ratio[, !(names(nef_ratio) %in% columns_to_remove)]

data_matrix <- as.matrix(PCA_filter_df)  # Convert dataframe to matrix

pca_result <- prcomp(data_matrix, scale. = TRUE)

# View PCA results
summary(pca_result)  # Summary of the PCA result was applied)

# Calculate the percentage of variance explained by each component
variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100

# Create a bar plot for the scree plot
par(mfrow = c(1,1))
barplot(variance_explained,
        names.arg = paste0("PC", seq(1, length(variance_explained))),
        xlab = "Principal Component",
        ylab = "Percentage of Variance Explained",
        ylim = c(1, 100),
        main = "Scree Plot with Percentage of Variance Explained",
        col = "skyblue")
axis(2, at = seq(0, 100, by = 10), labels = seq(0, 100, by = 10))
