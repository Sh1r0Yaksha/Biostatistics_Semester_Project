dev.off()
rm(list = ls())



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

serinc5 <- "ENSG00000164300"

################ Anushka Khobragade 20049 ################

# Importing nef data
nef_ratio <- read.csv("infectivity_with_rpm.csv", header = TRUE)
View(nef_ratio)

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

################ Akshat Singh 20031 ################

# Hypothesis testing
corr_spearman <- cor.test(nef_ratio$infectivity_ratio, rpm_serinc5 , method = "spearman")
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

proteins_below_threshold <- proteins_below_threshold[order(proteins_below_threshold$corr, decreasing = TRUE), ]
View(proteins_below_threshold)

plot_correlation <- function(number, df, method) {
  iterations <- 0
  i <- 1
  while (iterations < number) {
    protein_id <- df$id[i]
    protein_name <- df$hgnc_symbol[i]
    rpm <- nef_ratio[[paste0("RPM_", protein_id)]]
    shapiro <- shapiro.test(rpm)
    method_here <- ifelse(shapiro$p.value < 0.05, "spearman", "pearson")

    if ((method == "pearson" && method_here != method) || protein_name == "" || is.na(protein_name)) {
      i <- i + 1
      next
    }

    corr <- cor.test(nef_ratio$infectivity_ratio, rpm, method = method, exact = FALSE)
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

pos_normally_df <- proteins_below_threshold[order(proteins_below_threshold$corr, decreasing = TRUE), ]
par(mfrow = c(num_rows, num_cols))
plot_correlation(proteins_num, pos_normally_df, "pearson")
mtext("Top 5 Positive Correlation with normally distributed data",
      side = 3, line = -1.5, outer = TRUE, cex = 1.5, font = 2)

neg_normally_df <- proteins_below_threshold[order(proteins_below_threshold$corr), ]
par(mfrow = c(num_rows, num_cols))
plot_correlation(proteins_num, neg_normally_df, "pearson")
mtext("Top 5 Negative Correlation with normally distributed data",
      side = 3, line = -1.5, outer = TRUE, cex = 1.5, font = 2)

pos_non_normally_spearman_df <- proteins_below_threshold[order(proteins_below_threshold$corr_spearman, decreasing = TRUE), ]
par(mfrow = c(num_rows, num_cols))
plot_correlation(proteins_num, pos_non_normally_spearman_df, "spearman")
mtext("Top 5 Positive Correlation with non-normally distributed data",
      side = 3, line = -1.5, outer = TRUE, cex = 1.5, font = 2)

neg_non_normally_spearman_df <- proteins_below_threshold[order(proteins_below_threshold$corr_spearman), ]
par(mfrow = c(num_rows, num_cols))
plot_correlation(proteins_num, neg_non_normally_spearman_df, "spearman")
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



# Plot of Proteins

################ Sohum Ranade 20270 ################

# HIRA
HIRA <- "ENSG00000100084"
rpm_HIRA <- nef_ratio[, paste0("RPM_", HIRA)]

corr_spearman <- cor.test(nef_ratio$infectivity_ratio, rpm_HIRA , method = "spearman")

linear_model <- lm(rpm_HIRA ~ nef_ratio$infectivity_ratio)
summary_linear_model <- summary(linear_model)

r <- corr_spearman$estimate
p_value <- corr_spearman$p.value
r_squared <- summary_linear_model$r.squared

plot(nef_ratio$infectivity_ratio,
     rpm_HIRA,
     pch = 19,
     col = color_vector,
     xlab = "infectivity ratio",
     ylab = "RPM",
     main = "HIRA")
abline(linear_model)

text(4, max(rpm_HIRA),
     paste("r", ":", round(r, 3)), pos = 2)
text(4, max(rpm_HIRA - 1),
     paste("R" %p% supsc("2"), ":", round(r_squared, 4)), pos = 2)
text(4, max(rpm_HIRA - 2),
     "P < 0.001", pos = 2)


################ Gunashree Rathi 20123 ################

# VMP1
VMP1 <- "ENSG00000062716"
rpm_VMP1 <- nef_ratio[, paste0("RPM_", VMP1)]

corr_spearman <- cor.test(nef_ratio$infectivity_ratio, rpm_VMP1 , method = "spearman")

linear_model <- lm(rpm_VMP1 ~ nef_ratio$infectivity_ratio)
summary_linear_model <- summary(linear_model)

r <- corr_spearman$estimate
p_value <- corr_spearman$p.value
r_squared <- summary_linear_model$r.squared

plot(nef_ratio$infectivity_ratio,
     rpm_VMP1,
     pch = 19,
     col = color_vector,
     xlab = "infectivity ratio",
     ylab = "RPM",
     main = "VMP1")
abline(linear_model)

text(4, max(rpm_VMP1),
     paste("r", ":", round(r, 3)), pos = 2)
text(4, max(rpm_VMP1 - 10),
     paste("R" %p% supsc("2"), ":", round(r_squared, 4)), pos = 2)
text(4, max(rpm_VMP1 - 20),
     "P < 0.001", pos = 2)

# ATG13
ATG13 <- "ENSG00000175224"
rpm_ATG13 <- nef_ratio[, paste0("RPM_", ATG13)]

corr <- cor.test(nef_ratio$infectivity_ratio, rpm_ATG13 , method = "pearson")

linear_model <- lm(rpm_ATG13 ~ nef_ratio$infectivity_ratio)
summary_linear_model <- summary(linear_model)

r <- corr$estimate
p_value <- corr$p.value
r_squared <- summary_linear_model$r.squared

plot(nef_ratio$infectivity_ratio,
     rpm_ATG13,
     pch = 19,
     col = color_vector,
     xlab = "infectivity ratio",
     ylab = "RPM",
     main = "ATG13")
abline(linear_model)

text(4, max(rpm_ATG13),
     paste("r", ":", round(r, 3)), pos = 2)
text(4, max(rpm_ATG13 - 3),
     paste("R" %p% supsc("2"), ":", round(r_squared, 4)), pos = 2)
text(4, max(rpm_ATG13 - 6),
     "P < 0.001", pos = 2)

################ Manas Kulkarni 20166 ################

# PPP2R1B
PPP2R1B <- "ENSG00000137713"
rpm_PPP2R1B <- nef_ratio[, paste0("RPM_", PPP2R1B)]

corr_spearman <- cor.test(nef_ratio$infectivity_ratio, rpm_PPP2R1B , method = "spearman")

linear_model <- lm(rpm_PPP2R1B ~ nef_ratio$infectivity_ratio)
summary_linear_model <- summary(linear_model)

r <- corr_spearman$estimate
p_value <- corr_spearman$p.value
r_squared <- summary_linear_model$r.squared

plot(nef_ratio$infectivity_ratio,
     rpm_PPP2R1B,
     pch = 19,
     col = color_vector,
     xlab = "infectivity ratio",
     ylab = "RPM",
     main = "PPP2R1B")
abline(linear_model)

text(4, max(rpm_PPP2R1B),
     paste("r", ":", round(r, 3)), pos = 2)
text(4, max(rpm_PPP2R1B - 3),
     paste("R" %p% supsc("2"), ":", round(r_squared, 4)), pos = 2)
text(4, max(rpm_PPP2R1B - 6),
     "P < 0.001", pos = 2)

