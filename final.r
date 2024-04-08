#Clearing envirnment variables
rm(list = ls())

# Installing and importing libraries
if (rstudioapi::isAvailable()) {
  if (require("rstudioapi") != TRUE) {
    install.packages("rstudioapi")
  }else {
    library(rstudioapi)
  }
  wdir <- dirname(getActiveDocumentContext()$path)
}
library("biomaRt")
library(stringr)

# Setting Working Directory
setwd(paste0(wdir, "/Data"))


# Getting the start and end index of proteins
start_index <- 1
end_index <- 58302

# Importing nef data
nef_ratio <- read.delim("infect.txt", header = FALSE)


# Getting the unique IDs of all proteins and saving them in a csv file
get_unique_ids <- function() {
  all_ids <- c()

  for (file in nef_ratio$V1){
    f <- read.delim(paste0(file, ".htseq"), header = FALSE)
    f_clean <- head(f, -5)

    # Extract IDs from f_clean and append to all_ids vector
    ids_in_file <- unique(f_clean$V1)
    all_ids <- c(all_ids, ids_in_file)
  }

  unique(all_ids)
}

proteins <- data.frame(id = get_unique_ids())

for (i in seq(1, nrow(proteins))) {
  proteins$p_value[i] <- NA
  proteins$corr[i] <- NA
}

View(proteins)

write.csv(proteins, file = "proteins.csv", row.names = FALSE)

# Creating the proteins df from the previously made csv file
proteins <- read.csv(file = "proteins.csv", na.strings = c("NA", "N/A", ""))
View(proteins)

# Function to get the Reads per Millioin(RPM) of a given protein
get_rpm <- function(protein) {

  rpm <- c()

  for (file in nef_ratio$V1){
    f <- read.delim(paste0(file, ".htseq"), header = FALSE)
    f_clean <- head(f, -5)

    RPM <- (f_clean[f_clean$V1 == protein, ]$V2 / sum(f_clean$V2)) * 1000000
    rpm <- append(rpm, RPM)
  }
  rpm
}


# Function to calculate the correlation and p_value of a protein from the list of proteins
# with start and end index set previously
# This was done so that the code could be run many times with different indices
# on many laptops as there are a total of 58,302 unique proteins
get_correlation <- function(start, end) {
  for (i in seq(start, end)) {
    protein_id <- proteins$id[i]
    protein_p_value <- proteins$p_value[i]
    protein_corr <- proteins$corr[i]
    if (is.na(protein_p_value) && is.na(protein_corr)) {
      rpm_values <- get_rpm(protein_id)

      # Add RPM values to nef_ratio with a new column name based on the protein ID
      nef_ratio[[paste0("RPM_", toupper(protein_id))]] <<- rpm_values
      corr <- cor.test(nef_ratio$V2, rpm_values, method = "pearson")
      proteins[i, "corr"] <<- corr$estimate
      proteins[i, "p_value"] <<- corr$p.value
      cat("Protein", i, "out of", end, "done", "\n")
    }
  }
}

get_correlation(start_index, end_index)

View(proteins)

# Writing the updated proteins df with correlation and p_value to a csv
write.csv(proteins, file = "proteins.csv", row.names = FALSE)

# Filtering the proteins whose pearson correlation p-values gave value less than 0.05 and removing null values
proteins_filtered <- proteins[proteins$p_value < 0.05, ]
proteins_filtered_clean <- proteins_filtered[complete.cases(proteins_filtered), ]

View(proteins_filtered_clean)

# Using biomaRt library to find the gene hgnc symbol based on the ensemble id
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <-  proteins_filtered_clean$id

gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
              values = genes, mart= mart)

proteins_filtered_clean <- merge(proteins_filtered_clean, gene_IDs, by.x = "id", by.y = "ensembl_gene_id", all.x = TRUE)


View(proteins_filtered_clean)

# Sorting the proteins in ascending order of p-value to find the best correlated proteins
proteins_sorted <- proteins_filtered_clean[order(proteins_filtered_clean$p_value), ]

View(proteins_sorted)

# Writing the final dataframe to a csv
write.csv(proteins_sorted, file = "proteins_below_threshold.csv", row.names = FALSE)

# Saving the infectivity ratio with rpm in another csv file
for (i in seq(1, nrow(proteins_sorted))) {
  protein_id <- proteins_sorted$id[i]
  get_rpm(protein_id)
  cat("Protein: ", i, "out of ", nrow(proteins_sorted), "done\n")
}

write.csv(nef_ratio, file = "infectivity_with_rpm.csv", row.names = FALSE)

proteins_below_threshold <- read.csv("proteins_below_threshold.csv")
proteins_below_threshold <- subset(proteins_below_threshold, select = c(id, p_value, corr, hgnc_symbol))
proteins_below_threshold <- proteins_below_threshold[order(proteins_below_threshold$p_value), ]
infectivity_with_rpm <- read.csv("infectivity_with_rpm.csv")

for (i in 4:ncol(infectivity_with_rpm)) {
  col_name <- colnames(infectivity_with_rpm)[i]
  correlation_test <- cor.test(infectivity_with_rpm[[col_name]],
                               infectivity_with_rpm$infectivity_ratio,
                               method = "spearman",
                               exact = FALSE)
  id <- str_extract(col_name, "(?<=RPM_).*")  # Extract the ID using regex
  print(col_name)
  cat("corr estimate:", correlation_test$estimate, "of id:", id, "\n",
      "\n","number", i, "out of:", ncol(infectivity_with_rpm), "done\n")
  proteins_below_threshold$corr_spearman[proteins_below_threshold$id == id] <- correlation_test$estimate
  proteins_below_threshold$p_value_spearman[proteins_below_threshold$id == id] <- correlation_test$p.value
}

View(proteins_below_threshold)

write.csv(proteins_below_threshold, file = "proteins_below_threshold.csv", row.names = FALSE)
