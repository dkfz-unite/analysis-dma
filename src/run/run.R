library(minfi)
library(jsonlite)
source("helper.R")
source("da_analysis.R")

args = commandArgs(trailingOnly = TRUE)

inputFilePath <- args[1]
outputFilePath <- args[3]
optionsFilePath <- args[2]

options <- fromJSON(optionsFilePath)

metadata <- read.table(file = inputFilePath, header = T, sep = "\t", check.names = F)

# get m-values
m_values = get_m_values(metadata, options)

# do the differential expression analysis
results <- da_analysis(data_matrix=m_values,
            condition=metadata$condition,
            ref_category = tail(metadata$condition,1) # assign the second (i.e. last) condition as the reference category
            )

# Binding CpgId to results
results <- cbind(CpgId = rownames(results), results)

#refactored results
reduced_results = get_refactored_result(results)
rm(results)

# Get annotation results
annotated_results <- get_annotation_result(reduced_results)
rm(reduced_results)

# Save annotated results
write.table(annotated_results, outputFilePath, row.names = FALSE, quote = FALSE, sep = "\t")
rm(annotated_results)
gc()