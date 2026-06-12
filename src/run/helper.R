library(minfi)
args = commandArgs(trailingOnly = T)
analysis_package <- "IlluminaHumanMethylationEPICmanifest"
library(analysis_package, character.only = TRUE)
annotation_package <- "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
library(dplyr)


#' Get M-values
#' @param metadata A data frame containing the meta data, including the paths to the IDAT files.
#' @param opts A list containing options for preprocessing the IDAT files. The list should include:
#' @return A matrix of M-values obtained after preprocessing the IDAT files.
#' @return A matrix of M-values obtained after preprocessing the IDAT files.
#' default preprocessing type set is preprocessIllumina
get_m_values <- function(metadata, opts) {
    # Read IDAT files
    RGset <- read.metharray(metadata$path, extended = TRUE, force= TRUE)
    preprocess_method = opts$pp
    if (preprocess_method == "SWAN") {
        Mset <- preprocessSWAN(RGset)
    } else if (preprocess_method == "Quantile") {
        Mset <- preprocessQuantile(RGset)
    } else if (preprocess_method == "Noob") {
        Mset <- preprocessNoob(RGset)
    } else if (preprocess_method == "Raw") {
        Mset <- preprocessRaw(RGset)
    } else {
        Mset <- preprocessIllumina(RGset)
    } 
    # get m_ values
    m_values <- getM(Mset)
    return(m_values)
}

#' Refactor Results Table for Output
#' This function selects and rounds essential columns (`logFC`, `adj.P.Val`)
#' from the input results table. It is typically used after differential analysis
#' to prepare results for display or downstream processing.
#' @param results A data frame containing the results of a differential analysis,
#'                including columns `logFC` (log fold change) and `adj.P.Val`
#'                (adjusted p-value).
#' @return A data frame containing only the selected columns (`logFC`, `adj.P.Val`),
#'         with values rounded to six decimal places.
get_refactored_result <- function(results)
{
    # Select only essential columns
    essential_cols <- c("CpgId", "logFC", "adj.P.Val")
    print(results)
    print(essential_cols)
    results <- results[, essential_cols, drop = FALSE]
    results$adj.P.Val <- round(results$adj.P.Val, 4)
    results$logFC <- round(results$logFC, 4)

    # Step 1: Force proper numeric conversion
    results$logFc <- as.numeric(results$logFC)
    # Step 2: Remove non-finite values (NA, NaN, Inf, -Inf)
    results <- results[is.finite(results$logFC), ]
    # Step 3: Now safely create bins
    results$logFc_bin <- cut(results$logFC, breaks = 10)
    results <- results %>%
        group_by(logFc_bin) %>%
        sample_n(size = min(50000, n()), replace = FALSE) %>%
        ungroup()
    
    return (results);
}

#' Get annotation
#' This function runs annotation for differential analysis results
#' @param results A data frame containing the results of a differential analysis,
#' @return annotated results with additional information from the annotation package.
#' Steps:
#' 1. Load annotation package
#' 2. Merge DMA results with annotation
#' 3. Select essential columns to avoid processing of large file
#' 4. Enhancer is specific to 450/850K, so referring dynamically without specific types
#' 5. Return the annotated results
get_annotation_result<- function(results)
{
    anno <- getAnnotation(annotation_package)
    results <- merge(results, anno, by.x = "CpgId", by.y = "Name", all.x = TRUE)
    essential_cols <- c("CpgId", "logFC", "adj.P.Val","UCSC_RefGene_Name","Regulatory_Feature_Name")
    enhancer_cols <- grep("Enhancer", names(results), value = TRUE, ignore.case = TRUE)
    all_cols <- unique(c(essential_cols, enhancer_cols))
    missing_cols <- setdiff(all_cols, colnames(results))
    if(length(missing_cols) > 0) {
        stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
    }
    results <- results[, all_cols, drop = FALSE]
    results$adj.P.Val <- round(results$adj.P.Val, 4)
    results$logFC <- round(results$logFC, 4)
    
    # Step 1: Force proper numeric conversion
    results$logFc <- as.numeric(results$logFC)
    # Step 2: Remove non-finite values (NA, NaN, Inf, -Inf)
    results <- results[is.finite(results$logFC), ]
    # Step 3: Now safely create bins
    results$logFc_bin <- cut(results$logFC, breaks = 10)
    results <- as.data.frame(results)
    results <- results %>%
        group_by(logFc_bin) %>%
        sample_n(size = min(50000, n()), replace = FALSE) %>%
        ungroup()
    
    return(results);
}