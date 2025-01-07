#' Scan chromosomes by number of SNPs
#' 
#' This function will perform a scanning along the chromosomes 
#' taking into account ranges of SNPs of "window" length and "step"
#' number of SNPs betweent he starts of the windows.
#' 
#' @details
#' This function gets the dataframe obtained from the iDIG_loader function
#' and scans the genome to obtain the mean allelic value, the estandard
#' deviation, the mean position and the cluster obtained from the kmean
#' clustering algorithm for each window on "window" number of SNPs. The distance
#' from the start of one window and the next one is defined by the "step"
#' parameter. Only step equal or less than window are used. If the window
#' and/or step are vectors of more than one value, the combinations that
#' follow the previous condition mentioned are used.
#' 
#' Warning: the last window of each chromosome can be shorter than the given
#' window size if:
#' 
#' (total_size - window)/step - ceiling((total_size-window)/step) != 0
#' 
#' @param iDIG_input DataFrame obtained by with the iDIG_loader function.
#' This dataframe consist of one column with the label for the
#' chromosome ('CHR') and a column for the positions ('POS') followed
#' by as many columns as individuals analysed. Each row represents a SNP
#' and with values for the individuals from 0 to 2 (0, homozigote with the 
#' refference genome; 1 heterozigote; 2 homozigote for the alternative
#' allele)
#' 
#' @param window Numeric vector. The window is the number of SNPs that will
#' be used to obtain each descriptor value.
#' 
#' @param step Numeric vector. Distance in number of SNPs between two window
#' starts. If Step is lower than window, then an overlap between windows 
#' occur.
#' 
#' @returns
#' List of four dataframes:
#' 
#' - results: mean values of each window for each individual
#' 
#' - results_sd: standard deviation values for each window and individual.
#' 
#' - position: number of SNPs used to compute the mean values of each window.
#' 
#' - clusters: for each window, the individuals are clustered into the following
#' clusters using the kmean algorithm with 100 iterations: for centers equal to 
#' 3, this is three possible clusters, A, B and C from highest to lowest mean
#' window mean result; D and E if only two centers can be computed; and F if
#' no clustering can be computed (all individuals have equal genotype window)
#'
#' @export
#' 

scan_by_snps <- function(iDIG_input, window, step) {

    ##############
    # source("R/iDIG_loader.R")

    # gt_filename <- '../iDIG_Rpackage/2allele_tots_012_chr11.012'
    # indv_filename <- '../iDIG_Rpackage/2allele_tots_012_chr11.012.indv'
    # pos_filename <- '../iDIG_Rpackage/2allele_tots_012_chr11.012.pos'

    # idig_inputs <- iDIG_loader(gt_filename = gt_filename,
    #                         indv_filename = indv_filename,
    #                         pos_filename = pos_filename)

    # iDIG_input <- idig_inputs
    # window <- c(10,100,1000)
    # step <- c(10,100)
    ##############
    get_mean_window <- function(spot) {
        rows_selected <- spot:(spot + window - 1)
        rows_selected <- rows_selected[rows_selected <= nrow(iDIG_input_chr)]
        window_data <- iDIG_input_chr[rows_selected,,drop=FALSE]
        return(colMeans(window_data, na.rm = TRUE))
    }
    get_sd_window <- function(spot) {
        rows_selected <- spot:(spot + window - 1)
        rows_selected <- rows_selected[rows_selected <= nrow(iDIG_input_chr)]
        window_data <- iDIG_input_chr[rows_selected,,drop=FALSE]
        return(apply(window_data, 2, sd, na.rm = TRUE))
    }
    get_pos_window <- function(spot) {
        rows_selected <- spot:(spot + window - 1)
        rows_selected <- rows_selected[rows_selected <= nrow(iDIG_input_chr)]
        window_data <- iDIG_input_chr[rows_selected,,drop=FALSE]
        return(as.data.frame(t(colSums(!is.na(window_data)))))
    }
    get_cluster <- function(spot) {
        # for(spot in spots){
        rows_selected <- spot:(spot + window - 1)
        rows_selected <- rows_selected[rows_selected <= nrow(iDIG_input_chr)]
        window_data <- iDIG_input_chr[rows_selected,,drop=FALSE]

        # Perform k-means clustering of the mean values of each individual
        # we asume that there will be 3 clusters (populations)
        ind_means <- colMeans(window_data, na.rm = TRUE)
        # it assumes that the first column is the position
        
        kmeans_result <- try(kmeans(ind_means[ -1], centers = 3, iter.max = 100), silent = TRUE)
        # if only two clusters are found, try with 2, D and E
        if(inherits(kmeans_result, "try-error")){
            kmeans_result <- try(kmeans(ind_means[ -1], centers = 2, iter.max = 100), silent = TRUE)
            # if all individuals are in the same cluster name it F
            if(inherits(kmeans_result, "try-error")){
                clusters <- data.frame(cluster = rep("F", (length(ind_means)-1)))
            } else{
                clusters <- kmeans_result$cluster
                # rename cluster from highest to lowest mean
                # D for the highest mean, E for the lowest
                clusters <- factor(clusters, levels = order(kmeans_result$centers, decreasing = TRUE))
                levels(clusters) <- c("D", "E")
            }
        } else {
            clusters <- kmeans_result$cluster
            # rename cluster from highest to lowest mean
            # A for the highest mean, B for the second highest, C for the lowest
            clusters <- factor(clusters, levels = order(kmeans_result$centers, decreasing = TRUE))
            levels(clusters) <- c("A", "B", "C")
        }
        
        # following this population structure test its significance with the mean values
        # of the individuals
        # aovsum <- summary(aov(ind_means[ -1] ~ clusters))
        # pvalue <- aovsum[[1]][["Pr(>F)"]][1]
        output <- as.data.frame(t(ind_means))
        output[,-1] <- clusters
        # }
        # output$pvalue <- pvalue
        return(output)
    }

    # combine windows and steps to get all combinations in a data frame
    combinations <- expand.grid(window = window, step = step)

    # step can't be bigger than window
    combinations <- combinations[combinations$step <= combinations$window,]
    if(nrow(combinations) == 0){
        stop("Step can't be bigger than window")
    }

    chromosome <- iDIG_input[, 1]
    iDIG_input[iDIG_input == -1] <- NA
    data_cols <- iDIG_input[, -1]  # All columns except the first (chromosome)
    data_cols <- as.data.frame(apply(data_cols, 2, as.numeric))
    output <- list(results = data.frame(),
                   result_sd = data.frame(),
                   positions = data.frame(),
                   clusters = data.frame())

    for (chr in unique(chromosome)) {
    
        iDIG_input_chr <- data_cols[which(chromosome == chr),]

        total <- nrow(iDIG_input_chr)
        for(comb in seq_len(nrow(combinations))){
            window <- combinations[comb,1]
            step <- combinations[comb,2]
            spots <- seq(from = 1, to = (total - window + 1 + step), by = step) # I add the step to the total to make sure the last window is included even with NAs

            result <- as.data.frame(t(sapply(spots, get_mean_window)))

            result_sd <- as.data.frame(t(sapply(spots, get_sd_window)))

            positions <- as.data.frame(t(sapply(spots, get_pos_window)))

            clstr <- lapply(spots, get_cluster)
            clstr <- do.call(rbind, clstr)
          
            result[is.na(result)] <- NA


            result <- cbind(chromosome = chr, result, window = window, step = step)
            result_sd <- cbind(chromosome = chr, result_sd, window = window, step = step)
            positions <- cbind(chromosome = chr, positions, window = window, step = step)
            clstr <- cbind(chromosome = chr, clstr, window = window, step = step)
            output$results <- rbind(output$results, result)
            output$result_sd <- rbind(output$result_sd, result_sd)
            output$positions <- rbind(output$positions, positions)
            output$clusters <- rbind(output$clusters, clstr)
        }
    }
    return(output)
}