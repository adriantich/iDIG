#' Scan chromosomes by number of positions
#' 
#' This function will perform a scanning along the chromosomes
#' taking into account ranges of nucleotide positions of "window"
#' length and "step" number of positions betweent he starts of the
#' windows.
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
#' Warnings:
#' 
#' - If a window does not have any value, the mean and standard deviation
#' will be NA.
#' 
#' - the last window of each chromosome can be shorter than the given
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
#' @param window Numeric vector. The window is the number of nucleotides
#' positions that will be used to obtain each descriptor value.
#' 
#' @param step Numeric vector. Distance in number of nucleotides between
#' two window starts. If Step is lower than window, then an overlap between
#' windows occur.
#' 
#' @param min_size Numeric. Minimum number of positions in the chromosome between
#' the first SNP and the last SNP.
#' 
#' @returns
#' List of four dataframes:
#' 
#' - results: mean values of each window for each individual
#' 
#' - results_sd: standard deviation values for each window and individual.
#' 
#' - position: number of nucleotide positions used to compute the mean values
#' of each window.
#' 
#' - clusters: for each window, the individuals are clustered into the following
#' clusters using the kmean algorithm with 100 iterations: for centers equal to
#' 3, this is three possible clusters, A, B and C from highest to lowest mean
#' window mean result; D and E if only two centers can be computed; and F if
#' no clustering can be computed (all individuals have equal genotype window)
#'
#' @export
#' 

scan_by_pos <- function(iDIG_input, window, step, min_size = 1) {

    ##############
    # source("R/iDIG_loader.R")

    # gt_filename <- '../iDIG_Rpackage/2allele_tots_012_chr11.012'
    # indv_filename <- '../iDIG_Rpackage/2allele_tots_012_chr11.012.indv'
    # pos_filename <- '../iDIG_Rpackage/2allele_tots_012_chr11.012.pos'

    # idig_inputs <- iDIG_loader(gt_filename = gt_filename,
    #                         indv_filename = indv_filename,
    #                         pos_filename = pos_filename)

    # iDIG_input <- idig_inputs
    # window <- c(9,100,1000)
    # step <- c(7,100)
    ##############
    get_mean_window_pos <- function(spot) {
        rows_selected <- spot:(spot + window - 1)
        rows_selected <- rows_selected[rows_selected <= total]
        window_data <- iDIG_input_chr[iDIG_input_chr$POS %in% rows_selected,,drop=FALSE]
        return(colMeans(window_data, na.rm = TRUE))
    }
    get_sd_window_pos <- function(spot) {
        rows_selected <- spot:(spot + window - 1)
        rows_selected <- rows_selected[rows_selected <= total]
        window_data <- iDIG_input_chr[iDIG_input_chr$POS %in% rows_selected,,drop=FALSE]
        return(apply(window_data, 2, sd, na.rm = TRUE))
    }
    # get_pos_window_pos <- function(spot) {
    #     rows_selected <- spot:(spot + window - 1)
    #     rows_selected <- rows_selected[rows_selected <= total]
    #     window_data <- iDIG_input_chr[iDIG_input_chr$POS %in% rows_selected,,drop=FALSE]
    #     return(as.data.frame(t(colSums(!is.na(window_data)))))
    # }
    get_cluster_pos <- function(spot) {
        # for(spot in spots){
        rows_selected <- spot:(spot + window - 1)
        rows_selected <- rows_selected[rows_selected <= total]
        window_data <- iDIG_input_chr[iDIG_input_chr$POS %in% rows_selected,,drop=FALSE]
        
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

        # be sure that the positions are ordered
        iDIG_input_chr <- iDIG_input_chr[order(iDIG_input_chr$POS),]

        # check for the min size of the chromosome
        if(iDIG_input_chr$POS[nrow(iDIG_input_chr)] - iDIG_input_chr$POS[1] < min_size){
            warning("The chromosome ", chr, " is smaller than the minimum size")
            next
        }

        total <- iDIG_input_chr$POS[nrow(iDIG_input_chr)]
        first_pos <- iDIG_input_chr$POS[1]
        for(comb in seq_len(nrow(combinations))){
            window <- combinations[comb,1]
            step <- combinations[comb,2]

            # to get the start of each window I create a matrix 
            # with a vector from the first position to the last multiple
            # position of step that include the last window start.
            # The first row of the matrix will be the start of each
            # window.
            
            last_pos_matrix <- ceiling((total - first_pos+1- window)/step)*step+first_pos-1
            spot_starts <- matrix((first_pos:last_pos_matrix),
                                         nrow = step)
            spot_starts <- spot_starts[1,]

            # in a case of high window and low step, the position matrix
            # can be too large. In this case the process is done in chunks.
            chunks <- ceiling(window/step)
            first_spot <- 1
            if(length(spot_starts) < chunks){
                chunks <- length(spot_starts)
            }
            last_spot <- ceiling(length(spot_starts)/chunks)

            for (i in 1:chunks){
                spots <- spot_starts[first_spot:last_spot]
                positions_matrix <- matrix(rep(0:(window-1),last_spot-first_spot+1),
                                           nrow = last_spot-first_spot+1,
                                           byrow = TRUE)
                positions_matrix <- positions_matrix + spots


                positions_matrix[!positions_matrix %in% iDIG_input_chr$POS] <- NA
                # in order to have a better performance those windows with only one
                # position are excluded from the scanning but once I merge the
                # result with the excluded I need to sort them back.
                spots_1pos <- spots[rowSums(!is.na(positions_matrix))==1]
                spots <- spots[rowSums(!is.na(positions_matrix))>1]

                result_1pos <- iDIG_input_chr[iDIG_input_chr$POS %in% spots_1pos,]
                result <- as.data.frame(t(sapply(spots, get_mean_window_pos)))
                result <- rbind(result, result_1pos)
                result <- result[order(result$POS),]

                result_sd_1pos <- iDIG_input_chr[iDIG_input_chr$POS %in% spots_1pos,]
                result_sd_1pos[,-1] <- NA
                result_sd <- as.data.frame(t(sapply(spots, get_sd_window_pos)))
                result_sd <- rbind(result_sd, result_sd_1pos)
                result_sd <- result_sd[order(result_sd$POS),]

                positions <- cbind(data.frame(POS = rowMeans(positions_matrix, na.rm = TRUE)),
                                matrix(rowSums(!is.na(positions_matrix)),
                                        ncol = ncol(iDIG_input_chr)-1,
                                        nrow = nrow(positions_matrix)))
                colnames(positions) <- colnames(iDIG_input_chr)
                positions <- positions[!is.na(positions$POS),]

                clstr <- lapply(c(spots_1pos,spots), get_cluster_pos)
                clstr <- do.call(rbind, clstr)
                clstr <- clstr[order(clstr$POS),]
            
                result <- cbind(chromosome = chr, result, window = window, step = step)
                result_sd <- cbind(chromosome = chr, result_sd, window = window, step = step)
                positions <- cbind(chromosome = chr, positions, window = window, step = step)
                clstr <- cbind(chromosome = chr, clstr, window = window, step = step)
                output$results <- rbind(output$results, result)
                output$result_sd <- rbind(output$result_sd, result_sd)
                output$positions <- rbind(output$positions, positions)
                output$clusters <- rbind(output$clusters, clstr)

                first_spot <- last_spot + 1
                if(i == (chunks-1)){
                    last_spot <- length(spot_starts)
                } else {
                    last_spot <- last_spot + ceiling(length(spot_starts)/chunks)
                }
            }
        }
    }
    return(output)
}