####Detect_Inversions Function####

Detect_Inversions <- function(iDIG_input, window, step) {

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
        rows_selected <- rows_selected[rows_selected <= nrow(data_cols)]
        window_data <- data_cols[rows_selected,,drop=FALSE]
        return(colMeans(window_data, na.rm = TRUE))
    }
    get_sd_window <- function(spot) {
        rows_selected <- spot:(spot + window - 1)
        rows_selected <- rows_selected[rows_selected <= nrow(data_cols)]
        window_data <- data_cols[rows_selected,,drop=FALSE]
        return(apply(window_data, 2, sd, na.rm = TRUE))
    }
    get_pos_window <- function(spot) {
        rows_selected <- spot:(spot + window - 1)
        rows_selected <- rows_selected[rows_selected <= nrow(data_cols)]
        window_data <- data_cols[rows_selected,,drop=FALSE]
        return(as.data.frame(t(colSums(!is.na(window_data)))))
    }
    get_cluster <- function(spot) {
        # for(spot in spots){
        rows_selected <- spot:(spot + window - 1)
        rows_selected <- rows_selected[rows_selected <= nrow(data_cols)]
        window_data <- data_cols[rows_selected,,drop=FALSE]

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
    data_cols <- apply(data_cols, 2, as.numeric)
    output <- list(results = data.frame(),
                   result_sd = data.frame(),
                   positions = data.frame(),
                   clusters = data.frame())

    for (chr in unique(chromosome)) {
    
        iDIG_input_chr <- iDIG_input[which(chromosome == chr),]

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

window_step_plot <- function(iDIG_obj) {
    ##############
    # iDIG_obj <- output
    ##############
    library(ggplot2)
    library(reshape2)
    library(plotly)

    # from wide format to long format
    inv_long <- melt(iDIG_obj$results, id.vars = c("chromosome", "window", "step","POS"))

    p <- ggplot(inv_long,aes(x=POS,y = value)) +
        geom_line(aes( group = variable)) +
        facet_grid(window ~ step, labeller = labeller(window = label_both, step = label_both)) +
        labs(x = "Position", y = "Value") +
        theme(strip.text.x = element_text(size = 10, face = "bold"),
              strip.text.y = element_text(size = 10, face = "bold"))

    return(p)
    # ggsave("plot.png",p,width = 15, height = 10)

}

