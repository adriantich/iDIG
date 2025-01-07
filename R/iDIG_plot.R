#' Plot the results of the iDIG function
#'
#' 
#' @details
#' This function creates a plot.
#' 
#' @param iDIG_obj List of four dataframes. The output of the iDIG function.
#' 
#' @export
#' 
#' @returns
#' A ggplot object.


iDIG_plot <- function(iDIG_obj) {
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