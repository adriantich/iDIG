#' iDIG Loader
#'
#' This function helps to load the inputs to be processed by the iDIG package.
#' This function is specially designed to load tables with large amounts of
#' columns outperforming the processing time.
#'
#' @details
#' iDIG_loader is capable of optimally load the three main outputs required 
#' by the inversion detector. The output of this function will be a list
#' with three objects named as the following parameters:
#' 
#' - GT: these are the genotipes for each individual
#' 
#' - indv: labels for each individual
#' 
#' - pos: positions for all SNPs
#' 
#' The object list will contain only those objects that have been specified.d
#'
#' @param gt_filename Character vector with the names of the forward fastq or
#' fastq.gz files. Only needed for multiplexed libraries.
#'
#' @param indv_filename Character vector. Acronym for each sequencing library. This
#' acronym must be of 4 characters in capital letters. Do not mix up library and
#' experiment acronyms. The latter will be required in following steps. However
#' they can be the same.
#'
#' @param pos_filename Character string. Acronym for the experiment. This
#' acronym must be of 4 characters in capital letters. Do not mix up library and
#' experiment acronyms. However they can be the same.
#'
#' @export 
#' 
#' @examples
#' 
#' library(idig)
#' 
#' # define the path names for each file
#' 
#' gt_path <- "chr_biallele_005maf_5DP_Raw.012"
#' indv_path <- "chr_biallele_005maf_5DP_Raw.012.indv"
#' pos_path <- "chr_biallele_005maf_5DP_Raw.012.pos"
#' 
#' idig_inputs <- iDIG_loader(gt_filename = gt_path,
#'                            indv_filename = indv_path,
#'                            pos_filename = pos_path)


iDIG_loader <- function(gt_filename = NULL,
                        indv_filename = NULL,
                        pos_filename = NULL)