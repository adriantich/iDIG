#' iDIG Loader
#'
#' This function helps to load the inputs to be processed by the iDIG package.
#' This function is specially designed to load tables with large amounts of
#' columns outperforming the processing time.
#'
#' @details
#' iDIG_loader is capable of optimally load the three main outputs required 
#' by the inversion detector. The output of this function will be data frame
#' in which the first two columns will be the chromosome (CHR) and the position
#' (POS) of the SNPs. The rest of the columns will be the genotypes for each
#' individual. If no individual names are provided, the columns will be
#' named automatically following the order of the individuals in the file.
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
                        pos_filename = NULL) {

    # add the pos to the list if specified
    if (is.null(pos_filename)) {
        stop("pos_filename is NULL. pos_filename is required.")
    } else {
        if (!file.exists(pos_filename)) {
            stop("The file ", pos_filename, " does not exist.")
        }
        pos_table <- read.table(pos_filename, header = FALSE)
        colnames(pos_table) <- c("CHR", "POS")
    }

    # add the GT to the list if specified
    if (is.null(gt_filename)) {
        stop("gt_filename is NULL. gt_filename is required.")
    } else {
        if (!file.exists(gt_filename)) {
            stop("The file ", gt_filename, " does not exist.")
        }

        # load the file and transpose it
        gt_table <- importtransposedTSV(gt_filename)
        # add the indv to the list if specified
        if (is.null(indv_filename)) {
            message("indv_filename is NULL. The individuals won't have name")
            colnames(gt_table) <- gt_table[1,]
        } else {
            if (!file.exists(indv_filename)) {
                stop("The file ", indv_filename, " does not exist.")
            }
            names <- read.table(indv_filename, header = FALSE)
            colnames(gt_table) <- names[,1]
        }
        gt_table <- gt_table[-1,]
    }
    return(cbind(pos_table, gt_table))
}

library(Rcpp)

# Define the C++ function
cppFunction('
    #include <Rcpp.h>
    #include <fstream>
    #include <sstream>
    #include <vector>
    using namespace Rcpp;

    // [[Rcpp::export]]
    DataFrame importtransposedTSV(std::string input_file) {
        std::ifstream infile(input_file);
        if (!infile.is_open()) {
            stop("Could not open file");
        }

        std::vector<std::vector<std::string>> data;
        std::string line, cell;

        // Read the TSV file
        while (std::getline(infile, line)) {
            std::stringstream linestream(line);
            std::vector<std::string> row;
            while (std::getline(linestream, cell, \'\\t\')) {
                row.push_back(cell);
            }
            data.push_back(row);
        }

        if (data.empty()) {
            stop("No data read from file");
        }

        // Convert to DataFrame
        List df(data.size());
        for (size_t i = 0; i < data.size(); ++i) {
            df[i] = wrap(data[i]);
        }
        DataFrame result = DataFrame(df);
        result.attr("names") = R_NilValue; // Remove column names
        return result;
    }
')