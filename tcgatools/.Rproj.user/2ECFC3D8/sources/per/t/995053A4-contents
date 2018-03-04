##### TCGA TOOLS #####

#' Data for workable examples
#'
#' @name example.tcga.table
#' @docType data
NULL

#' @import graphics
#' @import stats
NULL

# ---- FUNCTION 1 ----

#' Extract indexes of columns of interest
#'
#' @description Extract the indexes of the columns with a given prefix (data type) and/or suffix (molecular class)
#'
#' @param tcga.data Dataframe of TCGA data
#' @param prefix Regular expression for a data type in the TCGA table (e.g. RNA, CNA)
#' @param suffix Regular expression for a gene defining a molecular class
#' @return Returns vector of indexes of columns with prefix or suffix of interest
#' @export
#'
#' @examples data("tcgaData")
#' @examples extract.cols.ix(example.tcga.table, prefix= 'RNA', suffix='PTEN')
#' @examples ## Use regular expressions to select multiple prefixes/suffixes
#' @examples extract.cols.ix(example.tcga.table, prefix='RNA|CNA')
extract.cols.ix <- function(tcga.data, prefix='.+', suffix='.+'){ # Regex for any
  cols.to.keep <- grep(colnames(tcga.data), pattern = paste0('^',prefix,'\\.',suffix, '$'))
  # Display error if the prefix/subtype is not found
  if (length(cols.to.keep) == 0) stop('There are no column names matching the specified pattern.')
  return(cols.to.keep)
  }


# ---- FUNCTION 2 ----

normalize <- function(expression.vector){
  norm.expression <- (expression.vector - stats::median(expression.vector, na.rm = TRUE))/stats::mad(expression.vector, na.rm=TRUE)
  return(norm.expression)
}

plot.normalization <- function(subset.tcga.data, norm.subset.tcga.data){
  # Boxplots to compare non-normalized and normalized data
  graphics::par(mfrow=c(1,2), mar=c(7.1,4.1,4.1,2.1))
  graphics::boxplot(subset.tcga.data, las=2, main = 'Before normalization')
  graphics::boxplot(norm.subset.tcga.data, las=2, main = 'After normalization')
}

#' Normalize selected columns
#'
#' @description Return a data frame with only selected columns normalized, based on median and mad
#'
#' @param tcga.data Dataframe of TCGA data
#' @param cols.to.normalize Index of columns to use for normalization
#' @param show.plot Boolean indicating wheather to show boxplots of distribution of in columns of interest before and after normalization
#' @return Returns TCGA dataframe with selected columns normalized
#' @export
#'
#' @examples ## Normalize columns containing RNA data
#' @examples data('tcgaData')
#' @examples columns.of.interest <- extract.cols.ix(example.tcga.table, prefix='RNA')
#' @examples norm.columns(example.tcga.table, cols.to.normalize = columns.of.interest)
norm.columns <- function(tcga.data, cols.to.normalize, show.plot=FALSE) {
  subset.tcga.data <- tcga.data[,cols.to.normalize]
  if (length(cols.to.normalize)==1) {
    norm.subset.tcga.data <- normalize(subset.tcga.data)
  } else{
    norm.subset.tcga.data <- apply(subset.tcga.data, MARGIN = 2, normalize)
  }
  # Show boxplots of normalized data before and after normalization
  if (show.plot) {
    plot.normalization(subset.tcga.data, norm.subset.tcga.data)
  }
  tcga.data[,cols.to.normalize] <- norm.subset.tcga.data
  return(tcga.data)
}

# ---- FUNCTION 3 ----

extract.metadata <- function(tcga.data) {
  metadata.cols <- grep(colnames(tcga.data), pattern = "\\.", invert = TRUE)
  return(tcga.data[,metadata.cols])
}

#' Extract columns of a subtype
#'
#' @description Extracts columns of TCGA dataframe relative to a given data subtype (RNA, RNAZ, PPA)
#'
#' @param tcga.data Dataframe of TCGA data
#' @param subtype Regular expression for data subtype (RNA, RNAZ, PPA)
#' @param with_metadata Boolean indicating wheather to include metadata columns in the output dataframe
#' @return Returns TCGA dataframe with selected columns
#' @export
#'
#'
#' @examples ## Normalize columns containing RNA data
#' @examples data('tcgaData')
#' @examples extract.subtype(example.tcga.table, 'RNA')
#' @examples ## Use pattern matching to select multiple subtypes
#' @examples extract.subtype(example.tcga.table, 'RNA|CNA')
extract.subtype <- function(tcga.data, subtype, with_metadata=TRUE) {
  subtype.cols <- extract.cols.ix(tcga.data, prefix = subtype)
  sample.col <- grep('sample', colnames(tcga.data))
  if (with_metadata==TRUE) {
    metadata <- extract.metadata(tcga.data)
    subset <- merge(tcga.data[,c(sample.col,subtype.cols)], metadata, by='sample')
  } else {
    subset <- tcga.data[,c(sample.col,subtype.cols)]
  }
  return(subset)
  }







