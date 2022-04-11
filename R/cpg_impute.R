#' Gestational DNAm age estimation using different DNA methylation clocks. IMPORTANT: This file includes the adapted function "DNAmGA".
#' It only completes the imputation and skips the other steps (that will be carried out differently). 
#' The file function cpg_impute carries out the imputation for all GA clocks and returns the object "cpg.imp" that will be used for further processing.
#' This adaptation has been made so that larger data sets can be handeled
#' @param x data.frame (Individual in columns, CpGs in rows, CpG names in first colum - i.e. Horvath's format), matrix (individuals in columns and Cpgs in rows having CpG names in the rownames), ExpressionSet or GenomicRatioSet.
#' @param toBetas Should data be transformed to beta values? Default is FALSE. If TRUE, it implies data are M values.
#' @param fastImp Is fast imputation performed if necessary? (see details). Default is FALSE
#' @param normalize Is Horvath's normalization performed? By default is FALSE
#' @param age individual's chronological age. Required to compute gestational age difference output
#' @param cell.count Are cell counts estimated? Default is TRUE.
#' @param cell.count.reference Used when 'cell.count' is TRUE. Default is "blood gse35069 complete". See 'meffil::meffil.list.cell.count.references()' for possible values.
#' @param min.perc Indicates the minimum conicidence percentage required between CpGs in or dataframee x and CpGs in clock coefficients to perform the calculation. If min.prec is too low, the estimated gestational DNAm age can be poor
#' @param ... Other arguments to be passed through impute package
#'
#' @details Imputation is performed when having missing data.
#'          Fast imputation is performed by ...
#'          what about imputing only when CpGs for the clock are missing?
#'
#' @examples
#' TestDataset[1:5, ]
#' ga.test <- DNAmGA(TestDataset)
#' 
#' @return the estimated gestational DNAm age
#' 
#' @import impute dplyr tidyverse tibble
#' @importFrom minfi getBeta
#' @importFrom Biobase featureNames exprs
#' @importFrom minfi getBeta
#' 
#' @export

cpg_impute <- function(x, toBetas = FALSE,
                   fastImp = FALSE,
                   normalize = FALSE,
                   age,
                   cell.count = TRUE,
                   cell.count.reference = "andrews and bakulski cord blood",
                   min.perc = 0.8,
                   ...) {
  if (inherits(x, "data.frame")) {
    cpgs.names <- as.character(x[, 1, drop = TRUE])
    if (length(grep("cg", cpgs.names)) == 0) {
      stop("First column should contain CpG names")
    }
    cpgs <- t(as.matrix(x[, -1]))
    colnames(cpgs) <- cpgs.names
  }
  else if (inherits(x, "matrix")) {
    cpgs <- t(x)
  }
  else if (inherits(x, "ExpressionSet")) {
    cpgs <- t(Biobase::exprs(x))
  }
  else if (inherits(x, "GenomicRatioSet")) {
    cpgs <- t(minfi::getBeta(x))
  }
  else {
    stop("x must be a data.frame, matrix, 'GenomicRatioSet' or an 'ExpressionSet' object")
  }
  
  if (toBetas) {
    toBeta <- function(m) {
      2^m / (2^m + 1)
    }
    cpgs <- toBeta(cpgs)
  }
  
  if (any(cpgs < -0.1 | cpgs > 1.1, na.rm = TRUE)) {
    stop("Data seems to do not be beta values. Check your data or set 'toBetas=TRUE'")
  }
  
  cpgs.all <- c(
    as.character(coefKnightGA$CpGmarker[-1]),
    as.character(coefBohlin$CpGmarker[-1]),
    as.character(coefMayneGA$CpGmarker[-1]),
    as.character(coefLeeGA$CpGmarker[-1]),
    as.character(coefEPIC$CpGmarker[-1])
  )
  
  if (any(!cpgs.all %in% colnames(cpgs))) {
    warning("CpGs in all Gestational Age clocks are not present in your data. Try 'checkClocksGA' function
            to find the missing CpGs of each method.")
  }
  
  cpgs.in <- intersect(cpgs.all, colnames(cpgs))
  
  miss <- apply(cpgs[, cpgs.in], 2, function(x) any(is.na(x)))
  
  if (any(miss)) {
    if (fastImp) {
      cat(paste("Imputing missing data of", sum(miss), "CpGs .... \n"))
      mm <- apply(cpgs[, cpgs.in], 2, median, na.rm = TRUE)
      cpgs.imp <- sweep(cpgs[, cpgs.in], 2,
                        STATS = mm,
                        FUN = function(x, s) ifelse(is.na(x), s, x)
      )
    }
    else {
      quiet <- function(x) {
        sink(tempfile())
        on.exit(sink())
        invisible(force(x))
      }
      cat(paste("Imputing missing data of the entire matrix .... \n"))
      cpgs.imp <- quiet(t(impute.knn(t(cpgs), ...)$data))
    }
    cat("Data imputed. Starting DNAm clock estimation ... \n")
  }
  else {
    cpgs.imp <- cpgs
  }
  cpgs.imp <-cpgs.imp
}
  

  
  