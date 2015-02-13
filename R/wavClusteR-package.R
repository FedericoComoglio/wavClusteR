

#' Components of the non-parametric mixture moodel fitted on Ago2 PAR-CLIP data
#' 
#' The non-parametric mixture model was fit on the entire Ago2 public available
#' PAR-CLIP dataset (Kishore et al.) using the \code{fitMixtureModel} function.
#' 
#' 
#' @name model
#' @docType data
#' @usage data(model) 
#' model
#' @format List of 5 items containing the estimated mixing coefficients and
#' model densities. See the help page of the \link{fitMixtureModel} function
#' for a detailed description of the output.
#' @seealso \code{\link{fitMixtureModel}}
#' @references Kishore et al. A quantitative analysis of CLIP methods for
#' identifying binding sites of RNA-binding proteins. Nat Methods (2011) vol. 8
#' (7) pp. 559-64
#' \url{http://www.nature.com/nmeth/journal/v8/n7/full/nmeth.1608.html}
#' @keywords datasets
NULL

#'
#' A comprehensive pipeline for the analysis of PAR-CLIP data. PAR-CLIP-induced #'transitions are first discriminated from sequencing errors, SNPs and additional non-#'experimental sources by a non-parametric mixture model. The protein binding sites #'(clusters) are then resolved at high resolution and cluster statistics are estimated #'using a rigorous Bayesian framework. Post-processing of the results, data export for #'UCSC genome browser visualization and motif search analysis are provided. In addition, #'the package allows to integrate RNA-Seq data to estimate the False Discovery Rate of #'cluster detection.  Key functions support parallel multicore computing. Note: while #'wavClusteR was designed for PAR-CLIP data analysis, it can be applied to the analysis of #'other NGS data obtained from experimental procedures that induce nucleotide #'substitutions (e.g. BisSeq).#' 
#'
#' \tabular{ll}{ Package: \tab wavClusteR\cr Type: \tab Package\cr Version:
#' \tab 1.99.4\cr Date: \tab 2015-13-02\cr License: \tab GPL-2\cr }
#' 
#' @name wavClusteR-package
#' @aliases wavClusteR-package wavClusteR
#' @docType package
#' @author Federico Comoglio and Cem Sievers
#' 
#' Epigenomics Group
#' 
#' Department of Biosystems Science and Engineering (D-BSSE)
#' 
#' ETH Zurich, Switzerland
#' 
#' Maintainer: Federico Comoglio
#' 
#' \email{federico.comoglio@@bsse.ethz.ch}
#'
#' @references Sievers C, Schlumpf T, Sawarkar R, Comoglio F and Paro R. (2012) Mixture
#' models and wavelet transforms reveal high confidence RNA-protein interaction
#' sites in MOV10 PAR-CLIP data, Nucleic Acids Res. 40(20):e160. doi:
#' 10.1093/nar/gks697
#' 
#' Comoglio F, Sievers C and Paro R (2015) Sensitive and highly resolved identification
#' of RNA-protein interaction sites in PAR-CLIP data, BMC Bioinformatics 16, 32.
#'
#' @keywords package
NULL



