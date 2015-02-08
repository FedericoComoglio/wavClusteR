

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





#' An R package for the analysis of Next-Generation Sequencing data obtained
#' from PAR-CLIP or other substitution-inducing experimental procedures.
#' 
#' Infer PAR-CLIP induced transitions and discriminate them from sequencing
#' error, SNPs, contaminants and additional non-experimental causes, using a
#' non-parametric mixture model. wavClusteR resolves cluster boundaries at high
#' resolution and provides robust estimation of cluster statistics. In
#' addition, the package allows to integrate RNA-Seq data to estimate FDR over
#' the entire range of relative substitution frequencies. Furthermore, the
#' package provides post-processing of results and functions to export results
#' for UCSC genome browser visualization and motif search analysis. Key
#' functions support parallel multicore computing. While wavClusteR was
#' designed for PAR-CLIP data analysis, it can be applied to the analysis of
#' other Next-Generation Sequencing data obtained from substitution inducing
#' experimental procedures (e.g. BisSeq)
#' 
#' \tabular{ll}{ Package: \tab wavClusteR\cr Type: \tab Package\cr Version:
#' \tab 1.99.3\cr Date: \tab 2015-15-01\cr License: \tab GPL-2\cr }
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



