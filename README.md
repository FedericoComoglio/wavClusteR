wavClusteR
==========

The new wavClusteR package (March 2014), version 1.99.2, in preview by BioC

This R package for PAR-CLIP data analysis can be used to:

+ Infer PAR-CLIP induced transitions and discriminate them from sequencing error, SNPs, contaminants and additional non-experimental causes, using a non-parametric mixture model. 
+ Resolve cluster boundaries at high resolution and provide robust estimation of cluster statistics
+ Integrate RNA-Seq data to estimate False Discovery Rates over the entire range of relative substitution frequencies. 
+ Post-processing of results and functions to export results for UCSC genome browser visualization and motif search analysis. 

Key functions *support parallel multicore computing*. While wavClusteR was designed for PAR-CLIP data analysis, it *can be applied to the analysis of other Next-Generation Sequencing data obtained from substitution inducing experimental procedures* (e.g. BisSeq). Please see the package vignette for getting started with wavClusteR.

Should you experience any issue using wavClusteR for your analysis, please do not hesitate to contact me at federico.comoglio@bsse.ethz.ch

If you are using wavClusteR for your analysis, we kindly ask you to cite the following work:

+ Sievers et al (2012) Mixture models and wavelet transforms reveal high confidence RNA-protein interaction sites in MOV10 PAR-CLIP data, Nucleic Acids Research, 40(20), e160. doi 10.1093/nar/gks697
+ Comoglio et al (2014) wavClusteR: an R package for PAR-CLIP data analysis, submitted