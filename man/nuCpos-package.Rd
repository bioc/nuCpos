\name{nuCpos-package}
\alias{nuCpos-package}
\alias{nuCpos}
\docType{package}
\title{
    An R package for nucleosome positioning prediction
}
\description{
    \pkg{nuCpos}, a derivative of \pkg{NuPoP}, 
    is an R package for prediction of nucleosome positions. 
    \pkg{nuCpos} calculates local and whole nucleosomal 
    histone binding affinity (HBA) scores for 
    a given 147-bp sequence. Note: This package was 
    designed to demonstrate the use of chemical maps in 
    prediction. As the parental package \pkg{NuPoP} now 
    provides chemical-map-based prediction, 
    the function for dHMM-based prediction was removed 
    from this package. Please refer to Xi et al. (2010) and 
    Wang et al. (2008) for technical details of 
    \pkg{NuPoP}. nuCpos continues to provide functions 
    for HBA calculation. The models are based on 
    chemical maps of nucleosomes from 
    budding yeast (Brogaard et al. (2012)), 
    fission yeast (Moyle-Heyrman et al. (2012)), 
    or mouse embryonic stem cells (Voong et al. (2016)). 
}
\details{
\tabular{ll}{
    Package: \tab nuCpos\cr
    Type: \tab Package\cr
    Version: \tab 1.17.3\cr
    Date: \tab 2023-02-15\cr
    License: \tab  GPL-2\cr
}

\code{\link{HBA}}: R function for calculation of the histone 
    binding affinity score of a whole nucleosome.

\code{\link{localHBA}}: R function for calculation of 
    the local histone binding affinity.

}
\author{
    Hiroaki Kato and Takeshi Urano

    Maintainer: Hiroaki Kato<hkato@med.shimane-u.ac.jp>
}
\references{
\enumerate{
    \item
    Wang JP, Fondufe-Mittendorf Y, Xi L, Tsai GF, Segal E 
    and Widom J (2008). Preferentially quantized linker DNA 
    lengths in \emph{Saccharomyces cerevisiae}. 
    \emph{PLoS Computational Biology}, 4(9):e1000175.

    \item    Xi L, Fondufe-Mittendorf Y, Xia L, Flatow J, Widom J 
    and Wang JP (2010). Predicting nucleosome positioning 
    using a duration hidden markov model. 
    \emph{BMC Bioinformatics}, 11:346

    \item
    Brogaard K, Xi L, and Widom J (2012). 
    A map of nucleosome positions in yeast at base-pair 
    resolution. 
    \emph{Nature}, 486(7404):496-501.

    \item
    Moyle-Heyrman G, Zaichuk T, Xi L, Zhang Q, 
    Uhlenbeck OC, Holmgren R, Widom J and Wang JP (2013). 
    Chemical map of \emph{Schizosaccharomyces pombe} reveals 
    species-specific features in nucleosome positioning. 
    \emph{Proc. Natl. Acad. Sci. U. S. A.}, 110(50):20158-63.

    \item
    Ichikawa Y, Morohoshi K, Nishimura Y, Kurumizaka H and 
    Shimizu M (2014). Telomeric repeats act as 
    nucleosome-disfavouring sequences in vivo. 
    \emph{Nucleic Acids Res.}, 42(3):1541-1552.

    \item
    Voong LN, Xi L, Sebeson AC, Xiong B, Wang JP and 
    Wang X (2016). Insights into Nucleosome Organization 
    in Mouse Embryonic Stem Cells through Chemical Mapping. 
    \emph{Cell}, 167(6):1555-1570.

    \item
    Fuse T, Katsumata K, Morohoshi K, Mukai Y, Ichikawa Y, 
    Kurumizaka H, Yanagida A, Urano T, Kato H, and 
    Shimizu M (2017). Parallel mapping with site-directed 
    hydroxyl radicals and micrococcal nuclease reveals 
    structural features of positioned nucleosomes in vivo. 
    \emph{Plos One}, 12(10):e0186974. 
}
}

\examples{
load(system.file("extdata","inseq.RData",package="nuCpos"))
HBA(inseq, species = "sc")
localHBA(inseq, species = "sc")
}
