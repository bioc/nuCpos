HBA  <- function(inseq, species = "mm", silent = FALSE){

    if(silent == FALSE) message("species: ", species, "\n")
    inseqList <- getWC(inseq = inseq, silent = silent)
    w <- inseqList$w
    c <- inseqList$c

    if(species == "sc") profiles <- nature11142_s2.147.HBA.prof
    if(species == "sp") profiles <- sd01.147.HBA.prof
    if(species == "mm") profiles <- chem.mm9.HBA.prof
    freqN4 <- profiles$freqN4
    freqL4 <- profiles$freqL4
    tranN4 <- profiles$tranN4
    tranL4 <- profiles$tranL4

    t <- 74
    z <- 147
    asc <- freqN4[w[t-73], w[t-72], w[t-71], w[t-70]] /
            freqL4[w[t-73], w[t-72], w[t-71], w[t-70]] *
            freqN4[c[z-t-72], c[z-t-71], c[z-t-70], c[z-t-69]] /
            freqL4[c[z-t-72], c[z-t-71], c[z-t-70], c[z-t-69]]
    for(i in 5:147){
        asc <- asc * 
                tranN4[i, w[t-78+i], w[t-77+i], w[t-76+i], w[t-75+i], w[t-74+i]] /
                tranL4[w[t-78+i], w[t-77+i], w[t-76+i], w[t-75+i], w[t-74+i]] *
                tranN4[i, c[z-t-77+i], c[z-t-76+i], c[z-t-75+i], c[z-t-74+i], c[z-t-73+i]] /
                tranL4[c[z-t-77+i], c[z-t-76+i], c[z-t-75+i], c[z-t-74+i], c[z-t-73+i]]
    }
    logasc <- log(asc)
    names(logasc) <- "HBA"
    return(logasc)

}
