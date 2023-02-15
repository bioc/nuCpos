HBA  <- function(inseq, species = "mm", silent = FALSE){

    if(silent == FALSE) message("species: ", species, "\n")
    # if(!is(inseq)[1] == "character"){
    #     if(is(inseq)[1] == "DNAString"){
    if(is(inseq, "character")[1] == FALSE){
        if(is(inseq, 'DNAString')[1] == TRUE){
            if(requireNamespace("Biostrings", quietly = TRUE)){
                inseq <- as.character(inseq)
                if(silent == FALSE){
                message(
                "The class of inseq was changed from DNAString to character")
                }
            }else{
            message("DNAString cannot be changed to a character string")
            out <- NA; names(out) <- "HBA"; return(out)
            }
        }else{
            message("The class of inseq must be DNAString or character")
            out <- NA; names(out) <- "HBA"; return(out)
        }
    }
    if(nchar(inseq) != 147){
        message("Length of inseq: ", nchar(inseq), "bp")
        message("The length of inseq must be 147 bp")
        out <- NA; names(out) <- "HBA"; return(out)
    }

    if(species == "sc") profiles <- nature11142_s2.147.HBA.prof
    if(species == "sp") profiles <- sd01.147.HBA.prof
    if(species == "mm") profiles <- chem.mm9.HBA.prof
    freqN4 <- profiles$freqN4
    freqL4 <- profiles$freqL4
    tranN4 <- profiles$tranN4
    tranL4 <- profiles$tranL4

    inseqW <- strsplit(inseq, split = "")[[1]]
    inseqW[inseqW == "A"] <- 1
    inseqW[inseqW == "C"] <- 2
    inseqW[inseqW == "G"] <- 3
    inseqW[inseqW == "T"] <- 4
    w <- as.integer(inseqW)

    inseqRC <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(inseq)))
    inseqC <- strsplit(inseqRC, split = "")[[1]]
    inseqC[inseqC == "A"] <- 1
    inseqC[inseqC == "C"] <- 2
    inseqC[inseqC == "G"] <- 3
    inseqC[inseqC == "T"] <- 4
    c <- as.integer(inseqC)

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
