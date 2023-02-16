getWC <- function(inseq, silent = FALSE){

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

    out <- list()
    out$w <- w
    out$c <- c
    return(out)
}
