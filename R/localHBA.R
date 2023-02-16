localHBA <- function(inseq, species = "mm", silent = FALSE){

    if(silent == FALSE) message("species: ", species, "\n")
    inseqList <- getWC(inseq = inseq, silent = silent)
    w <- inseqList$w
    c <- inseqList$c

    if(species == "sc") profiles <- nature11142_s2.147.LHBA.prof
    if(species == "sp") profiles <- sd01.147.LHBA.prof
    if(species == "mm") profiles <- chem.mm9.LHBA.prof

    out <- list()
    out$lHBA_A <- getLHBA(fragment = "A", w = w, c = c, profiles = profiles)
    out$lHBA_B <- getLHBA(fragment = "B", w = w, c = c, profiles = profiles)
    out$lHBA_C <- getLHBA(fragment = "C", w = w, c = c, profiles = profiles)
    out$lHBA_D <- getLHBA(fragment = "D", w = w, c = c, profiles = profiles)
    out$lHBA_E <- getLHBA(fragment = "E", w = w, c = c, profiles = profiles)
    out$lHBA_F <- getLHBA(fragment = "F", w = w, c = c, profiles = profiles)
    out$lHBA_G <- getLHBA(fragment = "G", w = w, c = c, profiles = profiles)
    out$lHBA_H <- getLHBA(fragment = "H", w = w, c = c, profiles = profiles)
    out$lHBA_I <- getLHBA(fragment = "I", w = w, c = c, profiles = profiles)
    out$lHBA_J <- getLHBA(fragment = "J", w = w, c = c, profiles = profiles)
    out$lHBA_K <- getLHBA(fragment = "K", w = w, c = c, profiles = profiles)
    out$lHBA_L <- getLHBA(fragment = "L", w = w, c = c, profiles = profiles)
    out$lHBA_M <- getLHBA(fragment = "M", w = w, c = c, profiles = profiles)

    return(unlist(out))
}
