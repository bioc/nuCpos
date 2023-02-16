getLHBA <- function(fragment = "A", w, c, profiles = profiles){

    freqN4SA <- profiles$freqN4SA
    freqN4SB <- profiles$freqN4SB
    freqN4SC <- profiles$freqN4SC
    freqN4SD <- profiles$freqN4SD
    freqN4SE <- profiles$freqN4SE
    freqN4SF <- profiles$freqN4SF
    freqN4SG <- profiles$freqN4SG
    freqN4SH <- profiles$freqN4SH
    freqN4SI <- profiles$freqN4SI
    freqN4SJ <- profiles$freqN4SJ
    freqN4SK <- profiles$freqN4SK
    freqN4SL <- profiles$freqN4SL
    freqN4SM <- profiles$freqN4SM
    tranN4 <- profiles$tranN4
    tranL4 <- profiles$tranL4
    freqL4 <- profiles$freqL4

    if(fragment == "A"){
        fragmentY <- "M"; w_start <- 1; c_start <- 127
        iw_start <- 5; iw_end <- 21; ic_start <- 131; ic_end <- 147
    }
    if(fragment == "C"){
        fragmentY <- "K"; w_start <- 22; c_start <- 106
        iw_start <- 26; iw_end <- 42; ic_start <- 110; ic_end <- 126
    }
    if(fragment == "E"){
        fragmentY <- "I"; w_start <- 43; c_start <- 85
        iw_start <- 47; iw_end <- 63; ic_start <- 89; ic_end <- 105
    }
    if(fragment == "G"){
        fragmentY <- "G"; w_start <- 64; c_start <- 64
        i_start <- 68; i_end <- 84
    }
    if(fragment == "I"){
        fragmentY <- "E"; w_start <- 85; c_start <- 43
        iw_start <- 89; iw_end <- 105; ic_start <- 47; ic_end <- 63
    }
    if(fragment == "K"){
        fragmentY <- "C"; w_start <- 106; c_start <- 22
        iw_start <- 110; iw_end <- 126; ic_start <- 26; ic_end <- 42
    }
    if(fragment == "M"){
        fragmentY <- "A"; w_start <- 127; c_start <- 1
        iw_start <- 131; iw_end <- 147; ic_start <- 5; ic_end <- 21
    }
    if(fragment == "B"){
        fragmentY <- "L"; w_start <- 12; c_start <- 117
        iw_start <- 16; iw_end <- 31; ic_start <- 121; ic_end <- 136
    }
    if(fragment == "D"){
        fragmentY <- "J"; w_start <- 33; c_start <- 96
        iw_start <- 37; iw_end <- 52; ic_start <- 100; ic_end <- 115
    }
    if(fragment == "F"){
        fragmentY <- "H"; w_start <- 54; c_start <- 75
        iw_start <- 58; iw_end <- 73; ic_start <- 79; ic_end <- 94
    }
    if(fragment == "H"){
        fragmentY <- "F"; w_start <- 75; c_start <- 54
        iw_start <- 79; iw_end <- 94; ic_start <- 58; ic_end <- 73
    }
    if(fragment == "J"){
        fragmentY <- "D"; w_start <- 96; c_start <- 33
        iw_start <- 100; iw_end <- 115; ic_start <- 37; ic_end <- 52
    }
    if(fragment == "L"){
        fragmentY <- "B"; w_start <- 117; c_start <- 12
        iw_start <- 121; iw_end <- 136; ic_start <- 16; ic_end <- 31
    }

    freqN4SxName <- paste("freqN4S", fragment, sep = "")
    freqN4Sx <- get(freqN4SxName)

    freqN4SyName <- paste("freqN4S", fragmentY, sep = "")
    freqN4Sy <- get(freqN4SyName)

    if(fragment != "G"){
        ascSx <- freqN4Sx[w[w_start], w[w_start+1], w[w_start+2], w[w_start+3]] /
                 freqL4[w[w_start], w[w_start+1], w[w_start+2], w[w_start+3]] *
                 freqN4Sy[c[c_start], c[c_start+1], c[c_start+2], c[c_start+3]] /
                 freqL4[c[c_start], c[c_start+1], c[c_start+2], c[c_start+3]]
        for(i in iw_start:iw_end){
            ascSx <- ascSx * 
                      tranN4[i, w[i-4], w[i-3], w[i-2], w[i-1], w[i]] /
                      tranL4[w[i-4], w[i-3], w[i-2], w[i-1], w[i]]
        }
        for(i in ic_start:ic_end){
            ascSx <- ascSx * 
                      tranN4[i, c[i-4], c[i-3], c[i-2], c[i-1], c[i]] /
                      tranL4[c[i-4], c[i-3], c[i-2], c[i-1], c[i]]
        }
    }

    if(fragment == "G"){
        ascSx <- freqN4Sx[w[w_start], w[w_start+1], w[w_start+2], w[w_start+3]] /
                 freqL4[w[w_start], w[w_start+1], w[w_start+2], w[w_start+3]] *
                 freqN4Sy[c[c_start], c[c_start+1], c[c_start+2], c[c_start+3]] /
                 freqL4[c[c_start], c[c_start+1], c[c_start+2], c[c_start+3]]
        for(i in i_start:i_end){
            ascSx <- ascSx * 
                      tranN4[i, w[i-4], w[i-3], w[i-2], w[i-1], w[i]] /
                      tranL4[w[i-4], w[i-3], w[i-2], w[i-1], w[i]] *
                      tranN4[i, c[i-4], c[i-3], c[i-2], c[i-1], c[i]] /
                      tranL4[c[i-4], c[i-3], c[i-2], c[i-1], c[i]]
        }
    }

    return(log(ascSx))
}
