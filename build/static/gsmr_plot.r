# ************************************************** #
#              Read exposure and outcome             #
# ************************************************** #
read_gsmr_trait = function(file_con) {
    expo_str = scan(file_con, nlines=1, quiet=TRUE, what="");
    outcome_str = scan(file_con, nlines=1, quiet=TRUE, what="");
    strbuf = scan(file_con, nlines=1, quiet=TRUE, what="");
    return(list(expo_str=expo_str, outcome_str=outcome_str))
}

# ************************************************** #
#                  Read GSMR result                  #
# ************************************************** #
read_gsmr_result = function(file_con) {
    expo_str = outcome_str = bxy = bxy_se = bxy_pval = bxy_m = c()
    while(1) {
        strbuf = scan(file_con, nlines=1, quiet=TRUE, what="");
        if(strbuf[1] == "#gsmr_end") break;
        if(strbuf[1] == "Exposure") next;
        expo_str = c(expo_str, strbuf[1]);
        outcome_str = c(outcome_str, strbuf[2]);
        bxy = c(bxy, as.numeric(strbuf[3]));
        bxy_se = c(bxy_se, as.numeric(strbuf[4]));
        bxy_pval = c(bxy_pval, as.numeric(strbuf[5]));
        bxy_m = c(bxy_m, as.numeric(strbuf[6]));
    }
    return(cbind(expo_str, outcome_str, bxy, bxy_se, bxy_pval, bxy_m))
}

# ************************************************** #
#                  Read SNP effects                  #
# ************************************************** #
read_snp_effect = function(file_con) {
    snp_effect = c()
    while(1) {
        strbuf = scan(file_con, nlines=1, quiet=TRUE, what="");
        if(strbuf[1] == "#effect_end") break;
        snp_effect = rbind(snp_effect, strbuf);
    }
    return(snp_effect)
}

# ************************************************** #
#                  Read SNP instruments              #
# ************************************************** #
read_snp_instru = function(file_con, snplist, nexpo, noutcome) {
    nrow = length(snplist); ncol = nexpo+noutcome
    snp_instru = matrix(NA, nrow, ncol)
    while(1) {
        strbuf = scan(file_con, nlines=1, quiet=TRUE, what="");
        if(strbuf[1] == "#marker_end") break;
        expo_indx = as.numeric(strbuf[1]); outcome_indx = as.numeric(strbuf[2]);
        forward_flag = TRUE;
        if(expo_indx < outcome_indx) {
            outcome_indx = outcome_indx - nexpo
        } else {
            expo_indx = expo_indx - nexpo
            forward_flag = FALSE;
        }
        snpbuf = scan(file_con, nlines=1, quiet=TRUE, what="");
        snp_indx = match(snpbuf, snplist)
        posbuf = rep(0, nrow); posbuf[snp_indx] = 1;
        indxbuf = expo_indx
        if(!forward_flag) indxbuf = indxbuf + nexpo
        if(length(which(!is.na(snp_instru[,indxbuf])))==0) {
            snp_instru[,indxbuf] = posbuf;
        } else {
            snp_instru[,indxbuf] = paste(snp_instru[,indxbuf], posbuf, sep="")
        }
    }
    return(snp_instru)
}

# ************************************************** #
#          Read output by GCTA-GSMR for plot         #
# ************************************************** #
read_gsmr_data = function(gsmr_effect_file) {
    trait_flag = gsmr_flag = marker_flag = effect_flag = FALSE;
    file_con = file(gsmr_effect_file, "r")
    while(1) {
        strbuf = scan(file_con, nlines=1, quiet=TRUE, what="");
        if(strbuf == "#trait_begin") {
            # Read the exposures and outcomes
            resbuf = read_gsmr_trait(file_con);
            expo_str = resbuf$expo_str; 
            outcome_str = resbuf$outcome_str;
            pheno_str = c(expo_str, outcome_str);
            nexpo = length(expo_str); noutcome = length(outcome_str)
            trait_flag = TRUE;
        } else if(strbuf == "#gsmr_begin") {
            # Read the GSMR result
            bxy_result = read_gsmr_result(file_con);
            colnames(bxy_result) = c("Exposure", "Outcome", "bxy", "se", "p", "n_snps")
            gsmr_flag = TRUE;
        } else if(strbuf == "#effect_begin") {
            # Read the summary statistics
            snp_effect = read_snp_effect(file_con);
            snplist = as.character(snp_effect[,1])
            effect_flag = TRUE;
        } else if(strbuf == "#marker_begin") {
            # Read the SNPs
            snp_instru = read_snp_instru(file_con, snplist, nexpo, noutcome);
            snp_effect = cbind(snp_effect[,1], snp_instru, snp_effect[,-1])
            marker_flag = TRUE;
        }
        if(trait_flag==T & gsmr_flag==T & marker_flag==T & effect_flag==T) break;
    }
    return(list(pheno=c(nexpo, noutcome, pheno_str), bxy_result=bxy_result, snp_effect = snp_effect))
}

# ************************************************** #
#         Display summary of the gsmr data           #
# ************************************************** #
gsmr_summary = function(gsmr_data) {
    message("\n## Exposure and outcome")
    pheno_str = gsmr_data$pheno[c(-1,-2)]
    # exposure
    nexpo = as.numeric(gsmr_data$pheno[1]);
    noutcome = as.numeric(gsmr_data$pheno[2]);
    logger_m = paste(nexpo, "exposure(s):");
    logger_m = paste(logger_m, gsmr_data$pheno[3])
    if(nexpo > 1) {
        for(i in 2 : nexpo) 
            logger_m = paste(logger_m, gsmr_data$pheno[i+2], sep=", ")
    } 
    message(logger_m)
    # outcome
    logger_m = paste(noutcome, "outcome(s):");
    logger_m = paste(logger_m, gsmr_data$pheno[3+nexpo])
    if(noutcome > 1) {
        for(i in 2 : noutcome) 
            logger_m = paste(logger_m, gsmr_data$pheno[i+2+nexpo], sep=", ")
    } 
    message(logger_m)

    message("\n## GSMR result")
    m_bxy_rst = data.frame(gsmr_data$bxy_result)
    print(m_bxy_rst)
}

# ************************************************** #
#                  Plot bzy vs bzx                   #
# ************************************************** #
plot_snp_effect = function(expo_str, outcome_str, bxy, bzx, bzx_se, bzy, bzy_se, effect_col=colors()[75]) {
    vals = c(bzx-bzx_se, bzx+bzx_se)
    xmin = min(vals); xmax = max(vals)
    vals = c(bzy-bzy_se, bzy+bzy_se)
    ymin = min(vals); ymax = max(vals)
    plot(bzx, bzy, pch=20, cex=0.8, bty="n", cex.axis=1.1, cex.lab=1.2,
         col=effect_col, xlim=c(xmin, xmax), ylim=c(ymin, ymax),
         xlab=substitute(paste(trait, " (", italic(b[zx]), ")", sep=""), list(trait=expo_str)),
         ylab=substitute(paste(trait, " (", italic(b[zy]), ")", sep=""), list(trait=outcome_str)))
    abline(0, bxy, lwd=1.5, lty=2, col="dim grey")
    ## Standard errors
    nsnps = length(bzx)
    for( i in 1:nsnps ) {
        # x axis
        xstart = bzx[i] - bzx_se[i]; xend = bzx[i] + bzx_se[i]
        ystart = bzy[i]; yend = bzy[i]
        segments(xstart, ystart, xend, yend, lwd=1.5, col=effect_col)
        # y axis
        xstart = bzx[i]; xend = bzx[i] 
        ystart = bzy[i] - bzy_se[i]; yend = bzy[i] + bzy_se[i]
        segments(xstart, ystart, xend, yend, lwd=1.5, col=effect_col)
    }
}

# ************************************************** #
#                  Effect size plot                  #
# ************************************************** #
# expo_str, exposure
# outcome_str, outcome
# effect_col, plotting colour
plot_gsmr_effect = function(gsmr_data, expo_str, outcome_str, effect_col=colors()[75]) {
    # index of SNP instruments
    pheno_str = as.character(gsmr_data$pheno[c(-1,-2)])
    nexpo = as.numeric(gsmr_data$pheno[1])
    noutcome = as.numeric(gsmr_data$pheno[2])
    expo_indx = match(expo_str, pheno_str)
    if(is.na(expo_indx)) stop("\"", expo_str, "\" is not found.")
    outcome_indx = match(outcome_str, pheno_str)
    if(is.na(outcome_indx)) stop("\"", outcome_str, "\" is not found.")
    forward_flag = TRUE;
    if(expo_indx > outcome_indx) forward_flag = FALSE;
    if(forward_flag) {
        outcome_indx = outcome_indx - nexpo;
    } else {
        expo_indx = expo_indx - nexpo;
    }
    indxbuf = expo_indx + 1
    if(!forward_flag) indxbuf = indxbuf + nexpo
    strbuf = as.character(substr(gsmr_data$snp_effect[,indxbuf], outcome_indx, outcome_indx))
    snpindx = which(strbuf=="1")
    if(length(snpindx) < 3) stop("Not enough SNPs retained for plot.")
    # bxy
    indxbuf = which(gsmr_data$bxy_result[,1]==expo_str & gsmr_data$bxy_result[,2]==outcome_str)
    bxy = as.numeric(gsmr_data$bxy_result[indxbuf, 3])
    # SNP effects
    if(forward_flag) {
        indxbuf1 = 1 + nexpo + noutcome + 3 + expo_indx
        indxbuf2 = 1 + nexpo + noutcome + 3 + nexpo*2 + outcome_indx
    } else {
        indxbuf1 = 1 + nexpo + noutcome + 3 + nexpo*2 + expo_indx
        indxbuf2 = 1 + nexpo + noutcome + 3 + expo_indx
    }
    bzx = as.numeric(gsmr_data$snp_effect[snpindx,indxbuf1]); indxbuf1 = indxbuf1 + 1;
    bzx_se = as.numeric(gsmr_data$snp_effect[snpindx,indxbuf1]);
    bzy = as.numeric(gsmr_data$snp_effect[snpindx,indxbuf2]); indxbuf2 = indxbuf2 + 1;
    bzy_se = as.numeric(gsmr_data$snp_effect[snpindx,indxbuf2]);
    # plot
    plot_snp_effect(expo_str, outcome_str, bxy, bzx, bzx_se, bzy, bzy_se, effect_col)
}