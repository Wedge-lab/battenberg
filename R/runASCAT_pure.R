runASCAT_cell_line=function (lrr, baf, lrrsegmented, bafsegmented, chromosomes, 
    dist_choice, distancepng = NA, copynumberprofilespng = NA, 
    nonroundedprofilepng = NA, cnaStatusFile = "copynumber_solution_status.txt", 
    gamma = 0.55, allow100percent, reliabilityFile = NA, min.ploidy = 1.6, 
    max.ploidy = 4.8, min.rho = 0.1, max.rho = 1, min.goodness = 63, 
    uninformative_BAF_threshold = 0.51, chr.names) 
{
    ch = chromosomes
    b = bafsegmented
    r = lrrsegmented[names(bafsegmented)]
    dist_min_psi = max(min.ploidy - 0.6, 0)
    dist_max_psi = max.ploidy + 0.6
    dist_min_rho = max(min.rho - 0.03, 0.05)
    dist_max_rho = max.rho + 0.03
    s = ASCAT::make_segments(r, b)
    dist_matrix_info <- create_distance_matrix(s, dist_choice, 
        gamma, uninformative_BAF_threshold = uninformative_BAF_threshold, 
        min_psi = dist_min_psi, max_psi = dist_max_psi, min_rho = dist_min_rho, 
        max_rho = dist_max_rho)
    d = dist_matrix_info$distance_matrix
    minimise = dist_matrix_info$minimise
    TheoretMaxdist = sum(rep(0.25, dim(s)[1]) * s[, "length"], 
        na.rm = T)
    if (!(minimise)) {
        d = -d
    }
    nropt = 0
    localmin = NULL
    optima = list()
    for (i in 4:(dim(d)[1] - 3)) {
        for (j in 4:(dim(d)[2] - 3)) {
            m = d[i, j]
            seld = d[(i - 3):(i + 3), (j - 3):(j + 3)]
            seld[4, 4] = max(seld)
            if (min(seld) > m) {
                psi = as.numeric(rownames(d)[i])
                rho = as.numeric(colnames(d)[j])
                nA = (rho - 1 - (s[, "b"] - 1) * 2^(s[, "r"]/gamma) * 
                  ((1 - rho) * 2 + rho * psi))/rho
                nB = (rho - 1 + s[, "b"] * 2^(s[, "r"]/gamma) * 
                  ((1 - rho) * 2 + rho * psi))/rho
                ploidy = sum((nA + nB) * s[, "length"])/sum(s[, 
                  "length"])
                ploidy_opt1 = ploidy
                percentzero = (sum((round(nA) == 0) * s[, "length"]) + 
                  sum((round(nB) == 0) * s[, "length"]))/sum(s[, 
                  "length"])
                perczeroAbb = (sum((round(nA) == 0) * s[, "length"] * 
                  ifelse(s[, "b"] == 0.5, 0, 1)) + sum((round(nB) == 
                  0) * s[, "length"] * ifelse(s[, "b"] == 0.5, 
                  0, 1)))/sum(s[, "length"] * ifelse(s[, "b"] == 
                  0.5, 0, 1))
                if (is.na(perczeroAbb)) {
                  perczeroAbb = 0
                }
                if (minimise) {
                  goodnessOfFit = (1 - m/TheoretMaxdist) * 100
                }
                else {
                  goodnessOfFit = -m/TheoretMaxdist * 100
                }
                print(paste("ploidy=", ploidy, ",rho=", rho, 
                  ",goodness=", goodnessOfFit, ",percentzero=", 
                  percentzero, ", perczerAbb=", perczeroAbb, 
                  sep = ""))
                if (ploidy >= min.ploidy & ploidy <= max.ploidy & 
                  rho >= min.rho & goodnessOfFit >= min.goodness & 
                  (percentzero > 0.01 | perczeroAbb > 0.1)) {
                  nropt = nropt + 1
                  optima[[nropt]] = c(m, i, j, ploidy, goodnessOfFit)
                  localmin[nropt] = m
                }
            }
        }
    }
    if (allow100percent & nropt == 0) {
        cold = which(as.numeric(colnames(d)) > 1)
        d[, cold] = 1e+20
        for (i in 4:(dim(d)[1] - 3)) {
            for (j in 4:(dim(d)[2] - 3)) {
                m = d[i, j]
                seld = d[(i - 3):(i + 3), (j - 3):(j + 3)]
                seld[4, 4] = max(seld)
                if (min(seld) > m) {
                  psi = as.numeric(rownames(d)[i])
                  rho = as.numeric(colnames(d)[j])
                  nA = (rho - 1 - (s[, "b"] - 1) * 2^(s[, "r"]/gamma) * 
                    ((1 - rho) * 2 + rho * psi))/rho
                  nB = (rho - 1 + s[, "b"] * 2^(s[, "r"]/gamma) * 
                    ((1 - rho) * 2 + rho * psi))/rho
                  ploidy = sum((nA + nB) * s[, "length"])/sum(s[, 
                    "length"])
                  percentzero = (sum((round(nA) == 0) * s[, "length"]) + 
                    sum((round(nB) == 0) * s[, "length"]))/sum(s[, 
                    "length"])
                  perczeroAbb = (sum((round(nA) == 0) * s[, "length"] * 
                    ifelse(s[, "b"] == 0.5, 0, 1)) + sum((round(nB) == 
                    0) * s[, "length"] * ifelse(s[, "b"] == 0.5, 
                    0, 1)))/sum(s[, "length"] * ifelse(s[, "b"] == 
                    0.5, 0, 1))
                  if (is.na(perczeroAbb)) {
                    perczeroAbb = 0
                  }
                  if (minimise) {
                    goodnessOfFit = (1 - m/TheoretMaxdist) * 
                      100
                  }
                  else {
                    goodnessOfFit = -m/TheoretMaxdist * 100
                  }
                  if (ploidy > min.ploidy & ploidy < max.ploidy & 
                    rho >= min.rho & goodnessOfFit >= min.goodness) {
                    nropt = nropt + 1
                    optima[[nropt]] = c(m, i, j, ploidy, goodnessOfFit)
                    localmin[nropt] = m
                  }
                }
            }
        }
    }
    psi_opt1_plot = vector(mode = "numeric")
    rho_opt1_plot = vector(mode = "numeric")
    if (nropt > 0) {
        write.table(paste(nropt, " copy number solutions found", 
            sep = ""), file = cnaStatusFile, quote = F, col.names = F, 
            row.names = F)
        optlim = sort(localmin)[1]
        for (i in 1:length(optima)) {
            if (optima[[i]][1] == optlim) {
                psi_opt1 = as.numeric(rownames(d)[optima[[i]][2]])
                rho_opt1 = as.numeric(colnames(d)[optima[[i]][3]])
                if (rho_opt1 > 1) {
                  rho_opt1 = 1
                }
                ploidy_opt1 = optima[[i]][4]
                goodnessOfFit_opt1 = optima[[i]][5]
                psi_opt1_plot = c(psi_opt1_plot, psi_opt1)
                rho_opt1_plot = c(rho_opt1_plot, rho_opt1)
            }
        }
    }
    else {
        write.table(paste("no copy number solutions found", sep = ""), 
            file = cnaStatusFile, quote = F, col.names = F, row.names = F)
        print("No suitable copy number solution found")
        psi = NA
        ploidy = NA
        rho = NA
    }
    if (nropt > 0) {
        rho = rho_opt1
        psi = psi_opt1
        ploidy = ploidy_opt1
        nAfull = (rho - 1 - (b - 1) * 2^(r/gamma) * ((1 - rho) * 
            2 + rho * psi))/rho
        nBfull = (rho - 1 + b * 2^(r/gamma) * ((1 - rho) * 2 + 
            rho * psi))/rho
        nA = pmax(round(nAfull), 0)
        nB = pmax(round(nBfull), 0)
        rBacktransform = gamma * log((rho * (nA + nB) + (1 - 
            rho) * 2)/((1 - rho) * 2 + rho * psi), 2)
        bBacktransform = (1 - rho + rho * nB)/(2 - 2 * rho + 
            rho * (nA + nB))
        rConf = ifelse(abs(rBacktransform) > 0.15, pmin(100, 
            pmax(0, 100 * (1 - abs(rBacktransform - r)/abs(r)))), 
            NA)
        bConf = ifelse(bBacktransform != 0.5, pmin(100, pmax(0, 
            ifelse(b == 0.5, 100, 100 * (1 - abs(bBacktransform - 
                b)/abs(b - 0.5))))), NA)
        if (!is.na(reliabilityFile)) {
            write.table(data.frame(segmentedBAF = b, backTransformedBAF = bBacktransform, 
                confidenceBAF = bConf, segmentedR = r, backTransformedR = rBacktransform, 
                confidenceR = rConf, nA = nA, nB = nB, nAfull = nAfull, 
                nBfull = nBfull), reliabilityFile, sep = ",", 
                row.names = F)
        }
        confidence = ifelse(is.na(rConf), bConf, ifelse(is.na(bConf), 
            rConf, (rConf + bConf)/2))
        if (!is.na(copynumberprofilespng)) {
            png(filename = copynumberprofilespng, width = 2000, 
                height = 500, res = 200)
        }
        ASCAT::ascat.plotAscatProfile(n1all = nA, n2all = nB, 
            heteroprobes = TRUE, ploidy = ploidy_opt1, rho = rho_opt1, 
            goodnessOfFit = goodnessOfFit_opt1, nonaberrant = FALSE, 
            ch = ch, lrr = lrr, bafsegmented = bafsegmented, 
            chrs = chr.names)
        if (!is.na(copynumberprofilespng)) {
            dev.off()
        }
        if (!is.na(nonroundedprofilepng)) {
            png(filename = nonroundedprofilepng, width = 2000, 
                height = 500, res = 200)
        }
        ASCAT::ascat.plotNonRounded(ploidy = ploidy_opt1, rho = rho_opt1, 
            goodnessOfFit = goodnessOfFit_opt1, nonaberrant = FALSE, 
            nAfull = nAfull, nBfull = nBfull, bafsegmented = bafsegmented, 
            ch = ch, lrr = lrr, chrs = chr.names)
        if (!is.na(nonroundedprofilepng)) {
            dev.off()
        }
    }
    output_optimum_pair = list(psi = psi, rho = rho, ploidy = ploidy)
    return(output_optimum_pair)
}
