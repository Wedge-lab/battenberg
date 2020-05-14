fit.copy.number_cell_line = function (samplename, outputfile.prefix, inputfile.baf.segmented, 
    inputfile.baf, inputfile.logr, dist_choice, ascat_dist_choice, 
    min.ploidy = 1.6, max.ploidy = 4.8, min.rho = 0.1, max.rho = 1, 
    min.goodness = 63, uninformative_BAF_threshold = 0.51, gamma_param = 1, 
    use_preset_rho_psi = F, preset_rho = NA, preset_psi = NA, 
    read_depth = 30) 
{
    source(paste0(Ref_files_dir, "runASCAT_pure.R"))
    environment(runASCAT_pure) <- environment(runASCAT)
    assert.file.exists(inputfile.baf.segmented)
    assert.file.exists(inputfile.baf)
    assert.file.exists(inputfile.logr)
    if ((max.ploidy - min.ploidy) < 0.05) {
        stop(paste("Supplied ploidy range must be larger than 0.05: ", 
            min.ploidy, "-", max.ploidy, sep = ""))
    }
    if ((max.rho - min.rho) < 0.01) {
        stop(paste("Supplied rho range must be larger than 0.01: ", 
            min.rho, "-", max.rho, sep = ""))
    }
    segmented.BAF.data = read.table(inputfile.baf.segmented, 
        header = T, stringsAsFactors = F)
    raw.BAF.data = as.data.frame(read_table_generic(inputfile.baf))
    raw.logR.data = as.data.frame(read_table_generic(inputfile.logr))
    identifiers = paste(segmented.BAF.data[, 1], segmented.BAF.data[, 
        2], sep = "_")
    dups = which(duplicated(identifiers))
    if (length(dups) > 0) {
        segmented.BAF.data = segmented.BAF.data[-dups, ]
        identifiers = identifiers[-dups]
    }
    rownames(segmented.BAF.data) = identifiers
    raw.BAF.data = raw.BAF.data[!is.na(raw.BAF.data[, 3]), ]
    raw.logR.data = raw.logR.data[!is.na(raw.logR.data[, 3]), 
        ]
    BAF.data = NULL
    logR.data = NULL
    segmented.logR.data = NULL
    matched.segmented.BAF.data = NULL
    chr.names = unique(segmented.BAF.data[, 1])
    for (chr in chr.names) {
        chr.BAF.data = raw.BAF.data[raw.BAF.data$Chromosome == 
            chr, ]
        if (nrow(chr.BAF.data) == 0) {
            next
        }
        chr.segmented.BAF.data = segmented.BAF.data[segmented.BAF.data[, 
            1] == chr, ]
        indices = match(chr.segmented.BAF.data[, 2], chr.BAF.data$Position)
        if (sum(is.na(indices)) == length(indices) | length(indices) == 
            0) {
            next
        }
        chr.segmented.BAF.data = chr.segmented.BAF.data[!is.na(indices), 
            ]
        matched.segmented.BAF.data = rbind(matched.segmented.BAF.data, 
            chr.segmented.BAF.data)
        BAF.data = rbind(BAF.data, chr.BAF.data[indices[!is.na(indices)], 
            ])
        chr.logR.data = raw.logR.data[raw.logR.data$Chromosome == 
            chr, ]
        indices = match(chr.segmented.BAF.data[, 2], chr.logR.data$Position)
        logR.data = rbind(logR.data, chr.logR.data[indices[!is.na(indices)], 
            ])
        chr.segmented.logR.data = chr.logR.data[indices[!is.na(indices)], 
            ]
        segs = rle(chr.segmented.BAF.data[, 5])$lengths
        cum.segs = c(0, cumsum(segs))
        for (s in 1:length(segs)) {
            chr.segmented.logR.data[(cum.segs[s] + 1):cum.segs[s + 
                1], 3] = mean(chr.segmented.logR.data[(cum.segs[s] + 
                1):cum.segs[s + 1], 3], na.rm = T)
        }
        segmented.logR.data = rbind(segmented.logR.data, chr.segmented.logR.data)
    }
    names(matched.segmented.BAF.data)[5] = samplename
    selection = c()
    for (chrom in chr.names) {
        matched.segmented.BAF.data.chr = matched.segmented.BAF.data[matched.segmented.BAF.data[, 
            1] == chrom, ]
        logR.data.chr = logR.data[logR.data[, 1] == chrom, ]
        selection = c(selection, matched.segmented.BAF.data.chr[, 
            2] %in% logR.data.chr[, 2])
    }
    matched.segmented.BAF.data = matched.segmented.BAF.data[selection, 
        ]
    segmented.logR.data = segmented.logR.data[selection, ]
    row.names(segmented.logR.data) = row.names(matched.segmented.BAF.data)
    row.names(logR.data) = row.names(matched.segmented.BAF.data)
    write.table(segmented.logR.data, paste(samplename, ".logRsegmented.txt", 
        sep = ""), sep = "\t", quote = F, col.names = F, row.names = F)
    segBAF = 1 - matched.segmented.BAF.data[, 5]
    segLogR = segmented.logR.data[, 3]
    logR = logR.data[, 3]
    names(segBAF) = rownames(matched.segmented.BAF.data)
    names(segLogR) = rownames(matched.segmented.BAF.data)
    names(logR) = rownames(matched.segmented.BAF.data)
    chr.segs = NULL
    for (ch in 1:length(chr.names)) {
        chr.segs[[ch]] = which(logR.data[, 1] == chr.names[ch])
    }
    if (use_preset_rho_psi) {
        ascat_optimum_pair = list(rho = preset_rho, psi = preset_psi, 
            ploidy = preset_psi)
    }
    else {
        distance.outfile = paste(outputfile.prefix, "distance.png", 
            sep = "", collapse = "")
        copynumberprofile.outfile = paste(outputfile.prefix, 
            "copynumberprofile.png", sep = "", collapse = "")
        nonroundedprofile.outfile = paste(outputfile.prefix, 
            "nonroundedprofile.png", sep = "", collapse = "")
        cnaStatusFile = paste(outputfile.prefix, "copynumber_solution_status.txt", 
            sep = "", collapse = "")
        ascat_optimum_pair = runASCAT_pure(logR, 1 - BAF.data[, 3], 
            segLogR, segBAF, chr.segs, ascat_dist_choice, distance.outfile, 
            copynumberprofile.outfile, nonroundedprofile.outfile, 
            cnaStatusFile = cnaStatusFile, gamma = gamma_param, 
            allow100percent = T, reliabilityFile = NA, min.ploidy = min.ploidy, 
            max.ploidy = max.ploidy, min.rho = min.rho, max.rho = max.rho, 
            min.goodness, chr.names = chr.names)
    }
    distance.outfile = paste(outputfile.prefix, "second_distance.png", 
        sep = "", collapse = "")
    copynumberprofile.outfile = paste(outputfile.prefix, "second_copynumberprofile.png", 
        sep = "", collapse = "")
    nonroundedprofile.outfile = paste(outputfile.prefix, "second_nonroundedprofile.png", 
        sep = "", collapse = "")
    out = run_clonal_ASCAT(logR, 1 - BAF.data[, 3], segLogR, 
        segBAF, chr.segs, matched.segmented.BAF.data, ascat_optimum_pair, 
        dist_choice, distance.outfile, copynumberprofile.outfile, 
        nonroundedprofile.outfile, gamma_param = gamma_param, 
        read_depth, uninformative_BAF_threshold, allow100percent = T, 
        reliabilityFile = NA, psi_min_initial = min.ploidy, psi_max_initial = max.ploidy, 
        rho_min_initial = min.rho, rho_max_initial = max.rho, 
        chr.names = chr.names)
    ascat_optimum_pair_fraction_of_genome = out$output_optimum_pair_without_ref
    ascat_optimum_pair_ref_seg = out$output_optimum_pair
    is.ref.better = out$is.ref.better
    rho_psi_output = data.frame(rho = c(ascat_optimum_pair$rho, 
        ascat_optimum_pair_fraction_of_genome$rho, ascat_optimum_pair_ref_seg$rho), 
        psi = c(ascat_optimum_pair$psi, ascat_optimum_pair_fraction_of_genome$psi, 
            ascat_optimum_pair_ref_seg$psi), ploidy = c(ascat_optimum_pair$ploidy, 
            ascat_optimum_pair_fraction_of_genome$ploidy, ascat_optimum_pair_ref_seg$ploidy), 
        distance = c(NA, out$distance_without_ref, out$distance), 
        is.best = c(NA, !is.ref.better, is.ref.better), row.names = c("ASCAT", 
            "FRAC_GENOME", "REF_SEG"))
    write.table(rho_psi_output, paste(outputfile.prefix, "rho_and_psi.txt", 
        sep = ""), quote = F, sep = "\t")
}
