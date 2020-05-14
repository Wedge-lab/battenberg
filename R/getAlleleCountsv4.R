getAlleleCountsv4=function (bam.file, output.file, g1000.loci, min.base.qual = 20, 
    min.map.qual = 35, allelecounter.exe = "alleleCounterv4") 
{
    cmd = paste(allelecounter.exe, "-b", bam.file, "-l", g1000.loci, 
        "-o", output.file, "-m", min.base.qual, "-q", min.map.qual, "-d")
    system(cmd, wait = T)
}

