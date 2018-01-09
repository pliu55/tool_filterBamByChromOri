#!/bin/env Rscript

library(data.table)
suppressMessages(library(GenomicAlignments))

main <- function() {
    argv = commandArgs(trailingOnly=T)

    if ( length(argv) != 4 ) {
        stop('filterBamByChromOri.R  [input] [output] [chrom] [strand]\n')
    }
    finbam  = argv[1] ## input Bam file name
    foutbam = argv[2] ## output Bam file name
    chrom   = argv[3] ## 'chr10' ..
    strand  = argv[4] ## 'plus' or 'minus'

    filterBamByChromOri(finbam, foutbam, chrom, strand)
}


filterBamByChromOri <- function(fuserbam, fchromoribam, chrom, strand) {
    max_yield_size = 600000000
    max_chrom_len = 536870912

    fr1ststrand2mate2flag = list(
        'plus' = list(
            '1stmate' = scanBamFlag( isProperPair     = T,
                                     hasUnmappedMate  = F,
                                     isMinusStrand    = T,
                                     isFirstMateRead  = T,
                                     isSecondMateRead = F  ),

            '2ndmate' = scanBamFlag( isProperPair     = T,
                                     hasUnmappedMate  = F,
                                     isMinusStrand    = F,
                                     isFirstMateRead  = F,
                                     isSecondMateRead = T  ) ),

        'minus' = list(
            '1stmate' = scanBamFlag( isProperPair     = T,
                                     hasUnmappedMate  = F,
                                     isMinusStrand    = F,
                                     isFirstMateRead  = T,
                                     isSecondMateRead = F  ),

            '2ndmate' = scanBamFlag( isProperPair     = T,
                                     hasUnmappedMate  = F,
                                     isMinusStrand    = T,
                                     isFirstMateRead  = F,
                                     isSecondMateRead = T  ) ) )

    flag1st = fr1ststrand2mate2flag[[strand]][['1stmate']]

    grs = GRanges(paste0(chrom, ':1-', max_chrom_len))
    mate1st = scanBam(fuserbam, param=ScanBamParam(flag=flag1st, tag=c('HI'),
                                                   what=c('qname'), which=grs))
    seldt = data.table()
    if ( length(mate1st[[1]]$qname) > 0 ) {
        seldt = data.table( qname = mate1st[[1]]$qname,
                            HI    = mate1st[[1]]$tag$HI )
        seldt[, qname_HI := paste0(qname, '_', HI)]
    }

    ## scanBam first for qname flag, and HI, save them and filterBam
    filter_func = function(x) {
        dt = data.table(as.data.frame(x))
        dt[, qname_HI := paste0(qname, '_', HI)]
        dt[, to_select := ifelse(qname_HI %in% seldt$qname_HI, T, F)]
        return(dt$to_select)
    }

    filter_rules = FilterRules(list(tmp=filter_func))
    inbam = BamFile(fuserbam, yieldSize=max_yield_size)
    filterBam(inbam, fchromoribam, filter=filter_rules,
              param=ScanBamParam(what=c('qname'), tag=c('HI'), which=grs))
}

system.time( main() )
