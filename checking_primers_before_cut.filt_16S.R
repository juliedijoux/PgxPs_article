library("dada2")
library("ShortRead")
library("Biostrings")

# Path to fastq.gz files
path <- "/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S"
list.files(path)

# Distinction between foward and reverse samples
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))

# Identification of primers
FWD <- "CCTACGGGNGGCWGCAG" # 341F
REV <- "GACTACHVGGGTATCTAATCC" # 805R

# Checking the presence and orientation of primers in sequences
allOrients <- function(primer) {
    require(Biostrings)
    dna <- DNAString(primer)
    orients <- c(Foward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna), RevComp = Biostrings::reverseComplement(dna))
    return(sapply(orients, toString))
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

# Occurrence of primers in sequences
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}
table_bef_cut <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs), FWD.ReverseReads = sapply(FWD.orients,
    primerHits, fn = fnRs), REV.ForwardReads = sapply(REV.orients, primerHits,
    fn = fnFs), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs))

write.table(table_bef_cut, file = "/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_tab/table_before_cut.filt_16S.csv", sep = ";", col.names = NA)