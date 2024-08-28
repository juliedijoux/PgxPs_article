library("dada2")
library("ShortRead")
library("Biostrings")

# Chemin d'accès aux fichier fastq.gz
path <- "/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/ITS"
list.files(path)

# Distinction entre les échantillons foward et reverse
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))

# Identification des amorces

FWD <- "CTTGGTCATTTAGAGGAAGTAA" # ITS 1F
REV <- "TGTGTTCTTCATCGATG" # ITS 2R

# Vérification de la présence et de l'orientation des primers dans les séquences
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

# Nombre de présence des primers
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}
table_av_cut <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs), FWD.ReverseReads = sapply(FWD.orients,
    primerHits, fn = fnRs), REV.ForwardReads = sapply(REV.orients, primerHits,
    fn = fnFs), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs))

write.table(table_av_cut, file = "/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/ITS_tab/table_before_cut.filt_ITS.csv", sep = ";", col.names = NA)

