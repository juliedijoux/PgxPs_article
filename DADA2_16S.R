library(anyLib)
packagesNeeded_dada2 <- c("dada2", "ShortRead", "Biostrings")
anyLib(packagesNeeded_dada2, force = FALSE)

#setwd("/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S/cut_fastq")

########## Part 1 ##########

path <- "/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S/cut_fastq" # variable that tells us where to retrieve the fastq
list.files(path) # viewing fastq

fnFs <- sort(list.files(path, pattern = "_1.fastq")) # list of fastq sense
fnRs <- sort(list.files(path, pattern = "_2.fastq")) # list of fastq antisense

sample.name <- sapply(strsplit(fnFs, "_cut_"), '[', 1) # variable to extract the names of the samples. We assume that the names of fastq files have the format SAMPLENAME_XXX.fastq
sample.name 

fnFs <- file.path(path, fnFs) # specifies the path of the variable
fnRs <- file.path(path, fnRs) # specifies the path of the variable

#svg(file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_figures/plotquali_sample1_sens.svg") # Saves the next plot in svg format
#plotQualityProfile(fnFs[1]) # display of the Qscore associated with each nucleotide in the 1st sense sample
#dev.off() # to delete the plot
#svg(file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_figures/plotquali_sample1_antisens.svg") # Saves the next plot in svg format
#plotQualityProfile(fnRs[1]) # display of the Qscore associated with each nucleotide in the 1st antisense sample
#dev.off() # to delete the plot

pdf(file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_figures/plotquali_agregate_sens.pdf") # Saves the next plot in pdf format
plotQualityProfile(fnFs, aggregate = TRUE) # visualisation of the Qscore associated with each nucleotide for all sense samples
dev.off() # to delete the plot
pdf(file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_figures/plotquali_agregate_antisens.pdf") # Saves the next plot in pdf format
plotQualityProfile(fnRs, aggregate = TRUE) # visualisation of the Qscore associated with each nucleotide for all antisense samples
dev.off() # to delete the plot

save(list = ls(), file = "/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_env/env_entier_16S_dada2_part1.RData") # save the entire environment in the chosen directory

########## Part 2 ##########

#load("/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_env/env_entier_16S_dada2_part1.RData") # Loading the R environment

filt_path <- file.path(path, "filtered_pairedend") # Creation of the filtered_pairedend folder and the filtFs and filtRs objects for storing filtered sequences

filtFs <- file.path(filt_path, paste0(sample.name, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.name, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncQ = 10, # truncQ = minimum score Q index; on the 1st instance of a quality score less than or equal to truncQ, the sequence is truncated
                     truncLen = c(227,205), # truncLen = définit à quelle longueur les séquences vont être tronquées ; les séquences les plus courtes que la longueur sont éliminées.
                     #trimLeft = c(18,20), # trimLeft = defines the length to be cut on the 5' side of the sequences; this allows the primers to be removed if this has not been done beforehand
                     maxEE = 2, # maxEE = the maximum number of ‘expected errors’ allowed in a reading. This filter is based on the Q score index. The higher the number, the less stringent the filter
                     multithread=TRUE) # to use several cores

# Checking and deleting empty FASTQ files
emptyFiles <- c()

for (file in c(fnFs, fnRs)) {
  if (file.size(file) == 0) {
    cat("Removing empty file:", file, "\n")
    file.remove(file)
    emptyFiles <- c(emptyFiles, file)
  }
}

cat("Empty files:", emptyFiles, "\n")

filtFs <- file.path(filt_path, paste0(sample.name, "_F_filt.fastq.gz")) # redefine the paths in case any files have been deleted
filtRs <- file.path(filt_path, paste0(sample.name, "_R_filt.fastq.gz")) # redefine the paths in case any files have been deleted

pourc <- cbind((out[,2]/out[,1])*100) # % of seq filtered/seq unfiltered
pourc_disc <- cbind(out, pourc) # combines out and pourcfiltre
pourc_disc
write.table(pourc_disc, file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_tab/pourc_seqTrimmed_sur_seq_noTrimmed_16S.csv", sep=";")

suivi_filterAndTrim <- (mean(out[,2])/mean(out[,1]))*100 # indicates the average % of sequences that have passed the filter parameters
write.table(pourc_disc, file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_tab/suivi_filterAndTrim_16S.csv", sep=";", col.names = NA)

save(list = ls(), file = "/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_env/env_entier_16S_dada2_part2.RData") # save the entire environment in the chosen directory

########## Part 3 ##########

#load("/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_env/env_entier_16S_dada2_part2.RData") # Loading the R environment

###Error rate###
errF <- learnErrors(filtFs) # Estimation of the sense sequencing error rate
errR <- learnErrors(filtRs) # Estimation of the antisense sequencing error rate

pdf(file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_figures/plotErrors_sens_16S.pdf")
plotErrors(errF, nominalQ = TRUE) # display of sense sequencing error rate
dev.off()
pdf(file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_figures/plotErrors_antisens_16S.pdf")
plotErrors(errR, nominalQ = TRUE) # display of antisense sequencing error rate
dev.off()

###Dereplication###
derepFs <- derepFastq(filtFs) # the identical sense sequences will be grouped together into single sequences to which abundances will be assigned
names(derepFs) <- sample.name # dereplicated sense sequences take on the name of the samples from which they originate

derepRs <- derepFastq(filtRs) # the identical antisense sequences will be grouped together into single sequences to which abundances will be assigned
names(derepRs) <- sample.name # dereplicated antisense sequences take on the name of the samples from which they originate

###Sample Inference###
dadaFs <- dada(derepFs,
               err = errF,
               multithread=TRUE,
               pool = TRUE)
dadaRs <- dada(derepRs,
               err = errR,
               multithread=TRUE,
               pool = TRUE)

dadaFs[[1]] # displays data relating to sense sample no. 1 (in alphabetical order, sample.name)
dadaRs[[1]] # displays data relating to antisense sample no. 1 (in alphabetical order, sample.name)

save(list = ls(), file = "/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_env/env_entier_16S_dada2_part3.RData") # save the entire environment in the chosen directory

########## Part 4 ##########

#load("/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_env/env_entier_16S_dada2_part3.RData")

###Merging Sequences###
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, 
                      minOverlap = 12, # minOverlap = defines the minimum size of overlap between the sense and antisense strands for their fusion to be accepted; sequences that do not fuse are eliminated
                      maxMismatch = 0) # maxMismatch = defines the maximum number of nucleotide incompatibilities in the overlap

head(mergers[[1]]) # visualisation of the merger
max(mergers[[1]]$nmatch) # Size of largest overlap
min(mergers[[1]]$nmatch) # Size of smallest overlap

###ASV table###
seqtab <- makeSequenceTable(mergers) # to store ASVs in the seqtab object
dim(seqtab)

seqtab[,1] # indicates the number of times the sequence of the first ASV is found in each sample
seqtab
save(seqtab, file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_tab/seqtab_16S.rdata")

pdf(file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_figures/hist_ASVsLength_16S.pdf")
hist(nchar(getSequences(seqtab)),
     xlab="Size",
     ylab="Frequency",
     main = "ASVs length",
     xlim=c(0,500),
     ylim=c(0,8000)) # visualise the length of ASVs
dev.off()

###Dechimerization###
seqtab.nochim <- removeBimeraDenovo(seqtab, 
                                    method = "pooled", 
                                    multithread = TRUE,
                                    verbose = TRUE) 
save(seqtab.nochim, file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_tab/seqtab.nochim_16S.rdata")
save(seqtab.nochim, file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/biostat_analysis/0_DADA2_output/seqtab.nochim_16S.rdata")

suivi_dechim <- round((sum(seqtab.nochim)/sum(seqtab)*100),2) # % of total number of non-chemical sequences / total number of sequences
write.table(suivi_dechim, file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_tab/suivi_dechim_16S.csv", sep=";", col.names = NA)

pdf(file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_figures/hist_NonChimericASVsLength_16S.pdf")
hist(nchar(getSequences(seqtab.nochim)),
     xlab="Size",
     ylab="Frequency",
     main = "Non-chimeric ASVs length",
     xlim=c(0,500),
     ylim=c(0,8000)) # Length of non-chimeric sequences
dev.off()

seqtab.nochim.bin <- ifelse(seqtab.nochim>0,1,0) # transforms the occurrences of ASVs into presence/absence; this makes it possible to quantify the number of ASVs per sample
save(seqtab.nochim.bin, file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_tab/seqtab.nochim.bin_16S.rdata")

save(list = ls(), file = "/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_env/env_entier_16S_dada2_part4.RData") # save the entire environment in the chosen directory

########## Part 5 ##########

#load("/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_env/env_entier_16S_dada2_part4.RData")

getN <- function(x) sum(getUniques(x))

track <- data.frame(Input=as.numeric(out[,1]), # Raw sequences
                    Filtered=as.numeric(out[,2]), # Filtered sequences
                    "Filt//In"=as.numeric(round(((out[,2]/out[,1])*100),2)),# % filtered/raw sequences
                    Merge = as.numeric(sapply(mergers, getN)), # sequences that have merged
                    "Mer//In"=as.numeric(round(((sapply(mergers, getN)/out[,1])*100),2)),# % merged/raw sequences
                    Nonchim = as.numeric(rowSums(seqtab.nochim)),# Non-chimeric                      
                    "Nonchim//In"=as.numeric(round(((rowSums(seqtab.nochim)/out[,1])*100),2)),# % non-chimeric/raw sequences
                    ASV = as.numeric(rowSums(seqtab.nochim.bin))) # Number of ASVs per sample

rownames(track) <- sample.name # Row names

head(track)
write.table(track, file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_tab/table_suivi_16S.csv", sep=";", col.names = NA)

save(list = ls(), file = "/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_env/env_entier_16S_dada2_part5.RData") # save the entire environment in the chosen directory

########## Part 6 ##########

#load("/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_env/env_entier_16S_dada2_part5.RData")
setwd("/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/DB_taxo")

taxotab <- assignTaxonomy(seqtab.nochim,
                          refFasta = "silva_nr99_v138.1_train_set.fa.gz",
                          minBoot = 50, # default = 50. Minimum Bootsrap (represents the confidence level for assignment to a taxonomic rank)
                          multithread=TRUE)

taxotab <- addSpecies(taxotab, "silva_species_assignment_v138.1.fa.gz")

save(taxotab, file = "/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_tab/taxo_tab_16S.rdata")
save(taxotab, file = "/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/biostat_analysis/0_DADA2_output/taxo_tab_16S.rdata")

save(list = ls(), file = "/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_env/env_entier_16S_dada2_part6.RData") # save the entire environment in the chosen directory