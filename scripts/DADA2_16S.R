library(anyLib)
packagesNeeded_dada2 <- c("dada2", "ShortRead", "Biostrings")
anyLib(packagesNeeded_dada2, force = FALSE)

#setwd("/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S/cut_fastq")

########## Partie 1 ##########

path <- "/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S/cut_fastq" # variable qui nous indique où récupérer les fastq
list.files(path) # visualisation des fastq

fnFs <- sort(list.files(path, pattern = "_1.fastq")) # liste des fastq sens
fnRs <- sort(list.files(path, pattern = "_2.fastq")) # liste des fastq antisens

sample.name <- sapply(strsplit(fnFs, "_cut_"), '[', 1) # variable pour extraire le nom des échantillons. On part du principe que les noms des fichiers fastq ont un format SAMPLENAME_XXX.fastq.
sample.name 

fnFs <- file.path(path, fnFs) # précise le chemin d'accés de la variable
fnRs <- file.path(path, fnRs) # précise le chemin d'accés de la variable

#svg(file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_figures/plotquali_sample1_sens.svg") # Permet de sauvegarder le prochain plot sous format svg
#plotQualityProfile(fnFs[1]) # visualisation du Qscore associé à chaque nucléotides du 1er échantillon sens
#dev.off() # pour supprimer le plot
#svg(file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_figures/plotquali_sample1_antisens.svg") # Permet de sauvegarder le prochain plot sous format svg
#plotQualityProfile(fnRs[1]) # visualisation du Qscore associé à chaque nucléotides du 1er échantillon antisens
#dev.off() # pour supprimer le plot

pdf(file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_figures/plotquali_agregate_sens.pdf") # Permet de sauvegarder le prochain plot sous format pdf
plotQualityProfile(fnFs, aggregate = TRUE) # visualisation agrégée du Qscore associé à chaque nucléotides sens
dev.off() # pour supprimer le plot
pdf(file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_figures/plotquali_agregate_antisens.pdf") # Permet de sauvegarder le prochain plot sous format pdf
plotQualityProfile(fnRs, aggregate = TRUE) # visualisation agrégée du Qscore associé à chaque nucléotides antisens
dev.off() # pour supprimer le plot

save(list = ls(), file = "/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_env/env_entier_16S_dada2_part1.RData") # sauvegarde de tout l'environnement dans le répertoire choisi

########## Partie 2 ##########

#load("/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_env/env_entier_16S_dada2_part1.RData") # Chargement de l'environnement R.

filt_path <- file.path(path, "filtered_pairedend") # Création du dossier filtered_pairedend et des objects filtFs et filtRs pour stocker les séquences filtrées

filtFs <- file.path(filt_path, paste0(sample.name, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.name, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncQ = 10, # truncQ = indice Q score minimal ; à la 1ère instance d'un score de qualité inférieur ou égal à truncQ, la séquence est tronquée.
                     truncLen = c(227,205), # truncLen = définit à quelle longueur les séquences vont être tronquées ; les séquences les plus courtes que la longueur sont éliminées.
                     #trimLeft = c(18,20), # trimLeft = définit la longueur que l'on coupe du côté 5' des séquences ; cela permet d'enlever les amorces si ça n'a pas été fait préalablement (ici ITS3KYO2 = 18bp et ITS4 = 20bp).
                     maxEE = 2, # maxEE = le nombre maximum "d'erreur attendues" autorisées dans une lecture. Ce filtre se base sur l'indice Q score. Plus on augment le chiffre, moins on est strict.
                     multithread=TRUE) # pour utiliser plusieurs coeurs.
# D'autres option peuvent être ajouter, voir la rebrique "help".

# Vérification et suppression des fichiers FASTQ vides
emptyFiles <- c()

for (file in c(fnFs, fnRs)) {
  if (file.size(file) == 0) {
    cat("Removing empty file:", file, "\n")
    file.remove(file)
    emptyFiles <- c(emptyFiles, file)
  }
}

cat("Empty files:", emptyFiles, "\n")

filtFs <- file.path(filt_path, paste0(sample.name, "_F_filt.fastq.gz")) # On redéfinit les paths au cas où des fichiers auraient été supprimés
filtRs <- file.path(filt_path, paste0(sample.name, "_R_filt.fastq.gz")) # On redéfinit les paths au cas où des fichiers auraient été supprimés

pourc <- cbind((out[,2]/out[,1])*100) # % de seq filtrées/seq non-filtrées
pourc_disc <- cbind(out, pourc) # combines out and pourcfiltre
pourc_disc
write.table(pourc_disc, file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_tab/pourc_seqTrimmed_sur_seq_noTrimmed_16S.csv", sep=";")

suivi_filterAndTrim <- (mean(out[,2])/mean(out[,1]))*100 # nous indique le % moyen de séquences ayant passées les paramètres de filtrage !
write.table(pourc_disc, file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_tab/suivi_filterAndTrim_16S.csv", sep=";", col.names = NA)

save(list = ls(), file = "/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_env/env_entier_16S_dada2_part2.RData") # sauvegarde de tout l'environnement dans le répertoire choisi

########## Partie 3 ##########

#load("/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_env/env_entier_16S_dada2_part2.RData") # Chargement de l'environnement R.

###Taux d'Erreur###
errF <- learnErrors(filtFs) # Estimation du taux d'erreur de séquençage des sens
errR <- learnErrors(filtRs) # Estimation du taux d'erreur de séquençage des antisens

pdf(file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_figures/plotErrors_sens_16S.pdf")
plotErrors(errF, nominalQ = TRUE) # visualisation du taux d'erreur de séquençage des sens
dev.off()
pdf(file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_figures/plotErrors_antisens_16S.pdf")
plotErrors(errR, nominalQ = TRUE) # visualisation du taux d'erreur de séquençage des antisens
dev.off()

###Déréplication###
derepFs <- derepFastq(filtFs) # les séquences sens identiques vont être regroupées en séquences uniques auxquelles sont attribuées des abondances
names(derepFs) <- sample.name # les séquences sens dérépliquées prennent le nom des échantillons d'où elles proviennent

derepRs <- derepFastq(filtRs) # les séquences antisens identiques vont être regroupées en séquences uniques auxquelles sont attribuées des abondances
names(derepRs) <- sample.name # les séquences antisens dérépliquées prennent le nom des échantillons d'où elles proviennent

###Inférence des Echantillons###
dadaFs <- dada(derepFs,
               err = errF,
               multithread=TRUE,
               pool = TRUE)
dadaRs <- dada(derepRs,
               err = errR,
               multithread=TRUE,
               pool = TRUE)

dadaFs[[1]] # permet de visualiser les données relatives à l'échantillon sens n°1 (dans l'ordre alphabétique, sample.name).
dadaRs[[1]] # permet de visualiser les données relatives à l'échantillon antisens n°1 (dans l'ordre alphabétique, sample.name).

save(list = ls(), file = "/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_env/env_entier_16S_dada2_part3.RData") # sauvegarde de tout l'environnement dans le répertoire choisi

########## Partie 4 ##########

#load("/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_env/env_entier_16S_dada2_part3.RData")

###Fusion des Séquences###
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, 
                      minOverlap = 12, # minOverlap = définit la taille minimale du chevauchement des brins sens et anti-sens pour que leur fusion soit acceptée ; les séquences qui ne fusionnent pas sont éliminées.
                      maxMismatch = 0) # maxMismatch = définit le nombre maximal d'incompatibilité nucléotidique dans le chevauchement.
# D'autres paramètres peuvent également être modifiés, ils sont accessibles à la page d'aide de la fonction: ?mergePairs. Par exemple, si returnRejects = TRUE, les paires qui ont été rejetées en raison de discordances dans la région de chevauchement sont conservés dans la sortie.

head(mergers[[1]]) # visualisation de la fusion
max(mergers[[1]]$nmatch) # Taille du plus grand chevauchement (overlap)
min(mergers[[1]]$nmatch) # Taille du plus petit chevauchement

###Tableau des ASV###
seqtab <- makeSequenceTable(mergers) #pour stocker les ASVs dans l'objet seqtab
dim(seqtab)

seqtab[,1] #renseigne le nombre de fois qu'on retrouve la séquence du premier ASV dans chaque échantillon !
seqtab
save(seqtab, file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_tab/seqtab_16S.rdata")

pdf(file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_figures/hist_ASVsLength_16S.pdf")
hist(nchar(getSequences(seqtab)),
     xlab="Size",
     ylab="Frequency",
     main = "ASVs length",
     xlim=c(200,500),
     ylim=c(0,8000)) #visualiser la longueur des ASVs
dev.off()

###Déchimérisation###
seqtab.nochim <- removeBimeraDenovo(seqtab, 
                                    method = "pooled", 
                                    multithread = TRUE,
                                    verbose = TRUE) 
save(seqtab.nochim, file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_tab/seqtab.nochim_16S.rdata")

suivi_dechim <- round((sum(seqtab.nochim)/sum(seqtab)*100),2) # Pourcentage du nombre total de séquences non-chimériques / nombre total de séquences
write.table(suivi_dechim, file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_tab/suivi_dechim_16S.csv", sep=";", col.names = NA)

pdf(file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_figures/hist_NonChimericASVsLength_16S.pdf")
hist(nchar(getSequences(seqtab.nochim)),
     xlab="Size",
     ylab="Frequency",
     main = "Non-chimeric ASVs length",
     xlim=c(200,500),
     ylim=c(0,8000)) # Longueur de séquences non-chimériques
dev.off()

seqtab.nochim.bin <- ifelse(seqtab.nochim>0,1,0) # transforme les occurences des ASVs en présence/absence. Ceci permet de quantifier le nombre d'ASVs par échantillons.
save(seqtab.nochim.bin, file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_tab/seqtab.nochim.bin_16S.rdata")

save(list = ls(), file = "/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_env/env_entier_16S_dada2_part4.RData") # sauvegarde de tout l'environnement dans le répertoire choisi

########## Partie 5 ##########

#load("/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_env/env_entier_16S_dada2_part4.RData")

getN <- function(x) sum(getUniques(x))

track <- data.frame(Input=as.numeric(out[,1]), # Séquences brutes 
                    Filtered=as.numeric(out[,2]), # Séquences filtres
                    "Filt//In"=as.numeric(round(((out[,2]/out[,1])*100),2)),# % (Filtrées/Brutes)
                    Merge = as.numeric(sapply(mergers, getN)), # Mergés 
                    "Mer//In"=as.numeric(round(((sapply(mergers, getN)/out[,1])*100),2)),# % (Mergés/Brutes)
                    Nonchim = as.numeric(rowSums(seqtab.nochim)),# Non-chimériques                       
                    "Nonchim//In"=as.numeric(round(((rowSums(seqtab.nochim)/out[,1])*100),2)),# % (Non-chimériques/Brutes)
                    ASV = as.numeric(rowSums(seqtab.nochim.bin))) # Nombre d'ASVs par échantillons 

rownames(track) <- sample.name # Noms des lignes

head(track)
write.table(track, file="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_tab/table_suivi_16S.csv", sep=";", col.names = NA)

save(list = ls(), file = "/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_env/env_entier_16S_dada2_part5.RData") # sauvegarde de tout l'environnement dans le répertoire choisi

########## Partie 6 ##########

#load("/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_env/env_entier_16S_dada2_part5.RData")
setwd("/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/DB_taxo")

taxotab <- assignTaxonomy(seqtab.nochim,
                          refFasta = "silva_nr99_v138.1_train_set.fa.gz",
                          minBoot = 50, # Par défaut = 50. Bootsrap minimum (représente le niveau de confiance pour l'assignation à un rang taxonomique).
                          multithread=TRUE)

taxotab <- addSpecies(taxotab, "silva_species_assignment_v138.1.fa.gz")

save(taxotab, file = "/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_tab/taxo_tab_16S.rdata")

save(list = ls(), file = "/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_env/env_entier_16S_dada2_part6.RData") # sauvegarde de tout l'environnement dans le répertoire choisi