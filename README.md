# Requirements
## Download Raw Data
The raw sequencing reads used in the study are available in the NCBI SRA database under **`BioProject PRJNA1157462`**.

Once downloaded and placed in a folder named **`raw_data`**, the data integrity can be verified with the file **`HN00186490_328samples_md5sum_DownloadLink.txt`** using the following lines of code in Anaconda3:
```bash
awk '{print $3 " " $1}' HN00186490_328samples_md5sum_DownloadLink.txt | grep -v File > md5sum.txt
cat md5sum.txt
md5sum -c md5sum.txt
```

## Workspace
Download all files in this directory and place them as follows:
```plaintext
~/
├── raw_data/
│   ├── sample1_F.fastq.gz
│   ├── sample1_R.fastq.gz
│   └── ...
└── work/
    ├── DB_taxo/
    │   ├── sh_general_release_dynamic_s_all_25.07.2023.fasta
    │   ├── silva_nr99_v138.1_train_set.fa.gz
    │   └── silva_species_assignment_v138.1.fa.gz
    └── scripts/
        ├── checking_primers_after_cut.filt_16S.R
        ├── checking_primers_after_cut.filt_ITS.R
        ├── checking_primers_before_cut.filt_16S.R
        ├── checking_primers_before_cut.filt_ITS.R
        ├── DADA2_16S.R
        ├── DADA2_ITS.R
        ├── pipeline_DADA2_16S.sh
        ├── pipeline_DADA2_ITS.sh
        └── tri_fichier.py 
```
Then, run the script **`tri_fichier.py`** to sort the raw data and complete the organization of the workspace.

# Pipeline
## Requirements
The following files allow you to set up the different Conda environments used to run the pipeline:
- **r420**
- **py27**
- **cutadapt**

## Quality Control
The **`checking_primers_before_cut.filt_16S.R`** and **`checking_primers_before_cut.filt_ITS.R`** scripts allow the tracking of primers in raw data, whereas the **`checking_primers_after_cut.filt_16S.R`** and **`checking_primers_after_cut.filt_ITS.R`** scripts allow for the validation of primer removal step in filtered reads.

The quality controls were performed using the *[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ "A quality control tool for high throughput sequence data")* and the summary quality reports were edited using the *[MultiQC](https://seqera.io/multiqc/ "MultiQC: summarize analysis results for multiple tools and samples in a single report")* packages.

## Removal of 5' and 3' Primers With Cutadapt
To remove primers, the *[cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html "Cutadapt removes adapter sequences from high-throughput sequencing reads")* package were used.

## DADA2 Workflows
The **`DADA2_16S.R`** and **`DADA2_ITS.R`** scripts were adapted from the [DADA2 workflow](https://benjjneb.github.io/dada2/index.html "DADA2: Fast and accurate sample inference from amplicon data with single-nucleotide resolution") of Benjamin Callahan.

## Taxonomic Assignment
The formatted **[SILVA database version 138.1](https://zenodo.org/records/4587955 "Silva 138.1 prokaryotic SSU taxonomic training data formatted for DADA2")** were used for the bacterial library. While the updated **[UNITE all eukaryote version 9.0](https://doi.plutof.ut.ee/doi/10.15156/BIO/2938070 "UNITE general FASTA release for eukaryotes 2")** were used for the fungal library.
- **SILVA database version 138.1**: **`silva_nr99_v138.1_train_set.fa.gz`** and **`silva_species_assignment_v138.1.fa.gz`**
- **UNITE all eukaryote version 9.0**: **`sh_general_release_dynamic_s_all_25.07.2023.fasta`**

## Run Pipeline
The entire pipeline can be run using the **`pipeline_DADA2_16S.sh`** and **`pipeline_DADA2_ITS.sh`** scripts:
```bash
sbatch ~/work/sripts/pipeline_DADA2_16S.sh
sbatch ~/work/sripts/pipeline_DADA2_ITS.sh
```
# References
Andrews, S., Krueger, F., Segonds-Pichon, A., Biggins, L., Krueger, C., & Wingett, S. (2012). FastQC. Babraham, UK.

Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047–3048.

Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10–12. Available from: https://journal.embnet.org/index.php/embnetjournal/article/view/200

Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. Nature Methods, 13(7), 581–583.

McLaren, M. R., & Callahan, B. J. (2021). Silva 138.1 prokaryotic SSU taxonomic training data formatted for DADA2 [Data set]. Zenodo. Available from: https://doi.org/10.5281/zenodo.4587955

Abarenkov, K., Zirk, A., Piirmann, T., Pöhönen, R., Ivanov, F., Nilsson, R. H., et al. (2023). UNITE general FASTA release for eukaryotes 2 [Data set]. UNITE Community. Available from: https://dx.doi.org/10.15156/BIO/2938070
