# Requirements
The raw sequencing reads used in the study are available in the NCBI SRA database under **`BioProject PRJNA1157462`**.

Once downloaded and placed in a folder named **`raw_data`**, the data integrity can be verified with the file **`HN00186490_328samples_md5sum_DownloadLink.txt`** using the following lines of code in Anaconda3:
```bash
awk '{print $3 " " $1}' HN00186490_328samples_md5sum_DownloadLink.txt | grep -v File > md5sum.txt
cat md5sum.txt
md5sum -c md5sum.txt
```
Download and place all the scripts in the folder **`~/work/scripts`**. Run the script **`tri_fichier.py`** to sort the raw data and organize the workspace.

# Pipeline
## Requirements
The following files allow you to set up the different Conda environments used to run the pipeline:
- **r420**
- **py27**
- **cutadapt**

## Quality control
The **`checking_primers_before_cut.filt_16S.R`** and **`checking_primers_before_cut.filt_ITS.R`** scripts allow the tracking of primers in raw data, whereas the **`checking_primers_after_cut.filt_16S.R`** and **`checking_primers_after_cut.filt_ITS.R`** scripts allow for the validation of primer removal step in filtered reads.

The quality controls were performed using the *[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ "A quality control tool for high throughput sequence data")* and the summary quality reports were edited using the *[MultiQC](https://seqera.io/multiqc/ "MultiQC: summarize analysis results for multiple tools and samples in a single report")* packages.

## Removal of 5' and 3' primers with cutadapt
To remove primers, the *[cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html "Cutadapt removes adapter sequences from high-throughput sequencing reads")* package were used.

## DADA2 workflow
The **`DADA2_16S.R`** and **`DADA2_ITS.R`** scripts were adapted from the [DADA2 workflow](https://benjjneb.github.io/dada2/index.html "DADA2: Fast and accurate sample inference from amplicon data with single-nucleotide resolution") of Benjamin Callahan.

## Run pipeline
```bash
sbatch ~/work/sripts/pipeline_DADA2_16S.sh
sbatch ~/work/sripts/pipeline_DADA2_ITS.sh
```
# References

