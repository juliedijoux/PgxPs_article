#!/bin/bash
#SBATCH -c 12
#SBATCH --mem 64G
#SBATCH --job-name=DADA2_ITS
#SBATCH --output=/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/log/pipeline_DADA2_ITS.log

hostname
date

module load anaconda3

# Checking the presence of primers

conda activate r420 
cd /home/unc/jdijoux/NGS_analysis_Pg.Ps/work/ITS

Rscript /home/unc/jdijoux/NGS_analysis_Pg.Ps/work/scripts/checking_primers_before_cut.filt_ITS.R

conda deactivate

# Quality control of raw data

conda activate py27

cd /home/unc/jdijoux/NGS_analysis_Pg.Ps/work/ITS
mkdir fastqc_outfiles
fastqc -o fastqc_outfiles *.fastq.gz

cd fastqc_outfiles
multiqc . -o MULTIQC

conda deactivate

# Removal of 5' and 3' primers with cutadapt

conda activate cutadapt

cd
cd /home/unc/jdijoux/NGS_analysis_Pg.Ps/work/ITS

# CTTGGTCATTTAGAGGAAGTAA FWD (5'-3') ITS-1F revcomp: TTACTTCCTCTAAATGACCAAG
# TGTGTTCTTCATCGATG REV (5'-3') ITS-2R revcomp: CATCGATGAAGAACACA

for f in *_2_[1,2].fastq.gz

do
    n=${f%%_2_[1,2].fastq.gz}
    cutadapt -g CTTGGTCATTTAGAGGAAGTAA -a CATCGATGAAGAACACA \
    -G TGTGTTCTTCATCGATG -A TTACTTCCTCTAAATGACCAAG -n 2 \
    -poly-a -m 20 --discard-untrimmed \
    --output ${n}_cut_1.fastq.gz -p ${n}_cut_2.fastq.gz \
    ${n}_2_1.fastq.gz ${n}_2_2.fastq.gz
done

mkdir cut_fastq
find . -type f -name "*_cut_*" -exec mv '{}' cut_fastq/ \;

conda deactivate

# DADA2 

conda activate r420

cd
cd /home/unc/jdijoux/NGS_analysis_Pg.Ps/work/ITS/cut_fastq
gunzip -k *.fastq.gz

mkdir cut_compressed
find . -type f -name "*.fastq.gz" -exec mv '{}' cut_compressed/ \;

Rscript /home/unc/jdijoux/NGS_analysis_Pg.Ps/work/scripts/DADA2_ITS.R

# Checking primer cut

cd
cd /home/unc/jdijoux/NGS_analysis_Pg.Ps/work/ITS/cut_fastq/filtered_pairedend

Rscript /home/unc/jdijoux/NGS_analysis_Pg.Ps/work/scripts/checking_primers_after_cut.filt_ITS.R

conda deactivate

# Checking the quality of filtered sequences with cutadapt and DADA2

conda activate py27

cd
cd /home/unc/jdijoux/NGS_analysis_Pg.Ps/work/ITS/cut_fastq/filtered_pairedend
mkdir fastqc_outfiles

fastqc -o fastqc_outfiles *.fastq.gz
cd fastqc_outfiles
multiqc . -o MULTIQC

conda deactivate

date