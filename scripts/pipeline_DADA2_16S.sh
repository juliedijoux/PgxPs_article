#!/bin/bash
#SBATCH -c 12
#SBATCH --mem 64G
#SBATCH --job-name=DADA2_16S
#SBATCH --output=/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/log/pipeline_DADA2_16S.log

hostname
date

module load anaconda3

# Checking the presence of primers

conda activate r420 
cd /home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S

Rscript /home/unc/jdijoux/NGS_analysis_Pg.Ps/work/scripts/checking_primers_before_cut.filt_16S.R

conda deactivate

# Quality control of raw data

conda activate py27

cd /home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S
mkdir fastqc_outfiles
fastqc -o fastqc_outfiles *.fastq.gz

cd fastqc_outfiles
multiqc . -o MULTIQC

conda deactivate

# Removal of 5' and 3' primers with cutadapt

conda activate cutadapt

cd
cd /home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S

# CCTACGGGNGGCWGCAG FWD (5'-3') 341F revcomp: CTGCWGCCNCCCGTAGG
# GACTACHVGGGTATCTAATCC REV (5'-3') 805R revcomp:GGATTAGATACCCBDGTAGTC

for f in *_1_[1,2].fastq.gz

do
    n=${f%%_1_[1,2].fastq.gz}
    cutadapt -g CCTACGGGNGGCWGCAG -a GGATTAGATACCCBDGTAGTC \
    -G GACTACHVGGGTATCTAATCC -A CTGCWGCCNCCCGTAGG -n 2 \
    -m 20 --discard-untrimmed \
    -j 0 \
    --output ${n}_cut_1.fastq.gz -p ${n}_cut_2.fastq.gz \
    ${n}_1_1.fastq.gz ${n}_1_2.fastq.gz
done

mkdir cut_fastq
find . -type f -name "*_cut_*" -exec mv '{}' cut_fastq/ \;

conda deactivate

# DADA2

conda activate r420

cd
cd /home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S/cut_fastq
gunzip -k *.fastq.gz

mkdir cut_compressed
find . -type f -name "*.fastq.gz" -exec mv '{}' cut_compressed/ \;

Rscript /home/unc/jdijoux/NGS_analysis_Pg.Ps/work/scripts/DADA2_16S.R

# Checking primer cut

cd
cd /home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S/cut_fastq/filtered_pairedend

Rscript /home/unc/jdijoux/NGS_analysis_Pg.Ps/work/scripts/checking_primers_after_cut.filt_16S.R

conda deactivate

# Checking the quality of filtered sequences with cutadapt and DADA2

conda activate py27

cd
cd /home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S/cut_fastq/filtered_pairedend
mkdir fastqc_outfiles

fastqc -o fastqc_outfiles *.fastq.gz

cd fastqc_outfiles
multiqc . -o MULTIQC

conda deactivate
date