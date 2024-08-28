#!/bin/bash
#SBATCH -c 12
#SBATCH --mem 64G
#SBATCH --job-name=PICRUSt2
#SBATCH --output=/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/log/PICRUSt2_16S.log

hostname
date

module load anaconda3
conda activate PLS_2

cd /home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_PICRUSt2

# All ASVs of all species
python /home/unc/jdijoux/.conda/envs/PLS_2/pkgs/picrust2-2.5.2-pyhdfd78af_0/python-scripts/picrust2_pipeline.py -s all.sp_all.asv.fasta -i all.sp_all.asv.txt -o all.sp_all.asv --verbose
python predict_metagenomes.py -i output/metagenome_contributions.tsv -o fichier_de_sortie_metagenomes.tsv

# Pg u
python /home/unc/jdijoux/.conda/envs/PLS_2/pkgs/picrust2-2.5.2-pyhdfd78af_0/python-scripts/picrust2_pipeline.py -s Pg_u.fasta -i Pg_u.txt -o Pg_u --verbose
python predict_metagenomes.py -i output/metagenome_contributions.tsv -o fichier_de_sortie_metagenomes.tsv

# Pg r
python /home/unc/jdijoux/.conda/envs/PLS_2/pkgs/picrust2-2.5.2-pyhdfd78af_0/python-scripts/picrust2_pipeline.py -s Pg_r.fasta -i Pg_r.txt -o Pg_r --verbose
python predict_metagenomes.py -i output/metagenome_contributions.tsv -o fichier_de_sortie_metagenomes.tsv

# Ps u
python /home/unc/jdijoux/.conda/envs/PLS_2/pkgs/picrust2-2.5.2-pyhdfd78af_0/python-scripts/picrust2_pipeline.py -s Ps_u.fasta -i Ps_u.txt -o Ps_u --verbose
python predict_metagenomes.py -i output/metagenome_contributions.tsv -o fichier_de_sortie_metagenomes.tsv

# Ps r
python /home/unc/jdijoux/.conda/envs/PLS_2/pkgs/picrust2-2.5.2-pyhdfd78af_0/python-scripts/picrust2_pipeline.py -s Ps_r.fasta -i Ps_r.txt -o Ps_r --verbose
python predict_metagenomes.py -i output/metagenome_contributions.tsv -o fichier_de_sortie_metagenomes.tsv

date