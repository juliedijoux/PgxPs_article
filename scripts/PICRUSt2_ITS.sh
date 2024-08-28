#!/bin/bash
#SBATCH -c 3
#SBATCH --mem 64G
#SBATCH --job-name=PICRUSt2
#SBATCH --output=/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/log/ITS_PICRUSt2.log

hostname
date

module load anaconda3
conda activate PLS_2

cd /home/unc/jdijoux/NGS_analysis_Pg.Ps/work/ITS_PICRUSt2

# Place amplicon sequence variants into reference phylogeny
place_seqs.py -s all.sp_all.asv.fasta \
              -o placed_seqs.tre \
              --ref_dir /home/unc/jdijoux/.conda/envs/PLS_2/lib/python3.8/site-packages/picrust2/default_files/fungi/fungi_ITS \
              -p 3 \
              --intermediate placement_working

# Run hidden-state prediction to get ITS copy numbers and E.C. number per predicted genome
hsp.py --observed_trait_table /home/unc/jdijoux/.conda/envs/PLS_2/lib/python3.8/site-packages/picrust2/default_files/fungi/ec_ITS_counts.txt.gz \
       -t placed_seqs.tre \
       -o EC_predicted.tsv.gz \
       -p 3

hsp.py --observed_trait_table /home/unc/jdijoux/.conda/envs/PLS_2/lib/python3.8/site-packages/picrust2/default_files/fungi/ITS_counts.txt.gz \
       -t placed_seqs.tre \
       -o ITS_counts_predicted.tsv.gz \
       -p 3

# Predict E.C. abundances in sequencing samples (adjusts gene family abundances by ITS sequence abundance)
metagenome_pipeline.py -i all.sp_all.asv.txt \
                       -m ITS_counts_predicted.tsv.gz \
                       -f EC_predicted.tsv.gz \
                       -o EC_metagenome_out

# Infer MetaCyc pathway abundances and coverages based on predicted E.C. number abundances
pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz \
                    -o pathways_out \
                    --map /home/unc/jdijoux/.conda/envs/PLS_2/lib/python3.8/site-packages/picrust2/default_files/pathway_mapfiles/metacyc_path2rxn_struc_filt_fungi.txt \
                    --intermediate pathways_working \
                    -p 3 \
                    --no_regroup

date