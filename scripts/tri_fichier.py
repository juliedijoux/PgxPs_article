import os
import shutil
from glob import glob
from tqdm import tqdm

in16S=glob("/home/unc/jdijoux/NGS_analysis_Pg.Ps/raw_data/*_1_*.fastq.gz")
inITS=glob("/home/unc/jdijoux/NGS_analysis_Pg.Ps/raw_data/*_2_*.fastq.gz")

out16S="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S/"
if not os.path.exists(out16S):
    os.makedirs(out16S)

for f in tqdm(in16S):
    fname=f.split("/")[-1]
    shutil.copyfile(f,out16S+fname)

outITS="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/ITS/"
if not os.path.exists(outITS):
    os.makedirs(outITS)

for f in tqdm(inITS):
    fname=f.split("/")[-1]
    shutil.copyfile(f,outITS+fname)

out16S_env="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_env/"
if not os.path.exists(out16S_env):
    os.makedirs(out16S_env)

outITS_env="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/ITS_env/"
if not os.path.exists(outITS_env):
    os.makedirs(outITS_env)

out16S_figures="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_figures/"
if not os.path.exists(out16S_figures):
    os.makedirs(out16S_figures)

outITS_figures="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/ITS_figures/"
if not os.path.exists(outITS_figures):
    os.makedirs(outITS_figures)

out16S_PICRUSt2="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_PICRUSt2/"
if not os.path.exists(out16S_PICRUSt2):
    os.makedirs(out16S_PICRUSt2)

outITS_PICRUSt2="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/ITS_PICRUSt2/"
if not os.path.exists(outITS_PICRUSt2):
    os.makedirs(outITS_PICRUSt2)

out16S_tab="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/16S_tab/"
if not os.path.exists(out16S_tab):
    os.makedirs(out16S_tab)

outITS_tab="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/ITS_tab/"
if not os.path.exists(outITS_tab):
    os.makedirs(outITS_tab)

outbio_stat="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/biostat_analysis/"
if not os.path.exists(outbio_stat):
    os.makedirs(outbio_stat)

outDADA2_output="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/biostat_analysis/0_DADA2_output/"
if not os.path.exists(outDADA2_output):
    os.makedirs(outDADA2_output)

out_log="/home/unc/jdijoux/NGS_analysis_Pg.Ps/work/log/"
if not os.path.exists(out_log):
    os.makedirs(out_log)