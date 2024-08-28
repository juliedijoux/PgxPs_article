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