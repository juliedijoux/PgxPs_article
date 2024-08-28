# How to check raw data integrity

cd /home/unc/jdijoux/NGS_analysis_Pg.Ps/raw_data

awk '{print $3 " " $1}' HN00186490_328samples_md5sum_DownloadLink.txt | grep -v File > md5sum.txt
cat md5sum.txt
md5sum -c md5sum.txt