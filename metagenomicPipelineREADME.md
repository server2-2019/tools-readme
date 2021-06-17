# tools-readme
Record readme of often used software

2020.08.19
bioBakery tools for meta'omic profiling
Huttenhower lab
ref: https://github.com/biobakery/biobakery/wiki
CD-HIT version 4.8.1

gene relative abundance
# /work/workspace/zhurj/script/5_wgspro/16_sumtaxo

参考网站
# 宏基因组教程Metagenomics Tutorial (HUMAnN2): https://blog.csdn.net/woodcorpse/article/details/80413379
# 代码放送 | 宏基因组实战： https://www.dxy.cn/bbs/newweb/pc/post/39249887

source activate humann
# https://anaconda.org/bioconda/bbmap

conda install -c bioconda bbmap

source activate humann
# https://anaconda.org/bioconda/kneaddata
# http://huttenhower.sph.harvard.edu/kneaddata
# https://bitbucket.org/biobakery/physalia_2018/wiki/Lab%20%235:%20KneadData#markdown-header-kneaddata-outputs
conda install -c bioconda kneaddata

metagenomics analysis ref: https://en.wikipedia.org/wiki/Metagenomics#Sequence_pre-filtering
1. /work/workspace/zhurj/script/5_wgspro/1_filter
assembly free:
2. metaphlan3 ( MetaPhlAn version 3.0.1 (25 Jun 2020) )
# env from zhengjiao： source activate /work/workspace/liangzj/program/Miniconda3/envs/metaphlan3
source activate humann
# humann humann v3.0.0.alpha.4， installed metaphlan3
ref: https://huttenhower.sph.harvard.edu/metaphlan3/
MetaPhlAn relies on unique clade-specific marker genes identified from ~17,000 reference genomes (~13,500 bacterial and archaeal, ~3,500 viral, and ~110 eukaryotic), allowing: 
```
cd /work/workspace/zhurj/project/1_metadata/metapipe/clean
metaphlan MNC00230.clean.1.fq.gz,MNC00230.clean.2.fq.gz --bowtie2out /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00230/metagenome.bowtie2.bz2 --input_type fastq -t rel_ab_w_read_stats --nproc 36 --sample_id_key MNC00230 -o /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00230/profile_1.txt
rm /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00230/metagenome.bowtie2.bz2 -f
```
cd /work/workspace/zhurj/project/1_metadata/metapipe/clean
metaphlan MNC00230.clean.1.fq.gz,MNC00230.clean.2.fq.gz --bowtie2out /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00230/metagenome.bowtie2.bz2 --input_type fastq --nproc 36 --sample_id_key MNC00230 -o /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00230/profile_1.txt
rm /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00230/metagenome.bowtie2.bz2 -f

metaphlan /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00230/metagenome.bowtie2.bz2 --input_type bowtie2out --tax_lev 'p' --nproc 36 --sample_id_key MNC00230 -o /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00230/phylum.txt

metaphlan /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00230/metagenome.bowtie2.bz2 --input_type bowtie2out --tax_lev 'g' --nproc 36 --sample_id_key MNC00230 -o /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00230/genus.txt

metaphlan /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00230/metagenome.bowtie2.bz2 --input_type bowtie2out --tax_lev 's' --nproc 36 --sample_id_key MNC00230 -o /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00230/species.txt

metaphlan /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00231/metagenome.bowtie2.bz2 --input_type bowtie2out --tax_lev 'p' --nproc 36 --sample_id_key MNC00231 -o /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00231/phylum.txt

metaphlan /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00231/metagenome.bowtie2.bz2 --input_type bowtie2out --tax_lev 'g' --nproc 36 --sample_id_key MNC00231 -o /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00231/genus.txt

metaphlan /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00231/metagenome.bowtie2.bz2 --input_type bowtie2out --tax_lev 's' --nproc 36 --sample_id_key MNC00231 -o /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00231/species.txt

metaphlan /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00230/metagenome.bowtie2.bz2 --input_type bowtie2out --nproc 36 --sample_id_key MNC00230 -o /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00230/profile_2.txt

metaphlan /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00231/metagenome.bowtie2.bz2 --input_type bowtie2out --tax_lev 'a' --nproc 36 -t rel_ab_w_read_stats --sample_id_key MNC00231 -o /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00231/profile.txt

MNC00231metaphlan.sh
```
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
metaphlan /work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00231.clean.1.fq.gz,/work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00231.clean.2.fq.gz --bowtie2out /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00231/metagenome.bowtie2.bz2 --input_type fastq --nproc 36 --sample_id_key MNC00231 -o /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00231/profile_1.txt
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

```
time srun -o MNC00231metaphlan.out -e MNC00231metaphlan.err -N 1 -c 20 -p slurm128 bash MNC00231metaphlan.sh &

source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
# 所有的profile 文件，建立一个软连接放到一个文件夹下， 再进行合并
merge_metaphlan_tables.py /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00231/profile_1.txt /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00230/profile_2.txt > metaphlan2_merged_sam2.tsv
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann
ln -s /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00231/profile_1.txt /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/profile/MNC00231.txt
ln -s /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00230/profile_2.txt /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/profile/MNC00230.txt
merge_metaphlan_tables.py  /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/profile/*.txt > /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/merge/metaphlan2_merged_softlink.tsv

ln -s /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00230/phylum.txt /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/profile/phylum/MNC00230.txt
ln -s /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00231/phylum.txt /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/profile/phylum/MNC00231.txt
ln -s /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00230/genus.txt /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/profile/genus/MNC00230.txt
ln -s /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00231/genus.txt /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/profile/genus/MNC00231.txt
ln -s /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00230/species.txt /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/profile/species/MNC00230.txt
ln -s /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00231/species.txt /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/profile/species/MNC00231.txt

merge_metaphlan_tables.py  /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/profile/phylum/*.txt > /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/merge/metaphlan3_merged_phylum.tsv
merge_metaphlan_tables.py  /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/profile/genus/*.txt > /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/merge/metaphlan3_merged_genus.tsv
merge_metaphlan_tables.py  /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/profile/species/*.txt > /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/merge/metaphlan3_merged_species.tsv

# profile with extimated read counts
ln -s  /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00230/profile.txt /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/abunread/MNC00230.txt
ln -s  /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00231/profile.txt /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/abunread/MNC00231.txt

merge_metaphlan_tables.py  /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/abunread/*.txt > /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/merge/metaphlan3_merged_read.tsv

alpha diversity analysis
https://www.drive5.com/usearch/manual/cmd_alpha_div.html
# Calculate alpha diversity using metaphlan2 results: https://www.programmersought.com/article/76941074682/
download microbiome_helper-master.zip to /work/workspace/zhurj/software from  https://github.com/LangilleLab/microbiome_helper
cd /work/workspace/zhurj/software
unzip microbiome_helper-master.zip
perl /work/workspace/zhurj/software/microbiome_helper-master/metaphlan_to_stamp.pl /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/merge/taxo_all.txt > /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/merge/taxo_all.spf
去除 /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/merge/taxo_all.spf 中的重复行
# 微生物多样性组间差异分析神器-STAMP： https://mp.weixin.qq.com/s/f02mPAXCobPQijkIis-c2w



usearch -alpha_div /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/merge/taxo_all_read.txt -output /work/workspace/zhurj/project/1_metadata/metapipe/alphaDiv/alpha.txt
usearch -alpha_div otutable.txt -output gini.txt -metrics gini_simpson
usearch -alpha_div otutable.txt -output alpha.txt -metrics chao1,berger_parker

beta diversity analysis
https://www.drive5.com/usearch/manual/cmd_beta_div.html

************
linux VERSION, CAN NOT WORK PROPERLY SINCE 
************
STAMP is a software package for analyzing taxonomic or metabolic profiles that promotes: https://beikolab.cs.dal.ca/software/STAMP
Quick installation instructions for STAMP: https://beikolab.cs.dal.ca/software/Quick_installation_instructions_for_STAMP
# https://anaconda.org/bioconda/stamp
source activate /work/workspace/zhurj/bin/miniconda3/envs/py2
conda install -c bioconda stamp
v2.1.3
also download version for windows

STAMP use windows version do the community differential analysis
lefse
lefse input file

```
Condition   NC_0    NC_0    NC_0    NC_0    NC_0    NC_0    NC_0    HFC_0   HFC_0   HFC_0   HFC_0   HFC_0   HFC_0   HFC_0   HFD_0   HFD_0   HFD_    0   HFD_0   HFD_0   HFD_0   HFD_0   HFP_0   HFP_0   HFP_0   HFP_0   HFP_0   HFP_0   HFP_0   HFP_0   HFDP_0  HFDP_0  HFDP_0  HFDP_0  HFDP_0  HFDP    _0  HFDP_0  HFDP_0  HFDP_0  NC_1    NC_1    NC_1    NC_1    NC_1    NC_1    NC_1    NC_1    NC_1    HFC_1   HFC_1   HFC_1   HFC_1   HFC_1   HFC_    1   HFC_1   HFC_1   HFC_1   HFD_1   HFD_1   HFD_1   HFD_1   HFD_1   HFD_1   HFD_1   HFD_1   HFD_1   HFP_1   HFP_1   HFP_1   HFP_1   HFP_1   HFP_    1   HFP_1   HFP_1   HFP_1   HFPD_1  HFPD_1  HFPD_1  HFPD_1  HFPD_1  HFPD_1  HFPD_1  HFPD_1  HFPD_1  NC_6    NC_6    NC_6    NC_6    NC_6    NC_6        NC_6    NC_6    NC_6    HFC_6   HFC_6   HFC_6   HFC_6   HFC_6   HFC_6   HFC_6   HFC_6   HFC_6   HFD_6   HFD_6   HFD_6   HFD_6   HFD_6   HFD_    6   HFD_6   HFD_6   HFD_6   HFP_6   HFP_6   HFP_6   HFP_6   HFP_6   HFP_6   HFP_6   HFP_6   HFPD_6  HFPD_6  HFPD_6  HFPD_6  HFPD_6  HFPD_6  HFPD    _6  HFPD_6  HFPD_6
Name    0-1 0-2 0-3 0-4 0-5 0-7 0-9 0-10    0-11    0-12    0-14    0-15    0-16    0-17    0-20    0-21    0-23    0-24    0-25    0-26    0-27        0-28    0-29    0-30    0-32    0-33    0-34    0-35    0-36    0-37    0-38    0-39    0-40    0-41    0-42    0-43    0-44    0-45    1-1     1-2 1-3 1-4 1-5 1-6 1-7 1-8 1-9 1-10    1-11    1-12    1-13    1-14    1-15    1-16    1-17    1-18    1-19    1-20    1-21    1-22    1-23        1-24    1-25    1-26    1-27    1-28    1-29    1-30    1-31    1-32    1-33    1-34    1-35    1-36    1-37    1-38    1-39    1-40    1-41        1-42    1-43    1-44    1-45    6-1 6-2 6-3 6-4 6-5 6-6 6-7 6-8 6-9 6-10    6-11    6-12    6-13    6-14    6-15    6-16    6-17    6-18    6-19        6-20    6-21    6-22    6-23    6-24    6-25    6-26    6-27    6-28    6-29    6-30    6-31    6-32    6-33    6-35    6-36    6-37    6-38        6-39    6-40    6-41    6-42    6-43    6-44    6-45
index   MNA00001    MNA00002    MNA00003    MNA00004    MNA00005    MNA00006    MNA00007    MNA00008    MNA00009    MNA00010    MNA00011    MNA0    0012    MNA00013    MNA00014    MNA00015    MNA00016    MNA00017    MNA00018    MNA00019    MNA00020    MNA00021    MNA00022    MNA00023    MNA0    0024    MNA00025    MNA00026    MNA00027    MNA00028    MNA00029    MNA00030    MNA00031    MNA00032    MNA00033    MNA00034    MNA00035    MNA0    0036    MNA00037    MNA00038    MNA00039    MNA00040    MNA00041    MNA00042    MNA00043    MNA00044    MNA00045    MNA00046    MNA00047    MNA0    0048    MNA00049    MNA00050    MNA00051    MNA00052    MNA00053    MNA00054    MNA00055    MNA00056    MNA00057    MNA00058    MNA00059    MNA0    0060    MNA00061    MNA00062    MNA00063    MNA00064    MNA00065    MNA00066    MNA00067    MNA00068    MNA00069    MNA00070    MNA00071    MNA0    0072    MNA00073    MNA00074    MNA00075    MNA00076    MNA00077    MNA00078    MNA00079    MNA00080    MNA00081    MNA00082    MNA00083    MNA0    0084    MNA00085    MNA00086    MNA00087    MNA00088    MNA00089    MNA00090    MNA00091    MNA00092    MNA00093    MNA00094    MNA00095    MNA0    0096    MNA00097    MNA00098    MNA00099    MNA00100    MNA00101    MNA00102    MNA00103    MNA00104    MNA00105    MNA00106    MNA00107    MNA0    0108    MNA00109    MNA00110    MNA00111    MNA00112    MNA00113    MNA00114    MNA00115    MNA00116    MNA00117    MNA00118    MNA00119    MNA0    0120    MNA00121    MNA00122    MNA00123    MNA00124    MNA00125    MNA00126    MNA00127
k__Bacteria_sUN  0   0   0   0   0   0   0   0   8   0   0   3   0   0   0   0   0   0   0   0   7   0   0   0   0   0   0   0   0   0   0   0       0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  13   0   0   0   0   0   0   0   0   0   0   0   0  41   0   0   0       0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   2   0   0  10  15  12   0   4  10  20   0   4   5   4   0       0   0   3   0  12   0   0   0   0   0   4   0   0   0   0  14   0   0   5   0  14   3   5
```




# https://beikolab.cs.dal.ca/software/STAMP
# read1, read2 as input for humann
# merge read1 + read2 into one file
# ref: https://github.com/biobakery/humann

```
# test
# 分析参考 https://github.com/biobakery/humann#guides-to-humann2-utility-scripts
humann -i /work/workspace/zhurj/bin/miniconda3/envs/humann/lib/python3.7/site-packages/humann/tests/data/demo.fastq -o /work/workspace/zhurj/project/1_metadata/metapipe/humann3/test1 
humann -i /work/workspace/zhurj/bin/miniconda3/envs/humann/lib/python3.7/site-packages/humann/tests/data/demo.fastq -o /work/workspace/zhurj/project/1_metadata/metapipe/humann3/test1
```
only run read1 on humann3
humann -i /work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00230.clean.1.fq.gz -o /work/workspace/zhurj/project/1_metadata/metapipe/humann3/MNC00230 --threads 36 


time srun -o phumann3.out -e phumann3.err -N 2 -c 20 -p slurm256 bash phumann3.sh &
```
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel --xapply -j 2 'humann --threads 20 --input {1} --output /work/workspace/zhurj/project/1_metadata/metapipe/humann3/humann3_out/{2}' :::: /work/workspace/zhurj/project/1_metadata/metapipe/clean/humann_in_file :::: /work/workspace/zhurj/project/1_metadata/metapipe/clean/humann_in_id
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

```

clumpify.sh -Xmx20g dedupe dupedist=12000 in1=/work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00230.clean.1.fq.gz in2=/work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00230.clean.2.fq.gz out1=/work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00230.dedup.1.fq.gz out2=/work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00230.dedup.2.fq.gz 

clumpify.sh -Xmx20g dedupe dupedist=12000 in1=/work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00231.clean.1.fq.gz in2=/work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00231.clean.2.fq.gz out1=/work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00231.dedup.1.fq.gz out2=/work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00231.dedup.2.fq.gz 


source activate python3.6
1. input create
perl /work/workspace/zhurj/script/5_wgspro/1_filter/mnid2rawdata.pl -i /work/workspace/zhurj/project/1_metadata/metapipe/input/ifile -t 0 -r metagenome -o /work/workspace/zhurj/project/1_metadata/metapipe/input -n test

#/work/rawdata/run/guangzhou/novogene/2020/07/20200714/run00057/rawdata
raw data QC
perl /work/workspace/zhurj/script/3_swgspro/strainscreen/dataqc.pl -i /work/workspace/zhurj/project/1_metadata/metapipe/input/test_r1r2_2c.txt -o /work/workspace/zhurj/project/1_metadata/metapipe/qc --thread 36 --basecutoff 2.7  --cutoffq20 85 --cutoffq30 80
srun -o qc_test.out -e qc_test.err -N 1 -c 20 -p slurm256 bash qc_test.sh &

2. data filter
perl /work/workspace/zhurj/script/5_wgspro/1_filter/datafilter.pl -i /work/workspace/zhurj/project/1_metadata/metapipe/input/test_r1r2_2c.txt -o /work/workspace/zhurj/project/1_metadata/metapipe
srun -o mnc00230.out -e mnc00230.err -N 1 -c 20 -p slurm256 bash mnc00230.sh &
srun -o x02.out -e x02.err -N 1 -c 20 -p slurm128 bash x02_20200110 &

kneaddata
download referen:
human_genome : bowtie2 = http://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_hg37_and_human_contamination_Bowtie2_v0.1.tar.gz
human_genome : bmtagger = http://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_BMTagger_v0.1.tar.gz
human_transcriptome : bowtie2 = http://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_hg38_transcriptome_Bowtie2_v0.1.tar.gz
ribosomal_RNA : bowtie2 = http://hut-i tenhower.sph.harvard.edu/kneadData_databases/SILVA_128_LSUParc_SSUParc_ribosomal_RNA_v0.1.tar.gz
mouse_C57BL : bowtie2 = http://huttenhower.sph.harvard.edu/kneadData_databases/mouse_C57BL_6NJ_Bowtie2_v0.1.tar.gz

time kneaddata --input /work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00231.clean.1.fq.gz --input /work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00231.clean.2.fq.gz --output /work/workspace/zhurj/project/1_metadata/metapipe/kneaddata/MNC00230 -v \
-db /work/workspace/zhurj/reference/human/kneaddata/20200921  \
--bypass-trim --run-trf \
-t 36 --bowtie2-options "--very-sensitive --dovetail" --remove-intermediate-output  --output-prefix MNC00231
103mins

kneaddata_read_count_table --input /work/workspace/zhurj/project/1_metadata/metapipe/kneaddata/MNC00230 --output /work/workspace/zhurj/project/1_metadata/metapipe/kneaddata/kneaddata_read_counts.txt
kneaddata_read_count_table --input /work/workspace/zhurj/project/1_metadata/metapipe/kneaddata/MNC00231 --output /work/workspace/zhurj/project/1_metadata/metapipe/kneaddata/MNC00231/kneaddata_read_counts.txt



#parallel run kneaddata
cd /work/workspace/zhurj/project/1_metadata/metapipe/p
content in pkneaddata.sh
```
parallel -j 2 --link 'kneaddata -i {1} -i {2} --output /work/workspace/zhurj/project/1_metadata/metapipe/kneaddata/MNC00230/parallel -v \
-db /work/workspace/zhurj/reference/human/kneaddata/20200921  \
--bypass-trim --run-trf \
-t 20 -p 20 --bowtie2-options "--very-sensitive --dovetail" --remove-intermediate-output' \
::: /work/workspace/zhurj/project/1_metadata/metapipe/clean/*.clean.1.fq.gz ::: /work/workspace/zhurj/project/1_metadata/metapipe/clean/*.clean.2.fq.gz
```
time srun -o pkeaddata.out -e pkeaddata.err -N 2 -c 20 -p slurm256 bash pkneaddata.sh &

cd /work/workspace/zhurj/project/1_metadata/metapipe/kneaddata/parallel

kneaddata_read_count_table --input /work/workspace/zhurj/project/1_metadata/metapipe/kneaddata/parallel --output /work/workspace/zhurj/project/1_metadata/metapipe/kneaddata/parallel/kneaddata_read_counts.txt

合并双端
cd 
zcat MNC00230.clean.1.fq.gz MNC00230.clean.2.fq.gz | pigz - > MNC00230.merge.fq.gz
zcat MNC00231.clean.1.fq.gz MNC00231.clean.2.fq.gz | pigz - > MNC00231.merge.fq.gz

parallel -j 9 'humann2 --threads 9 --input {} --output humann2_out/{/.}' ::: cat_reads/*fastq


3. contig assemble
cd /work/workspace/zhurj/project/1_metadata/metapipe/p
srun -o mnc00230.out -e mnc00230.err -N 1 -c 20 -p slurm256 bash mnc00230.sh
content in mnc00230.sh
```
python /work/workspace/moonis/bin/genoassemble/soft/SPAdes-3.13.1-Linux/bin/spades.py --meta -1 /work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00230.clean.1.fq.gz -2 /work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00230.clean.2.fq.gz -t 36 -o /work/workspace/zhurj/project/1_metadata/metapipe/spades/MNC00230
```
```

4. bwa map
ln -s /work/workspace/zhurj/project/1_metadata/metapipe/spades/MNC00230/contigs.fasta /work/workspace/zhurj/project/1_metadata/metapipe/bwa/db/ref.fa
bwa index /work/workspace/zhurj/project/1_metadata/metapipe/bwa/db/ref.fa
bwa mem -t 40 /work/workspace/zhurj/project/1_metadata/metapipe/bwa/db/ref.fa /work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00230.clean.1.fq.gz /work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00230.clean.2.fq.gz > /work/workspace/zhurj/project/1_metadata/metapipe/bwa/MNC00230.sam
samtools view -@ 40 -S -b /work/workspace/zhurj/project/1_metadata/metapipe/bwa/MNC00230.sam -o /work/workspace/zhurj/project/1_metadata/metapipe/bwa/MNC00230.bam
samtools sort -@ 40 /work/workspace/zhurj/project/1_metadata/metapipe/bwa/MNC00230.bam -o /work/workspace/zhurj/project/1_metadata/metapipe/bwa/MNC00230.sorted.bam

5. MetaBAT2 binning
ref: https://bitbucket.org/berkeleylab/metabat/src/master/
```
runMetaBat.sh -o /work/workspace/zhurj/project/1_metadata/metapipe/metabat /work/workspace/zhurj/project/1_metadata/metapipe/bwa/db/ref.fa /work/workspace/zhurj/project/1_metadata/metapipe/bwa/MNC00230.sorted.bam # does not work
```

jgi_summarize_bam_contig_depths --outputDepth /work/workspace/zhurj/project/1_metadata/metapipe/metabat/depth.txt  /work/workspace/zhurj/project/1_metadata/metapipe/bwa/MNC00230.sorted.bam
--minContigLength   arg  The mimimum length of contig to include for mapping and shredding
--minContigDepth    arg  The minimum depth along contig at which to break the contig
--outputGC          arg  The file to print the gc coverage histogram
--gcWindow          arg  The sliding window size for GC calculations

metabat2 -i /work/workspace/zhurj/project/1_metadata/metapipe/bwa/db/ref.fa -a /work/workspace/zhurj/project/1_metadata/metapipe/metabat/depth.txt -o /work/workspace/zhurj/project/1_metadata/metapipe/metabat/bin/MNHC00230 -m 2000 
6. checkm check MAG
CheckM v1.1.1
checkm lineage_wf -f /work/workspace/zhurj/project/1_metadata/metapipe/checkm/MNC00230/MNC00230 --tab_table -x fa -t 36 /work/workspace/zhurj/project/1_metadata/metapipe/metabat/bin /work/workspace/zhurj/project/1_metadata/metapipe/checkm/MNC00230

7. scaffold filter
Scaftigs（>=500bp）
source activate python3.6
cd /work/workspace/zhurj/project/1_metadata/metapipe/spades/MNC00230
seqtk seq -L 500 scaffolds.fasta > scaffolds_500.fasta

7. gene prediction
GeneMark.hmm version 3.38
INSTALL
cp /work/workspace/zhurj/software/gm_key_64 ~/.gm_key
cd /work/workspace/zhurj/software/MetaGeneMark_linux_64/mgm
#gmhmmp -a -A /work/workspace/zhurj/project/1_metadata/metapipe/metagenemark/MNC00230/MNC00230.faa -d -D /work/workspace/zhurj/project/1_metadata/metapipe/metagenemark/MNC00230/MNC00230.fna -m MetaGeneMark_v1.mod /work/workspace/zhurj/project/1_metadata/metapipe/spades/MNC00230/scaffolds_500.fasta  -o /work/workspace/zhurj/project/1_metadata/metapipe/metagenemark/MNC00230/MNC00230.lst
gmhmmp -a -A /work/workspace/zhurj/project/1_metadata/metapipe/metagenemark/MNC00230/MNC00230.faa -d -D /work/workspace/zhurj/project/1_metadata/metapipe/metagenemark/MNC00230/MNC00230.fna -m MetaGeneMark_v1.mod /work/workspace/zhurj/project/1_metadata/metapipe/spades/MNC00230/scaffolds_500.fasta  -o /work/workspace/zhurj/project/1_metadata/metapipe/metagenemark/MNC00230/MNC00230.lst

cd /work/workspace/zhurj/project/1_metadata/metapipe/metagenemark/MNC00230/
source activate python3.6
# filtering seq with length < 100
seqtk seq -L 100 MNC00230.fna > MNC00230_100.fna
# change long gene title to >gene_1, gene_2, ...
awk -F "|" '{if( $1~/^>/) print $1;else print $0}' MNC00230_100.fna > MNC00230_100_adj.fna
seqtk seq -L 100 MNC00230.faa > MNC00230_100.faa
awk -F "|" '{if( $1~/^>/) print $1;else print $0}' MNC00230_100.faa > MNC00230_100_adj.faa


8. 基因聚类
CD-HIT version 4.8.1 
ref: https://www.jianshu.com/p/d4e4009f0dc0
ref: http://weizhongli-lab.org/cd-hit/
cd /work/workspace/zhurj/project/1_metadata/metapipe/metagenemark/MNC00230
cd-hit-est -i /work/workspace/zhurj/project/1_metadata/metapipe/metagenemark/MNC00230/MNC00230_100_adj.fna -o /work/workspace/zhurj/project/1_metadata/metapipe/cdhit/MNC00230/cdhit.fna -c 0.95 -G 0 -aS 0.9 -g 1 -d 0 -T 36 
cd-hit -i /work/workspace/zhurj/project/1_metadata/metapipe/metagenemark/MNC00230/MNC00230_100_adj.faa -o /work/workspace/zhurj/project/1_metadata/metapipe/cdhit/MNC00230/cdhit.faa -c 0.95 -G 0 -aS 0.9 -g 1 -d 0 -T 36 
protein size
organism	median protein length amino acides
H.sapien	375
D.melanogaster	373
C.elegans	344
S.cerevisiae	379
A.thaliana	356
5 eukaryotes(above)	361
67 bacterial	267
15 archaeal	247
===========================================================
单个粪便样品，宏基因组测序，预测基因数目范围？
CD-HIT 去冗余后基因数目范围
===========================================================

9. gene abundance
# calculate gene with start codon, stop codon
```
ref: https://digitalinsights.qiagen.com/files/user_manuals/MetaGeneMark_User_Manual.pdf
“Genetic code” parameter block
For prokaryotic and phage metagenomes select “Genetic code
11”. This option defines use of ATG, GTG and TTG as starts
codons and TAA, TGA and TAG as stop codons.

For eukaryotic or viral metagenomes or eukaryotic metatranscriptomes select “Genetic code 1”. This option allows ATG
start codon only and TAA, TGA and TAG as stop codons. 
```
http://emboss.sourceforge.net/apps/cvs/emboss/apps/getorf.html

with start codon & stop codon 
cat cdhit.fna | grep -v -P "^>" | grep -P "^ATG|^GTG|^TTG" | grep -P "TAA$|TGA$|TAG$" | wc -l
cat cdhit.fna | grep -v -P "^>" | grep -P "^ATG|^GTG|^TTG" | grep -v -P "TAA$|TGA$|TAG$" | wc -l
cat cdhit.fna | grep -v -P "^>" | grep -v -P "^ATG|^GTG|^TTG" | grep -P "TAA$|TGA$|TAG$" | wc -l
cat cdhit.fna | grep -v -P "^>" | grep -v -P "^ATG|^GTG|^TTG|TAA$|TGA$|TAG$" | wc -l

cd /work/workspace/zhurj/project/1_metadata/metapipe/cdhit/MNC00230
bowtie2-build  cdhit.fna cdhit
bowtie2 -x /work/workspace/zhurj/project/1_metadata/metapipe/cdhit/MNC00230/cdhit \
-1 /work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00230.clean.1.fq.gz \
-2 /work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00230.clean.2.fq.gz \
--end-to-end --sensitive -I 200 -X 400 \
-p 36  --no-sq --no-unal --omit-sec-seq \
-S /work/workspace/zhurj/project/1_metadata/metapipe/genecatalog/MNC00230/bowtie2/MNC00230.sam

···
113 read paired，read reverse strand，mate reverse strand，first in pair
129 read paired，second in pair
*** 137 read paired，mate unmapped，second in pair
145 read paired，read reverse strand，second in pair
147 read paired，read mapped in proper pair，read reverse strand，second in pair
*** 153 read paired，mate unmapped，read reverse strand，second in pair
161 read paired，mate reverse strand，second in pair
163 read paired，read mapped in proper pair，mate reverse strand，second in pair
177 read paired，read reverse strand，mate reverse strand，second in pair
65 read paired，first in pair
*** 73 read paired，mate unmapped，first in pair
81 read paired，read reverse strand，first in pair
83 read paired，read mapped in proper pair，read reverse strand，first in pair
*** 89 read paired，mate unmapped，read reverse strand，first in pair
97 read paired，mate reverse strand， first in pair
99 read paired，read mapped in proper pair，mate reverse strand，first in pair
···

cat file.fa | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'
计算基因，及基因长度
infoseq -auto -nocolumns -delimiter ',' -only -noheading -name -length -sformat pearson /work/workspace/zhurj/project/1_metadata/metapipe/cdhit/MNC00230/cdhit.fna > /work/workspace/zhurj/project/1_metadata/metapipe/cdhit/MNC00230/cdhit_namesize.txt

1. 取 sequenceid + ref_geneid > 1 (paired reads mapped to the same gene ) ， command： uniq -d
2. 统计每个geneid 有多少个read map上 awk -F "," '{print $2}' | sort | uniq -c
cat /work/workspace/zhurj/project/1_metadata/metapipe/genecatalog/MNC00230/bowtie2/MNC00230.sam | grep -v -P "^@" | awk -F "\t" '{print $1","$3}' | uniq -d | awk -F "," '{print $2}' | sort | uniq -c > /work/workspace/zhurj/project/1_metadata/metapipe/genecatalog/MNC00230/bowtie2/preMappGene.txt

3. 挑选 >2 个reads mapped 的gen
cat /work/workspace/zhurj/project/1_metadata/metapipe/genecatalog/MNC00230/bowtie2/preMappGene.txt | awk '{ if ( $(NF-1) > 2 ) print $NF"\t"$(NF-1)}' > /work/workspace/zhurj/project/1_metadata/metapipe/genecatalog/MNC00230/bowtie2/mapGene.txt


/work/workspace/zhurj/project/1_metadata/metapipe/cdhit/MNC00230/cdhit_namesize.txt

#  合并数组： https://blog.csdn.net/qq_41853758/article/details/83280104
source activate python3.6
python /work/workspace/zhurj/script/5_wgspro/16_sumtaxo/test.py /work/workspace/zhurj/project/1_metadata/metapipe/genecatalog/MNC00230/bowtie2/mapGene.txt /work/workspace/zhurj/project/1_metadata/metapipe/cdhit/MNC00230/cdhit_namesize.txt

zcat test.1.fq.gz | awk '{ if ($1 ~ /^@/) print $1;else print $0}' | pigz - >  testadj.1.fq.gz
zcat test.2.fq.gz | awk '{ if ($1 ~ /^@/) print $1;else print $0}' | pigz - >  testadj.2.fq.gz

source activate humann
# kneaddata
# https://blog.csdn.net/woodcorpse/article/details/80413379
kneaddata --input /work/workspace/zhurj/project/1_metadata/metapipe/clean/test.1.fq.gz --input /work/workspace/zhurj/project/1_metadata/metapipe/clean/test.2.fq.gz --output /work/workspace/zhurj/project/1_metadata/metapipe/kneaddata/MNC00230 --bypass-trim --run-trf -db /work/workspace/zhurj/reference/human/kneaddata/20200921 --output-prefix MNC00230 -t 36  -p 2 --remove-intermediate-output

kneaddata --input /work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00230.clean.1.fq.gz --input /work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00230.clean.2.fq.gz --output /work/workspace/zhurj/project/1_metadata/metapipe/kneaddata/MNC00230 -v \
-db /work/workspace/zhurj/reference/human/kneaddata/20200921  \
--trimmomatic  /work/workspace/zhurj/bin/miniconda3/envs/humann/share/trimmomatic-0.39-1 --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" \
-t 36 --bowtie2-options "--very-sensitive --dovetail" --remove-intermediate-output  --output-prefix MNC00230


kneaddata --input /work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00230.clean.1.fq.gz --input /work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00230.clean.2.fq.gz --output /work/workspace/zhurj/project/1_metadata/metapipe/kneaddata/MNC00230 -v \
-db /work/workspace/zhurj/reference/human/kneaddata/20200921  \
--bypass-trim --run-trf \
-t 36 --bowtie2-options "--very-sensitive --dovetail" --remove-intermediate-output  --output-prefix MNC00230


kneaddata_read_count_table --input /work/workspace/zhurj/project/1_metadata/metapipe/kneaddata/MNC00230 --output /work/workspace/zhurj/project/1_metadata/metapipe/kneaddata/kneaddata_read_counts.txt

# sra-stat --quick SRR077487

pigz -dc /work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00231.clean.1.fq.gz | awk 'NR%4==2{c++; l+=length($0)} END{ print "Number of reads: "c; print "Number of bases in reads: "l }'


reformat.sh in1=MNC00230.clean.1.fq.gz in2=MNC00230.clean.2.fq.gz qchist=file.qchist.txt

200000
2.7 
85
80

perl /work/workspace/zhurj/script/3_swgspro/strainscreen/dataqc.pl -i /work/workspace/zhurj/project/1_metadata/metapipe/input/test_r1r2_2c.txt -o /work/workspace/zhurj/project/1_metadata/metapipe/qc --thread 36 --basecutoff 5.4  --cutoffq20 85 --cutoffq30 80

srun -o qc_test.out -e qc_test.err -N 1 -c 20 -p slurm256 bash qc_test.sh &

```
parallel -j 1 --link 'kneaddata -i {1} -i {2} --output /work/workspace/zhurj/project/1_metadata/metapipe/kneaddata/MNC00230/parallel -v \
-db /work/workspace/zhurj/reference/human/kneaddata/20200921  \
--bypass-trim --run-trf \
-t 36 -p 36 --bowtie2-options "--very-sensitive --dovetail" --remove-intermediate-output' \
::: /work/workspace/zhurj/project/1_metadata/metapipe/clean/*.clean.1.fq.gz ::: /work/workspace/zhurj/project/1_metadata/metapipe/clean/*.clean.2.fq.gz
```
parallel -j 1 --link 'kneaddata -i {1} -i {2} -o kneaddata_out/ \
-db /home/shared/bowtiedb/GRCh38_PhiX --trimmomatic 
/usr/local/prg/Trimmomatic-0.36/ \
-t 4 --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" \
--bowtie2-options "--very-sensitive --dovetail" --remove-intermediate-output' \
::: raw_data/*_R1.fastq ::: raw_data/*_R2.fastq

关联分析
REF: Correlation matrix : A quick start guide to analyze, format and visualize a correlation matrix using R software
URL: http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
对于cor()函数中use参数的认识（cor in R）， https://www.douban.com/note/725123309/
source activate /work/workspace/zhurj/bin/miniconda3/envs/R3.6
# 


#-------------------------------------------------------------------------
rawdata extraction
1. input create
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
perl /work/workspace/zhurj/script/5_wgspro/1_filter/mnid2rawdata.pl -i /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/mnid -t 0 -r metagenome -o /work/workspace/zhurj/project/1_metadata/metapipe20200921/input -n sam11

2. raw data QC
perl /work/workspace/zhurj/script/2_metapro/metapipe/dataqc.pl -i /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/sam11_r1r2_2c.txt -o /work/workspace/zhurj/project/1_metadata/metapipe20200921/QC --thread 36 --basecutoff 2.7  --cutoffq20 85 --cutoffq30 80
time srun -o qc_sam11.out -e qc_sam11.err -N 1 -c 20 -p slurm128 bash qc_sam11.sh &

3. data filter
# automatically generate folder: clean
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
perl /work/workspace/zhurj/script/5_wgspro/1_filter/datafilter.pl -i /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/sam11_r1r2_2c.txt -o /work/workspace/zhurj/project/1_metadata/metapipe20200921
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
cd /work/workspace/zhurj/project/1_metadata/metapipe20200921/p
time srun -o datafilter_sam11.out -e datafilter_sam11.err -N 1 -c 20 -p slurm256 bash datafilter_sam11.sh &

4. remove host DNA
# ref： Kneaddata FASTQ header problem： https://forum.biobakery.org/t/kneaddata-fastq-header-problem/482/3
#parallel run kneaddata
cd /work/workspace/zhurj/project/1_metadata/metapipe20200921/clean
find `pwd` | grep 1.fq.gz | sort > ../input/knead_sam11_r1
find `pwd` | grep 2.fq.gz | sort > ../input/knead_sam11_r2
find `pwd` | grep 2.fq.gz | awk -F "[/.]" '{print $(NF-4)}' | sort > ../input/knead_sam11_prefix

cd /work/workspace/zhurj/project/1_metadata/metapipe20200921/p
content in pkneaddata_sam11.sh
```
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel -j 4 --link 'kneaddata -i {1} -i {2} --output /work/workspace/zhurj/project/1_metadata/metapipe20200921/kneaddata/parallel -v \
-db /work/workspace/zhurj/reference/human/kneaddata/20200921  \
--bypass-trim --run-trf --output-prefix {3} \
-t 20 -p 20 --bowtie2-options "--very-sensitive --dovetail" --remove-intermediate-output ' \
:::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_r1 :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_r2 :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_prefix
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann


```
should add para "--reorder "

time srun -o pkneaddata_sam11.out -e pkneaddata_sam11.err -N 2 -c 20 -p slurm256 bash pkneaddata_sam11.sh &
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
kneaddata_read_count_table --input /work/workspace/zhurj/project/1_metadata/metapipe20200921/kneaddata/parallel --output /work/workspace/zhurj/project/1_metadata/metapipe20200921/kneaddata/kneaddata_read_counts_sam11.txt

kneaddataMNC00233.sh
```
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
kneaddata -i /work/workspace/zhurj/project/1_metadata/metapipe20200921/clean/MNC00233.clean.1.fq.gz \
-i /work/workspace/zhurj/project/1_metadata/metapipe20200921/clean/MNC00233.clean.2.fq.gz \
--output /work/workspace/zhurj/project/1_metadata/metapipe20200921/kneaddata/parallel -v \
-db /work/workspace/zhurj/reference/human/kneaddata/20200921  \
--bypass-trim --run-trf --output-prefix MNC00233 \
-t 20 -p 20 --bowtie2-options "--very-sensitive --dovetail" --remove-intermediate-output
kneaddata -i /work/workspace/zhurj/project/1_metadata/metapipe20200921/clean/MNC00231.clean.1.fq.gz \
-i /work/workspace/zhurj/project/1_metadata/metapipe20200921/clean/MNC00231.clean.2.fq.gz \
--output /work/workspace/zhurj/project/1_metadata/metapipe20200921/kneaddata/parallel -v \
-db /work/workspace/zhurj/reference/human/kneaddata/20200921  \
--bypass-trim --run-trf --output-prefix MNC00231 \
-t 20 -p 20 --bowtie2-options "--very-sensitive --dovetail" --remove-intermediate-output
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

```
cd /work/workspace/zhurj/project/1_metadata/metapipe20200921/p
time srun -o kneaddataMNC00233.out -e kneaddataMNC00233.err -N 1 -c 20 -p slurm256 bash kneaddataMNC00233.sh &

cd /work/workspace/zhurj/project/1_metadata/metapipe20200921/kneaddata/parallel
sed -i 's/\#0\/1//g' MNC00233.repeats.removed.1.fastq
sed -i 's/\#0\/2//g' MNC00233.repeats.removed.2.fastq


5. remove duplication(先去重，再去宿主)
# cluspify.sh 
# clumpify.sh -Xmx20g dedupe dupedist=12000 in1=/work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00230.clean.1.fq.gz in2=/work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00230.clean.2.fq.gz out1=/work/workspace/zhurj/project/1_metadata/metapipe20200921/test/MNC00230.dedup.1.fq.gz out2=/work/workspace/zhurj/project/1_metadata/metapipe20200921/test/MNC00230.dedup.2.fq.gz 

```
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel -j 2 --link 'clumpify.sh -Xmx80g dedupe dupedist=12000 in1={1} in2={2} \
out1=/work/workspace/zhurj/project/1_metadata/metapipe20200921/dedup/{3}.1.fq.gz out2=/work/workspace/zhurj/project/1_metadata/metapipe20200921/dedup/{3}.2.fq.gz' \
:::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_r1 :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_r2 :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_prefix
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

```
cd /work/workspace/zhurj/project/1_metadata/metapipe20200921/p
time srun -o pclumpify_sam11.out -e pclumpify_sam11.err -N 2 -c 20 -p slurm256 bash pclumpify_sam11.sh &


6. metaphlan
test
metaphlan /work/workspace/zhurj/project/1_metadata/metapipe20200921/kneaddata/parallel/MNC00233.repeats.removed.1.fastq,/work/workspace/zhurj/project/1_metadata/metapipe20200921/kneaddata/parallel/MNC00233.repeats.removed.2.fastq --bowtie2out /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/test/metagenome.bowtie2.bz2 \
--input_type fastq --nproc 36 \
--sample_id_key test -o /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/test/profile.txt

pmetaphlan3_sam11.sh
```
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel -j 1 'mkdir /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/{} ' :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_prefix
parallel -j 1 --link 'metaphlan {1},{2} --bowtie2out /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/{3}/metagenome.bowtie2.bz2 \
--input_type fastq --nproc 36 \
--sample_id_key {3} -o /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/{3}/profile.txt ' \
:::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_r1 :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_r2 :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_prefix
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

```
cd /work/workspace/zhurj/project/1_metadata/metapipe20200921/p
time srun -o pmetaphlan3_sam11.out -e pmetaphlan3_sam11.err -N 1 -c 20 -p slurm128 bash pmetaphlan3_sam11.sh &


metaphlan /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00231/metagenome.bowtie2.bz2 --input_type bowtie2out --tax_lev 's' --nproc 36 --sample_id_key MNC00231 -o /work/workspace/zhurj/project/1_metadata/metapipe/metaphlan3/MNC00231/species.txt
pmetaphlan3_sam11_species.sh
```
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel -j 1 --link 'metaphlan /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/{1}/metagenome.bowtie2.bz2 \
--input_type bowtie2out --tax_lev 's' --nproc 36 \
--sample_id_key {1} -o /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/{1}/species.txt ' \
:::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_prefix
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

```
cd /work/workspace/zhurj/project/1_metadata/metapipe20200921/p
time srun -o pmetaphlan3_sam11_species.out -e pmetaphlan3_sam11_species.err -N 1 -c 20 -p slurm128 bash pmetaphlan3_sam11_species.sh &

mkdir /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/species
parallel -j 1 --link 'ln -s /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/{1}/species.txt /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/species/{1}.txt' \
:::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_prefix

merge_metaphlan_tables.py  /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/species/*.txt > /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/merge/species_merge_sam11.tsv
#perl /work/workspace/zhurj/software/microbiome_helper-master/metaphlan_to_stamp.pl /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/merge/species_merge_sam11.tsv >/work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/merge/species_merge_sam11.spf
sed -i '1d' /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/merge/species_merge_sam11.tsv
# sed -i 'nd' filename 删除第n行
# sed -i '$d' filename 删除最后一行
cat /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/merge/species_merge_sam11.tsv | awk '{$2=null;print $0}' > /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/merge/speciea_adj.txt
sed -i -r 's/\s+/\t/g' /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/merge/speciea_adj.txt

#STAMP input /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/merge/speciea_adj.txt

# Lefse in
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/merge/speciea_adj.txt /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/merge/meta /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/metaphlan3/lefse_in
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

```
source activate /work/workspace/zhurj/bin/miniconda3/envs/py2
lefse-format_input.py /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/metaphlan3/lefse_in /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/metaphlan3/lefse_format.in -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/metaphlan3/lefse_format.in /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/metaphlan3/lefse_run.res -a 0.05 -w 0.05 -l 2 -y 0
#所有级别的物种输出的结果
# lefse-plot_res.py /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/metaphlan3/lefse_run.res  /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/metaphlan3/lefse_run.png --feature_font_size 8 --width 10 --format png
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/py2

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py /work/workspace/zhurj/project/1_metadata/mouse22/metaphlan3/merge/genus_merge_adj_20201010.txt /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group /work/workspace/zhurj/project/1_metadata/mouse22/lefse/metaphlan3/genus/lefse_in
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
# genus difference
source activate /work/workspace/zhurj/bin/miniconda3/envs/py2
lefse-format_input.py /work/workspace/zhurj/project/1_metadata/mouse22/lefse/metaphlan3/genus/lefse_in /work/workspace/zhurj/project/1_metadata/mouse22/lefse/metaphlan3/genus/lefse_format.in -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/zhurj/project/1_metadata/mouse22/lefse/metaphlan3/genus/lefse_format.in /work/workspace/zhurj/project/1_metadata/mouse22/lefse/metaphlan3/genus/lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/py2
```


'''
cd /work/workspace/zhurj/project/1_metadata/p/p
time srun -o sgb_lefse.out -e sgb_lefse.err -N 2 -c 20 -p slurm128 bash sgb_lefse.sh &

parallel --link 'mkdir -p {1}/{2}_dg ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/lefse/SGB ::: phylum genus species
source activate python3.6
python /work/workspace/zhurj/script/2_metapro/metapipe/orgSGBdg.py /work/workspace/zhurj/project/1_metadata/mouse22/SGB/taxonomy_in \
/work/workspace/zhurj/project/1_metadata/mouse22/SGB/taxonomy_species_dg_out species

python /work/workspace/zhurj/script/2_metapro/metapipe/orgSGBdg.py /work/workspace/zhurj/project/1_metadata/mouse22/SGB/taxonomy_in \
/work/workspace/zhurj/project/1_metadata/mouse22/SGB/taxonomy_genus_dg_out genus

python /work/workspace/zhurj/script/2_metapro/metapipe/orgSGBdg.py /work/workspace/zhurj/project/1_metadata/mouse22/SGB/taxonomy_in \
/work/workspace/zhurj/project/1_metadata/mouse22/SGB/taxonomy_phylum_dg_out phylum

python /work/workspace/zhurj/script/2_metapro/metapipe/orgSGBdg.py /work/workspace/zhurj/project/1_metadata/mouse22/SGB/taxonomy_in \
/work/workspace/zhurj/project/1_metadata/mouse22/SGB/taxonomy_all_dg_out all

python /work/workspace/zhurj/script/2_metapro/metapipe/orgSGBdg_norm.py /work/workspace/zhurj/project/1_metadata/mouse22/SGB/taxonomy_in /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup

parallel --link 'python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py {1}/SGB/taxonomy_{2}_dg_out {1}/input/meta_group {1}/lefse/SGB/{2}_dg/lefse_in' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: phylum genus species 
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
# genus difference

source activate /work/workspace/zhurj/bin/miniconda3/envs/py2
parallel --link -j 1 'lefse-format_input.py {1}/lefse/SGB/{2}_dg/lefse_in {1}/lefse/SGB/{2}_dg/lefse_format.in -c 1 -u 2 -f r -o 1000000 && \
run_lefse.py {1}/lefse/SGB/{2}_dg/lefse_format.in {1}/lefse/SGB/{2}_dg/lefse_run.res -a 0.01 -w 0.01 -l 2 -y 1 ' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: phylum genus species
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/py2

source activate python3.6
parallel --link 'python /work/workspace/zhurj/lib/python/script/lefse_fitler.py {1}/lefse/SGB/{2}_dg/lefse_run.res {1}/lefse/SGB/{2}_dg/filter_lefse_run.res ' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: phylum genus species
source deactivate python3.6


STAMP 不需要这个格式
import pandas as pd
df = pd.read_table("/work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/merge/speciea_adj.txt",sep='\s+',index_col=0,engine='python')
df2 = df.stack()
df3 = df2.unstack(0)
df3.to_csv("/work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/merge/speciea_adj_format.txt",sep='\t')

# for rarefaction
import pandas as pd
df = pd.read_table("/work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/merge/speciea_adj.txt",sep='\s+',index_col=0,engine='python')
df2 = df.stack()
df3 = df2.unstack(0)
df4 = df3*10000 
df5 = df4.astype(int)
df5.to_csv("/work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/merge/speciea_rarefaction_in.txt",sep='\t')

rarefaction
source activate /work/workspace/zhurj/bin/miniconda3/envs/R3.6
# rarefy: Rarefaction Species Richness https://rdrr.io/rforge/vegan/man/rarefy.html
library(vegan)
df <- read.table("/work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/merge/speciea_rarefaction_in.txt",sep='\t',header=TRUE,row.names=1)
S <- specnumber(df)
(raremax <- min(rowSums(df)))
Srare <- rarefy(df,raremax)
plot(S,Srare, xlab="Observed No. of Species", ylab="Rarefied No. of Species")
abline(0,1)
rarecurve(df,step = 200, sample = raremax, col = "blue", cex = 0.6)

rarecurve(df,step = 200, sample = raremax, col = "blue", cex = 0.6,xlim=c(0, 20000))

7. humann
parallel --xapply -j 2 'zcat /work/workspace/zhurj/project/1_metadata/metapipe20200921/clean/{1}.clean.1.fq.gz /work/workspace/zhurj/project/1_metadata/metapipe20200921/clean/{1}.clean.2.fq.gz | /work/workspace/zhurj/bin/pigz - > /work/workspace/zhurj/project/1_metadata/metapipe20200921/clean/{1}.merge.fq.gz' :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_prefix

```
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel --xapply -j 2 'humann --threads 20 --input /work/workspace/zhurj/project/1_metadata/metapipe20200921/clean/{1}.merge.fq.gz --output /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3/humann3_out/{1}' :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_prefix
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

```
cd /work/workspace/zhurj/project/1_metadata/metapipe20200921/p
time srun -o phumann3_sam11.out -e phumann3_sam11.err -N 1 -c 20 -p slurm256 bash phumann3_sam11.sh &

source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel --xapply  'humann_join_tables -s --input {2}/humann3_out --file_name {1} \
--output {2}/humann3_final_out/humann_{1}.tsv' ::: pathcoverage pathabundance genefamilies ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3

parallel --xapply  'humann_renorm_table --input {2}/humann3_final_out/humann_{1}.tsv --units relab \
--output {2}/humann3_final_out/{1}_relab.tsv' ::: pathcoverage pathabundance genefamilies ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3

humann_regroup_table -i humann_genefamilies.tsv -g uniref90_ko --output uniref90_ko_known.tsv
humann_renorm_table -i humann_genefamilies.tsv -u relab --output genefamilies_relab.tsv
humann_renorm_table -i humann_genefamilies.tsv -u relab --output genefamilies_relab.tsv

#humann_join_tables -s --input /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3/humann3_out/ --file_name pathcoverage --output humann_final_out/humann_pathcoverage.tsv
#humann_join_tables -s --input /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3/humann3_out/ --file_name pathabundance --output humann_final_out/humann_pathcoverage.tsv
#humann_join_tables -s --input /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3/humann3_out/ --file_name genefamilies --output humann_final_out/humann_genefamilies.tsv
humann_join_tables -s --input /work/workspace/zhurj/project/1_metadata/mouse22/humann3/ --file_name metacyc_genefamilies --output /work/workspace/zhurj/project/1_metadata/mouse22/humann3/merge/uniref90_genefamilies.tsv

# Gene Family   MNC00223.merge_Abundance-RPKs
# Pathway       MNC00223.merge_Coverage 
# Gene Family   MNC00223.merge_Abundance-RPKs

cd /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3/humann3_final_out
sed -i -r 's/\#\s+//;s/\.merge_Abundance-RPKs//g;s/Gene Family/GeneFamily/' /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3/humann3_final_out/humann_genefamilies.tsv
sed -i -r 's/\#\s+//;s/\.merge_Abundance//g' /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3/humann3_final_out/humann_pathabundance.tsv
sed -i -r 's/\#\s+//;s/\.merge_Coverage//g' /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3/humann3_final_out/humann_pathcoverage.tsv

# Lefse in
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3/humann3_final_out/humann_pathabundance.tsv /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/merge/meta /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/humann3/pathabundance
python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3/humann3_final_out/humann_genefamilies.tsv /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/merge/meta /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/humann3/genefamilies
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

```
source activate /work/workspace/zhurj/bin/miniconda3/envs/py2
lefse-format_input.py /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/humann3/pathabundance /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/humann3/pathabundanc_format.in -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/humann3/pathabundanc_format.in /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/humann3/pathabundanc.res -a 0.05 -w 0.05 -l 2 -y 1
#所有级别的物种输出的结果
# lefse-plot_res.py /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/metaphlan3/lefse_run.res  /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/metaphlan3/lefse_run.png --feature_font_size 8 --width 10 --format png

lefse-format_input.py /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/humann3/genefamilies /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/humann3/genefamilies_format.in -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/humann3/genefamilies_format.in /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/humann3/genefamilies.res -a 0.05 -w 0.05 -l 2 -y 1

source deactivate /work/workspace/zhurj/bin/miniconda3/envs/py2



# kegg
#/work/workspace/liangyj/script/tools/kegg
parallel 'mkdir -p {1}/{2}' ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3/humann3_kegg :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_prefix
parallel --link 'mv {1}/{2}/uniref90_ko_known_pathabundance.tsv {1}/{2}/{2}_uniref90_ko_known_pathabundance.tsv' ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3/humann3_kegg :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_prefix

```
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel --link -j 1 'humann_regroup_table -i {1}/humann3_out/{2}/{2}.merge_genefamilies.tsv -g uniref90_ko --output {1}/humann3_kegg/{2}/{2}_uniref90_ko_known.tsv'  ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3 :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_prefix

parallel --link -j 1 'humann -i {1}/humann3_kegg/{2}/uniref90_ko_known.tsv --pathways-database /work/workspace/zhurj/database/humann/current/kegg/keggc --output {1}/humann3_kegg --threads 36 '  ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3 

parallel -j 1 --link 'humann_join_tables -s --input {1}/humann3_kegg/ --file_name pathabundance --output {1}/humann3_final_out/kegg_pathabundance.tsv'  ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3 :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_prefix
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

```
cd /work/workspace/zhurj/project/1_metadata/metapipe20200921/p
time srun -o pkegg_sam11.out -e pkegg_sam11.err -N 1 -c 20 -p slurm256 bash pkegg_sam11.sh &

sed -i -r 's/_uniref90_ko_known_pathabundance//g;s/# Pathway/Pathway/' kegg_pathabundance.tsv

# Lefse in
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3/humann3_final_out/kegg_pathabundance.tsv /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/merge/meta /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/humann3/kegg_pathabundance
source activate /work/workspace/zhurj/bin/miniconda3/envs/py2
lefse-format_input.py /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/humann3/kegg_pathabundance /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/humann3/kegg_pathabundanc_format.in -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/humann3/kegg_pathabundanc_format.in /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/humann3/kegg_pathabundanc.res -a 0.05 -w 0.05 -l 2 -y 1
echo "" > diff_kegg.txt
cat kegg_pathabundanc.res | awk '{if ($NF != "-" && $2 >= 2 ) print $0}' >> diff_kegg.txt
cat kegg_pathabundanc.res | awk '{if ($NF != "-" && $2 <= -2 ) print $0}' >> diff_kegg.txt

source activate python3.6
parallel -j 1 --link 'python {2}/add_kegg_description.py {1}/diff_kegg.txt {1}/diff_kegg_desc.txt ' \
::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/lefse/humann3 ::: /work/workspace/zhurj/script/2_metapro/metapipe

test
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel --link -j 1 'humann --threads 20 --input {2}/kneaddata/parallel/{1}.repeats.removed.1.fastq --output {2}/humann3/test --output-basename {1}_metacyc && \
humann_regroup_table -i {2}/humann3/test/{1}_metacyc_genefamilies.tsv -g uniref90_go --output {2}/humann3/test/{1}_uniref90_go.tsv && \
humann_regroup_table -i {2}/humann3/test/{1}_metacyc_genefamilies.tsv -g uniref90_level4ec --output {2}/humann3/test/{1}_uniref90_level4ec.tsv && \
humann_regroup_table -i {2}/humann3/test/{1}_metacyc_genefamilies.tsv -g uniref90_pfam --output {2}/humann3/test/{1}_uniref90_pfam.tsv && \
humann_regroup_table -i {2}/humann3/test/{1}_metacyc_genefamilies.tsv -g uniref90_eggnog --output {2}/humann3/test/{1}_uniref90_eggnog.tsv && \
humann_regroup_table -i {2}/humann3/test/{1}_metacyc_genefamilies.tsv -g uniref90_ko --output {2}/humann3/test/{1}_uniref90_ko.tsv && \
humann -i {2}/humann3/test/{1}_uniref90_ko.tsv --pathways-database {1} --output {2}/humann3/test --threads 36 --output-basename {1}_kegg \
' ::: MNC00233 ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921 ::: /work/workspace/zhurj/database/humann/current/kegg/keggc
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann
cd /work/workspace/zhurj/project/1_metadata/metapipe20200921/p
time srun -o test_humann.out -e test_humann.err -N 1 -c 20 -p slurm256 -w mnclient01 bash test_humann.sh &

time humann -i /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3/test/MNC00233_uniref90_ko.tsv \
--pathways-database /work/workspace/zhurj/database/humann/current/kegg/keggc \
--output /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3/test \
--threads 36 --output-basename MNC00233_kegg &

test --------------------------------------
parallel --link -j 1 ' \
humann_regroup_table -i {2}/humann3/{1}/{1}_metacyc_genefamilies.tsv -g uniref90_go --output {2}/humann3/{1}/{1}_uniref90_go.tsv && \
humann_regroup_table -i {2}/humann3/{1}/{1}_metacyc_genefamilies.tsv -g uniref90_level4ec --output {2}/humann3/{1}/{1}_uniref90_level4ec.tsv && \
humann_regroup_table -i {2}/humann3/{1}/{1}_metacyc_genefamilies.tsv -g uniref90_pfam --output {2}/humann3/{1}/{1}_uniref90_pfam.tsv && \
humann_regroup_table -i {2}/humann3/{1}/{1}_metacyc_genefamilies.tsv -g uniref90_eggnog --output {2}/humann3/{1}/{1}_uniref90_eggnog.tsv && \
humann_regroup_table -i {2}/humann3/{1}/{1}_metacyc_genefamilies.tsv -g uniref90_ko --output {2}/humann3/{1}/{1}_uniref90_ko.tsv && \
parallel --link -j 1 ' \
humann -i {2}/humann3/{1}/{1}_uniref90_ko.tsv --pathways-database {3} --output {2}/humann3/{1} --threads 36 --output-basename {1}_kegg \
' ::: MNA00233 ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: /work/workspace/zhurj/database/humann/current/kegg/keggc

humannmapping_sam11.sh
```
MNC00230.merge_genefamilies.tsv
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann

parallel --link -j 2 'humann --threads 20 --input {2}/kneaddata/parallel/{1}.repeats.removed.1.fastq --output {2}/humann3/test --output-basename {1}_metacyc ' \
:::: MNC00233 ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921

parallel --link -j 2 'humann --threads 20 --input {2}/kneaddata/{1}/{1}.fastq --output {2}/humann3/{1} --output-basename {1}_metacyc ' \
:::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22 ::: /work/workspace/zhurj/project/1_metadata/mouse22

parallel --link -j 1 'humann_regroup_table -i {1}/humann3_out/{2}/{2}.merge_genefamilies.tsv -g uniref90_go --output {1}/humann3_kegg/{2}/{2}_uniref90_go.tsv'  ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3 :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_prefix
parallel --link -j 1 'humann_regroup_table -i {1}/humann3_out/{2}/{2}.merge_genefamilies.tsv -g uniref90_level4ec --output {1}/humann3_kegg/{2}/{2}_uniref90_level4ec.tsv'  ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3 :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_prefix
parallel --link -j 1 'humann_regroup_table -i {1}/humann3_out/{2}/{2}.merge_genefamilies.tsv -g uniref90_pfam --output {1}/humann3_kegg/{2}/{2}_uniref90_pfam.tsv'  ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3 :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_prefix
parallel --link -j 1 'humann_regroup_table -i {1}/humann3_out/{2}/{2}.merge_genefamilies.tsv -g uniref90_eggnog --output {1}/humann3_kegg/{2}/{2}_uniref90_eggnog.tsv'  ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3 :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_prefix
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

```
cd /work/workspace/zhurj/project/1_metadata/metapipe20200921/p
time srun -o humannmapping_sam11.out -e humannmapping_sam11.err -N 1 -c 20 -p slurm256 bash humannmapping_sam11.sh &

parallel -j 1 --link 'humann_join_tables -s --input {1}/humann3_kegg/ --file_name pathabundance --output {1}/humann3_final_out/kegg_pathabundance.tsv'  ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3 :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_prefix

source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel --xapply  'humann_renorm_table --input {2}/humann3_final_out/humann_{1}.tsv --units relab \
--output {2}/humann3_final_out/{1}_relab.tsv' ::: pathcoverage pathabundance genefamilies ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3

source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann


ref: 
/work/workspace/zhurj/database/humann/current/kegg/keggc
```
pathway_pipeline.py -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/EC/pred_metagenome_strat.tsv \
  -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/metacyc \
  -p 20

  add_descriptions.py -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/metacyc/path_abun_unstrat.tsv \
  -m METACYC \
  -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/metacyc/path_abun_unstrat_descrip.tsv
 
  pathway_pipeline.py -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/KO/pred_metagenome_strat.tsv \
  -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/kegg_module \
  --no_regroup \
  --map /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/ref/pathway_mapfiles/KEGG_modules_to_KO.tsv \
  -p 20

  add_descriptions.py -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/kegg_module/path_abun_unstrat.tsv \
  --custom_map_table /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/ref/description_mapfiles/KEGG_modules_info_adj.tsv \
  -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/kegg_module/path_abun_unstrat_descrip.tsv

  pathway_pipeline.py -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/KO/pred_metagenome_strat.tsv \
  -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/kegg_path/all \
  --no_regroup \
  --map /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/ref/pathway_mapfiles/KEGG_pathways_to_KO.tsv \
  -p 20

  add_descriptions.py -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/kegg_path/all/path_abun_unstrat.tsv \
  --custom_map_table /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/ref/description_mapfiles/KEGG_pathways_info.tsv.gz \
  -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/kegg_path/all/path_abun_unstrat_descrip.tsv

  pathway_pipeline.py -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/KO/pred_metagenome_strat.tsv \
  -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/kegg_path/level2 \
  --no_regroup \
  --map /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/ref/kegg/kegg_path_ko_l2.txt \
  -p 20

  add_descriptions.py -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/kegg_path/level2/path_abun_unstrat.tsv \
  --custom_map_table /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/ref/kegg/kegg_path_description.txt \
  -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/kegg_path/level2/path_abun_unstrat_descrip.tsv

  pathway_pipeline.py -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/KO/pred_metagenome_strat.tsv \
  -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/kegg_path/level3 \
  --no_regroup \
  --map /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/ref/kegg/kegg_path_ko_l3.txt \
  -p 20

  add_descriptions.py -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/kegg_path/level3/path_abun_unstrat.tsv \
  --custom_map_table /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/ref/kegg/kegg_path_description.txt \
  -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/kegg_path/level3/path_abun_unstrat_descrip.tsv
```

source activate humann3
humann2_infer_taxonomy -i genefamilies.tsv -r uniref90


#-------------------------------------------------------------------------
组装
1. contig assemble
pspedes_sam11.sh
```
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --xapply -j 1 'mkdir -p /work/workspace/zhurj/project/1_metadata/metapipe20200921/spades/{3} ' :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_prefix
parallel --xapply -j 1 'python /work/workspace/moonis/bin/genoassemble/soft/SPAdes-3.13.1-Linux/bin/spades.py --meta -1 {1} -2 {2} -t 36 -o /work/workspace/zhurj/project/1_metadata/metapipe20200921/spades/{3}' :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_r1 :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_r2 :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_prefix
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

```
cd /work/workspace/zhurj/project/1_metadata/metapipe20200921/p
time srun -o pspedes_sam11.out -e pspedes_sam11.err -N 1 -c 20 -p slurm256 bash pspedes_sam11.sh &


2. bin
#bwa map
```
/work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel -j 1 -k 'mkdir -p {1}/{2} ' ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/metabat :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_prefix
parallel -j 1 -k 'mkdir -p {1}/{2} ' ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/bwa :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_prefix
parallel -j 1 -k 'bwa index {1}/spades/{2}/contigs.fasta && \
bwa mem -t {3} {1}/spades/{2}/contigs.fasta {4} {5} > {1}/bwa/{2}/{2}.sam && \
samtools view -@ {5} -S -b {1}/bwa/{2}/{2}.sam -o {1}/bwa/{2}/{2}.bam && \
samtools sort -@ {5} {1}/bwa/{2}/{2}.bam -o {1}/bwa/{2}/{2}.sorted.bam && \
jgi_summarize_bam_contig_depths --outputDepth {1}/metabat/{2}/depth.txt  {1}/bwa/{2}/{2}.sorted.bam && \
metabat2 -i {1}/spades/{2}/contigs.fasta -a {1}/metabat/{2}/depth.txt -o {1}/metabat/{2}/{2} -m 2000' \
::: /work/workspace/zhurj/project/1_metadata/metapipe20200921 :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_prefix ::: 40 :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_r1 :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_r2
/work/workspace/zhurj/bin/miniconda3/envs/python3.6

```
cd /work/workspace/zhurj/project/1_metadata/metapipe20200921/p
time srun -o metabat_sam11.out -e metabat_sam11.err -N 1 -c 20 -p slurm256 bash metabat_sam11.sh &


metabattest.sh
test
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel -j 1 -k 'mkdir -p {1}/{2} ' ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/metabat :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_prefix
parallel -j 1 -k 'mkdir -p {1}/{2} ' ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/bwa :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/knead_sam11_prefix
parallel -j 1 -k 'bwa index {1}/spades/{2}/contigs.fasta && \
bwa mem -t {3} {1}/spades/{2}/contigs.fasta {4} {5} > {1}/bwa/{2}/{2}.sam && \
samtools view -@ {5} -S -b {1}/bwa/{2}/{2}.sam -o {1}/bwa/{2}/{2}.bam && \
samtools sort -@ {5} {1}/bwa/{2}/{2}.bam -o {1}/bwa/{2}/{2}.sorted.bam && \
jgi_summarize_bam_contig_depths --outputDepth {1}/metabat/{2}/depth.txt  {1}/bwa/{2}/{2}.sorted.bam && \
metabat2 -i {1}/spades/{2}/contigs.fasta -a {1}/metabat/{2}/depth.txt -o {1}/metabat/{2}/MNHC00230 -m 2000' \
::: /work/workspace/zhurj/project/1_metadata/metapipe20200921 ::: MNC00223 ::: 40 ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/clean/MNC00223.clean.1.fq.gz ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/clean/MNC00223.clean.2.fq.gz
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
cd /work/workspace/zhurj/project/1_metadata/metapipe20200921/p
time srun -o metabattest.out -e metabattest.err -N 1 -c 20 -p slurm256 bash metabattest.sh &


```
jgi_summarize_bam_contig_depths --outputDepth /work/workspace/zhurj/project/1_metadata/metapipe/metabat/depth.txt  /work/workspace/zhurj/project/1_metadata/metapipe/bwa/MNC00230.sorted.bam
--minContigLength   arg  The mimimum length of contig to include for mapping and shredding
--minContigDepth    arg  The minimum depth along contig at which to break the contig
--outputGC          arg  The file to print the gc coverage histogram
--gcWindow          arg  The sliding window size for GC calculations
```

python /work/workspace/zhurj/script/2_metapro/metapipe/metapipe_ref.py /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid /work/workspace/zhurj/project/1_metadata/mouse22
cd /work/workspace/zhurj/project/1_metadata/metatest/p

sed -n '20,29p' > kneaddata.sh

parallel 'time srun -o {}/kneaddata.out -e {}/kneaddata.err -N 1 -c 20 -p slurm256 bash {}/kneaddata.sh & ' ::: /work/workspace/zhurj/project/1_metadata/metatest/p 


parallel echo  ::: /work/workspace/zhurj/project/1_metadata/metatest/metaphlan3 ::: all phylum genus species
humann_regroup_table -i MNC00233.merge_genefamilies.tsv -g uniref90_ko --output  MNC00233_uniref90_ko.tsv
humann -i MNC00233_uniref90_ko.tsv --pathways-database /work/workspace/zhurj/database/humann/current/kegg/keggc --output /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3/humann3_out/MNC00233 --threads 20 --output-basename MNC00233_kegg

parallel -j -1 --link 'metaphlan {1}/{2}/metagenome.bowtie2.bz2 --input_type bowtie2out --tax_lev 'g' --nproc 40 --sample_id_key {2} -o {1}/{2}/genus.txt '  ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3 :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/mnid
parallel -j 1 --link 'ln -s {1}/{2}/genus.txt {3}/{2}.txt ' ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3 :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/mnid ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/genus
parallel -j 1 --link ' merge_metaphlan_tables.py  {1}/{2}/*.txt > {1}/merge/{2}_merge.txt '    ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3 ::: genus
parallel -j 1 --link ' merge_metaphlan_tables.py  {1}/{2}/*.txt > {1}/merge/{2}_merge.txt '    ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3 ::: species

parallel -j -1 --link 'num=$(cat {1}/{2}/profile.txt | wc -l) && echo -e "{2}\t$num" ' \
::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3 :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/mnid

sed -n '31,57p' metapipe_ref.sh > metaphlan.sh
parallel 'time srun -o {}/metaphlan.out -e {}/metaphlan.err -N 1 -c 20 -p slurm256 bash {}/metaphlan.sh & ' ::: /work/workspace/zhurj/project/1_metadata/metatest/p
sed -n '49,64p' metapipe_ref.sh > metalefse.sh
parallel 'time srun -o {}/metalefse.out -e {}/metalefse.err -N 1 -c 20 -p slurm256 bash {}/metalefse.sh & ' ::: /work/workspace/zhurj/project/1_metadata/metatest/p


parallel -j 1 --link ' merge_metaphlan_tables.py  {1}/{2}/*.txt > {1}/merge/{2}_merge.txt '    ::: /work/workspace/zhurj/project/1_metadata/metatest/metaphlan3 ::: all phylum genus species
sed -n '67,77p' metapipe_ref.sh > humann.sh
parallel 'time srun -o {}/humann.out -e {}/humann.err -N 1 -c 20 -p slurm256 bash {}/humann.sh & ' ::: /work/workspace/zhurj/project/1_metadata/metatest/p

parallel -j 2 --link 'humann -i {1}/{2}/{2}_uniref90_ko.tsv --pathways-database /work/workspace/zhurj/database/humann/current/kegg/keggc --output {1}/kegg --threads 20 '   ::: /work/workspace/zhurj/project/1_metadata/metatest/humann3 :::: /work/workspace/zhurj/project/1_metadata/metatest/input/m2_kneaddata_id

parallel --link -j 1 'humann -i {1}/humann3_kegg/{2}/uniref90_ko_known.tsv --pathways-database /work/workspace/zhurj/database/humann/current/kegg/keggc --output {1}/humann3_kegg --threads 36 '  ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/humann3 

parallel --link -j 1 'metaphlan {1}/MNH04863.1.fq.gz,{2}/MNH04863.1.fq.gz --bowtie2out {2}/metagenome.bowtie2.bz2 \
--input_type fastq --nproc 40 \
--sample_id_key MNH04863 -o {2}/profile.txt' ::: /work/rawdata/fastq/genome/MNH/MNH048/MNH04863/E1/L1/S1 ::: /work/workspace/zhurj/project/2_swgs/MNH04863/metaphlan3_20200929
time srun -o metaphlan3.out -e metaphlan3.err -N 1 -c 20 -p slurm256 bash metaphlan3.sh &

merge_metaphlan_tables.py  /work/workspace/zhurj/project/1_metadata/mouse22/metaphlan3/genus/*.txt > /work/workspace/zhurj/project/1_metadata/mouse22/metaphlan3/merge/genus_test_merge.txt 

source activate /work/workspace/zhurj/bin/miniconda3/envs/kraken2
parallel --link -j -1 'kraken2 --db /work/workspace/zhurj/reference/kraken2/20200929/current --threads 36  \
--report {1}/MNH04863.kraken2.report --output {1}/MNH04863.kraken2.output --paired {2}/MNH04863.1.fq.gz {2}/MNH04863.2.fq.gz' \
::: /work/workspace/zhurj/project/2_swgs/MNH04863/kraken2 ::: /work/rawdata/fastq/genome/MNH/MNH048/MNH04863/E1/L1/S1
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/kraken2
time srun -o kraken.out -e kraken.err -N 1 -c 20 -p slurm256 bash kraken.sh &
# 第二步，使用Braken校正
# 第三步，将Braken的report格式转换成--use-mpa-style格式
source activate /work/workspace/zhurj/bin/miniconda3/envs/kraken2
parallel 'bracken -d /work/workspace/zhurj/reference/kraken2/20200929/current -i {1}/MNH04863.kraken2.report -o {1}/MNH04863.S.bracken -w {1}/MNH04863.S.bracken.report -r 150 -l S ' ::: /work/workspace/zhurj/project/2_swgs/MNH04863/kraken2
parallel 'kreport2mpa.py -r {1}/MNH04863.S.bracken.report -o {1}/MNH04863.new.report ' ::: /work/workspace/zhurj/project/2_swgs/MNH04863/kraken2
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/kraken2

time srun -o bracken.out -e bracken.err -N 1 -c 20 -p slurm256 bash bracken.sh &
#--------------------------------------------------------------------------------------------------------------
python /work/workspace/zhurj/script/2_metapro/metapipe/metapipe_ref.py /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22 /work/workspace/zhurj/project/1_metadata/mouse22
cd /work/workspace/zhurj/project/1_metadata/mouse22/p
time srun -o metapipe.out -e metapipe.err -N 1 -c 20 -p slurm256 -w mnclient01 bash metapipe.sh &
拆分分析流程
qc+filter： qcfilter.sh
sed -n '1,18p' metapipe_ref.sh > qcfilter.sh
kneaddata： kneaddata.sh
sed -n '21,31p' metapipe_ref.sh > kneaddata.sh
metaphlan: metaphlan3.sh
sed -n '34,45p' metapipe_ref.sh > metaphlan3-1.sh
sed -n '47,55p' metapipe_ref.sh > metaphlan3-2.sh
sed -n '57,62p' metapipe_ref.sh > metaphlan3-3.sh
humann: humann3.sh
sed -n '66,76p' metapipe_ref.sh > humann3.sh


cd /work/workspace/zhurj/project/1_metadata/mouse22/p
time srun -o qcfilter.out -e qcfilter.err -N 1 -c 20 -p slurm256 bash qcfilter.sh &
time srun -o kneaddata.out -e kneaddata.err -N 1 -c 20 -p slurm256 bash kneaddata.sh &
time srun -o metaphlan3.out -e metaphlan3.err -N 2 -c 20 -p slurm256 bash metaphlan3-1.sh &
time srun -o metaphlan3.out -e metaphlan3.err -N 1 -c 20 -p slurm256 bash metaphlan3-2.sh &
time srun -o metaphlan3.out -e metaphlan3.err -N 1 -c 20 -p slurm256 bash metaphlan3-3.sh &
time srun -o humann3.out -e humann3.err -N 1 -c 20 -p slurm128 -w mnclient04 bash humann3.sh &

# 两个样品
python /work/workspace/zhurj/script/2_metapro/metapipe/metapipe_ref.py /work/workspace/zhurj/project/1_metadata/mouse2/input/mnid /work/workspace/zhurj/project/1_metadata/mouse2
拆分分析流程
qc+filter： qcfilter.sh
sed -n '1,18p' metapipe_ref.sh > qcfilter.sh
kneaddata： kneaddata.sh
sed -n '21,31p' metapipe_ref.sh > kneaddata.sh

/work/workspace/zhurj/project/1_metadata/mouse2/p
time srun -o qcfilter.out -e qcfilter.err -N 1 -c 20 -p slurm256 bash qcfilter.sh &
time srun -o kneaddata.out -e kneaddata.err -N 1 -c 20 -p slurm256 bash kneaddata.sh &

# remove host reads with bowtie
http://www.metagenomics.wiki/tools/short-read/remove-host-sequences
cd /work/workspace/zhurj/reference/mouse/GRCm38

parallel --link -j 1 'mkdir -p {2}/{1}' :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22 \
::: /work/workspace/zhurj/project/1_metadata/mouse22/bwarmhost

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
bwa index -a bwtsw /work/workspace/zhurj/reference/mouse/GRCm38/grcm38.fa
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
cd /work/workspace/zhurj/reference/mouse/GRCm38
time srun -o index.out -e index.err -N 1 -c 20 -p slurm256 bash index.sh &


source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link -j 2 'bwa mem {4} {1}/{2}.1.fq.gz {1}/{2}.2.fq.gz -t 20 -o {3}/{2}/{2}.sam && \
samtools view -@ 20 -bS {3}/{2}/{2}.sam -o {3}/{2}/{2}.bam && \
samtools view -b -F 2 -@ 20 {3}/{2}/{2}.bam -o {3}/{2}/{2}_PEunmapped.bam && \
samtools sort -n -@ 20 {3}/{2}/{2}_PEunmapped.bam -o {3}/{2}/{2}_PEunmapped_sorted.bam && \
bedtools bamtofastq -i {3}/{2}/{2}_PEunmapped_sorted.bam -fq {3}/{2}/{2}_r1.fastq -fq2 {3}/{2}/{2}_r2.fastq && \
rm {3}/{2}/*.sam {3}/{2}/*.bam -f ' \
::: /work/workspace/zhurj/project/1_metadata/mouse22/clean :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22 \
::: /work/workspace/zhurj/project/1_metadata/mouse22/bwarmhost ::: /work/workspace/zhurj/reference/mouse/GRCm38/grcm38.fa
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

cd /work/workspace/zhurj/project/1_metadata/mouse22/p
time srun -o bwarmhost.out -e bwarmhost.err -N 1 -c 20 -p slurm256 bash bwarmhost.sh &

bwarmhost_MNA00253.sh
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link -j 1 'bwa mem {4} {1}/{2}.1.fq.gz {1}/{2}.2.fq.gz -t 20 -o {3}/{2}/{2}.sam && \
samtools view -@ 20 -bS {3}/{2}/{2}.sam -o {3}/{2}/{2}.bam && \
samtools view -b -F 2 -@ 20 {3}/{2}/{2}.bam -o {3}/{2}/{2}_PEunmapped.bam && \
samtools sort -n -@ 20 {3}/{2}/{2}_PEunmapped.bam -o {3}/{2}/{2}_PEunmapped_sorted.bam && \
bedtools bamtofastq -i {3}/{2}/{2}_PEunmapped_sorted.bam -fq {3}/{2}/{2}_r1.fastq -fq2 {3}/{2}/{2}_r2.fastq && \
rm {3}/{2}/*.sam {3}/{2}/*.bam -f ' \
::: /work/workspace/zhurj/project/1_metadata/mouse22/clean ::: MNA00253 \
::: /work/workspace/zhurj/project/1_metadata/mouse22/bwarmhost ::: /work/workspace/zhurj/reference/mouse/GRCm38/grcm38.fa
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

cd /work/workspace/zhurj/project/1_metadata/mouse22/p
time srun -o bwarmhost_MNA00253.out -e bwarmhost_MNA00253.err -N 1 -c 20 -p slurm256 -w mnclient01 bash bwarmhost_MNA00253.sh &


# 提取map到GCF_000020225.1(Akkermansia muciniphila)的read

index reference: /work/workspace/zhurj/reference/NCBI/genome/GCA_000020225/bwa/GCA_000020225.fna.gz

parallel --link -j 2 'mkdir -p {2}/{1}' :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22 \
::: /work/workspace/zhurj/project/1_metadata/mouse22/Akk

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link -j 1 'bwa mem {4} {1}/{2}.1.fq.gz {1}/{2}.2.fq.gz -t 40 -o {3}/{2}/{2}.sam && \
samtools view -@ 40 -bS {3}/{2}/{2}.sam -o {3}/{2}/{2}.bam && \
samtools view -b -f 2 -@ 40 {3}/{2}/{2}.bam -o {3}/{2}/{2}_PEmapped.bam && \
samtools sort -n -@ 40 {3}/{2}/{2}_PEmapped.bam -o {3}/{2}/{2}_PEmapped_sorted.bam && \
bedtools bamtofastq -i {3}/{2}/{2}_PEmapped_sorted.bam -fq {3}/{2}/{2}_mapped.r1.fastq -fq2 {3}/{2}/{2}_mapped.r2.fastq && \
rm {3}/{2}/*.sam {3}/{2}/*.bam -f ' \
::: /work/workspace/zhurj/project/1_metadata/mouse22/clean :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22 \
::: /work/workspace/zhurj/project/1_metadata/mouse22/Akk ::: /work/workspace/zhurj/reference/NCBI/genome/GCA_000020225/bwa/GCA_000020225.fna.gz
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

cd /work/workspace/zhurj/project/1_metadata/mouse22/p
time srun -o Akkmap.out -e Akkmap.err -N 1 -c 20 -p slurm256 -w mnclient01 bash Akkmap.sh &


# 提取map到MNH04863的read
cd /work/workspace/zhurj/reference/genome/MNH04863
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
bwa index -a is /work/workspace/zhurj/reference/genome/MNH04863/mno863.fa
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
cd /work/workspace/zhurj/reference/genome/MNH04863
time srun -o index.out -e index.err -N 1 -c 20 -p slurm128 bash index.sh &

parallel --link -j 2 'mkdir -p {2}/{1}' :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22 \
::: /work/workspace/zhurj/project/1_metadata/mouse22/MNH04863
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link -j 1 'bwa mem {4} {1}/{2}.1.fq.gz {1}/{2}.2.fq.gz -t 40 -o {3}/{2}/{2}.sam && \
samtools view -@ 40 -bS {3}/{2}/{2}.sam -o {3}/{2}/{2}.bam && \
samtools view -b -f 2 -@ 40 {3}/{2}/{2}.bam -o {3}/{2}/{2}_PEmapped.bam && \
samtools sort -n -@ 40 {3}/{2}/{2}_PEmapped.bam -o {3}/{2}/{2}_PEmapped_sorted.bam && \
bedtools bamtofastq -i {3}/{2}/{2}_PEmapped_sorted.bam -fq {3}/{2}/{2}_mapped.r1.fastq -fq2 {3}/{2}/{2}_mapped.r2.fastq && \
rm {3}/{2}/*.sam {3}/{2}/*.bam -f ' \
::: /work/workspace/zhurj/project/1_metadata/mouse22/clean :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22 \
::: /work/workspace/zhurj/project/1_metadata/mouse22/MNH04863 ::: /work/workspace/zhurj/reference/genome/MNH04863/mno863.fa
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

cd /work/workspace/zhurj/project/1_metadata/mouse22/p
time srun -o mno863map.out -e mno863map.err -N 2 -c 20 -p slurm128 bash mno863map.sh &

# mapped MNH04863 from swgs MNH04863
 {2}/MNH04863.1.fq.gz {2}/MNH04863.2.fq.gz' \
::: /work/workspace/zhurj/project/2_swgs/MNH04863/kraken2 ::: /work/rawdata/fastq/genome/MNH/MNH048/MNH04863/E1/L1/S1
mkdir -p /work/workspace/zhurj/project/1_metadata/mouse22/MNH04863/MNH04863

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link -j 1 'bwa mem {4} {1}/{2}.1.fq.gz {1}/{2}.2.fq.gz -t 40 -o {3}/{2}/{2}.sam && \
samtools view -@ 40 -bS {3}/{2}/{2}.sam -o {3}/{2}/{2}.bam && \
samtools view -b -f 2 -@ 40 {3}/{2}/{2}.bam -o {3}/{2}/{2}_PEmapped.bam && \
samtools sort -n -@ 40 {3}/{2}/{2}_PEmapped.bam -o {3}/{2}/{2}_PEmapped_sorted.bam && \
bedtools bamtofastq -i {3}/{2}/{2}_PEmapped_sorted.bam -fq {3}/{2}/{2}_mapped.r1.fastq -fq2 {3}/{2}/{2}_mapped.r2.fastq && \
rm {3}/{2}/*.sam {3}/{2}/*.bam -f ' \
::: /work/rawdata/fastq/genome/MNH/MNH048/MNH04863/E1/L1/S1 ::: MNH04863 \
::: /work/workspace/zhurj/project/1_metadata/mouse22/MNH04863 ::: /work/workspace/zhurj/reference/genome/MNH04863/mno863.fa
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
cd /work/workspace/zhurj/project/1_metadata/mouse22/p
time srun -o mno863mapself.out -e mno863mapself.err -N 1 -c 20 -p slurm256 -w mnclient01 bash mno863mapself.sh &

mapped read count
rm /work/workspace/zhurj/project/1_metadata/mouse22/MNH04863/stat/mapMNH04863_sum.txt -f 
parallel -j 1 --link 'rm {1}/stat/mapMNH04863_sum.txt -f && num=$(cat {1}/{2}/{2}_mapped.r1.fastq | grep @ | wc -l ) && echo -e "{2}\t$num" >> {1}/stat/mapMNH04863_sum.txt ' \
::: /work/workspace/zhurj/project/1_metadata/mouse22/MNH04863 :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22

rm /work/workspace/zhurj/project/1_metadata/mouse22/MNH04863/stat/clean_sum.txt -f 
parallel -j 1 --link 'num=$(zcat {1}/{2}.1.fq.gz | grep @ | wc -l ) && echo -e "{2}\t$num" >> {3}/stat/clean_sum.txt ' \
::: /work/workspace/zhurj/project/1_metadata/mouse22/clean :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22 ::: /work/workspace/zhurj/project/1_metadata/mouse22/MNH04863

rm /work/workspace/zhurj/project/1_metadata/mouse22/MNH04863/stat/dehost_sum.txt -f 
parallel -j 1 --link 'num=$(cat {1}/{2}/{2}_r1.fastq | grep @ | wc -l ) && echo -e "{2}\t$num" >> {3}/stat/dehost_sum.txt ' \
::: /work/workspace/zhurj/project/1_metadata/mouse22/bwarmhost :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/dehost ::: /work/workspace/zhurj/project/1_metadata/mouse22/MNH04863  

rm /work/workspace/zhurj/project/1_metadata/mouse22/MNH04863/stat/dehost_249_254.txt -f 
parallel -j 1 --link 'num=$(cat {1}/{2}/{2}_r1.fastq | grep @ | wc -l ) && echo -e "{2}\t$num" >> {3}/stat/dehost_249_254.txt ' \
::: /work/workspace/zhurj/project/1_metadata/mouse22/bwarmhost :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/dehost_249_254 ::: /work/workspace/zhurj/project/1_metadata/mouse22/MNH04863 


python /work/workspace/zhurj/script/2_metapro/metapipe/metapipe_ref_all.py /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22 /work/workspace/zhurj/project/1_metadata/mouse22
/work/workspace/zhurj/project/1_metadata/mouse22/p
sed -n '87,95p' metapipe_ref.sh > metabat.sh
sed -n '98,114p' metapipe_ref.sh > prokka-cdhit.sh
#time srun -o metabat.out -e metabat.err -N 2 -c 20 -p slurm256 --dependency=afterok:257987 bash metabat.sh &
#time srun -o prokka-cdhit.out -e prokka-cdhit.err -N 2 -c 20 -p slurm128 --dependency=afterok:257987 bash prokka-cdhit.sh &
time srun -o metabat.out -e metabat.err -N 2 -c 20 -p slurm256 bash metabat.sh &
time srun -o prokka-cdhit.out -e prokka-cdhit.err -N 1 -c 20 -p slurm128 -w mnclient04 bash prokka-cdhit.sh &

cpanm XML::Simple -l /work/workspace/zhurj/bin/miniconda3/lib/site_perl/5.26.2 --mirror http://mirrors.aliyun.com/CPAN

1. contig assemble
spades.sh
删除中间文件只保留每个样本的 contigs.fasta scaffolds.fasta
parallel --link -j 1 'mkdir -p {1}/spades/{2} ' ::: /work/workspace/zhurj/project/1_metadata/mouse22 :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22
mkdir -p /work/workspace/zhurj/project/1_metadata/mouse22/spades/tmp
```
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --xapply -j 2 ' \
python /work/workspace/moonis/bin/genoassemble/soft/SPAdes-3.13.1-Linux/bin/spades.py --meta -1 {1}/bwarmhost/{2}/{2}_r1.fastq \
-2 {1}/bwarmhost/{2}/{2}_r2.fastq -t 20 -o {1}/spades/{2}  && \
rm {1}/spades/tmp/* -rf && \
mv {1}/spades/{2}/contigs.fasta  {1}/spades/tmp && \
mv {1}/spades/{2}/scaffolds.fasta  {1}/spades/tmp && \
rm {1}/spades/{2}/* -rf && \
mv {1}/spades/tmp/* {1}/spades/{2} ' \
::: /work/workspace/zhurj/project/1_metadata/mouse22 :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

```
cd /work/workspace/zhurj/project/1_metadata/mouse22/p
time srun -o spades.out -e spades.err -N 2 -c 20 -p slurm256 bash spades.sh &

MNA00251, MNA00253 spades 第一次失败，从新run - 2020.10.06

spadesmetabat_MNA00251.sh
'''
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
/work/workspace/moonis/bin/genoassemble/soft/SPAdes-3.13.1-Linux/bin/spades.py --meta \
-1 /work/workspace/zhurj/project/1_metadata/mouse22/bwarmhost/MNA00251/MNA00251_r1.fastq \
-2 /work/workspace/zhurj/project/1_metadata/mouse22/bwarmhost/MNA00251/MNA00251_r2.fastq -t 36 \
-o /work/workspace/zhurj/project/1_metadata/mouse22/spades/MNA00251

parallel -j 1 --link 'bwa index {1}/spades/{2}/contigs.fasta && \
bwa mem -t {3} {1}/spades/{2}/contigs.fasta {1}/bwarmhost/{2}/{2}_r1.fastq {1}/bwarmhost/{2}/{2}_r2.fastq > {1}/bwa/{2}/{2}.sam && \
samtools view -@ {3} -S -b {1}/bwa/{2}/{2}.sam -o {1}/bwa/{2}/{2}.bam && \
samtools sort -@ {3} {1}/bwa/{2}/{2}.bam -o {1}/bwa/{2}/{2}.sorted.bam && \
jgi_summarize_bam_contig_depths --outputDepth {1}/metabat/{2}/depth.txt  {1}/bwa/{2}/{2}.sorted.bam && \
metabat2 -i {1}/spades/{2}/contigs.fasta -a {1}/metabat/{2}/depth.txt -o {1}/metabat/{2}/{2} -m 2000' \
::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: MNA00251 ::: 36
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
'''

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
seqtk seq -L 500 /work/workspace/zhurj/project/1_metadata/mouse22/spades/MNA00251/scaffolds.fasta > /work/workspace/zhurj/project/1_metadata/mouse22/spades/MNA00251/scaffolds_500.fasta
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
prokkaMNA00251.sh
'''
source activate /work/workspace/zhurj/bin/miniconda3/envs/prokka
parallel -j 1 --link 'prokka --outdir {1}/prokkaSlurm/{2} -cpus 20 --prefix {2} --force --metagenome --compliant --locustag {2} {1}/spades/{2}/scaffolds_500.fasta' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: MNA00251
echo `date +"%Y-%m-%d %H:%M:%S"`
echo "Process8: gene prediction finished "
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/prokka

'''
prokka242_253clean.sh
'''
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
seqtk seq -L 500 /work/workspace/zhurj/project/1_metadata/mouse22/spades/MNA00253_clean/scaffolds.fasta > /work/workspace/zhurj/project/1_metadata/mouse22/spades/MNA00253_clean/scaffolds_500.fasta
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

source activate /work/workspace/zhurj/bin/miniconda3/envs/prokka
parallel -j 1 --link 'prokka --outdir {1}/prokkaSlurm/{2} -cpus 20 --prefix {2} --force --metagenome --compliant --locustag {2} {1}/spades/MNA00253_clean/scaffolds_500.fasta' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: MNA00253
parallel -j 1 --link 'prokka --outdir {1}/prokkaSlurm/{2} -cpus 20 --prefix {2} --force --metagenome --compliant --locustag {2} {1}/spades/{2}/scaffolds_500.fasta' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: MNA00242
echo `date +"%Y-%m-%d %H:%M:%S"`
echo "Process8: gene prediction finished "
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/prokka

cd /work/workspace/zhurj/project/1_metadata/mouse22/p
time srun -o prokka242_253clean.out -e prokka242_253clean.err -N 1 -c 20 -p slurm256 -w mnclient02 bash prokka242_253clean.sh &


spadesmetabat_MNA00253_dehost.sh

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
/work/workspace/moonis/bin/genoassemble/soft/SPAdes-3.13.1-Linux/bin/spades.py --meta \
-1 /work/workspace/zhurj/project/1_metadata/mouse22/bwarmhost/MNA00253/MNA00253_r1.fastq \
-2 /work/workspace/zhurj/project/1_metadata/mouse22/bwarmhost/MNA00253/MNA00253_r2.fastq -t 36 \
-o /work/workspace/zhurj/project/1_metadata/mouse22/spades/MNA00253
parallel -j 1 --link 'bwa index {1}/spades/{2}/contigs.fasta && \
bwa mem -t {3} {1}/spades/{2}/contigs.fasta {1}/clean/{2}.1.fq.gz {1}/clean/{2}.1.fq.gz > {1}/bwa/{2}/{2}.sam && \
samtools view -@ {3} -S -b {1}/bwa/{2}/{2}.sam -o {1}/bwa/{2}/{2}.bam && \
samtools sort -@ {3} {1}/bwa/{2}/{2}.bam -o {1}/bwa/{2}/{2}.sorted.bam && \
jgi_summarize_bam_contig_depths --outputDepth {1}/metabat/{2}/depth.txt  {1}/bwa/{2}/{2}.sorted.bam && \
metabat2 -i {1}/spades/{2}/contigs.fasta -a {1}/metabat/{2}/depth.txt -o {1}/metabat/{2}/{2} -m 2000' \
::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: MNA00253  ::: 36
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6


spadesmetabat_MNA00253.sh

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
/work/workspace/moonis/bin/genoassemble/soft/SPAdes-3.13.1-Linux/bin/spades.py --meta \
-1 /work/workspace/zhurj/project/1_metadata/mouse22/clean/MNA00253.1.fq.gz \
-2 /work/workspace/zhurj/project/1_metadata/mouse22/clean/MNA00253.2.fq.gz -t 36 \
-o /work/workspace/zhurj/project/1_metadata/mouse22/spades/MNA00253
parallel -j 1 --link 'bwa index {1}/spades/{2}/contigs.fasta && \
bwa mem -t {3} {1}/spades/{2}/contigs.fasta {1}/clean/{2}.1.fq.gz {1}/clean/{2}.1.fq.gz > {1}/bwa/{2}/{2}.sam && \
samtools view -@ {3} -S -b {1}/bwa/{2}/{2}.sam -o {1}/bwa/{2}/{2}.bam && \
samtools sort -@ {3} {1}/bwa/{2}/{2}.bam -o {1}/bwa/{2}/{2}.sorted.bam && \
jgi_summarize_bam_contig_depths --outputDepth {1}/metabat/{2}/depth.txt  {1}/bwa/{2}/{2}.sorted.bam && \
metabat2 -i {1}/spades/{2}/contigs.fasta -a {1}/metabat/{2}/depth.txt -o {1}/metabat/{2}/{2} -m 2000' \
::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: MNA00253  ::: 36
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

time srun -o spadesmetabat_MNA00253_dehost.out -e spadesmetabat_MNA00253_dehost.err -N 1 -c 20 -p slurm256 -w mnclient02 bash spadesmetabat_MNA00253_dehost.sh &
time srun -o spadesmetabat_MNA00253.out -e spadesmetabat_MNA00253.err -N 1 -c 20 -p slurm256 -w mnclient02 bash spadesmetabat_MNA00253.sh &
time srun -o spadesmetabat_MNA00251.out -e spadesmetabat_MNA00251.err -N 1 -c 20 -p slurm256 -w mnclient02 bash spadesmetabat_MNA00251.sh &
time srun -o prokkaMNA00251.out -e prokkaMNA00251.err -N 1 -c 20 -p slurm256 -w mnclient02 bash prokkaMNA00251.sh &



2. bin
metabat.sh
parallel -j 1 'mkdir -p {1}/{2}/{3} ' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: metabat bwa :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22
paralle -j 1 'find {1}/spades | grep "/contigs.fasta" | awk -F "/" \'{print $(NF-1)}\' > {1}/input/spades_contigs ' ::: /work/workspace/zhurj/project/1_metadata/mouse22
paralle -j 1 'find {1}/spades | grep "/scaffolds.fasta" | awk -F "/" \'{print $(NF-1)}\' > {1}/input/spades_scaffolds ' ::: /work/workspace/zhurj/project/1_metadata/mouse22
```
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel -j 2 --link 'bwa index {1}/spades/{2}/contigs.fasta && \
bwa mem -t {3} {1}/spades/{2}/contigs.fasta {1}/bwarmhost/{2}/{2}_r1.fastq {1}/bwarmhost/{2}/{2}_r2.fastq > {1}/bwa/{2}/{2}.sam && \
samtools view -@ {3} -S -b {1}/bwa/{2}/{2}.sam -o {1}/bwa/{2}/{2}.bam && \
samtools sort -@ {3} {1}/bwa/{2}/{2}.bam -o {1}/bwa/{2}/{2}.sorted.bam && \
jgi_summarize_bam_contig_depths --outputDepth {1}/metabat/{2}/depth.txt  {1}/bwa/{2}/{2}.sorted.bam && \
metabat2 -i {1}/spades/{2}/contigs.fasta -a {1}/metabat/{2}/depth.txt -o {1}/metabat/{2}/{2} -m 2000' \
::: /work/workspace/zhurj/project/1_metadata/mouse22 :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22 ::: 20
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

```

3. gene prediction
Scaftigs（>=500bp）
source activate python3.6
cd /work/workspace/zhurj/project/1_metadata/metapipe/spades/MNC00230
parallel -j 1 'mkdir -p {1}/prokka/{2} ' ::: /work/workspace/zhurj/project/1_metadata/mouse22 :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel -j 2 --link ' \
seqtk seq -L 500 {1}/spades/{2}/scaffolds.fasta > {1}/spades/{2}/scaffolds_500.fasta && \
/work/workspace/zhurj/bin/miniconda3/envs/prokka/bin/prokka --outdir {1}/prokka/{2} -cpus 20 --prefix {2} --force \
--metagenome --compliant --locustag {2} {1}/spades/{2}/scaffolds_500.fasta \
' ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921 ::: MNC00233 
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
--centre X --compliant

cd /work/workspace/zhurj/project/1_metadata/metapipe20200921/p
time srun -o test_prokka.out -e test_prokka.err -N 1 -c 20 -p slurm256 bash test_prokka.sh &


+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MetaGeneMark can not be run on slurm, because this setting "cp /work/workspace/zhurj/software/gm_key_64 ~/.gm_key"
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

7. gene prediction
GeneMark.hmm version 3.38
INSTALL
cp /work/workspace/zhurj/software/gm_key_64 ~/.gm_key
cd /work/workspace/zhurj/software/MetaGeneMark_linux_64/mgm
#gmhmmp -a -A /work/workspace/zhurj/project/1_metadata/metapipe/metagenemark/MNC00230/MNC00230.faa -d -D /work/workspace/zhurj/project/1_metadata/metapipe/metagenemark/MNC00230/MNC00230.fna -m MetaGeneMark_v1.mod /work/workspace/zhurj/project/1_metadata/metapipe/spades/MNC00230/scaffolds_500.fasta  -o /work/workspace/zhurj/project/1_metadata/metapipe/metagenemark/MNC00230/MNC00230.lst

/work/workspace/zhurj/software/MetaGeneMark_linux_64/mgm/gmhmmp -a -A /work/workspace/zhurj/project/1_metadata/metapipe/metagenemark/MNC00230/MNC00230.faa -d -D /work/workspace/zhurj/project/1_metadata/metapipe/metagenemark/MNC00230/MNC00230.fna -m MetaGeneMark_v1.mod /work/workspace/zhurj/project/1_metadata/metapipe/spades/MNC00230/scaffolds_500.fasta  -o /work/workspace/zhurj/project/1_metadata/metapipe/metagenemark/MNC00230/MNC00230.lst

cd /work/workspace/zhurj/project/1_metadata/metapipe/metagenemark/MNC00230/
source activate python3.6
# filtering seq with length < 100
seqtk seq -L 100 MNC00230.fna > MNC00230_100.fna
# change long gene title to >gene_1, gene_2, ...
awk -F "|" '{if( $1~/^>/) print $1;else print $0}' MNC00230_100.fna > MNC00230_100_adj.fna
seqtk seq -L 100 MNC00230.faa > MNC00230_100.faa
awk -F "|" '{if( $1~/^>/) print $1;else print $0}' MNC00230_100.faa > MNC00230_100_adj.faa


8. 基因聚类
CD-HIT version 4.8.1 
ref: https://www.jianshu.com/p/d4e4009f0dc0
ref: http://weizhongli-lab.org/cd-hit/
cd /work/workspace/zhurj/project/1_metadata/metapipe/metagenemark/MNC00230
cd-hit-est -i /work/workspace/zhurj/project/1_metadata/metapipe/metagenemark/MNC00230/MNC00230_100_adj.fna -o /work/workspace/zhurj/project/1_metadata/metapipe/cdhit/MNC00230/cdhit.fna -c 0.95 -G 0 -aS 0.9 -g 1 -d 0 -T 36 
cd-hit -i /work/workspace/zhurj/project/1_metadata/metapipe/metagenemark/MNC00230/MNC00230_100_adj.faa -o /work/workspace/zhurj/project/1_metadata/metapipe/cdhit/MNC00230/cdhit.faa -c 0.95 -G 0 -aS 0.9 -g 1 -d 0 -T 36 
protein size
organism	median protein length amino acides
H.sapien	375
D.melanogaster	373
C.elegans	344
S.cerevisiae	379
A.thaliana	356
5 eukaryotes(above)	361
67 bacterial	267
15 archaeal	247
===========================================================
单个粪便样品，宏基因组测序，预测基因数目范围？
CD-HIT 去冗余后基因数目范围
===========================================================

test
cd-hit-est -i /work/workspace/zhurj/project/1_metadata/metapipe20200921/prokka/MNC00233/MNC00233.ffn -o /work/workspace/zhurj/project/1_metadata/metapipe20200921/cdhit/cdhit.fna -c 0.95 -G 0 -aS 0.9 -g 1 -d 0 -T 36


time srun -o prokkaSlurm.out -e prokkaSlurm.err -N 2 -c 20 -p slurm256 bash prokkaSlurm.sh &

cd /work/workspace/zhurj/project/1_metadata/mouse22/diff/data
sed -i 's/_unclassified/_uncl./g' bar_genus_in 
sed -i 's/_unclassified/_uncl./g' bar_species_in
cd /work/workspace/zhurj/project/1_metadata/mouse22/diff
cd /work/workspace/zhurj/project/1_metadata/mouse22/diff && python /work/workspace/zhurj/lib/python/script/lefse2boxplot.py data/bar_genus_in  input/meta graph_in/genus_format_in
cd /work/workspace/zhurj/project/1_metadata/mouse22/diff && python /work/workspace/zhurj/lib/python/script/lefse2boxplot.py data/bar_phylum_in  input/meta graph_in/phylum_format_in
cd /work/workspace/zhurj/project/1_metadata/mouse22/diff && python /work/workspace/zhurj/lib/python/script/lefse2boxplot.py data/bar_species_in  input/meta graph_in/species_format_in

axid   kindom     phylum          class         order           family           genus         species
816     Bacteria   Bacteroidetes   Bacteroidia   Bacteroidales   Bacteroidaceae   Bacteroides
Species Sample  Abundance   Group
p__Firmicutes   X1  0.38509871  NIG
p__Firmicutes   X2  0.462668997 NIG
p__Firmicutes   X3  0.713981572 NIG
p__Firmicutes   X4  0.708546589 NIG

library('ggplot2')
library('RColorBrewer')
library('ggsci')
D <- read.table(file="/work/workspace/zhurj/project/1_metadata/mouse22/diff/graph_in/phylum_format_in", head=T, sep='\t')
png("/work/workspace/zhurj/project/1_metadata/mouse22/diff/png/phylum_20201010.png", width=15, height=8, units = "in", res=100);
a<-D[grep('p__Firmicutes',D$Species),]
b<-a[order(a$Abundance,decreasing = T),]
#D$Sample<-factor(D$SDample,levels = b$Sample)
#D$Sample<-as.character(D$Sample)
D$Sample<-as.name(D$Sample)
D$Group<-factor(D$Group,levels = c("NCD","HFD","MNO863"))
mypal = c("#E41A1C","#4DAF4A","#377EB8","#A65628","#FF7F00","#FFFF33")
mycolor = colorRampPalette(mypal)
ggplot(D) +
  geom_bar(aes(x=Sample, y=Abundance, group=Sample, fill=factor(Species)),stat="identity", width=1, color= "black") +
  #scale_x_discrete(breaks=NULL,expand=c(0.01,0)) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_continuous(expand = c(0.00,0.00)) +
  facet_grid(~Group, scales = "free", space = "free_x", switch = "x") +
  theme_classic() +
  theme(
      axis.title=element_text( family="Helvetica", face="bold", colour="black",size=16),
      axis.text=element_text(family="Helvetica",colour="black",size=12),
      axis.text.x=element_text(angle=90,hjust=1, vjust=0.5,size=10),
      strip.text=element_text(family="Helvetica",colour="black",angle=0,hjust=0.5,vjust=0.5,size=14),
      #strip.text.x=element_blank(),
      #strip.background=element_blank(),
      #legend.position=c(0.90,0.80),
      legend.title=element_text(family="Helvetica",colour="black",size=16),
      legend.text=element_text(family="Helvetica", colour="black",size=14)
  ) +
  guides(fill = guide_legend("Phylum", ncol=1)) +
  scale_fill_manual(values = mycolor(9), breaks = rev(as.vector(unique(D$Species)))) +
  xlab('') +
  ylab('Abundance')
  dev.off()

library('ggplot2')
library('RColorBrewer')
library('ggsci')
D <- read.table(file="/work/workspace/zhurj/project/1_metadata/mouse22/diff/graph_in/genus_format_in", head=T, sep='\t')
png("/work/workspace/zhurj/project/1_metadata/mouse22/diff/png/genus_20201010.png", width=15, height=8, units = "in", res=100);
a<-D[grep('g__Lachnospiraceae_uncl.',D$Species),]
b<-a[order(a$Abundance,decreasing = T),]
#D$Sample<-factor(D$SDample,levels = b$Sample)
D$Sample<-as.character(D$Sample)
D$Group<-factor(D$Group,levels = c("NCD","HFD","MNO863"))
mypal = c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black",
"gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F",
"gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1",
"steelblue4", "darkturquoise", "green1", "yellow4", "yellow3",
"darkorange4", "brown")
mycolor = colorRampPalette(mypal)
ggplot(D) +
  geom_bar(aes(x=Sample, y=Abundance, group=Sample, fill=factor(Species)),stat="identity", width=1, color= "black") +
  #scale_x_discrete(breaks=NULL,expand=c(0.01,0)) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_continuous(expand = c(0.00,0.00)) +
  facet_grid(~Group, scales = "free", space = "free_x", switch = "x") +
  theme_classic() +
  theme(
      axis.title=element_text( family="Helvetica", face="bold", colour="black",size=16),
      axis.text=element_text(family="Helvetica",colour="black",size=12),
      axis.text.x=element_text(angle=90,hjust=1, vjust=0.5,size=10),
      strip.text=element_text(family="Helvetica",colour="black",angle=0,hjust=0.5,vjust=0.5,size=14),
      #strip.text.x=element_blank(),
      #strip.background=element_blank(),
      #legend.position=c(0.90,0.80),
      legend.title=element_text(family="Helvetica",colour="black",size=16),
      legend.text=element_text(family="Helvetica", colour="black",size=14)
  ) +
  guides(fill = guide_legend("Phylum", ncol=1)) +
  scale_fill_manual(values = mycolor(21), breaks = rev(as.vector(unique(D$Species)))) +
  xlab('') +
  ylab('Abundance')
  dev.off()

library('ggplot2')
library('RColorBrewer')
library('ggsci')
D <- read.table(file="/work/workspace/zhurj/project/1_metadata/mouse22/diff/graph_in/species_format_in", head=T, sep='\t')
png("/work/workspace/zhurj/project/1_metadata/mouse22/diff/png/species_20201010.png", width=15, height=8, units = "in", res=100);
a<-D[grep('s__Lachnospiraceae_bacterium_28_4',D$Species),]
b<-a[order(a$Abundance,decreasing = T),]
#D$Sample<-factor(D$SDample,levels = b$Sample)
D$Sample<-as.character(D$Sample)
D$Group<-factor(D$Group,levels = c("NCD","HFD","MNO863"))
mypal = c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black",
"gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F",
"gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1",
"steelblue4", "darkturquoise", "green1", "yellow4", "yellow3",
"darkorange4", "brown")
mycolor = colorRampPalette(mypal)
ggplot(D) +
  geom_bar(aes(x=Sample, y=Abundance, group=Sample, fill=factor(Species)),stat="identity", width=1, color= "black") +
  #scale_x_discrete(breaks=NULL,expand=c(0.01,0)) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_continuous(expand = c(0.00,0.00)) +
  facet_grid(~Group, scales = "free", space = "free_x", switch = "x") +
  theme_classic() +
  theme(
      axis.title=element_text( family="Helvetica", face="bold", colour="black",size=16),
      axis.text=element_text(family="Helvetica",colour="black",size=12),
      axis.text.x=element_text(angle=90,hjust=1, vjust=0.5,size=10),
      strip.text=element_text(family="Helvetica",colour="black",angle=0,hjust=0.5,vjust=0.5,size=14),
      #strip.text.x=element_blank(),
      #strip.background=element_blank(),
      #legend.position=c(0.90,0.80),
      legend.title=element_text(family="Helvetica",colour="black",size=16),
      legend.text=element_text(family="Helvetica", colour="black",size=14)
  ) +
  guides(fill = guide_legend("Phylum", ncol=1)) +
  scale_fill_manual(values = mycolor(21), breaks = rev(as.vector(unique(D$Species)))) +
  xlab('') +
  ylab('Abundance')
  dev.off()

diff
cd /work/workspace/zhurj/project/1_metadata/mouse22/diff/data && sed -i 's/_unclassified/_uncl./g' dif_genus_in 
cd /work/workspace/zhurj/project/1_metadata/mouse22/diff/data && sed -i 's/_unclassified/_uncl./g' dif_species_in
cd /work/workspace/zhurj/project/1_metadata/mouse22/diff && python /work/workspace/zhurj/lib/python/script/lefse2boxplot.py data/dif_genus_in  input/meta graph_in/dif_genus_format_in
cd /work/workspace/zhurj/project/1_metadata/mouse22/diff && python /work/workspace/zhurj/lib/python/script/lefse2boxplot.py data/dif_phylum_in  input/meta graph_in/dif_phylum_format_in
cd /work/workspace/zhurj/project/1_metadata/mouse22/diff && python /work/workspace/zhurj/lib/python/script/lefse2boxplot.py data/dif_species_in  input/meta graph_in/dif_species_format_in

library(ggplot2)
library('RColorBrewer')
library('ggsci')
library('ggpubr')

#library('ggbiplot')
#library(reshape2)
#library(dplyr)
# stat_compare_means(comparisons=my_comparisions, method = "wilcox.test", aes(label=..p.adj..)) +
# label = "p.signif", aes(label=..p.adj..)
# method = "t.test", method = "wilcox.test"

Data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/diff/graph_in/dif_phylum_format_in",header=T)
my_comparisions <- list(c("NCD", "HFD"),c("NCD", "MNO863"),c("HFD", "MNO863"))
p <- ggboxplot(Data, x = "Group", y = "Abundance", color = "Group", palette = "npg",  facet.by = "Species", ylim = c(0, 150)) + 
ggtitle("Boxplot of phylum")+
theme(legend.position="none",plot.title = element_text(face="bold.italic", hjust=0.5),)+
scale_y_continuous(name="Relative Abundance", breaks=seq(0,150,20)) +
scale_x_discrete(name="") +
stat_compare_means(comparisons=my_comparisions, method = "wilcox.test", label = "p.signif") +
guides(fill=FALSE)
ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/diff/png/dif_phylum_20201010.jpg", p, width = 8, height =6)


Data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/diff/graph_in/dif_genus_format_in",header=T)
my_comparisions <- list(c("NCD", "HFD"),c("NCD", "MNO863"),c("HFD", "MNO863"))
p <- ggboxplot(Data, x = "Group", y = "Abundance", color = "Group", palette = "npg",  facet.by = "Species", ylim = c(0, 150)) + 
ggtitle("Boxplot of genus")+
theme(legend.position="none",plot.title = element_text(face="bold.italic", hjust=0.5),)+
scale_y_continuous(name="Relative Abundance", breaks=seq(0,150,20)) +
scale_x_discrete(name="") +
stat_compare_means(comparisons=my_comparisions, method = "wilcox.test", label = "p.signif") +
guides(fill=FALSE)
ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/diff/png/dif_genus_20201010.jpg", p, width = 10, height =10)

Data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/diff/graph_in/dif_species_format_in",header=T)
my_comparisions <- list(c("NCD", "HFD"),c("NCD", "MNO863"),c("HFD", "MNO863"))
p <- ggboxplot(Data, x = "Group", y = "Abundance", color = "Group", palette = "npg",  facet.by = "Species", ylim = c(0, 150)) + 
ggtitle("Boxplot of species")+
theme(legend.position="none",plot.title = element_text(face="bold.italic", hjust=0.5),)+
scale_y_continuous(name="Relative Abundance", breaks=seq(0,150,20)) +
scale_x_discrete(name="") +
stat_compare_means(comparisons=my_comparisions, method = "wilcox.test", label = "p.signif") +
guides(fill=FALSE)
ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/diff/png/dif_species_20201010.jpg", p, width = 15, height =12)

# http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
library(corrplot)
dev.off()
pdf(file = "/work/workspace/zhurj/project/1_metadata/mouse22/diff/png/corr_item36.pdf",width = 400, height = 350)
Data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/diff/graph_in/corr_col36_in",header=T,row.names=1)
res <- cor(Data)
corrplot(res, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
dev.off()


corrplot(correlations_history, method = "number", type = "lower", 
title = "Regional Factor Correlation Matrix over history", 
mar = c(0,0,1,0), number.cex = 0.5, number.digits = 2)

library("PerformanceAnalytics")


pdf(file = "/work/workspace/zhurj/project/1_metadata/mouse22/diff/png/corr_pearson_item36.pdf",width = 100, height = 100)
Data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/diff/graph_in/corr_col36_simple",header=T,row.names=1)
chart.Correlation(Data, histogram=TRUE, method="pearson", pch=19)
dev.off()

pdf(file = "/work/workspace/zhurj/project/1_metadata/mouse22/diff/png/corr_spearson_item36.pdf",width = 100, height = 100)
Data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/diff/graph_in/corr_col36_simple",header=T,row.names=1)
chart.Correlation(Data, histogram=TRUE, method="spearman", pch=19)
dev.off()

#  hist.col="#00FA9A", 
library("psych")
pdf(file = "/work/workspace/zhurj/project/1_metadata/mouse22/diff/png/corr_item36.pdf",width = 100, height = 100)
Data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/diff/graph_in/corr_col36_simple",header=T,row.names=1)
pairs.panels(Data, 
             hist.col="red", 
             show.points=TRUE, 
             stars=TRUE, 
             gap=0.05, 
             pch=".", 
             ellipses=FALSE, 
             scale=TRUE,
             jiggle=TRUE,
             factor=2,
             main="Correlation", 
             col="#ADFF2F", 
             pty="m", 
             font=2)
dev.off()

## rarefaction curve
library(vegan)
pdf(file = "/work/workspace/zhurj/project/1_metadata/mouse22/diversity/pic/rarefaction_20201011.pdf")
df <- read.table("/work/workspace/zhurj/project/1_metadata/mouse22/diversity/input/rarefaction_in",sep='\t',header=TRUE,row.names=1)
S <- specnumber(df)
(raremax <- min(rowSums(df)))
Srare <- rarefy(df,raremax)
plot(S,Srare, xlab="Observed No. of Species", ylab="Rarefied No. of Species")
abline(0,1)
rarecurve(df,step = 200, sample = raremax, col = "blue", cex = 0.6, xlim=c(0, 20000))
dev.off()
rarecurve(df,step = 200, sample = raremax, col = "blue", cex = 0.6)

## humann gene filter
source activate humann
parallel --link -j 1 'humann_join_tables -s --input {1} --file_name _{2} --output {1}/merge/{2}_merge.txt'   ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3 ::: uniref90_ko uniref90_go uniref90_level4ec uniref90_pfam uniref90_eggnog metacyc_pathabundance kegg_pathabundance metacyc_genefamilies
source deactivate python3.6

python /work/workspace/zhurj/lib/python/script/humann_gene_filter.py /work/workspace/zhurj/project/1_metadata/mouse22/humann3/merge/uniref90_genefamilies.tsv /work/workspace/zhurj/project/1_metadata/mouse22/humann3/merge_filter/filter_genefamilies.tsv

cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3 && python /work/workspace/zhurj/lib/python/script/humann_gene_filter.py merge/uniref90_genefamilies.tsv merge_filter/filter_genefamilies.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3 && python /work/workspace/zhurj/lib/python/script/humann_others_filter.py merge/kegg_pathabundance_merge.txt merge_filter/filter_kegg.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3 && python /work/workspace/zhurj/lib/python/script/humann_others_filter.py merge/metacyc_pathabundance_merge.txt merge_filter/filter_metacyc.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3 && python /work/workspace/zhurj/lib/python/script/humann_others_filter.py merge/uniref90_eggnog_merge.txt merge_filter/filter_eggnog.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3 && python /work/workspace/zhurj/lib/python/script/humann_others_filter.py merge/uniref90_go_merge.txt merge_filter/filter_go.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3 && python /work/workspace/zhurj/lib/python/script/humann_others_filter.py merge/uniref90_ko_merge.txt merge_filter/filter_ko.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3 && python /work/workspace/zhurj/lib/python/script/humann_others_filter.py merge/uniref90_level4ec_merge.txt merge_filter/filter_level4ec.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3 && python /work/workspace/zhurj/lib/python/script/humann_others_filter.py merge/uniref90_pfam_merge.txt merge_filter/filter_pfam.tsv

source activate humann
# relab: relative abundance [relab]
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3 && humann_renorm_table --units relab --input merge_filter/filter_genefamilies.tsv --output filter_relab/relab_genefamilies.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3 && humann_renorm_table --units relab --input merge_filter/filter_kegg.tsv --output filter_relab/relab_kegg.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3 && humann_renorm_table --units relab --input merge_filter/filter_metacyc.tsv --output filter_relab/relab_metacyc.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3 && humann_renorm_table --units relab --input merge_filter/filter_eggnog.tsv --output filter_relab/relab_eggnog.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3 && humann_renorm_table --units relab --input merge_filter/filter_go.tsv --output filter_relab/relab_go.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3 && humann_renorm_table --units relab --input merge_filter/filter_ko.tsv --output filter_relab/relab_ko.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3 && humann_renorm_table --units relab --input merge_filter/filter_level4ec.tsv --output filter_relab/relab_level4ec.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3 && humann_renorm_table --units relab --input merge_filter/filter_pfam.tsv --output filter_relab/relab_pfam.tsv

# copies per million [cpm]
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3 && humann_renorm_table --units cpm --input merge_filter/filter_genefamilies.tsv --output filter_cpm/cpm_genefamilies.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3 && humann_renorm_table --units cpm --input merge_filter/filter_kegg.tsv --output filter_cpm/cpm_kegg.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3 && humann_renorm_table --units cpm --input merge_filter/filter_metacyc.tsv --output filter_cpm/cpm_metacyc.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3 && humann_renorm_table --units cpm --input merge_filter/filter_eggnog.tsv --output filter_cpm/cpm_eggnog.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3 && humann_renorm_table --units cpm --input merge_filter/filter_go.tsv --output filter_cpm/cpm_go.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3 && humann_renorm_table --units cpm --input merge_filter/filter_ko.tsv --output filter_cpm/cpm_ko.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3 && humann_renorm_table --units cpm --input merge_filter/filter_level4ec.tsv --output filter_cpm/cpm_level4ec.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3 && humann_renorm_table --units cpm --input merge_filter/filter_pfam.tsv --output filter_cpm/cpm_pfam.tsv
source deactivate humann

2020-10-29
# humann 差异ko及
source activate python3.6
python /work/workspace/zhurj/lib/python/script/humann_gene_filter.py /work/workspace/zhurj/project/1_metadata/mouse22/humann3/merge/genefamilies_merge.txt /work/workspace/zhurj/project/1_metadata/mouse22/humann3/merge_filter/filter_genefamilies.tsv
parallel --link 'python /work/workspace/zhurj/lib/python/script/humann_others_filter.py {1}/merge/{2}_merge.txt {1}/merge_filter/filter_{2}.tsv ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3 ::: uniref50_ko uniref50_go uniref50_eggnog uniref50_level4ec uniref50_pfam
source deactivate python3.6

source activate humann
parallel --link 'humann_renorm_table --units cpm --input {1}/merge_filter/filter_{2}.tsv --output {1}/merge_cpm/cpm_{2}.tsv' \
::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3 ::: genefamilies uniref50_eggnog uniref50_go uniref50_ko uniref50_level4ec uniref50_pfam
source deactivate humann


parallel --link 'mkdir -p {1}/humann3/{2} ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/lefse ::: genefamilies kegg metacyc eggnog go ko level4ec pfam

cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3/merge_filter && sed -i 's/_metacyc_Abundance-RPKs//g' filter_eggnog.tsv && sed -i 's/# //g' filter_eggnog.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3/merge_filter && sed -i 's/_metacyc_Abundance-RPKs//g' filter_genefamilies.tsv && sed -i 's/# //g' filter_genefamilies.tsv

cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3/merge_filter && sed -i 's/_metacyc_Abundance-RPKs//g' filter_go.tsv && sed -i 's/# //g' filter_go.tsv 
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3/merge_filter && sed -i 's/_kegg_Abundance//g' filter_kegg.tsv && sed -i 's/# //g' filter_kegg.tsv 
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3/merge_filter && sed -i 's/_metacyc_Abundance-RPKs//g' filter_ko.tsv && sed -i 's/# //g' filter_ko.tsv 
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3/merge_filter && sed -i 's/_metacyc_Abundance-RPKs//g' filter_level4ec.tsv && sed -i 's/# //g' filter_level4ec.tsv 
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3/merge_filter && sed -i 's/_metacyc_Abundance//g' filter_metacyc.tsv && sed -i 's/# //g' filter_metacyc.tsv 
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3/merge_filter && sed -i 's/_metacyc_Abundance-RPKs//g' filter_pfam.tsv && sed -i 's/# //g' filter_pfam.tsv 

cd /work/workspace/zhurj/project/1_metadata/mouse22/compare/GD/humann3
source activate humann
parallel --link -j 1 'humann_join_tables -s --input {1} --file_name _{2} --output {1}/merge/{2}_merge.txt'   ::: /work/workspace/zhurj/project/1_metadata/mouse22/compare/GD/humann3 ::: uniref90_ko uniref90_go uniref90_level4ec uniref90_pfam uniref90_eggnog metacyc_pathabundance kegg_pathabundance genefamilies
parallel --link -j 1 'humann_join_tables -s --input {1} --file_name _{2} --output {1}/merge/{2}_merge.txt'   ::: /work/workspace/zhurj/project/1_metadata/mouse22/compare/GD/humann3 ::: kegg_pathcoverage metacyc_pathcoverage 
source deactivate humann
source activate python3.6
python /work/workspace/zhurj/lib/python/script/humann_gene_filter.py merge/genefamilies_merge.txt merge_filter/filter_genefamilies.tsv
python /work/workspace/zhurj/lib/python/script/humann_others_filter.py merge/kegg_pathabundance_merge.txt merge_filter/filter_kegg.tsv
python /work/workspace/zhurj/lib/python/script/humann_others_filter.py merge/metacyc_pathabundance_merge.txt merge_filter/filter_metacyc.tsv
python /work/workspace/zhurj/lib/python/script/humann_others_filter.py merge/uniref90_eggnog_merge.txt merge_filter/filter_eggnog.tsv
python /work/workspace/zhurj/lib/python/script/humann_others_filter.py merge/uniref90_go_merge.txt merge_filter/filter_go.tsv
python /work/workspace/zhurj/lib/python/script/humann_others_filter.py merge/uniref90_ko_merge.txt merge_filter/filter_ko.tsv
python /work/workspace/zhurj/lib/python/script/humann_others_filter.py merge/uniref90_level4ec_merge.txt merge_filter/filter_level4ec.tsv
python /work/workspace/zhurj/lib/python/script/humann_others_filter.py merge/uniref90_pfam_merge.txt merge_filter/filter_pfam.tsv

cd /work/workspace/zhurj/project/1_metadata/mouse22/compare/GD/humann3/merge_filter
sed -i 's/_metacyc_Abundance-RPKs//g' filter_eggnog.tsv && sed -i 's/# //g' filter_eggnog.tsv
sed -i 's/_metacyc_Abundance-RPKs//g' filter_genefamilies.tsv && sed -i 's/# //g' filter_genefamilies.tsv
sed -i 's/_metacyc_Abundance-RPKs//g' filter_go.tsv && sed -i 's/# //g' filter_go.tsv 
sed -i 's/_kegg_Abundance//g' filter_kegg.tsv && sed -i 's/# //g' filter_kegg.tsv 
sed -i 's/_metacyc_Abundance-RPKs//g' filter_ko.tsv && sed -i 's/# //g' filter_ko.tsv 
sed -i 's/_metacyc_Abundance-RPKs//g' filter_level4ec.tsv && sed -i 's/# //g' filter_level4ec.tsv 
sed -i 's/_metacyc_Abundance//g' filter_metacyc.tsv && sed -i 's/# //g' filter_metacyc.tsv 
sed -i 's/_metacyc_Abundance-RPKs//g' filter_pfam.tsv && sed -i 's/# //g' filter_pfam.tsv 
source deactivate python3.6
parallel --link -j 1 'humann_join_tables -s --input {1} --file_name _{2} --output {1}/merge/{2}_merge.txt'  ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref90 ::: metacyc_pathcoverage kegg_pathcoverage

# remove unmapped genes
2020-10-28
uniref50
humann_uniref50.sh
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel --link -j 2 'humann --threads 20 --input {2}/kneaddata/{1}/{1}.fastq --output {2}/humann3/{1} --output-basename {1}_metacyc --protein-database {3} && \
' ::: MNA00239 ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: /work/workspace/zhurj/database/humann/201901/uniref50
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann
cd /work/workspace/zhurj/project/1_metadata/p/p
time srun -o humann_uniref50.out -e humann_uniref50.err -N 1 -c 20 -p slurm256 -w mnclient01 bash humann_uniref50.sh &

parallel 'mkdir -p {1}/{2} ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3 ::: MNA00239 MNA00240 MNA00241 MNA00242 MNA00243 MNA00244 MNA00245 MNA00246 MNA00247 MNA00248 MNA00249 MNA00250 MNA00251 MNA00252 MNA00253 MNA00254 merge merge_cpm merge_filter merge_reb

humann_uniref50_1.sh
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel --link -j 1 'humann --threads 20 --input {2}/kneaddata/{1}/{1}.fastq --output {2}/humann3/{1} --output-basename {1}_metacyc --protein-database {3} \
' ::: MNA00239 MNA00240 MNA00241 MNA00242 ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: /work/workspace/zhurj/database/humann/201901/uniref50
parallel --link -j 1 'rm {1}/{2}/{2}_metacyc_humann_temp -rf '  \
::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3 ::: MNA00239 MNA00240 NNA00241 MNA00242
parallel -j 1 --link ' humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_ko --output  {1}/{2}/{2}_uniref50_ko.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_go --output  {1}/{2}/{2}_uniref50_go.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_level4ec --output  {1}/{2}/{2}_uniref50_level4ec.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_pfam --output  {1}/{2}/{2}_uniref50_pfam.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_eggnog --output  {1}/{2}/{2}_uniref50_eggnog.tsv \
 ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3 ::: MNA00239 MNA00240 NNA00241 MNA00242
 parallel -j 1 --link 'humann -i {1}/{2}/{2}_uniref50_ko.tsv --pathways-database {3} --output {1}/{2} --threads 20 --output-basename {2}_kegg' \
 ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3 ::: MNA00239 MNA00240 NNA00241 MNA00242 ::: /work/workspace/zhurj/database/humann/current/kegg/keggc
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

humann_uniref50_2.sh
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel --link -j 1 'humann --threads 20 --input {2}/kneaddata/{1}/{1}.fastq --output {2}/humann3/{1} --output-basename {1}_metacyc --protein-database {3} \
' ::: MNA00243 MNA00244 MNA00245 MNA00246 ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: /work/workspace/zhurj/database/humann/201901/uniref50
parallel --link -j 1 'rm {1}/{2}/{2}_metacyc_humann_temp -rf '  \
::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3 ::: MNA00243 MNA00244 MNA00245 MNA00246
parallel -j 1 --link ' humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_ko --output  {1}/{2}/{2}_uniref50_ko.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_go --output  {1}/{2}/{2}_uniref50_go.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_level4ec --output  {1}/{2}/{2}_uniref50_level4ec.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_pfam --output  {1}/{2}/{2}_uniref50_pfam.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_eggnog --output  {1}/{2}/{2}_uniref50_eggnog.tsv \
 ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3 ::: MNA00243 MNA00244 MNA00245 MNA00246
 parallel -j 1 --link 'humann -i {1}/{2}/{2}_uniref50_ko.tsv --pathways-database {3} --output {1}/{2} --threads 20 --output-basename {2}_kegg' \
 ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3 ::: MNA00243 MNA00244 MNA00245 MNA00246 ::: /work/workspace/zhurj/database/humann/current/kegg/keggc
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

humann_uniref50_3.sh
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel --link -j 1 'humann --threads 20 --input {2}/kneaddata/{1}/{1}.fastq --output {2}/humann3/{1} --output-basename {1}_metacyc --protein-database {3} \
' ::: MNA00247 MNA00248 MNA00249 MNA00250 ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: /work/workspace/zhurj/database/humann/201901/uniref50
parallel --link -j 1 'rm {1}/{2}/{2}_metacyc_humann_temp -rf '  \
::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3 ::: MNA00247 MNA00248 MNA00249 MNA00250
parallel -j 1 --link ' humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_ko --output  {1}/{2}/{2}_uniref50_ko.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_go --output  {1}/{2}/{2}_uniref50_go.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_level4ec --output  {1}/{2}/{2}_uniref50_level4ec.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_pfam --output  {1}/{2}/{2}_uniref50_pfam.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_eggnog --output  {1}/{2}/{2}_uniref50_eggnog.tsv \
 ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3 ::: MNA00247 MNA00248 MNA00249 MNA00250
 parallel -j 1 --link 'humann -i {1}/{2}/{2}_uniref50_ko.tsv --pathways-database {3} --output {1}/{2} --threads 20 --output-basename {2}_kegg' \
 ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3 ::: MNA00247 MNA00248 MNA00249 MNA00250 ::: /work/workspace/zhurj/database/humann/current/kegg/keggc
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

humann_uniref50_4.sh
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel --link -j 1 'humann --threads 20 --input {2}/kneaddata/{1}/{1}.fastq --output {2}/humann3/{1} --output-basename {1}_metacyc --protein-database {3} \
' ::: MNA00251 MNA00252 MNA00253 MNA00254 ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: /work/workspace/zhurj/database/humann/201901/uniref50
parallel --link -j 1 'rm {1}/{2}/{2}_metacyc_humann_temp -rf '  \
::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3 ::: MNA00251 MNA00252 MNA00253 MNA00254
parallel -j 1 --link ' humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_ko --output  {1}/{2}/{2}_uniref50_ko.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_go --output  {1}/{2}/{2}_uniref50_go.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_level4ec --output  {1}/{2}/{2}_uniref50_level4ec.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_pfam --output  {1}/{2}/{2}_uniref50_pfam.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_eggnog --output  {1}/{2}/{2}_uniref50_eggnog.tsv \
 ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3 ::: MNA00251 MNA00252 MNA00253 MNA00254
 parallel -j 1 --link 'humann -i {1}/{2}/{2}_uniref50_ko.tsv --pathways-database {3} --output {1}/{2} --threads 20 --output-basename {2}_kegg' \
 ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3 ::: MNA00251 MNA00252 MNA00253 MNA00254 ::: /work/workspace/zhurj/database/humann/current/kegg/keggc
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

humann_MNA00241.sh
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel --link -j 1 'humann --threads 20 --input {2}/kneaddata/{1}/{1}.fastq --output {2}/humann3/{1} --output-basename {1}_metacyc --protein-database {3} \
' ::: MNA00241 ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: /work/workspace/zhurj/database/humann/201901/uniref50
parallel --link -j 1 'rm {1}/{2}/{2}_metacyc_humann_temp -rf '  \
::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3 ::: MNA00241
parallel -j 1 --link ' humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_ko --output  {1}/{2}/{2}_uniref50_ko.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_go --output  {1}/{2}/{2}_uniref50_go.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_level4ec --output  {1}/{2}/{2}_uniref50_level4ec.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_pfam --output  {1}/{2}/{2}_uniref50_pfam.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_eggnog --output  {1}/{2}/{2}_uniref50_eggnog.tsv \
 ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3 ::: MNA00241
 parallel -j 1 --link 'humann -i {1}/{2}/{2}_uniref50_ko.tsv --pathways-database {3} --output {1}/{2} --threads 20 --output-basename {2}_kegg' \
 ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3 ::: MNA00241 ::: /work/workspace/zhurj/database/humann/current/kegg/keggc
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

cd /work/workspace/zhurj/project/1_metadata/p/p
time srun -o humann_uniref50_1.out -e humann_uniref50_1.err -N 1 -c 20 -p slurm128 -w mnclient03 bash humann_uniref50_1.sh &
time srun -o humann_uniref50_2.out -e humann_uniref50_2.err -N 1 -c 20 -p slurm128 -w mnclient04 bash humann_uniref50_2.sh &
time srun -o humann_uniref50_3.out -e humann_uniref50_3.err -N 1 -c 20 -p slurm256 -w mnclient01 bash humann_uniref50_3.sh &
time srun -o humann_uniref50_4.out -e humann_uniref50_4.err -N 1 -c 20 -p slurm256 -w mnclient02 bash humann_uniref50_4.sh &
time srun -o humann_MNA00241.out -e humann_MNA00241.err -N 1 -c 20 -p slurm256 -w mnclient02 bash humann_MNA00241.sh &

parallel --link -j 1 'humann_join_tables -s --input {1} --file_name _{2} --output {1}/merge/{2}_merge.txt'   ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3 ::: uniref50_ko uniref50_go uniref50_level4ec uniref50_pfam uniref50_eggnog metacyc_pathabundance kegg_pathabundance kegg_pathcoverage metacyc_pathcoverage genefamilies
parallel -j 1 --link 'humann_renorm_table {1}/merge/{2}_merge.txt -u relab --output {1}/merge/{2}_relab.txt ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3 ::: uniref50_ko uniref50_go uniref50_level4ec uniref50_pfam uniref50_eggnog  metacyc_pathabundance kegg_pathabundance genefamilies
 
# use uniref50_all
parallel 'mkdir -p {1}/{2} ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all ::: MNA00239 MNA00240 MNA00241 MNA00242 MNA00243 MNA00244 MNA00245 MNA00246 MNA00247 MNA00248 MNA00249 MNA00250 MNA00251 MNA00252 MNA00253 MNA00254 merge merge_cpm merge_filter merge_reb

humann_uniref50_all_1.sh
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel --link -j 1 'humann --threads 20 --input {2}/kneaddata/{1}/{1}.fastq --output {2}/humann3_uniref50_all/{1} --output-basename {1}_metacyc --protein-database {3} \
' ::: MNA00239 MNA00240 MNA00241 MNA00242 ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: /work/workspace/zhurj/database/humann/201901/uniref50_all
parallel --link -j 1 'rm {1}/{2}/{2}_metacyc_humann_temp -rf '  \
::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all ::: MNA00239 MNA00240 NNA00241 MNA00242
parallel -j 1 --link ' humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_ko --output  {1}/{2}/{2}_uniref50_ko.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_go --output  {1}/{2}/{2}_uniref50_go.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_level4ec --output  {1}/{2}/{2}_uniref50_level4ec.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_pfam --output  {1}/{2}/{2}_uniref50_pfam.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_eggnog --output  {1}/{2}/{2}_uniref50_eggnog.tsv \
 ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all ::: MNA00239 MNA00240 NNA00241 MNA00242
 parallel -j 1 --link 'humann -i {1}/{2}/{2}_uniref50_ko.tsv --pathways-database {3} --output {1}/{2} --threads 20 --output-basename {2}_kegg' \
 ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all ::: MNA00239 MNA00240 NNA00241 MNA00242 ::: /work/workspace/zhurj/database/humann/current/kegg/keggc
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

humann_uniref50_all_2.sh
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel --link -j 1 'humann --threads 20 --input {2}/kneaddata/{1}/{1}.fastq --output {2}/humann3_uniref50_all/{1} --output-basename {1}_metacyc --protein-database {3} \
' ::: MNA00243 MNA00244 MNA00245 MNA00246 ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: /work/workspace/zhurj/database/humann/201901/uniref50_all
parallel --link -j 1 'rm {1}/{2}/{2}_metacyc_humann_temp -rf '  \
::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all ::: MNA00243 MNA00244 MNA00245 MNA00246
parallel -j 1 --link ' humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_ko --output  {1}/{2}/{2}_uniref50_ko.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_go --output  {1}/{2}/{2}_uniref50_go.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_level4ec --output  {1}/{2}/{2}_uniref50_level4ec.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_pfam --output  {1}/{2}/{2}_uniref50_pfam.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_eggnog --output  {1}/{2}/{2}_uniref50_eggnog.tsv \
 ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all ::: MNA00243 MNA00244 MNA00245 MNA00246
 parallel -j 1 --link 'humann -i {1}/{2}/{2}_uniref50_ko.tsv --pathways-database {3} --output {1}/{2} --threads 20 --output-basename {2}_kegg' \
 ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all ::: MNA00243 MNA00244 MNA00245 MNA00246 ::: /work/workspace/zhurj/database/humann/current/kegg/keggc
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

humann_uniref50_all_3.sh
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel --link -j 1 'humann --threads 20 --input {2}/kneaddata/{1}/{1}.fastq --output {2}/humann3_uniref50_all/{1} --output-basename {1}_metacyc --protein-database {3} \
' ::: MNA00247 MNA00248 MNA00249 MNA00250 ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: /work/workspace/zhurj/database/humann/201901/uniref50_all
parallel --link -j 1 'rm {1}/{2}/{2}_metacyc_humann_temp -rf '  \
::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all ::: MNA00247 MNA00248 MNA00249 MNA00250
parallel -j 1 --link ' humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_ko --output  {1}/{2}/{2}_uniref50_ko.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_go --output  {1}/{2}/{2}_uniref50_go.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_level4ec --output  {1}/{2}/{2}_uniref50_level4ec.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_pfam --output  {1}/{2}/{2}_uniref50_pfam.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_eggnog --output  {1}/{2}/{2}_uniref50_eggnog.tsv \
 ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all ::: MNA00247 MNA00248 MNA00249 MNA00250
 parallel -j 1 --link 'humann -i {1}/{2}/{2}_uniref50_ko.tsv --pathways-database {3} --output {1}/{2} --threads 20 --output-basename {2}_kegg' \
 ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all ::: MNA00247 MNA00248 MNA00249 MNA00250 ::: /work/workspace/zhurj/database/humann/current/kegg/keggc
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

humann_uniref50_all_4.sh
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel --link -j 1 'humann --threads 20 --input {2}/kneaddata/{1}/{1}.fastq --output {2}/humann3_uniref50_all/{1} --output-basename {1}_metacyc --protein-database {3} \
' ::: MNA00251 MNA00252 MNA00253 MNA00254 ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: /work/workspace/zhurj/database/humann/201901/uniref50_all
parallel --link -j 1 'rm {1}/{2}/{2}_metacyc_humann_temp -rf '  \
::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all ::: MNA00251 MNA00252 MNA00253 MNA00254
parallel -j 1 --link ' humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_ko --output  {1}/{2}/{2}_uniref50_ko.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_go --output  {1}/{2}/{2}_uniref50_go.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_level4ec --output  {1}/{2}/{2}_uniref50_level4ec.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_pfam --output  {1}/{2}/{2}_uniref50_pfam.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_eggnog --output  {1}/{2}/{2}_uniref50_eggnog.tsv \
 ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all ::: MNA00251 MNA00252 MNA00253 MNA00254
 parallel -j 1 --link 'humann -i {1}/{2}/{2}_uniref50_ko.tsv --pathways-database {3} --output {1}/{2} --threads 20 --output-basename {2}_kegg' \
 ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all ::: MNA00251 MNA00252 MNA00253 MNA00254 ::: /work/workspace/zhurj/database/humann/current/kegg/keggc
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

cd /work/workspace/zhurj/project/1_metadata/p/p
time srun -o humann_uniref50_all_1.out -e humann_uniref50_all_1.err -N 1 -c 20 -p slurm128 -w mnclient03 bash humann_uniref50_all_1.sh &
time srun -o humann_uniref50_all_2.out -e humann_uniref50_all_2.err -N 1 -c 20 -p slurm128 -w mnclient04 bash humann_uniref50_all_2.sh &
time srun -o humann_uniref50_all_3.out -e humann_uniref50_all_3.err -N 1 -c 20 -p slurm256 -w mnclient01 bash humann_uniref50_all_3.sh &
time srun -o humann_uniref50_all_4.out -e humann_uniref50_all_4.err -N 1 -c 20 -p slurm256 -w mnclient02 bash humann_uniref50_all_4.sh &

source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel -j 1 --link ' humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_ko --output  {1}/{2}/{2}_uniref50_ko.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_go --output  {1}/{2}/{2}_uniref50_go.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_level4ec --output  {1}/{2}/{2}_uniref50_level4ec.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_pfam --output  {1}/{2}/{2}_uniref50_pfam.tsv && \
 humann_regroup_table -i {1}/{2}/{2}_metacyc_genefamilies.tsv -g uniref50_eggnog --output  {1}/{2}/{2}_uniref50_eggnog.tsv \
 ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all ::: MNA00241
 parallel -j 1 --link 'humann -i {1}/{2}/{2}_uniref50_ko.tsv --pathways-database {3} --output {1}/{2} --threads 20 --output-basename {2}_kegg' \
 ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all ::: MNA00241 ::: /work/workspace/zhurj/database/humann/current/kegg/keggc
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

parallel --link -j 1 'humann_join_tables -s --input {1} --file_name _{2} --output {1}/merge/{2}_merge.txt'   ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all ::: uniref50_ko uniref50_go uniref50_level4ec uniref50_pfam uniref50_eggnog metacyc_pathabundance kegg_pathabundance kegg_pathcoverage metacyc_pathcoverage genefamilies
parallel -j 1 --link 'humann_renorm_table -i {1}/merge/{2}_merge.txt -u cpm --output {1}/merge_cpm/{2}_cpm.txt ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all ::: uniref50_ko uniref50_go uniref50_level4ec uniref50_pfam uniref50_eggnog  metacyc_pathabundance kegg_pathabundance genefamilies
#parallel -j 1 --link 'humann_renorm_table -i {1}/merge/{2}_merge.txt -u relab --output {1}/merge/{2}_relab.txt ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all ::: uniref50_ko uniref50_go uniref50_level4ec uniref50_pfam uniref50_eggnog  metacyc_pathabundance kegg_pathabundance genefamilies

source activate python3.6
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all && python /work/workspace/zhurj/lib/python/script/humann_gene_filter.py merge_cpm/genefamilies_cpm.txt merge_filter/filter_genefamilies.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all && python /work/workspace/zhurj/lib/python/script/humann_others_filter.py merge_cpm/kegg_pathabundance_cpm.txt merge_filter/filter_kegg_pathabundance.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all && python /work/workspace/zhurj/lib/python/script/humann_others_filter.py merge_cpm/metacyc_pathabundance_cpm.txt merge_filter/filter_metacyc_pathabundance.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all && python /work/workspace/zhurj/lib/python/script/humann_others_filter.py merge_cpm/uniref50_eggnog_cpm.txt merge_filter/filter_eggnog.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all && python /work/workspace/zhurj/lib/python/script/humann_others_filter.py merge_cpm/uniref50_go_cpm.txt merge_filter/filter_go.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all && python /work/workspace/zhurj/lib/python/script/humann_others_filter.py merge_cpm/uniref50_ko_cpm.txt merge_filter/filter_ko.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all && python /work/workspace/zhurj/lib/python/script/humann_others_filter.py merge_cpm/uniref50_level4ec_cpm.txt merge_filter/filter_level4ec.tsv
cd /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all && python /work/workspace/zhurj/lib/python/script/humann_others_filter.py merge_cpm/uniref50_pfam_cpm.txt merge_filter/filter_pfam.tsv
source deactivate python3.6


# 2020-10-29
metacyc different ana
source activate /work/workspace/zhurj/bin/miniconda3/envs/py2
parallel --link 'lefse-format_input.py {1}/lefse/humann3_uniref50_ec/{2}/lefse_in {1}/lefse/humann3_uniref50_ec/{2}/lefse_format.in -c 1 -u 2 -f r -o 1000000 && \
run_lefse.py {1}/lefse/humann3_uniref50_ec/{2}/lefse_format.in {1}/lefse/humann3_uniref50_ec/{2}/lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1 \
' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: metacyc 
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/py2

source activate python3.6
parallel --link 'python /work/workspace/zhurj/lib/python/script/lefse_fitler.py {1}/lefse/humann3_uniref50_ec/{2}/lefse_run.res {1}/lefse/humann3_uniref50_ec/{2}/filter_lefse_run.res ' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: metacyc
source deactivate python3.6

kegg different ana
metacyc different ana
mkdir /work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3_uniref50_ec/kegg
source activate /work/workspace/zhurj/bin/miniconda3/envs/py2
parallel --link 'lefse-format_input.py {1}/lefse/humann3_uniref50_ec/{2}/lefse_in {1}/lefse/humann3_uniref50_ec/{2}/lefse_format.in -c 1 -u 2 -f r -o 1000000 && \
run_lefse.py {1}/lefse/humann3_uniref50_ec/{2}/lefse_format.in {1}/lefse/humann3_uniref50_ec/{2}/lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1 \
' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: kegg 
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/py2

source activate python3.6
parallel --link 'python /work/workspace/zhurj/lib/python/script/lefse_fitler.py {1}/lefse/humann3_uniref50_ec/{2}/lefse_run.res {1}/lefse/humann3_uniref50_ec/{2}/filter_lefse_run.res ' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: kegg
source deactivate python3.6

# Lefse in
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link 'python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py {1}/humann3/merge_filter/filter_{2}.tsv {1}/input/meta_group {1}/lefse/humann3_uniref50_ec/{2}/lefse_in' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: genefamilies kegg metacyc eggnog go ko level4ec pfam
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

# lefse
source activate /work/workspace/zhurj/bin/miniconda3/envs/py2
parallel --link 'lefse-format_input.py {1}/lefse/humann3/{2}/lefse_in {1}/lefse/humann3/{2}/lefse_format.in -c 1 -u 2 -f r -o 1000000 && \
run_lefse.py {1}/lefse/humann3/{2}/lefse_format.in {1}/lefse/humann3/{2}/lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1 \
' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: genefamilies kegg metacyc eggnog go ko level4ec pfam
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/py2

# lefse result filter
source activate python3.6
parallel --link 'python /work/workspace/zhurj/lib/python/script/lefse_fitler.py {1}/lefse/humann3/{2}/lefse_run.res {1}/lefse/humann3/{2}/filter_lefse_run.res ' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: kegg metacyc eggnog go ko level4ec pfam
source deactivate python3.6

# 2020-10-29


COG1961 2.63496553797   HFD     2.08595209301   0.000287113654735
ENOG4105YPD     2.67416601121   HFD     2.22160363873   0.00237446248379

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link 'python /work/workspace/zhurj/lib/python/script/humann_annotation.py {1}/lefse/humann3/{2}/filter_lefse_run.res \
/work/workspace/zhurj/database/humann/current/misc/map_{2}_name.txt.gz {1}/lefse/humann3/{2}/filter_lefse_run_desc.txt' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: pfam ko eggnog 
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link 'python /work/workspace/zhurj/lib/python/script/humann_kegg_annotation.py {1}/filter_lefse_run.res \
/work/workspace/zhurj/database/humann/current/kegg/kegg_path_description.txt {1}/filter_lefse_run_desc.txt ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/kegg
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link 'python /work/workspace/zhurj/lib/python/script/humann_ec_annotation.py {1}/filter_lefse_run.res \
/work/workspace/zhurj/database/humann/current/misc/map_ec_name.txt.gz {1}/filter_lefse_run_desc.txt ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/level4ec
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link 'python /work/workspace/zhurj/lib/python/script/humann_go_annotation.py {1}/filter_lefse_run.res \
/work/workspace/zhurj/database/humann/current/misc/map_go_name.txt.gz {1}/filter_lefse_run_desc.txt ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/go
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link 'python /work/workspace/zhurj/lib/python/script/humann2boxplot.py {1}/{2}/filter_lefse_run_desc.txt \
{1}/{2}/lefse_in {1}/{2}/graph_pre_in ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3 ::: pfam ko eggnog go kegg
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

parallel --link 'python /work/workspace/zhurj/lib/python/script/humann2boxplot.py {1}/{2}/filter_lefse_run_desc.txt \
{1}/{2}/lefse_in {1}/{2}/graph_pre_in ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3 ::: eggnog
parallel --link 'python /work/workspace/zhurj/lib/python/script/humann2boxplot.py {1}/{2}/filter_lefse_run_desc.txt \
{1}/{2}/lefse_in {1}/{2}/graph_pre_in ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3 ::: level4ec


parallel --link 'python /work/workspace/zhurj/lib/python/script/lefse2boxplot.py {1}/lefse/humann3/{2}/graph_pre_in {1}/input/meta_group {1}/lefse/humann3/{2}/graph_in ' \
::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: pfam ko eggnog go kegg metacyc level4ec


source activate /work/workspace/zhurj/bin/miniconda3/envs/R3.6
library(ggplot2)
library('RColorBrewer')
library('ggsci')
library('ggpubr')

# stat_compare_means(comparisons=my_comparisions, method = "wilcox.test", aes(label=..p.adj..)) +
# label = "p.signif", aes(label=..p.adj..)
# method = "t.test", method = "wilcox.test"
# kegg
Data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/kegg/graph_in",header=T,sep='\t')
my_comparisions <- list(c("NCD", "HFD"),c("NCD", "MNO863"),c("HFD", "MNO863"))
p <- ggboxplot(Data, x = "Group", y = "Abundance", color = "Group", palette = "npg",  facet.by = "Species", ylim = c(0, 16)) + 
ggtitle("Boxplot of KEGG")+
theme(legend.position="none",plot.title = element_text(face="bold.italic", hjust=0.5),)+
scale_y_continuous(name="Log2CPM(copy per million)", breaks=seq(0,16,2)) +
scale_x_discrete(name="") +
stat_compare_means(comparisons=my_comparisions, method = "wilcox.test", label = "p.signif") +
guides(fill=FALSE)
ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/pic/dif_kegg_wilcox_20201010.jpg", p, width = 10, height =8)

Data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/kegg/graph_in",header=T,sep='\t')
my_comparisions <- list(c("NCD", "HFD"),c("NCD", "MNO863"),c("HFD", "MNO863"))
p <- ggboxplot(Data, x = "Group", y = "Abundance", color = "Group", palette = "npg",  facet.by = "Species", ylim = c(0, 16)) + 
ggtitle("Boxplot of KEGG")+
theme(legend.position="none",plot.title = element_text(face="bold.italic", hjust=0.5),)+
scale_y_continuous(name="Log2CPM(copy per million)", breaks=seq(0,16,2)) +
scale_x_discrete(name="") +
stat_compare_means(comparisons=my_comparisions, method = "t.test", label = "p.signif") +
guides(fill=FALSE)
ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/pic/dif_kegg_ttest_20201010.jpg", p, width = 10, height =8)

#------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library('RColorBrewer')
library('ggsci')
library('ggpubr')

# pfam ko eggnog go kegg 
# stat_compare_means(comparisons=my_comparisions, method = "wilcox.test", aes(label=..p.adj..)) +
# label = "p.signif", aes(label=..p.adj..)
# method = "t.test", method = "wilcox.test"
# level4ec
Data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/level4ec/graph_in",header=T,sep='\t')
my_comparisions <- list(c("NCD", "HFD"),c("NCD", "MNO863"),c("HFD", "MNO863"))
p <- ggboxplot(Data, x = "Group", y = "Abundance", color = "Group", palette = "npg",  facet.by = "Species", ylim = c(8, 22)) + 
ggtitle("Boxplot of ec")+
theme(legend.position="none",plot.title = element_text(face="bold.italic", hjust=0.5),)+
scale_y_continuous(name="Log2CPM(copy per million)", breaks=seq(8,22,2)) +
scale_x_discrete(name="") +
stat_compare_means(comparisons=my_comparisions, method = "wilcox.test", label = "p.signif") +
guides(fill=FALSE)
ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/pic/dif_level4ec_wilcox_20201010.jpg", p, width = 6, height =6)

Data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/level4ec/graph_in",header=T,sep='\t')
my_comparisions <- list(c("NCD", "HFD"),c("NCD", "MNO863"),c("HFD", "MNO863"))
p <- ggboxplot(Data, x = "Group", y = "Abundance", color = "Group", palette = "npg",  facet.by = "Species", ylim = c(8, 22)) + 
ggtitle("Boxplot of ec")+
theme(legend.position="none",plot.title = element_text(face="bold.italic", hjust=0.5),)+
scale_y_continuous(name="Log2CPM(copy per million)", breaks=seq(8,22,2)) +
scale_x_discrete(name="") +
stat_compare_means(comparisons=my_comparisions, method = "t.test", label = "p.signif") +
guides(fill=FALSE)
ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/pic/dif_level4ec_ttest_20201010.jpg", p, width = 6, height =6)

#------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library('RColorBrewer')
library('ggsci')
library('ggpubr')

# pfam ko eggnog go kegg 
# stat_compare_means(comparisons=my_comparisions, method = "wilcox.test", aes(label=..p.adj..)) +
# label = "p.signif", aes(label=..p.adj..)
# method = "t.test", method = "wilcox.test"
# kegg
Data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/pfam/graph_in",header=T,sep='\t')
my_comparisions <- list(c("NCD", "HFD"),c("NCD", "MNO863"),c("HFD", "MNO863"))
p <- ggboxplot(Data, x = "Group", y = "Abundance", color = "Group", palette = "npg",  facet.by = "Species", ylim = c(0, 24)) + 
ggtitle("Boxplot of pfam")+
theme(legend.position="none",plot.title = element_text(face="bold.italic", hjust=0.5),)+
scale_y_continuous(name="Log2CPM(copy per million)", breaks=seq(0,24,2)) +
scale_x_discrete(name="") +
stat_compare_means(comparisons=my_comparisions, method = "wilcox.test", label = "p.signif") +
guides(fill=FALSE)
ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/pic/dif_pfam_wilcox_20201010.jpg", p, width = 30, height =25)

Data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/pfam/graph_in",header=T,sep='\t')
my_comparisions <- list(c("NCD", "HFD"),c("NCD", "MNO863"),c("HFD", "MNO863"))
p <- ggboxplot(Data, x = "Group", y = "Abundance", color = "Group", palette = "npg",  facet.by = "Species", ylim = c(0, 24)) + 
ggtitle("Boxplot of pfam")+
theme(legend.position="none",plot.title = element_text(face="bold.italic", hjust=0.5),)+
scale_y_continuous(name="Log2CPM(copy per million)", breaks=seq(0,24,2)) +
scale_x_discrete(name="") +
stat_compare_means(comparisons=my_comparisions, method = "t.test", label = "p.signif") +
guides(fill=FALSE)
ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/pic/dif_pfam_ttest_20201010.jpg", p, width = 30, height =25)


#------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library('RColorBrewer')
library('ggsci')
library('ggpubr')

# pfam ko eggnog go kegg 
# stat_compare_means(comparisons=my_comparisions, method = "wilcox.test", aes(label=..p.adj..)) +
# label = "p.signif", aes(label=..p.adj..)
# method = "t.test", method = "wilcox.test"
# kegg
Data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/eggnog/graph_in",header=T,sep='\t')
my_comparisions <- list(c("NCD", "HFD"),c("NCD", "MNO863"),c("HFD", "MNO863"))
p <- ggboxplot(Data, x = "Group", y = "Abundance", color = "Group", palette = "npg",  facet.by = "Species", ylim = c(0, 20)) + 
ggtitle("Boxplot of eggnog")+
theme(legend.position="none",plot.title = element_text(face="bold.italic", hjust=0.5),)+
scale_y_continuous(name="Log2CPM(copy per million)", breaks=seq(0,20,2)) +
scale_x_discrete(name="") +
stat_compare_means(comparisons=my_comparisions, method = "wilcox.test", label = "p.signif") +
guides(fill=FALSE)
ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/pic/dif_eggnog_wilcox_20201010.jpg", p, width = 4, height =4)

Data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/eggnog/graph_in",header=T,sep='\t')
my_comparisions <- list(c("NCD", "HFD"),c("NCD", "MNO863"),c("HFD", "MNO863"))
p <- ggboxplot(Data, x = "Group", y = "Abundance", color = "Group", palette = "npg",  facet.by = "Species", ylim = c(0, 20)) + 
ggtitle("Boxplot of eggnog")+
theme(legend.position="none",plot.title = element_text(face="bold.italic", hjust=0.5),)+
scale_y_continuous(name="Log2CPM(copy per million)", breaks=seq(0,20,2)) +
scale_x_discrete(name="") +
stat_compare_means(comparisons=my_comparisions, method = "t.test", label = "p.signif") +
guides(fill=FALSE)
ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/pic/dif_eggnog_ttest_20201010.jpg", p, width = 4, height =4)

#------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library('RColorBrewer')
library('ggsci')
library('ggpubr')

# pfam ko eggnog go kegg 
# stat_compare_means(comparisons=my_comparisions, method = "wilcox.test", aes(label=..p.adj..)) +
# label = "p.signif", aes(label=..p.adj..)
# method = "t.test", method = "wilcox.test"
# kegg
Data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/go/graph_in",header=T,sep='\t')
my_comparisions <- list(c("NCD", "HFD"),c("NCD", "MNO863"),c("HFD", "MNO863"))
p <- ggboxplot(Data, x = "Group", y = "Abundance", color = "Group", palette = "npg",  facet.by = "Species", ylim = c(0, 24)) + 
ggtitle("Boxplot of go")+
theme(legend.position="none",plot.title = element_text(face="bold.italic", hjust=0.5),)+
scale_y_continuous(name="Log2CPM(copy per million)", breaks=seq(0,24,2)) +
scale_x_discrete(name="") +
stat_compare_means(comparisons=my_comparisions, method = "wilcox.test", label = "p.signif") +
guides(fill=FALSE)
ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/pic/dif_go_wilcox_20201010.jpg", p, width = 45, height =40)

Data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/go/graph_in",header=T,sep='\t')
my_comparisions <- list(c("NCD", "HFD"),c("NCD", "MNO863"),c("HFD", "MNO863"))
p <- ggboxplot(Data, x = "Group", y = "Abundance", color = "Group", palette = "npg",  facet.by = "Species", ylim = c(0, 24)) + 
ggtitle("Boxplot of go")+
theme(legend.position="none",plot.title = element_text(face="bold.italic", hjust=0.5),)+
scale_y_continuous(name="Log2CPM(copy per million)", breaks=seq(0,24,2)) +
scale_x_discrete(name="") +
stat_compare_means(comparisons=my_comparisions, method = "t.test", label = "p.signif") +
guides(fill=FALSE)
ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/pic/dif_go_ttest_20201010.jpg", p, width = 45, height =40)

#------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library('RColorBrewer')
library('ggsci')
library('ggpubr')

# pfam ko eggnog go kegg 
# stat_compare_means(comparisons=my_comparisions, method = "wilcox.test", aes(label=..p.adj..)) +
# label = "p.signif", aes(label=..p.adj..)
# method = "t.test", method = "wilcox.test"
# kegg
Data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/metacyc/graph_in",header=T,sep='\t')
my_comparisions <- list(c("NCD", "HFD"),c("NCD", "MNO863"),c("HFD", "MNO863"))
p <- ggboxplot(Data, x = "Group", y = "Abundance", color = "Group", palette = "npg",  facet.by = "Species", ylim = c(8, 16)) + 
ggtitle("Boxplot of metacyc")+
theme(legend.position="none",plot.title = element_text(face="bold.italic", hjust=0.5),)+
scale_y_continuous(name="Log2CPM(copy per million)", breaks=seq(8,16,2)) +
scale_x_discrete(name="") +
stat_compare_means(comparisons=my_comparisions, method = "wilcox.test", label = "p.signif") +
guides(fill=FALSE)
ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/pic/dif_metacyc_wilcox_20201010.jpg", p, width = 4, height =4)

Data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/metacyc/graph_in",header=T,sep='\t')
my_comparisions <- list(c("NCD", "HFD"),c("NCD", "MNO863"),c("HFD", "MNO863"))
p <- ggboxplot(Data, x = "Group", y = "Abundance", color = "Group", palette = "npg",  facet.by = "Species", ylim = c(8, 16)) + 
ggtitle("Boxplot of metacyc")+
theme(legend.position="none",plot.title = element_text(face="bold.italic", hjust=0.5),)+
scale_y_continuous(name="Log2CPM(copy per million)", breaks=seq(8,16,2)) +
scale_x_discrete(name="") +
stat_compare_means(comparisons=my_comparisions, method = "t.test", label = "p.signif") +
guides(fill=FALSE)
ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3/pic/dif_metacyc_ttest_20201010.jpg", p, width = 4, height =4)

# --------------------------------------------------------------------------------------------------------------------------------
species PCA analysis
source activate R3.6
/work/workspace/zhurj/project/1_metadata/mouse22/metaphlan3/pic
library('ggbiplot')
data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/metaphlan3/merge/species_pca_in",sep='\t',header=TRUE,row.names=1)
data.pca <- prcomp(data, center = TRUE,scale. = TRUE)
#summary(data.pca)
#str(data.pca)
dtgroup <- c(rep("NCD",6),rep("HFD",8),rep("MNO863",8))
#ggbiplot(data.pca,ellipse=TRUE,  labels=rownames(data), groups=dtgroup)
p <- ggbiplot(data.pca,ellipse=TRUE,obs.scale = 1, var.scale = 1,var.axes=FALSE, labels=rownames(data), groups=dtgroup)
ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/metaphlan3/pic/pca_species34_20201010.jpg", p)

data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/metaphlan3/merge/genus_pca_in",sep='\t',header=TRUE,row.names=1)
data.pca <- prcomp(data, center = TRUE,scale. = TRUE)
dtgroup <- c(rep("NCD",6),rep("HFD",8),rep("MNO863",8))
#ggbiplot(data.pca,ellipse=TRUE,  labels=rownames(data), groups=dtgroup)
p <- ggbiplot(data.pca,ellipse=TRUE,obs.scale = 1, var.scale = 1,var.axes=FALSE, labels=rownames(data), groups=dtgroup)
ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/metaphlan3/pic/pca_genus23_20201010.jpg", p)


192.168.2.74
Principal Component Analysis in R: https://www.datacamp.com/community/tutorials/pca-analysis-r
library('ggbiplot')
data = read.table("/work/workspace/ruijuan/project/7_metapipe/mouse22/pca/species_pca_in",sep='\t',header=TRUE,row.names=1)
data.pca <- prcomp(data, center = TRUE,scale. = TRUE)
summary(data.pca)
str(data.pca)
ggbiplot(data.pca)
ggbiplot(data.pca, labels=rownames(data))
dtgroup <- c(rep("NCD",6),rep("HFD",8),rep("MNO863",8))
ggbiplot(data.pca,ellipse=TRUE,  labels=rownames(data), groups=dtgroup)
p <- ggbiplot(data.pca,ellipse=TRUE,obs.scale = 1, var.scale = 1,var.axes=FALSE, labels=rownames(data), groups=dtgroup)
ggsave("/work/workspace/ruijuan/project/7_metapipe/mouse22/pic/pca_species34_20201010.jpg", p)

library('ggbiplot')
data = read.table("/work/workspace/ruijuan/project/7_metapipe/mouse22/pca/genus_pca_in",sep='\t',header=TRUE,row.names=1)
data.pca <- prcomp(data, center = TRUE,scale. = TRUE)
dtgroup <- c(rep("NCD",6),rep("HFD",8),rep("MNO863",8))
#ggbiplot(data.pca,ellipse=TRUE,  labels=rownames(data), groups=dtgroup)
p <- ggbiplot(data.pca,ellipse=TRUE,obs.scale = 1, var.scale = 1,var.axes=FALSE, labels=rownames(data), groups=dtgroup)
ggsave("/work/workspace/ruijuan/project/7_metapipe/mouse22/pic/pca_genus23_20201010.jpg", p)


source activate /work/workspace/liangzj/program/Miniconda3/envs/metaphlan3
/work/workspace/liangzj/program/Miniconda3/envs/metaphlan3/bin/metaphlan /work/workspace/liangzj/project/metagenomics/Obesity/PRJEB37249/data/M0x10MCx3421.clean.rmHost.fastq.gz --input_type fastq --nproc 12 -o M0x10MCx3421.txt
source activate /work/workspace/liangzj/program/Miniconda3/envs/metaphlan3

metaphlan_onein.sh
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel --link -j 1 'metaphlan {1}/kneaddata/{2}/{2}.fastq --input_type fastq  --nproc 12 --sample_id_key {2} \
--bowtie2out {1}/metaphlan20201012/{2}/metagenome.bowtie2.bz2 -o {1}/metaphlan20201012/{2}/all.txt' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: MNA00233
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

cd /work/workspace/zhurj/project/1_metadata/metapipe20200921/p
time srun -o metaphlan_onein.out -e metaphlan_onein.err -N 1 -c 20 -p slurm128 -w mnclient03 bash metaphlan_onein.sh &

metaphlan_cleanone.sh
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel --link -j 1 'metaphlan {1}/clean/{2}.1.fq.gz,{1}/clean/{2}.2.fq.gz --input_type fastq  --nproc 12 --sample_id_key {2} \
--bowtie2out {1}/metaphlan20201012/{2}_clean/metagenome.bowtie2.bz2 -o {1}/metaphlan20201012/{2}_clean/all.txt' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: MNA00233
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

cd /work/workspace/zhurj/project/1_metadata/metapipe20200921/p
time srun -o metaphlan_cleanone.out -e metaphlan_cleanone.err -N 1 -c 20 -p slurm128 -w mnclient03 bash metaphlan_cleanone.sh &

metaphlanclean_MNA00233.sh
source activate humann
parallel -j -2 --link ' metaphlan {1}/{2}/metagenome.bowtie2.bz2 --input_type bowtie2out --tax_lev "p" --nproc 12 --sample_id_key {2} -o {1}/{2}/phylum.txt && metaphlan {1}/{2}/metagenome.bowtie2.bz2 --input_type bowtie2out --tax_lev "g" --nproc 12 --sample_id_key {2} -o {1}/{2}/genus.txt && metaphlan {1}/{2}/metagenome.bowtie2.bz2 --input_type bowtie2out --tax_lev "s" --nproc 12 --sample_id_key {2} -o {1}/{2}/species.txt '  ::: /work/workspace/zhurj/project/1_metadata/mouse22/metaphlan20201012 ::: MNA00233_clean
source deactivate humann

cd /work/workspace/zhurj/project/1_metadata/metapipe20200921/p
time srun -o metaphlanclean_MNA00233.out -e metaphlanclean_MNA00233.err -N 1 -c 20 -p slurm256 -w mnclient01 bash metaphlanclean_MNA00233.sh &


kraken2_MNA00233.sh
source activate /work/workspace/zhurj/bin/miniconda3/envs/kraken2
parallel --link -j -1 'kraken2 --db /work/workspace/zhurj/reference/kraken2/20200929/current --threads 36  \
--report  {1}/kraken2/{2}/kraken2.report --output {1}/kraken2/{2}/kraken2.output --paired {1}/bwarmhost/{2}/{2}_r1.fastq {1}/bwarmhost/{2}/{2}_r2.fastq && \
bracken -d /work/workspace/zhurj/reference/kraken2/20200929/current -i {1}/kraken2/{2}/kraken2.report -o {1}/kraken2/{2}/bracken -w {1}/kraken2/{2}/bracken.report -r 150 -l S && \
kreport2mpa.py -r {1}/kraken2/{2}/bracken.report -o {1}/kraken2/{2}/bracken.new.report' \
::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: MNA00233
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/kraken2

time srun -o kraken2_MNA00233.out -e kraken2_MNA00233.err -N 1 -c 20 -p slurm256 -w mnclient01 bash kraken2_MNA00233.sh &

kraken2_20201013.sh
parallel --link 'mkdir -p {1}/{2} ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/kraken2 :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22
parallel --link 'mkdir -p {1}/{2} ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/kraken2 ::: combine
parallel --link 'mkdir -p {1}/{2} ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/kraken2 ::: count percentage

source activate /work/workspace/zhurj/bin/miniconda3/envs/kraken2

parallel --link -j 1 'kraken2 --db /work/workspace/zhurj/reference/kraken2/20200929/current --threads 36  \
--report  {1}/kraken2/{2}/kraken2.report --output {1}/kraken2/{2}/kraken2.output --paired {1}/bwarmhost/{2}/{2}_r1.fastq {1}/bwarmhost/{2}/{2}_r2.fastq && \
bracken -d /work/workspace/zhurj/reference/kraken2/20200929/current -i {1}/kraken2/{2}/kraken2.report -o {1}/kraken2/{2}/bracken -w {1}/kraken2/{2}/bracken.report -r 150 -l S && \
kreport2mpa.py -r {1}/kraken2/{2}/bracken.report -o {1}/kraken2/{2}/bracken.new.report && \
ln -s {1}/kraken2/{2}/bracken.new.report {1}/kraken2/count/{2} && \
kreport2mpa.py --percentages -r {1}/kraken2/{2}/bracken.report -o {1}/kraken2/{2}/bracken.percentage.report  && \
ln -s {1}/kraken2/{2}/bracken.percentage.report {1}/kraken2/percentage/{2} ' \
::: /work/workspace/zhurj/project/1_metadata/mouse22 :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22

parallel --link -j 1 'cd {1}/count && \
combine_mpa.py -i MNA00233 MNA00234 MNA00235 MNA00236 MNA00237 MNA00238 MNA00239 MNA00240 MNA00241 MNA00242 MNA00243 MNA00244 MNA00245 MNA00246 MNA00247 MNA00248 MNA00249 MNA00250 MNA00251 MNA00252 MNA00253 MNA00254 -o {1}/combine/combine_count && \
cd {1}/percentage && \
combine_mpa.py -i MNA00233 MNA00234 MNA00235 MNA00236 MNA00237 MNA00238 MNA00239 MNA00240 MNA00241 MNA00242 MNA00243 MNA00244 MNA00245 MNA00246 MNA00247 MNA00248 MNA00249 MNA00250 MNA00251 MNA00252 MNA00253 MNA00254 -o {1}/combine/combine_percentage ' \
::: /work/workspace/zhurj/project/1_metadata/mouse22/kraken2

source deactivate /work/workspace/zhurj/bin/miniconda3/envs/kraken2

time srun -o kraken2_20201013.out -e kraken2_20201013.err -N 1 -c 20 -p slurm256 -w mnclient01 bash kraken2_20201013.sh &


test MAG mapped read

index reference: /work/workspace/zhurj/reference/NCBI/genome/GCA_000020225/bwa/GCA_000020225.fna.gz

parallel --link -j 2 'mkdir -p {2}/{1}' :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22 \
::: /work/workspace/zhurj/project/1_metadata/mouse22/Akk
bwa index 

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
bwa index /work/workspace/zhurj/project/1_metadata/mouse22/metabat/merge/ref233.fa
parallel --link -j 1 'bwa mem {3} {1}/bwarmhost/{2}/{2}_r1.fastq {1}/bwarmhost/{2}/{2}_r2.fastq -t 40 -o {1}/metabat/maptest/{2}.sam && \
samtools view -@ 40 -bS {1}/metabat/maptest/{2}.sam -o {1}/metabat/maptest/{2}.bam && \
samtools view -b -f 2 -@ 40 {1}/metabat/maptest/{2}.bam -o {1}/metabat/maptest/{2}_PEmapped.bam && \
samtools sort -n -@ 40 {1}/metabat/maptest/{2}_PEmapped.bam -o {1}/metabat/maptest/{2}_PEmapped_sorted.bam && \
bedtools bamtofastq -i {1}/metabat/maptest/{2}_PEmapped_sorted.bam -fq {1}/metabat/maptest/{2}_mapped.r1.fastq -fq2 {1}/metabat/maptest/{2}_mapped.r2.fastq && \
rm {1}/metabat/maptest/*.sam {1}/metabat/maptest/*.bam -f ' \
::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: MNA00233 ::: /work/workspace/zhurj/project/1_metadata/mouse22/metabat/merge/ref233.fa
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

cd /work/workspace/zhurj/project/1_metadata/mouse22/p
time srun -o metabatmapped.out -e metabatmapped.err -N 1 -c 20 -p slurm256 -w mnclient01 bash metabatmapped.sh &


metabatmap_20201014.sh
parallel -j 1 'mkdir -p {1}/{2}/{3}' ::: /work/workspace/zhurj/project/1_metadata/mouse22/metabat ::: map merge :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link -j 1 'cat {1}/metabat/{2}/*.fa > {1}/metabat/merge/{2}/ref.fa && \
bwa index {1}/metabat/merge/{2}/ref.fa && \
bwa mem {1}/metabat/merge/{2}/ref.fa {1}/bwarmhost/{2}/{2}_r1.fastq {1}/bwarmhost/{2}/{2}_r2.fastq -t 40 -o {1}/metabat/map/{2}/{2}.sam && \
samtools view -@ 40 -bS {1}/metabat/map/{2}/{2}.sam -o {1}/metabat/map/{2}/{2}.bam && \
samtools view -b -f 2 -@ 40 {1}/metabat/map/{2}/{2}.bam -o {1}/metabat/map/{2}/{2}_PEmapped.bam && \
samtools sort -n -@ 40 {1}/metabat/map/{2}/{2}_PEmapped.bam -o {1}/metabat/map/{2}/{2}_PEmapped_sorted.bam && \
bedtools bamtofastq -i {1}/metabat/map/{2}/{2}_PEmapped_sorted.bam -fq {1}/metabat/map/{2}/{2}_mapped.r1.fastq -fq2 {1}/metabat/map/{2}/{2}_mapped.r2.fastq && \
rm {1}/metabat/map/{2}/*.sam {1}/metabat/map/{2}/*.bam -f ' \
::: /work/workspace/zhurj/project/1_metadata/mouse22 :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

cd /work/workspace/zhurj/project/1_metadata/mouse22/p
time srun -o metabatmap_20201014.out -e metabatmap_20201014.err -N 1 -c 20 -p slurm256 -w mnclient01 bash metabatmap_20201014.sh &



parallel --link 'mkdir -p {1}/{2} ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/metaphlan3 ::: phylum_count genus_count species_count merge_count
metaphlan3_count.sh
'''
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel -j 1 --link ' metaphlan {1}/{2}/metagenome.bowtie2.bz2 -t "reads_map" --input_type bowtie2out --tax_lev "p" --nproc 12 --sample_id_key {2} -o {1}/{2}/phylum_count.txt '  \
::: /work/workspace/zhurj/project/1_metadata/mouse22/metaphlan3 :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/m2_kneaddata_id
parallel -j 1 --link 'ln -s {1}/{2}/phylum_count.txt {3}/{2}.txt ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/metaphlan3 :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/m2_kneaddata_id ::: /work/workspace/zhurj/project/1_metadata/mouse22/metaphlan3/phylum_count
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

'''
time srun -o metaphlan3_count.out -e metaphlan3_count.err -N 1 -c 20 -p slurm128 -w mnclient03 bash metaphlan3_count.sh &
map 到phylum，genus，species read 相同
cat MNA00235.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cd /work/workspace/zhurj/project/1_metadata/mouse22/metaphlan3/phylum_count
cat MNA00233.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNA00234.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNA00235.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNA00236.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNA00237.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNA00238.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNA00239.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNA00240.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNA00241.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNA00242.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNA00243.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNA00244.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNA00245.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNA00246.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNA00247.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNA00248.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNA00249.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNA00250.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNA00251.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNA00252.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNA00253.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNA00254.txt | awk -F "#" '{print $1}' | sort -u | wc -l

human faces
metaphlan3_human_count.sh
'''
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel -j 1 --link ' metaphlan {1}/{2}/metagenome.bowtie2.bz2 -t "reads_map" --input_type bowtie2out --tax_lev "p" --nproc 12 --sample_id_key {2} -o {1}/{2}/phylum_count.txt '  \
::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3 :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/mnid
parallel -j 1 --link 'ln -s {1}/{2}/phylum_count.txt {3}/{2}.txt ' ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3 :::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/input/mnid ::: /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/phylum_count
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

'''
cd /work/workspace/zhurj/project/1_metadata/metapipe20200921/p
time srun -o metaphlan3_human_count.out -e metaphlan3_human_count.err -N 1 -c 20 -p slurm128 -w mnclient03 bash metaphlan3_human_count.sh &
cd /work/workspace/zhurj/project/1_metadata/metapipe20200921/clean
zcat MNC00223.clean.2.fq.gz | grep @ | wc -l
zcat MNC00224.clean.2.fq.gz | grep @ | wc -l
zcat MNC00225.clean.2.fq.gz | grep @ | wc -l
zcat MNC00226.clean.2.fq.gz | grep @ | wc -l
zcat MNC00227.clean.2.fq.gz | grep @ | wc -l
zcat MNC00228.clean.2.fq.gz | grep @ | wc -l
zcat MNC00229.clean.2.fq.gz | grep @ | wc -l
zcat MNC00230.clean.2.fq.gz | grep @ | wc -l
zcat MNC00231.clean.2.fq.gz | grep @ | wc -l
zcat MNC00232.clean.2.fq.gz | grep @ | wc -l
zcat MNC00233.clean.2.fq.gz | grep @ | wc -l

cd /work/workspace/zhurj/project/1_metadata/metapipe20200921/metaphlan3/phylum_count
cat MNC00223.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNC00224.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNC00225.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNC00226.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNC00227.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNC00228.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNC00229.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNC00230.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNC00231.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNC00232.txt | awk -F "#" '{print $1}' | sort -u | wc -l
cat MNC00233.txt | awk -F "#" '{print $1}' | sort -u | wc -l


quast -l MNA00233.1 -t 16 -o /work/workspace/zhurj/project/1_metadata/mouse22/quast/MNA00233/mag1 /work/workspace/zhurj/project/1_metadata/mouse22/metabat/MNA00233/MNA00233.1.fa

cd /work/workspace/zhurj/project/1_metadata/mouse22/metabat
find `pwd` | grep fa$ | grep -v ref > input/sgbf_in
find `pwd` | grep fa$ | grep -v ref | awk -F "/" '{print $(NF-1)}'> input/sgbsam_in
find `pwd` | grep fa$ | grep -v ref | awk -F "[.]" '{print $(NF-1)}'> input/sgbindex_in

find `pwd` | grep fa$ | grep -v ref > input/sgbf_20201215_in

quast_1389.sh
parallel --link -j 1 'mkdir -p {1}/quast/{2}/{3}' \
::: /work/workspace/zhurj/project/1_metadata/mouse22 :::: /work/workspace/zhurj/project/1_metadata/mouse22/metabat/input/sgbsam_in :::: /work/workspace/zhurj/project/1_metadata/mouse22/metabat/input/sgbindex_in

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link -j 1 'quast -l {2}.{3} -t 36 -o {1}/quast/{2}/{3} {4} ' \
::: /work/workspace/zhurj/project/1_metadata/mouse22 :::: /work/workspace/zhurj/project/1_metadata/mouse22/metabat/input/sgbsam_in :::: /work/workspace/zhurj/project/1_metadata/mouse22/metabat/input/sgbindex_in :::: /work/workspace/zhurj/project/1_metadata/mouse22/metabat/input/sgbf_in
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6


cd /work/workspace/zhurj/project/1_metadata/mouse22/checkm/merge
cat *.txt > checkm
cat checkm | sort -u > checkm_uniq

checkm lineage_wf /work/workspace/zhurj/project/1_metadata/mouse22/metabat/test /work/workspace/zhurj/project/1_metadata/mouse22/checkm/test1 -x fa -t 36 --tab_table -f /work/workspace/zhurj/project/1_metadata/mouse22/checkm/test/checkm.txt

parallel --link -j 1 'mkdir -p {1}/checkm/{2} ' ::: /work/workspace/zhurj/project/1_metadata/mouse22 :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22
checkm_sam22.sh
'''
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link -j 1 'checkm lineage_wf {1}/metabat/{2} {1}/checkm/{2} -x fa -t 36 --tab_table -f {1}/checkm/{2}/checkm.txt ' \
::: /work/workspace/zhurj/project/1_metadata/mouse22  :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

'''
cd /work/workspace/zhurj/project/1_metadata/mouse22/p
time srun -o checkm_sam22.out -e checkm_sam22.err -N 1 -c 20 -p slurm256 -w mnclient01 bash checkm_sam22.sh &

# 2020-10-30
MNA00252, MNA00245 quast, checkm
cd /work/workspace/zhurj/project/1_metadata/mouse22/metabat
find `pwd` | grep fa$ | grep -v ref | grep -P "MNA00245|MNA00252" > input/sam2_sgbf_in
find `pwd` | grep fa$ | grep -v ref | grep -P "MNA00245|MNA00252" | awk -F "/" '{print $(NF-1)}'> input/sam2_sgbsam_in
find `pwd` | grep fa$ | grep -v ref | grep -P "MNA00245|MNA00252" | awk -F "[.]" '{print $(NF-1)}'> input/sam2_sgbindex_in

quast_sam2_134.sh
parallel --link -j 1 'mkdir -p {1}/quast/{2}/{3}' \
::: /work/workspace/zhurj/project/1_metadata/mouse22 :::: /work/workspace/zhurj/project/1_metadata/mouse22/metabat/input/sam2_sgbsam_in :::: /work/workspace/zhurj/project/1_metadata/mouse22/metabat/input/sam2_sgbindex_in

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link -j 1 'quast -l {2}.{3} -t 36 -o {1}/quast/{2}/{3} {4} ' \
::: /work/workspace/zhurj/project/1_metadata/mouse22 :::: /work/workspace/zhurj/project/1_metadata/mouse22/metabat/input/sam2_sgbsam_in :::: /work/workspace/zhurj/project/1_metadata/mouse22/metabat/input/sam2_sgbindex_in :::: /work/workspace/zhurj/project/1_metadata/mouse22/metabat/input/sam2_sgbf_in
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

cd /work/workspace/zhurj/project/1_metadata/p/p
time srun -o quast_sam2_134.out -e quast_sam2_134.err -N 1 -c 20 -p slurm128 -w mnclient03 bash quast_sam2_134.sh &
未在slurm上run

checkm lineage_wf /work/workspace/zhurj/project/1_metadata/mouse22/metabat/test /work/workspace/zhurj/project/1_metadata/mouse22/checkm/test1 -x fa -t 36 --tab_table -f /work/workspace/zhurj/project/1_metadata/mouse22/checkm/test/checkm.txt

#parallel --link -j 1 'mkdir -p {1}/checkm/{2} ' ::: /work/workspace/zhurj/project/1_metadata/mouse22 :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22
checkm_sam2.sh
'''
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link -j 1 'checkm lineage_wf {1}/metabat/{2} {1}/checkm/{2} -x fa -t 36 --tab_table -f {1}/checkm/{2}/checkm.txt ' \
::: /work/workspace/zhurj/project/1_metadata/mouse22  ::: MNA00245 MNA00252
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

'''



time srun -o prokka253clean.out -e prokka253clean.err -N 1 -c 20 -p slurm256 -w mnclient02 bash prokka253clean.sh &
time srun -o prokka242.out -e prokka242.err -N 1 -c 20 -p slurm128 -w mnclient03 bash prokka242.sh &
time srun -o checkm_sam8.out -e checkm_sam8.err -N 1 -c 20 -p slurm256 -w mnclient01 bash checkm_sam8.sh &

parallel --link -j 1 'ln -s {1}/{2}/checkm.txt {1}/merge/{2}.txt ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/checkm :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22

/work/workspace/zhurj/project/1_metadata/mouse22/quast
find `pwd` | grep transposed_report.tsv > report/file1389
find `pwd` | grep transposed_report.tsv | awk -F "/" '{print $(NF-2)}' > report/samf
find `pwd` | grep transposed_report.tsv | awk -F "/" '{print $(NF-1)}' > report/gindex


source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link -j 1 'ln -s {2} {1}/softlink/{3}_{4}.txt ' \
::: /work/workspace/zhurj/project/1_metadata/mouse22/quast :::: /work/workspace/zhurj/project/1_metadata/mouse22/quast/report/file1389 :::: /work/workspace/zhurj/project/1_metadata/mouse22/quast/report/samf :::: /work/workspace/zhurj/project/1_metadata/mouse22/quast/report/gindex
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

# ----------------------------------------------------------------------------------
cdhit
parallel --link 'ln -s {1}/{2}/{2}.faa {1}/faa/{2}.faa ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/prokka :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22
parallel --link 'ln -s {1}/{2}/{2}.ffn {1}/ffn/{2}.ffn ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/prokka :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22

parallel 'cat {1}/prokka/faa/*.faa > {1}/cdhit/in/merge.faa ' ::: /work/workspace/zhurj/project/1_metadata/mouse22
parallel 'cat {1}/prokka/ffn/*.ffn > {1}/cdhit/in/merge.ffn ' ::: /work/workspace/zhurj/project/1_metadata/mouse22

# remove length < 100nt items
source activate python3.6
cd /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/in
seqtk seq -L 100 merge.ffn > merge_100.ffn
seqtk seq -L 30 merge.faa > merge_30.faa

cdhit_ffn.sh
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link 'cd-hit-est -i {1}/in/merge_100.ffn -o {1}/out/cdhit.ffn -c 0.95 -G 0 -aS 0.9 -g 1 -d 0 -T 36 -M 0 ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/cdhit
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
cd  /work/workspace/zhurj/project/1_metadata/p/p
time srun -o cdhit_ffn.out -e cdhit_ffn.err -N 1 -c 20 -p slurm256 -w mnclient01 bash cdhit_ffn.sh &

cdhit_faa.sh
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link 'cd-hit -i {1}/in/merge_30.faa -o {1}/faa/cdhit.faa -c 0.95 -G 0 -aS 0.9 -g 1 -d 0 -T 36 -M 0 ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/cdhit
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

cd  /work/workspace/zhurj/project/1_metadata/p/p
time srun -o cdhit_ffn.out -e cdhit_ffn.err -N 1 -c 20 -p slurm256 -w mnclient01 bash cdhit_ffn.sh &
time srun -o cdhit_faa.out -e cdhit_faa.err -N 1 -c 20 -p slurm256 -w mnclient01 bash cdhit_faa.sh &

# Bowtie2使用方法与参数详细介绍 https://www.plob.org/article/4540.html
parallel --link 'mkdir -p {1}/{2}' ::: /work/workspace/zhurj/project/1_metadata/mouse22/genecatalog :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22
cdhit_bowtie2.sh
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel 'bowtie2-build --thread 36 {1}/cdhit.fna {1}/cdhit ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/out
parallel --link -j 1 'bowtie2 -x {1}/cdhit/out/cdhit -1 {1}/bwarmhost/{2}/{2}_r1.fastq -2 {1}/bwarmhost/{2}/{2}_r2.fastq \
--end-to-end --sensitive -I 200 -X 400 \
--threads 36  --reorder  --no-mixed --no-sq --no-unal \
-S {1}/genecatalog/{2}/{2}.sam \
' ::: /work/workspace/zhurj/project/1_metadata/mouse22 :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

cd  /work/workspace/zhurj/project/1_metadata/metapipe/p
time srun -o cdhit_bowtie2.out -e cdhit_bowtie2.err -N 1 -c 20 -p slurm256 -w mnclient01 bash cdhit_bowtie2.sh &
···
145 read paired，read reverse strand，second in pair
147 read paired，read mapped in proper pair，read reverse strand，second in pair
161 read paired，mate reverse strand，second in pair
163 read paired，read mapped in proper pair，mate reverse strand，second in pair
81 read paired，read reverse strand，first in pair
83 read paired，read mapped in proper pair，read reverse strand，first in pair
97 read paired，mate reverse strand， first in pair
99 read paired，read mapped in proper pair，mate reverse strand，first in pair
···

计算基因，及基因长度
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
genesize.sh
parallel --link 'infoseq -auto -nocolumns -delimiter ',' -only -noheading -name -length -sformat pearson {1}/cdhit.fna > {1}/cdhit_namesize.txt' ::: /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/out
source deactiate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
time srun -o genesize.out -e genesize.err -N 1 -c 20 -p slurm256 -w mnclient01 bash genesize.sh &

mapgenestat.sh

import pandas as pd
df = pd.read_table('test',skiprows=2,header=None,usecols=[0,2])
df.drop_duplicates(inplace = True)
dfcount = df.iloc[:,1].value_counts()
dfcount = dfcount [ dfcount > 2 ]
dfref = pd.read_table('ref',sep=',',header=None,index_col=0)
dfmerge = pd.merge(dfref, pd.DataFrame(dfcount), left_index=True, right_index=True)
dfmerge.columns=["size","count"]
dfmerge['cs'] = dfmerge['count'] / dfmerge['size']
dfmerge['csum'] = dfmerge['cs'].sum(axis=0)
dfmerge['abundance'] = dfmerge['cs'] / dfmerge['csum']
dfnew = dfmerge['abundance']
dfnew.to_csv('abundance.txt',sep='\t',header=None,index=True,encoding='utf-8')

esult = pd.merge(left, right, how='outer', on=['key1', 'key2'])

3. 挑选 >2 个reads 的基因，并计算每个基因相对丰度
python /work/workspace/zhurj/script/2_metapro/metapipe/geneabundance.py /work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/input/sam_in \
/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/out/cdhit_namesize.txt /work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out

time srun -o geneabundance.out -e geneabundance.err -N 1 -c 20 -p slurm128 -w mnclient03 bash geneabundance.sh &

物种鉴定-有参
python /work/workspace/zhurj/script/2_metapro/metapipe/orgSGB.py /work/workspace/zhurj/project/1_metadata/mouse22/SGB/taxonomy_in \
/work/workspace/zhurj/project/1_metadata/mouse22/SGB/taxonomy_species_out species

python /work/workspace/zhurj/script/2_metapro/metapipe/orgSGB.py /work/workspace/zhurj/project/1_metadata/mouse22/SGB/taxonomy_in \
/work/workspace/zhurj/project/1_metadata/mouse22/SGB/taxonomy_genus_out genus

python /work/workspace/zhurj/script/2_metapro/metapipe/orgSGB.py /work/workspace/zhurj/project/1_metadata/mouse22/SGB/taxonomy_in \
/work/workspace/zhurj/project/1_metadata/mouse22/SGB/taxonomy_phylum_out phylum


#-----------------------------------------------------------------------------------------------------------------------------------
物种统计及分析
(df > 0).astype(int).sum(axis=1)


#-----------------------------------------------------------------------------------------------------------------------------------
coabundance
方法一.Bonferroni
“最简单严厉的方法”
例如，如果检验1000次，我们就将阈值设定为5%/ 1000 = 0.00005；即使检验1000次，犯错误的概率还是保持在N×1000 = 5%。最终使得预期犯错误的次数不到1次，抹杀了一切假阳性的概率。
该方法虽然简单，但是检验过于严格，导致最后找不到显著表达的蛋白（假阴性）。
方法二.FalseDiscovery Rate
“比较温和的方法校正P值”
FDR（假阳性率）错误控制法是Benjamini于1995年提出的一种方法，基本原理是通过控制FDR值来决定P值的值域。相对Bonferroni来说，FDR用比较温和的方法对p值进行了校正。其试图在假阳性和假阴性间达到平衡，将假/真阳性比例控制到一定范围之内。例如，如果检验1000次，我们设定的阈值为0.05（5%），那么无论我们得到多少个差异蛋白，这些差异蛋白出现假阳性的概率保持在5%之内，这就叫FDR＜5%。
那么我们怎么从p value 来估算FDR呢，人们设计了几种不同的估算模型。其中使用最多的是Benjamini and Hochberg方法，简称BH法。虽然这个估算公式并不够完美，但是也能解决大部分的问题，主要还是简单好用！
FDR的计算方法
除了可以使用excel的BH计算方法外，对于较大的数据，我们推荐使用R命令p.adjust。
*************
SparCC的微生物网络构建示例:  https://blog.csdn.net/woodcorpse/article/details/106554536
微生物网络分析network: https://www.jianshu.com/p/e6735d65c5e9

Cytoscape
# 使用网络图展示Venn图集合及Cytoscape操作视频: https://blog.csdn.net/woodcorpse/article/details/106554605 #### 很好的网站
# cytoscape的cytohubba及MCODE插件寻找子网络hub基因： http://www.bio-info-trainee.com/5879.html
genus
cd /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/cytoscape/genus
source activate R3.6
R
library(plyr)
symbols_list = lapply(list.files(pattern = '*.csv'), function(f){
	data = read.csv(f, skip = 1, header = TRUE,  quote = '', check.names = F, row.names = NULL)
	symbols = as.character(data$Name[1:10])
	return(symbols)
})
#Reduce(intersect,symbols_list[2:9])
ldf <- ldply (symbols_list, data.frame)
colnames(ldf) = 'taxo'
freq_df <- as.data.frame(table(ldf['taxo']))
sort_df <- freq_df[order(-freq_df$Freq),]
write.table(sort_df,file='/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/cytoscape/genus/hub_genus.txt',sep='\t',row.names = TRUE,col.names = TRUE)


cd /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/cytoscape/species
source activate R3.6
R
library(plyr)
symbols_list = lapply(list.files(pattern = '*.csv'), function(f){
	data = read.csv(f, skip = 1, header = TRUE,  quote = '', check.names = F, row.names = NULL)
	symbols = as.character(data$Name[1:10])
	return(symbols)
})
#Reduce(intersect,symbols_list[2:9])
ldf <- ldply (symbols_list, data.frame)
colnames(ldf) = 'taxo'
freq_df <- as.data.frame(table(ldf['taxo']))
sort_df <- freq_df[order(-freq_df$Freq),]
write.table(sort_df,file='/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/cytoscape/species/hub_genus.txt',sep='\t',row.names = TRUE,col.names = TRUE)


source activate py2
SparCC.py in/
parallel --link 'SparCC.py {1}/in/otu_table.txt -i 20 --cor_file={1}/out/cor_phylum_sparcc.out' ::: /work/workspace/zhurj/project/1_metadata/mouse22/coabundance

SGB genus
cd /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data
source activate python3.6
import pandas as pd
df = pd.read_table('lefse_genus_in',sep='\t',header=0,index_col = 0)
ndf = df * 50000
nndf = ndf.astype(int)
nndf.to_csv('sparcc_genus_in',sep='\t',header=True,index=True)

species
cd /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data
source activate python3.6
import pandas as pd
df = pd.read_table('lefse_species_in',sep='\t',header=0,index_col = 0)
ndf = df * 50000
nndf = ndf.astype(int)
nndf.to_csv('sparcc_species_in',sep='\t',header=True,index=True)
source deactivate python3.6

parallel 'mkdir -p {1}/sparcc/{2}' ::: /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup ::: genus species
'''
source activate /work/workspace/zhurj/bin/miniconda3/envs/py2
parallel --link 'SparCC.py {1}/data/sparcc_genus_in -i 5 --cor_file={1}/sparcc/{2}/cor_sparcc.out.txt > {1}/sparcc/{2}/sparcc.log && \
MakeBootstraps.py {1}/data/sparcc_genus_in -n 100 -t {1}/sparcc/{2}/bootstrap_#.txt -p {1}/sparcc/{2}/pvals/ >> {1}/sparcc/{2}/sparcc.log && \
for n in {0..99}; do SparCC.py {1}/sparcc/{2}/pvals/bootstrap_${n}.txt -i 5 --cor_file={1}/sparcc/{2}/pvals/bootstrap_cor_${n}.txt >> {1}/sparcc/{2}/sparcc.log; done && \
PseudoPvals.py {1}/sparcc/{2}/cor_sparcc.out.txt {1}/sparcc/{2}/pvals/bootstrap_cor_#.txt 100 -o {1}/sparcc/{2}/pvals/pvals.two_sided.txt -t two_sided >> {1}/sparcc/{2}/sparcc.log && \
' ::: /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup ::: genus
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/py2
'''
cd /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/sparcc/genus
SparCC.py sparcc_genus_in -i 5 --cor_file=cor_sparcc.out.txt > sparcc.log
第 2 步，通过自举抽样获得随机数据集
#MakeBootstraps.py -h
MakeBootstraps.py sparcc_genus_in -n 100 -t bootstrap_#.txt -p pvals/ >> sparcc.log
#第 3 步，计算伪 p 值，作为评估相关性显著的依据
#首先通过循环语句批处理，获得各随机数据集中变量的相关矩阵（随机值的相关矩阵）
for n in {0..99}; do SparCC.py pvals/bootstrap_${n}.txt -i 5 --cor_file=pvals/bootstrap_cor_${n}.txt >> sparcc.log; done
#通过观测值的相关矩阵中系数（cor0），以及随机值的相关矩阵中系数（corN），考虑 |cor0|>|corN| 的频率，获得伪 p 值（我猜的应该是这样......）
#PseudoPvals.py -h
PseudoPvals.py cor_sparcc.out.txt pvals/bootstrap_cor_#.txt 100 -o pvals/pvals.two_sided.txt -t two_sided >> sparcc.log

R3.6
#观测值的相关矩阵
cor_sparcc <- read.delim('cor_sparcc.out.txt', row.names = 1, sep = '\t', check.names = FALSE)
#伪 p 值矩阵
pvals <- read.delim('pvals/pvals.two_sided.txt', row.names = 1, sep = '\t', check.names = FALSE)
#保留 |相关性|≥0.8 且 p<0.01的值
cor_sparcc[abs(cor_sparcc) < 0.3] <- 0
#pvals[pvals>=0.05] <- -1
pvals[pvals>=0.05] <- 0
pvals[pvals<0.05 & pvals>=0] <- 1
#pvals[pvals==-1] <- 0
 
#筛选后的邻接矩阵
adj <- as.matrix(cor_sparcc) * as.matrix(pvals)
diag(adj) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
write.table(data.frame(adj, check.names = FALSE), 'network.adj.txt', col.names = NA, sep = '\t', quote = FALSE)

##网络格式转换
library(igraph)
 
#输入数据，邻接矩阵
network_adj <- read.delim('network.adj.txt', row.names = 1, sep = '\t', check.names = FALSE)
head(neetwork_adj)[1:6]    #邻接矩阵类型的网络文件
 
#邻接矩阵 -> igraph 的邻接列表，获得含权的无向网络
g <- graph_from_adjacency_matrix(as.matrix(network_adj), mode = 'undirected', weighted = TRUE, diag = FALSE)
g    #igraph 的邻接列表
 
#这种转换模式下，默认的边权重代表了 sparcc 计算的相关性（存在负值）
#由于边权重通常为正值，因此最好取个绝对值，相关性重新复制一列作为记录
E(g)$sparcc <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)
 
#再转为其它类型的网络文件，例如
#再由 igraph 的邻接列表转换回邻接矩阵
adj_matrix <- as.matrix(get.adjacency(g, attr = 'sparcc'))
write.table(data.frame(adj_matrix, check.names = FALSE), 'network.adj_matrix.txt', col.names = NA, sep = '\t', quote = FALSE)
 
#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(g, 'network.graphml', format = 'graphml')
 
#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(g, 'network.gml', format = 'gml')
 
#边列表，也可以直接导入至 gephi 或 cytoscape 等网络可视化软件中进行编辑
edge <- data.frame(as_edgelist(g))
 
edge_list <- data.frame(
    source = edge[[1]],
    target = edge[[2]],
    weight = E(g)$weight,
    sparcc = E(g)$sparcc
)
head(edge_list)
 
write.table(edge_list, 'network.edge_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)
 
#节点属性列表，对应边列表，记录节点属性，例如
node_list <- data.frame(
    nodes_id = V(g)$name,    #节点名称
    degree = degree(g)    #节点度
)
head(node_list)
write.table(node_list, 'network.node_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)
/work/workspace/zhurj/software/gephi-0.9.2/bin/gephi

genus
cd /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data
source activate python3.6
import pandas as pd
df = pd.read_table('lefse_genus_in',sep='\t',header=0,index_col = 0)
ndf = df * 50000
nndf = ndf.astype(int)
nndf.to_csv('sparcc_genus_in',sep='\t',header=True,index=True)

species
parallel 'mkdir -p {1}/sparcc/{2}' ::: /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup ::: genus species
cd /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data
source activate python3.6
import pandas as pd
df = pd.read_table('lefse_species_in',sep='\t',header=0,index_col = 0)
ndf = df * 50000
nndf = ndf.astype(int)
nndf.to_csv('sparcc_species_in',sep='\t',header=True,index=True)
source deactivate python3.6
cp /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data/sparcc_species_in /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/sparcc/species
cd /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/sparcc/species
source activate py2
SparCC.py sparcc_species_in -i 5 --cor_file=cor_sparcc.out.txt > sparcc.log
#第 2 步，通过自举抽样获得随机数据集
#MakeBootstraps.py -h
MakeBootstraps.py sparcc_species_in -n 100 -t bootstrap_#.txt -p pvals/ >> sparcc.log
#第 3 步，计算伪 p 值，作为评估相关性显著的依据
#首先通过循环语句批处理，获得各随机数据集中变量的相关矩阵（随机值的相关矩阵）
for n in {0..99}; do SparCC.py pvals/bootstrap_${n}.txt -i 5 --cor_file=pvals/bootstrap_cor_${n}.txt >> sparcc.log; done
#通过观测值的相关矩阵中系数（cor0），以及随机值的相关矩阵中系数（corN），考虑 |cor0|>|corN| 的频率，获得伪 p 值（我猜的应该是这样......）
#PseudoPvals.py -h
PseudoPvals.py cor_sparcc.out.txt pvals/bootstrap_cor_#.txt 100 -o pvals/pvals.two_sided.txt -t two_sided >> sparcc.log

R3.6
#观测值的相关矩阵
library(igraph)
cor_sparcc <- read.delim('cor_sparcc.out.txt', row.names = 1, sep = '\t', check.names = FALSE)
#伪 p 值矩阵
pvals <- read.delim('pvals/pvals.two_sided.txt', row.names = 1, sep = '\t', check.names = FALSE)
#保留 |相关性|≥0.8 且 p<0.01的值
cor_sparcc[abs(cor_sparcc) <= 0.7] <- 0
#pvals[pvals>=0.05] <- -1
pvals[pvals>=0.05] <- 0
pvals[pvals<0.05 & pvals>=0] <- 1
#pvals[pvals==-1] <- 0
 
#筛选后的邻接矩阵
adj <- as.matrix(cor_sparcc) * as.matrix(pvals)
diag(adj) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
write.table(data.frame(adj, check.names = FALSE), 'network.adj.txt', col.names = NA, sep = '\t', quote = FALSE)
##网络格式转换
#输入数据，邻接矩阵
network_adj <- read.delim('network.adj.txt', row.names = 1, sep = '\t', check.names = FALSE)
head(neetwork_adj)[1:6]    #邻接矩阵类型的网络文件
#邻接矩阵 -> igraph 的邻接列表，获得含权的无向网络
g <- graph_from_adjacency_matrix(as.matrix(network_adj), mode = 'undirected', weighted = TRUE, diag = FALSE)
#g    #igraph 的邻接列表
#这种转换模式下，默认的边权重代表了 sparcc 计算的相关性（存在负值）
#由于边权重通常为正值，因此最好取个绝对值，相关性重新复制一列作为记录
E(g)$sparcc <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)

#再转为其它类型的网络文件，例如
#再由 igraph 的邻接列表转换回邻接矩阵
adj_matrix <- as.matrix(get.adjacency(g, attr = 'sparcc'))
write.table(data.frame(adj_matrix, check.names = FALSE), 'network.adj_matrix.txt', col.names = NA, sep = '\t', quote = FALSE)
#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(g, 'network.graphml', format = 'graphml')
#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(g, 'network.gml', format = 'gml')
#边列表，也可以直接导入至 gephi 或 cytoscape 等网络可视化软件中进行编辑
edge <- data.frame(as_edgelist(g))
 
edge_list <- data.frame(
    source = edge[[1]],
    target = edge[[2]],
    weight = E(g)$weight,
    sparcc = E(g)$sparcc
)
#head(edge_list)
write.table(edge_list, 'network.edge_list.txt', sep = '\t', row.names = FALSE, quote = FALSE) 
#节点属性列表，对应边列表，记录节点属性，例如
node_list <- data.frame(
    nodes_id = V(g)$name,    #节点名称
    degree = degree(g)    #节点度
)
#head(node_list)
write.table(node_list, 'network.node_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)
/work/workspace/zhurj/software/gephi-0.9.2/bin/gephi

测试：
****************************************************************************************
NETWORK
SPIEC-EASI的微生物网络构建示例: https://blog.csdn.net/woodcorpse/article/details/106554486
测试okay 环境R3.6
Co-occurrence网络图在R中的实现: https://blog.csdn.net/woodcorpse/article/details/78737867
zdk123/SpiecEasi: https://github.com/zdk123/SpiecEasi
微生物组统计和可视化——phyloseq入门: https://blog.csdn.net/woodcorpse/article/details/106554382
****************************************************************************************
library(SpiecEasi)
library(igraph)

#data(amgut1.filt)
#class(amgut1.filt)
#amgut1.filt[1:6,1:6]

df <- read.csv('merge_species_count', header=TRUE, sep='\t',row.names=1)
dm = t(data.matrix(df))
#method 提供了两种估计方法，即 glasso 和 mb，分别展示（二者间可能会差别较大）
se.gl.amgut <- spiec.easi(data = dm, method = 'glasso', lambda.min.ratio = 0.01,
    nlambda = 20, pulsar.params = list(rep.num = 50))
se.mb.amgut <- spiec.easi(dm, method = 'mb', lambda.min.ratio = 0.01,
    nlambda = 20, pulsar.params = list(rep.num=50))

adjacency_unweight <- data.frame(as.matrix(se.gl.amgut$refit$stars))
rownames(adjacency_unweight) <- colnames(dm)
colnames(adjacency_unweight) <- colnames(dm)

ig.gl <- adj2igraph(getRefit(se.gl.amgut), vertex.attr = list(label = colnames(dm)))
ig.mb <- adj2igraph(getRefit(se.mb.amgut), vertex.attr = list(label = colnames(dm)))

vsize <- rowMeans(clr(dm, 1)) + 6
am.coord <- layout.fruchterman.reingold(ig.mb)

par(mfrow = c(1, 2))
plot(ig.gl, layout = am.coord, vertex.size = vsize, vertex.label = NA, main = 'glasso')
plot(ig.mb, layout = am.coord, vertex.size = vsize, vertex.label = NA, main = 'MB')


library(SpiecEasi)
library(igraph)
sparcc.amgut <- sparcc(dm)
sparcc.graph <- abs(sparcc.amgut$Cor) >= 0.8
diag(sparcc.graph) <- 0
library(Matrix)
sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
ig.sparcc <- adj2igraph(sparcc.graph)


## set size of vertex proportional to clr-mean
vsize    <- rowMeans(clr(dm, 1))+6
am.coord <- layout.fruchterman.reingold(ig.mb)
par(mfrow=c(1,3))
#plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB")
#plot(ig.gl, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc")

#估计边的权重
library(Matrix)
 
secor <- cov2cor(getOptCov(se.gl.amgut))
sebeta <- symBeta(getOptBeta(se.mb.amgut), mode = 'maxabs')
elist.gl <- summary(triu(secor*getRefit(se.gl.amgut), k = 1))
elist.mb <- summary(sebeta)
 
head(elist.gl)
head(elist.mb)

#将权重合并到网络中，以 glasso 的结果为例
elist.gl <- data.frame(elist.gl)
names(elist.gl) <- c('source', 'target', 'weight')
elist.gl <- elist.gl[order(elist.gl$source, elist.gl$target), ]
E(ig.gl)$weight <- elist.gl$weight
 
#输出带权重的邻接矩阵，不再是 0-1 关系，而替换为具体数值（可表示两变量关联程度）
adjacency_weight <- as.matrix(get.adjacency(ig.gl, attr = 'weight'))
rownames(adjacency_weight) <- colnames(dm)
colnames(adjacency_weight) <- colnames(dm)
write.table(data.frame(adjacency_weight), 'adjacency_weight.glasso.txt', col.names = NA, sep = '\t', quote = FALSE)
 
#以及含权重的边列表
edge.gl <- data.frame(as_edgelist(ig.gl))
edge.gl <- data.frame(source = edge.gl[[1]], target = edge.gl[[2]], weight = E(ig.gl)$weight)
write.table(edge.gl, 'weight.glasso_edge.txt', sep = '\t', row.names = FALSE, quote = FALSE)
 
#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(ig.gl, 'weight.glasso.graphml', format = 'graphml')
 
#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(ig.gl, 'weight.glasso.gml', format = 'gml')

dim(matrix)
nrow(m)
ncol(m)
colnames(m)
rownames(m)

#------------------------------------------------------------------------------------------------------------------------------------
补测数据分析
parallel --link 'ln -s {1}/{2}/{2}_1.fq.gz {4}/{3}/{3}_1.fq.gz && ln -s {1}/{2}/{2}_2.fq.gz {4}/{3}/{3}_2.fq.gz  ' ::: /work/rawdata/run/guangzhou/novogene/2020/10/20201021/run00070/X101SC20073487-Z01-J001/1.rawdata ::: DCO-7_FDSW202266364-1r GCO-6_FDSW202266371-1r ::: MNA00245 MNA00252 ::: /work/workspace/zhurj/project/1_metadata/mouse22/rawdata

2. data filter
/work/workspace/zhurj/project/1_metadata/mouse22/rawdata/MNA00245/MNA00245_1.fq.gz	/work/workspace/zhurj/project/1_metadata/mouse22/rawdata/MNA00245/MNA00245_2.fq.gz
/work/workspace/zhurj/project/1_metadata/mouse22/rawdata/MNA00252/MNA00252_1.fq.gz	/work/workspace/zhurj/project/1_metadata/mouse22/rawdata/MNA00252/MNA00252_2.fq.gz
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
perl /work/workspace/zhurj/script/5_wgspro/1_filter/datafilter.pl -i /work/workspace/zhurj/project/1_metadata/mouse22/com/input/df_in -o /work/workspace/zhurj/project/1_metadata/mouse22/com
parallel --link -j 2 'bwa mem {4} {1}/{2}.clean.1.fq.gz {1}/{2}.clean.2.fq.gz -t 20 -o {3}/{2}/{2}.sam && \
samtools view -@ 20 -bS {3}/{2}/{2}.sam -o {3}/{2}/{2}.bam && \
samtools view -b -F 2 -@ 20 {3}/{2}/{2}.bam -o {3}/{2}/{2}_PEunmapped.bam && \
samtools sort -n -@ 20 {3}/{2}/{2}_PEunmapped.bam -o {3}/{2}/{2}_PEunmapped_sorted.bam && \
bedtools bamtofastq -i {3}/{2}/{2}_PEunmapped_sorted.bam -fq {3}/{2}/{2}_r1.fastq -fq2 {3}/{2}/{2}_r2.fastq && \
rm {3}/{2}/*.sam {3}/{2}/*.bam -f ' \
::: /work/workspace/zhurj/project/1_metadata/mouse22/com/clean ::: MNA00252 MNA00245 \
::: /work/workspace/zhurj/project/1_metadata/mouse22/com/bwarmhost ::: /work/workspace/zhurj/reference/mouse/GRCm38/grcm38.fa
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6


source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel -j 1 --link ' kneaddata  -i {1}/clean/{2}.1.fq.gz -i {1}/clean/{2}.2.fq.gz --output /work/workspace/zhurj/project/1_metadata/mouse22/com/kneaddata/{2} -db /work/workspace/zhurj/reference/mouse/kneaddata  --bypass-trim --run-trf -t 20 -p 20 --output-prefix  {3} --reorder  --bowtie2-options "--very-sensitive --dovetail" --remove-intermediate-output'  ::: /work/workspace/zhurj/project/1_metadata/mouse22/com ::: MNA00252 MNA00245
parallel -j 1 --link ' ln -s /work/workspace/zhurj/project/1_metadata/mouse22/com/kneaddata/{1}/{1}.log /work/workspace/zhurj/project/1_metadata/mouse22/com/kneaddata/log/{1}.log ' ::: MNA00252 MNA00245
parallel -j 2 --link ' cat {1}/{2}/{2}.repeats.removed.1.fastq {1}/{2}/{2}.repeats.removed.2.fastq > {1}/{2}/{2}.fastq ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/com/kneaddata ::: MNA00252 MNA00245
kneaddata_read_count_table --input /work/workspace/zhurj/project/1_metadata/mouse22/com/kneaddata/log --output /work/workspace/zhurj/project/1_metadata/mouse22/com/kneaddata/merge/kneaddata_read_counts.txt
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

merge kneaddata
parallel 'mkdir -p {1}/{2}/{3}' ::: /work/workspace/zhurj/project/1_metadata/mouse22/com/merge ::: metaphlan3 humann3 ::: MNA00245 MNA00252
parallel 'mkdir -p {1}/{2}/{3}' ::: /work/workspace/zhurj/project/1_metadata/mouse22/com/merge/metaphlan3 ::: phylum genus species all
parallel 'mkdir -p {1}/{2}/{3}' ::: /work/workspace/zhurj/project/1_metadata/mouse22/com/merge ::: humann3 ::: MNA00245 MNA00252 merge

parallel 'cat {1}/com/kneaddata/{2}/{2}.repeats.removed.{3}.fastq {1}/kneaddata/{2}/{2}.repeats.removed.{3}.fastq > {1}/com/merge/kneaddata/{2}/{2}.repeats.removed.{3}.fastq \
' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: MNA00245 MNA00252 ::: 1 2
parallel --link 'cat {1}/com/kneaddata/{2}/{2}.fastq {1}/kneaddata/{2}/{2}.fastq > {1}/com/merge/kneaddata/{2}/{2}.fastq \
' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: MNA00245 MNA00252

source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel -j 1 --link 'metaphlan {1}/kneaddata/{2}/{2}.repeats.removed.1.fastq,{1}/kneaddata/{2}/{2}.repeats.removed.2.fastq --bowtie2out {1}/metaphlan3/{2}/metagenome.bowtie2.bz2 \
--input_type fastq --nproc 20 --sample_id_key {2} -o {1}/metaphlan3/{2}/all.txt && \
metaphlan {1}/metaphlan3/{2}/metagenome.bowtie2.bz2 --input_type bowtie2out --tax_lev "p" --nproc 20 --sample_id_key {2} -o {1}/metaphlan3/{2}/phylum.txt && \
metaphlan {1}/metaphlan3/{2}/metagenome.bowtie2.bz2 --input_type bowtie2out --tax_lev "g" --nproc 20 --sample_id_key {2} -o {1}/metaphlan3/{2}/genus.txt && \
metaphlan {1}/metaphlan3/{2}/metagenome.bowtie2.bz2 --input_type bowtie2out --tax_lev "s" --nproc 12 --sample_id_key {2} -o {1}/metaphlan3/{2}/species.txt && \
ln -s {1}/metaphlan3/{2}/all.txt {1}/metaphlan3/all/all.txt && \
ln -s {1}/metaphlan3/{2}/phylum.txt {1}/metaphlan3/phylum/phylum.txt && \
ln -s {1}/metaphlan3/{2}/genus.txt {1}/metaphlan3/genus/genus.txt && \
ln -s {1}/metaphlan3/{2}/species.txt {1}/metaphlan3/species/species.txt \
' ::: /work/workspace/zhurj/project/1_metadata/mouse22/com/merge ::: MNA00245 MNA00252

parallel -j 1 --link 'humann --threads 36 --input {1}/kneaddata/{2}/{2}.fastq --output {1}/humann3/{2} --output-basename {2}_metacyc && \
rm {1}/humann3/{2}/{2}_metacyc_humann_temp -rf &&\
humann_regroup_table -i {1}/humann3/{2}/{2}_metacyc_genefamilies.tsv -g uniref90_ko --output  {1}/humann3/{2}/{2}_uniref90_ko.tsv && \
humann_regroup_table -i {1}/humann3/{2}/{2}_metacyc_genefamilies.tsv -g uniref90_go --output  {1}/humann3/{2}/{2}_uniref90_go.tsv && \
humann_regroup_table -i {1}/humann3/{2}/{2}_metacyc_genefamilies.tsv -g uniref90_level4ec --output  {1}/humann3/{2}/{2}_uniref90_level4ec.tsv && \
humann_regroup_table -i {1}/humann3/{2}/{2}_metacyc_genefamilies.tsv -g uniref90_pfam --output  {1}/humann3/{2}/{2}_uniref90_pfam.tsv && \
humann_regroup_table -i {1}/humann3/{2}/{2}_metacyc_genefamilies.tsv -g uniref90_eggnog --output  {1}/humann3/{2}/{2}_uniref90_eggnog.tsv \
' ::: /work/workspace/zhurj/project/1_metadata/mouse22/com/merge ::: MNA00245 MNA00252

parallel -j 1 --link 'humann -i {1}/{2}/{2}_uniref90_ko.tsv --pathways-database {3} --output {1}/{2} --threads 20 --output-basename {2}_kegg \
' ::: /work/workspace/zhurj/project/1_metadata/mouse22/com/merge/humann3 ::: MNA00245 MNA00252 ::: /work/workspace/zhurj/database/humann/current/kegg/keggc
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann


merge bwarmdata
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel 'mkdir -p {1}/merge/{2}/{3}' ::: /work/workspace/zhurj/project/1_metadata/mouse22/com ::: bwarmhost kneaddata ::: MNA00245 MNA00252
parallel 'mkdir -p {1}/merge/{2}/{3}' ::: /work/workspace/zhurj/project/1_metadata/mouse22/com ::: spades metabat ::: MNA00245 MNA00252
parallel 'mkdir -p {1}/merge/{2}/tmp' ::: /work/workspace/zhurj/project/1_metadata/mouse22/com ::: spades
parallel 'mkdir -p {1}/merge/{2}/{3}' ::: /work/workspace/zhurj/project/1_metadata/mouse22/com ::: bwa ::: MNA00245 MNA00252
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link 'cat {1}/com/bwarmhost/{2}/{2}_r1.fastq {1}/bwarmhost/{2}/{2}_r1.fastq > {1}/com/merge/bwarmhost/{2}/{2}_r1.fastq && \
cat {1}/com/bwarmhost/{2}/{2}_r2.fastq {1}/bwarmhost/{2}/{2}_r2.fastq > {1}/com/merge/bwarmhost/{2}/{2}_r2.fastq \
' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: MNA00245 MNA00252

parallel --link -j 1 'python /work/workspace/moonis/bin/genoassemble/soft/SPAdes-3.13.1-Linux/bin/spades.py --meta \
-1 {1}/bwarmhost/{2}/{2}_r1.fastq -2 {1}/bwarmhost/{2}/{2}_r2.fastq -o {1}/spades/{2} -t 36 && \
rm {1}/spades/tmp/* -rf && \
mv {1}/spades/{2}/contigs.fasta  {1}/spades/tmp && \
mv {1}/spades/{2}/scaffolds.fasta  {1}/spades/tmp && \
rm {1}/spades/{2}/* -rf && \
mv {1}/spades/tmp/* {1}/spades/{2} ' \
::: /work/workspace/zhurj/project/1_metadata/mouse22/com/merge :::  MNA00245 MNA00252

parallel -j 1 --link 'bwa index {1}/spades/{2}/contigs.fasta && \
bwa mem -t {3} {1}/spades/{2}/contigs.fasta {1}/bwarmhost/{2}/{2}_r1.fastq {1}/bwarmhost/{2}/{2}_r2.fastq > {1}/bwa/{2}/{2}.sam && \
samtools view -@ {3} -S -b {1}/bwa/{2}/{2}.sam -o {1}/bwa/{2}/{2}.bam && \
samtools sort -@ {3} {1}/bwa/{2}/{2}.bam -o {1}/bwa/{2}/{2}.sorted.bam && \
jgi_summarize_bam_contig_depths --outputDepth {1}/metabat/{2}/depth.txt  {1}/bwa/{2}/{2}.sorted.bam && \
metabat2 -i {1}/spades/{2}/contigs.fasta -a {1}/metabat/{2}/depth.txt -o {1}/metabat/{2}/{2} -m 2000 && \
rm {1}/bwa/{2}/{2}.sam {1}/bwa/{2}/{2}.bam -f ' \
::: /work/workspace/zhurj/project/1_metadata/mouse22/com/merge ::: MNA00245 MNA00252 ::: 36
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

cd /work/workspace/zhurj/project/1_metadata/p/p
time srun -o com2.out -e com2.err -N 1 -c 20 -p slurm256 -w mnclient01 bash com2.sh &
time srun -o dehost2.out -e dehost2.err -N 1 -c 20 -p slurm256 -w mnclient01 bash dehost2.sh &
time srun -o kneaddata2.out -e kneaddata2.err -N 1 -c 20 -p slurm256 -w mnclient02 bash kneaddata2.sh &
time srun -o spadesmetabatcom2.out -e spadesmetabatcom2.err -N 1 -c 20 -p slurm256 -w mnclient01 bash spadesmetabatcom2.sh &
time srun -o metaphlanhumanncom2.out -e metaphlanhumanncom2.err -N 1 -c 20 -p slurm256 -w mnclient02 bash metaphlanhumanncom2.sh &
#------------------------------------------------------------------------------------------------------------------------------------
MEGAN 物种鉴定
parallel --link 'mkdir -p {1}/megan/{2} ' ::: /work/workspace/zhurj/project/1_metadata/mouse22 :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22
parallel --link 'mkdir -p {1}/diamond/{2} ' ::: /work/workspace/zhurj/project/1_metadata/mouse22 :::: /work/workspace/zhurj/project/1_metadata/mouse22/input/mnid22
parallel --link 'mkdir -p {1}/diamond/tmp ' ::: /work/workspace/zhurj/project/1_metadata/mouse22
ln -s /work/database/ncbi/nr/2019.Oct12/microbe/sub_microbe.fa.gz /work/workspace/zhurj/reference/NCBI/nr/db/microbe/microbe.gz
ln -s /work/database/ncbi/nr/2019.Oct12/microbe/sub_microbe.fa.tab /work/workspace/zhurj/reference/NCBI/nr/db/microbe/microbe.tab
unpigz -k -p 20 /work/workspace/zhurj/reference/NCBI/nr/db/microbe/microbe.gz
diamond makedb --in /work/workspace/zhurj/reference/NCBI/nr/db/nr -d /work/workspace/zhurj/reference/NCBI/nr/db/nr -p 20

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
unpigz -k -p 20 /work/workspace/zhurj/reference/NCBI/nr/db/nr.gz
diamond makedb --in /work/workspace/zhurj/reference/NCBI/nr/db/nr -d /work/workspace/zhurj/reference/NCBI/nr/db/nr -p 20
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

megantest.sh
fastq to taxonomy
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link -j 1 'diamond blastx -c 1 --db {1} -t {2}/megan/{3}/tmp -p 20 -q {2}/bwarmhost/{3}/{3}_r1.fastq --daa {2}/diamond/{3}/{3}.1.daa && \
diamond blastx -c 1 --db {1} -t {2}/megan/{3}/tmp -p 20 -q {2}/bwarmhost/{3}/{3}_r2.fastq --daa {2}/diamond/{3}/{3}.2.daa \
'::: /work/workspace/zhurj/database/nr/diamond/nr_microbe.dmnd ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: MNA00233

#转化双端daa文件为MEGAN特有额rma文件。
parallel --link -j 1 ' {4}/daa2rma -i  {2}/diamond/{3}/{3}.1.daa  {2}/diamond/{3}/{3}.2.daa \
--paired -ms 50 -me 0.01 -top 50  -mdb {1}  -o {2}/megan/{3}/{3}.rma && \
{4}/rma2info -i {2}/megan/{3}/{3}.rma -r2c Taxonomy -n true --paths true --ranks true --list true -v >  {2}/megan/{3}/{3}_taxonomy.txt && \
{4}/rma2info -i {2}/megan/{3}/{3}.rma -r2c SEED -n true --paths true --ranks true --list true --listMore true -v -v >  {2}/megan/{3}/{3}_seed.txt && \
{4}/rma2info -i {2}/megan/{3}/{3}.rma -r2c EGGNOG -n true --paths true --ranks true --list true --listMore true -v -v >  {2}/megan/{3}/{3}_eggnog.txt && \
{4}/rma2info -i {2}/megan/{3}/{3}.rma -r2c INTERPRO2GO -n true --paths true --ranks true --list true --listMore true -v -v >  {2}/megan/{3}/{3}_interpro2go.txt \
' ::: /work/workspace/zhurj/database/megan/megan-map-Jul2020-2.db ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: MNA00233 ::: /work/workspace/zhurj/instal/megan/tools
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

faa to taxonomy
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link -j 1 'diamond blastp -c 1 --db {1} -t {2}/diamond/tmp -p 20 -q {2}/prokka/faa/{3}.faa --daa {2}/diamond/{3}/{3}.daa \
' ::: /work/workspace/zhurj/database/nr/diamond/nr_microbe.dmnd ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: MNA00233

parallel --link -j 1 ' {4}/daa2rma -i  {2}/diamond/{3}/{3}.daa \
--paired -ms 50 -me 0.01 -top 50  -mdb {1}  -o {2}/megan/{3}/{3}.rma && \
{4}/rma2info -i {2}/megan/{3}/{3}.rma -r2c Taxonomy -n true --paths true --ranks true --list true -v >  {2}/megan/{3}/{3}_taxonomy.txt && \
{4}/rma2info -i {2}/megan/{3}/{3}.rma -r2c SEED -n true --paths true --ranks true --list true --listMore true -v -v >  {2}/megan/{3}/{3}_seed.txt && \
{4}/rma2info -i {2}/megan/{3}/{3}.rma -r2c EGGNOG -n true --paths true --ranks true --list true --listMore true -v -v >  {2}/megan/{3}/{3}_eggnog.txt && \
{4}/rma2info -i {2}/megan/{3}/{3}.rma -r2c INTERPRO2GO -n true --paths true --ranks true --list true --listMore true -v -v >  {2}/megan/{3}/{3}_interpro2go.txt \
' ::: /work/workspace/zhurj/database/megan/megan-map-Jul2020-2.db ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: MNA00233 ::: /work/workspace/zhurj/instal/megan/tools

source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

time srun -o megantest.out -e megantest.err -N 1 -c 20 -p slurm256 -w mnclient02 bash megantest.sh &



#-------------------------------------------------------------------------------------
SGB 分析结果整理
PCA 分析
source activate R3.6
/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/pic
library('ggbiplot')
data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/species_pca_in",sep='\t',header=TRUE,row.names=1)
data.pca <- prcomp(data, center = TRUE,scale. = TRUE)
#summary(data.pca)
#str(data.pca)
dtgroup <- c(rep("HFD",8),rep("MNO863",8))
#ggbiplot(data.pca,ellipse=TRUE,  labels=rownames(data), groups=dtgroup)
p <- ggbiplot(data.pca,ellipse=TRUE,obs.scale = 1, var.scale = 1,var.axes=FALSE, groups=dtgroup) +
scale_color_discrete(name = '') +
theme(legend.direction = 'horizontal', legend.position = 'top')
ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/pic/pca_species160_20201026.jpg", p)

data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/genus_pca_in",sep='\t',header=TRUE,row.names=1)
data.pca <- prcomp(data, center = TRUE,scale. = TRUE)
dtgroup <- c(rep("HFD",8),rep("MNO863",8))
#ggbiplot(data.pca,ellipse=TRUE,  labels=rownames(data), groups=dtgroup)
p <- ggbiplot(data.pca,ellipse=TRUE,obs.scale = 1, var.scale = 1,var.axes=FALSE, groups=dtgroup) +
scale_color_discrete(name = '') +
theme(legend.direction = 'horizontal', legend.position = 'top')
ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/pic/pca_genus48_20201026.jpg", p)

物种组成bar图

python3.6
cd /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup
cd /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup
cd /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup && python /work/workspace/zhurj/lib/python/script/lefse2boxplot.py data/bar_genus_in  input/meta data/genus_format_in
cd /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup && python /work/workspace/zhurj/lib/python/script/lefse2boxplot.py data/bar_phylum_in  input/meta data/phylum_format_in
cd /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup && python /work/workspace/zhurj/lib/python/script/lefse2boxplot.py data/bar_species_in  input/meta data/species_format_in

R3.6
library('ggplot2')
library('RColorBrewer')
library('ggsci')
D <- read.table(file="/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data/phylum_format_in", head=T, sep='\t')
png("/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/pic/bar_phylum_20201026.png", width=15, height=8, units = "in", res=100);
a<-D[grep('p__Firmicutes',D$Species),]
b<-a[order(a$Abundance,decreasing = T),]
#D$Sample<-factor(D$SDample,levels = b$Sample)
D$Sample<-as.character(D$Sample)
D$Group<-factor(D$Group,levels = c("HFD","MNO863"))
mypal = c("#E41A1C","#4DAF4A","#377EB8","#A65628","#FF7F00","#FFFF33")
mycolor = colorRampPalette(mypal)
ggplot(D) +
  geom_bar(aes(x=Sample, y=Abundance, group=Sample, fill=factor(Species)),stat="identity", width=1, color= "black") +
  #scale_x_discrete(breaks=NULL,expand=c(0.01,0)) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_continuous(expand = c(0.00,0.00)) +
  facet_grid(~Group, scales = "free", space = "free_x", switch = "x") +
  theme_classic() +
  theme(
      axis.title=element_text( family="Helvetica", face="bold", colour="black",size=16),
      axis.text=element_text(family="Helvetica",colour="black",size=12),
      axis.text.x=element_text(angle=90,hjust=1, vjust=0.5,size=10),
      strip.text=element_text(family="Helvetica",colour="black",angle=0,hjust=0.5,vjust=0.5,size=14),
      #strip.text.x=element_blank(),
      #strip.background=element_blank(),
      #legend.position=c(0.90,0.80),
      legend.title=element_text(family="Helvetica",colour="black",size=16),
      legend.text=element_text(family="Helvetica", colour="black",size=14)
  ) +
  guides(fill = guide_legend("Phylum", ncol=1)) +
  scale_fill_manual(values = mycolor(9), breaks = rev(as.vector(unique(D$Species)))) +
  xlab('') +
  ylab('Relative Abundance')
  dev.off()

R3.6
library('ggplot2')
library('RColorBrewer')
library('ggsci')
D <- read.table(file="/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data/genus_format_in", head=T, sep='\t')
png("/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/pic/bar_genus_20201026.png", width=15, height=8, units = "in", res=100);
a<-D[grep('g__Pseudoflavonifractor',D$Species),]
b<-a[order(a$Abundance,decreasing = T),]
#D$Sample<-factor(D$SDample,levels = b$Sample)
D$Sample<-as.character(D$Sample)
D$Group<-factor(D$Group,levels = c("HFD","MNO863"))
mypal = c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black",
"gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F",
"gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1",
"steelblue4", "darkturquoise", "green1", "yellow4", "yellow3",
"darkorange4", "brown")
mycolor = colorRampPalette(mypal)
ggplot(D) +
  geom_bar(aes(x=Sample, y=Abundance, group=Sample, fill=factor(Species)),stat="identity", width=1, color= "black") +
  #scale_x_discrete(breaks=NULL,expand=c(0.01,0)) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_continuous(expand = c(0.00,0.00)) +
  facet_grid(~Group, scales = "free", space = "free_x", switch = "x") +
  theme_classic() +
  theme(
      axis.title=element_text( family="Helvetica", face="bold", colour="black",size=16),
      axis.text=element_text(family="Helvetica",colour="black",size=12),
      axis.text.x=element_text(angle=90,hjust=1, vjust=0.5,size=10),
      strip.text=element_text(family="Helvetica",colour="black",angle=0,hjust=0.5,vjust=0.5,size=14),
      #strip.text.x=element_blank(),
      #strip.background=element_blank(),
      #legend.position=c(0.90,0.80),
      legend.title=element_text(family="Helvetica",colour="black",size=16),
      legend.text=element_text(family="Helvetica", colour="black",size=14)
  ) +
  guides(fill = guide_legend("Phylum", ncol=1)) +
  scale_fill_manual(values = mycolor(21), breaks = rev(as.vector(unique(D$Species)))) +
  xlab('') +
  ylab('Relative Abundance')
  dev.off()

library('ggplot2')
library('RColorBrewer')
library('ggsci')
D <- read.table(file="/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data/species_format_in", head=T, sep='\t')
png("/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/pic/bar_species_20201026.png", width=15, height=12, units = "in", res=100);
a<-D[grep('s__Pseudoflavonifractor_sp',D$Species),]
b<-a[order(a$Abundance,decreasing = T),]
#D$Sample<-factor(D$SDample,levels = b$Sample)
D$Sample<-as.character(D$Sample)
D$Group<-factor(D$Group,levels = c("HFD","MNO863"))
mypal = c("#0048BA","#B0BF1A","#7CB9E8","#C0E8D5","#B284BE","#72A0C1","#EDEAE0","#F0F8FF","#C46210","#EFDECD",
"#E52B50","#9F2B68","#F19CBB","#AB274F","#D3212D","#3B7A57","#FFBF00","#FF7E00","#9966CC","#A4C639","#CD9575",
"#665D1E","#915C83","#841B2D","#FAEBD7","#008000","#8DB600","#FBCEB1","#00FFFF","#7FFFD4","#D0FF14","#4B5320",
"#8F9779","#E9D66B","#B2BEB5","#87A96B","#FF9966","#A52A2A","#FDEE00","#568203")
mycolor = colorRampPalette(mypal)
ggplot(D) + 
  geom_bar(aes(x=Sample, y=Abundance, group=Sample, fill=factor(Species)),stat="identity", width=1, color= "black") +
  #scale_x_discrete(breaks=NULL,expand=c(0.01,0)) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_continuous(expand = c(0.00,0.00)) +
  facet_grid(~Group, scales = "free", space = "free_x", switch = "x") +
  theme_classic() +
  theme(
      axis.title=element_text( family="Helvetica", face="bold", colour="black",size=16),
      axis.text=element_text(family="Helvetica",colour="black",size=12),
      axis.text.x=element_text(angle=90,hjust=1, vjust=0.5,size=10),
      strip.text=element_text(family="Helvetica",colour="black",angle=0,hjust=0.5,vjust=0.5,size=14),
      #strip.text.x=element_blank(),
      #strip.background=element_blank(),
      #legend.position=c(0.90,0.80),
      legend.title=element_text(family="Helvetica",colour="black",size=16),
      legend.text=element_text(family="Helvetica", colour="black",size=14)
  ) +
  guides(fill = guide_legend("Species", ncol=1)) +
  scale_fill_manual(values = mycolor(40), breaks = rev(as.vector(unique(D$Species)))) +
  xlab('') +
  ylab('Relative Abundance')
  dev.off()

# lefse 差异物种
parallel --link 'mkdir -p {1}/{2} ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/lefse/SGB/SGB_dg_adj ::: phylum genus species
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link 'python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py {1}/SGB/GDgroup/data/lefse_{2}_in {1}/SGB/GDgroup/input/meta {1}/lefse/SGB/SGB_dg_adj/{2}/lefse_in' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: phylum genus species 
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
# genus difference
sgb_lefse.sh
'''
source activate /work/workspace/zhurj/bin/miniconda3/envs/py2
parallel --link -j 1 'lefse-format_input.py {1}/lefse/SGB/SGB_dg_adj/{2}/lefse_in {1}/lefse/SGB/SGB_dg_adj/{2}/lefse_format.in -c 1 -u 2 -f r -o 1000000 && \
run_lefse.py {1}/lefse/SGB/SGB_dg_adj/{2}/lefse_format.in {1}/lefse/SGB/SGB_dg_adj/{2}/lefse_run.res -a 0.01 -w 0.01 -l 2 -y 1 ' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: phylum genus species
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/py2

source activate python3.6
parallel --link 'python /work/workspace/zhurj/lib/python/script/lefse_fitler.py {1}/lefse/SGB/SGB_dg_adj/{2}/lefse_run.res {1}/lefse/SGB/SGB_dg_adj/{2}/filter_lefse_run.res ' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: phylum genus species
source deactivate python3.6

source activate /work/workspace/zhurj/bin/miniconda3/envs/py2
parallel --link -j 1 'lefse-format_input.py {1}/lefse/SGB/SGB_dg_adj/{2}/lefse_in {1}/lefse/SGB/SGB_dg_adj/{2}/lefse_format.in -c 1 -u 2 -f r -o 1000000 && \
run_lefse.py {1}/lefse/SGB/SGB_dg_adj/{2}/lefse_format.in {1}/lefse/SGB/SGB_dg_adj/{2}/lefse_run_0.05.res -a 0.05 -w 0.05 -l 2 -y 1 ' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: phylum genus species
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/py2

source activate python3.6
parallel --link 'python /work/workspace/zhurj/lib/python/script/lefse_fitler.py {1}/lefse/SGB/SGB_dg_adj/{2}/lefse_run_0.05.res {1}/lefse/SGB/SGB_dg_adj/{2}/filter_lefse_rus_0.05.res ' ::: /work/workspace/zhurj/project/1_metadata/mouse22 ::: phylum genus species
source deactivate python3.6

R3.6
phylum
library(ggplot2)
library(dplyr)
df <- read.table("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/SGB/SGB_dg_adj/input/pic_phylum_in",header = T,sep="\t")
df <- df %>%
  mutate(
    kegg = factor(kegg, levels = kegg[order(value, decreasing = FALSE)]),
    label_y = ifelse(value < 0, 0.2, -0.2),
    label_hjust = ifelse(value < 0, 0, 1)
  )

my_plot <- ggplot(df, aes(x = kegg, y = value, fill = class)) +
  geom_bar(stat = "identity") +
  #geom_bar(stat = "identity", col = "black") +
  geom_text(aes(y = label_y, label = kegg, hjust = label_hjust)) +
  coord_flip() +
  scale_fill_manual(values = c(HFD = "#00AFBB", MNO863 = "#FC4E07")) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.direction='vertical',
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey80", linetype = "dashed"),
        panel.grid.minor.x = element_blank()) +
  scale_y_continuous(expression(log[10](italic("LDA score"))),
                     breaks = -5:5, limits = c(-5, 5))

ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/SGB/SGB_dg_adj/pic/hbar_phylum_20201026.jpg",my_plot,width = 6,height = 2)

# genus
library(ggplot2)
library(dplyr)

df <- read.table("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/SGB/SGB_dg_adj/input/pic_genus_in",header = T,sep="\t")
df <- df %>%
  mutate(
    kegg = factor(kegg, levels = kegg[order(value, decreasing = FALSE)]),
    label_y = ifelse(value < 0, 0.2, -0.2),
    label_hjust = ifelse(value < 0, 0, 1)
  )

my_plot <- ggplot(df, aes(x = kegg, y = value, fill = class)) +
  geom_bar(stat = "identity") +
  #geom_bar(stat = "identity", col = "black") +
  geom_text(aes(y = label_y, label = kegg, hjust = label_hjust)) +
  coord_flip() +
  scale_fill_manual(values = c(HFD = "#00AFBB", MNO863 = "#FC4E07")) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.direction='vertical',
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey80", linetype = "dashed"),
        panel.grid.minor.x = element_blank()) +
  scale_y_continuous(expression(log[10](italic("LDA score"))),
                     breaks = -5:5, limits = c(-5, 5))

ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/SGB/SGB_dg_adj/pic/hbar_genus_20201026.jpg",my_plot,width = 6,height = 8)

# species
library(ggplot2)
library(dplyr)

df <- read.table("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/SGB/SGB_dg_adj/input/pic_species_in",header = T,sep="\t")
df <- df %>%
  mutate(
    kegg = factor(kegg, levels = kegg[order(value, decreasing = FALSE)]),
    label_y = ifelse(value < 0, 0.2, -0.2),
    label_hjust = ifelse(value < 0, 0, 1)
  )

my_plot <- ggplot(df, aes(x = kegg, y = value, fill = class)) +
  geom_bar(stat = "identity") +
  #geom_bar(stat = "identity", col = "black") +
  geom_text(aes(y = label_y, label = kegg, hjust = label_hjust)) +
  coord_flip() +
  scale_fill_manual(values = c(HFD = "#00AFBB", MNO863 = "#FC4E07")) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.direction='vertical',
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey80", linetype = "dashed"),
        panel.grid.minor.x = element_blank()) +
  scale_y_continuous(expression(log[10](italic("LDA score"))),
                     breaks = -8:8, limits = c(-8, 8))

ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/SGB/SGB_dg_adj/pic/hbar_species_20201026.jpg",my_plot,width = 6,height = 14)

================
0.05
phylum
library(ggplot2)
library(dplyr)
df <- read.table("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/SGB/SGB_dg_adj/input/pic_phylum_in_0.05",header = T,sep="\t")
df <- df %>%
  mutate(
    kegg = factor(kegg, levels = kegg[order(value, decreasing = FALSE)]),
    label_y = ifelse(value < 0, 0.2, -0.2),
    label_hjust = ifelse(value < 0, 0, 1)
  )

my_plot <- ggplot(df, aes(x = kegg, y = value, fill = class)) +
  geom_bar(stat = "identity") +
  #geom_bar(stat = "identity", col = "black") +
  geom_text(aes(y = label_y, label = kegg, hjust = label_hjust)) +
  coord_flip() +
  scale_fill_manual(values = c(HFD = "#00AFBB", MNO863 = "#FC4E07")) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.direction='vertical',
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey80", linetype = "dashed"),
        panel.grid.minor.x = element_blank()) +
  scale_y_continuous(expression(log[10](italic("LDA score"))),
                     breaks = -5:5, limits = c(-5, 5))

ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/SGB/SGB_dg_adj/pic/hbar_0.05_phylum_20201027.jpg",my_plot,width = 6,height = 2)

# genus
library(ggplot2)
library(dplyr)

df <- read.table("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/SGB/SGB_dg_adj/input/pic_genus_in_0.05",header = T,sep="\t")
df <- df %>%
  mutate(
    kegg = factor(kegg, levels = kegg[order(value, decreasing = FALSE)]),
    label_y = ifelse(value < 0, 0.2, -0.2),
    label_hjust = ifelse(value < 0, 0, 1)
  )

my_plot <- ggplot(df, aes(x = kegg, y = value, fill = class)) +
  geom_bar(stat = "identity") +
  #geom_bar(stat = "identity", col = "black") +
  geom_text(aes(y = label_y, label = kegg, hjust = label_hjust)) +
  coord_flip() +
  scale_fill_manual(values = c(HFD = "#00AFBB", MNO863 = "#FC4E07")) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.direction='vertical',
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey80", linetype = "dashed"),
        panel.grid.minor.x = element_blank()) +
  scale_y_continuous(expression(log[10](italic("LDA score"))),
                     breaks = -5:5, limits = c(-5, 5))

ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/SGB/SGB_dg_adj/pic/hbar_0.05_genus_20201027.jpg",my_plot,width = 6,height = 8)

# species
library(ggplot2)
library(dplyr)

df <- read.table("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/SGB/SGB_dg_adj/input/pic_species_in_0.05",header = T,sep="\t")
df <- df %>%
  mutate(
    kegg = factor(kegg, levels = kegg[order(value, decreasing = FALSE)]),
    label_y = ifelse(value < 0, 0.2, -0.2),
    label_hjust = ifelse(value < 0, 0, 1)
  )

my_plot <- ggplot(df, aes(x = kegg, y = value, fill = class)) +
  geom_bar(stat = "identity") +
  #geom_bar(stat = "identity", col = "black") +
  geom_text(aes(y = label_y, label = kegg, hjust = label_hjust)) +
  coord_flip() +
  scale_fill_manual(values = c(HFD = "#00AFBB", MNO863 = "#FC4E07")) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.direction='vertical',
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey80", linetype = "dashed"),
        panel.grid.minor.x = element_blank()) +
  scale_y_continuous(expression(log[10](italic("LDA score"))),
                     breaks = -8:8, limits = c(-8, 8))

ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/SGB/SGB_dg_adj/pic/hbar_0.05_species_20201027.jpg",my_plot,width = 6,height = 14)

# metacyc_uniref50ec
library(ggplot2)
library(dplyr)

df <- read.table("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/SGB/SGB_dg_adj/input/pic_metacyc_uniref50ec_in",header = T,sep="\t")
df <- df %>%
  mutate(
    kegg = factor(kegg, levels = kegg[order(value, decreasing = FALSE)]),
    label_y = ifelse(value < 0, 0.2, -0.2),
    label_hjust = ifelse(value < 0, 0, 1)
  )

my_plot <- ggplot(df, aes(x = kegg, y = value, fill = class)) +
  geom_bar(stat = "identity") +
  #geom_bar(stat = "identity", col = "black") +
  geom_text(aes(y = label_y, label = kegg, hjust = label_hjust)) +
  coord_flip() +
  scale_fill_manual(values = c(HFD = "#00AFBB", MNO863 = "#FC4E07")) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = c(0,0.9),
        legend.justification = c(0,0.9),
        legend.direction='vertical',
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey80", linetype = "dashed"),
        panel.grid.minor.x = element_blank()) +
  scale_y_continuous(expression(log[10](italic("LDA score"))),
                     breaks = -5:5, limits = c(-5, 5))

ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/lefse/SGB/SGB_dg_adj/pic/hbar_0.05_metacycUniref50ec_20201029.jpg",my_plot,width = 20,height = 40)

差异item heatmap 图

一文详解如何用 R 语言绘制热图： https://www.leiphone.com/news/201706/IIshduMIhla9SQG6.html
R语言绘制热图——pheatmap： https://blog.csdn.net/sinat_38163598/article/details/72770404  ## 很好的网站
物种丰度聚类热图: https://www.jianshu.com/p/5f4d279c70e2
result： /work/workspace/zhurj/project/1_metadata/mouse22/heatmap

heatmap into creation
source activate python3.6
python /work/workspace/zhurj/script/2_metapro/metapipe/toheadmapin.py /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data/lefse_genus_in /work/workspace/zhurj/project/1_metadata/mouse22/lefse/SGB/SGB_dg_adj/genus/filter_lefse_rus_0.05.res /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data/heatmap_genus_in

python /work/workspace/zhurj/script/2_metapro/metapipe/toheadmapin.py /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data/lefse_species_in /work/workspace/zhurj/project/1_metadata/mouse22/lefse/SGB/SGB_dg_adj/species/filter_lefse_rus_0.05.res /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data/heatmap_species_in

2020-10-29
metacyc
cp  /work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3_uniref50_ec/metacyc/lefse_in /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data/lefse_uniref50ec_metacyc_in
cd /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data
title=$(cat heatmap_species_in | head -n 1) && sed -i "1i$title" lefse_uniref50ec_metacyc_in
sed -i '2,3d' lefse_uniref50ec_metacyc_in
python /work/workspace/zhurj/script/2_metapro/metapipe/toheadmapin.py /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data/lefse_uniref50ec_metacyc_in /work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3_uniref50_ec/metacyc/filter_lefse_run.res /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data/pre_uniref50ec_metacyc_in
# pre_uniref50ec_metacyc_in 添加pathway 名称
去除 untargetd rows from pre_heatmap_uniref50ec_metacyc_in
source activate python3.6
cd /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data
import pandas as pd
import numpy as np
df = pd.read_table('pre_heatmap_uniref50ec_metacyc_in',header=0,index_col=0)
ndf = df * 1000
nndf = ndf.apply(np.log1p)
nndf.to_csv('heatmap_uniref50ec_metacyc_in',index=True,sep='\t',header=True)

kegg
cp  /work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3_uniref50_ec/kegg/lefse_in /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data/lefse_uniref50ec_kegg_in
cd /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data
title=$(cat heatmap_species_in | head -n 1) && sed -i "1i$title" lefse_uniref50ec_kegg_in
sed -i '2,3d' lefse_uniref50ec_kegg_in
lefse_uniref50ec_kegg_in 名字中的|g__ 改为 _g__
.s__ 改为 _s__, .unclassified 改为 _unclassified

python /work/workspace/zhurj/script/2_metapro/metapipe/toheadmapin.py /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data/lefse_uniref50ec_kegg_in /work/workspace/zhurj/project/1_metadata/mouse22/lefse/humann3_uniref50_ec/kegg/adj_filter_lefse_run.res /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data/pre_uniref50ec_kegg_in
替换了pre_uniref50ec_kegg_in pathway的名字 br/path:pathway descrition_s
source activate python3.6
cd /work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data
import pandas as pd
import numpy as np
df = pd.read_table('pre_uniref50ec_kegg_in',header=0,index_col=0)
ndf = df * 1000
nndf = ndf.apply(np.log1p)
nndf.to_csv('heatmap_uniref50ec_kegg_in',index=True,sep='\t',header=True)

# output file : 

R3.6
genus
library(pheatmap)
infile <- '/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data/heatmap_genus_in'
ref_col <- '/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/input/heatmeta'
ref_row <- '/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data/lefse_genus_group'
raw_df <-read.table(infile,sep='\t',header=TRUE,row.names = 1,encoding='utf-8')
raw_mt <- data.matrix(raw_df)
ref_col_df <- read.table(ref_col,sep='\t',header=FALSE,row.names = 1, encoding='utf-8')
ref_row_df <- read.table(ref_row,sep='\t',header=FALSE,row.names = 1, encoding='utf-8')
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/lefse/SGB/SGB_dg_adj/pic/heatmap_genus_20201027.pdf'

pheatmap(raw_mt)
# #main可设置热图的标题，fontsize设置字体大小，filename可直接将热图存出，支持格式png, pdf, tiff, bmp, jpeg，并且可以通过width, height设置图片的大小；
# pheatmap(test, cellwidth = 15, cellheight = 12, main = "Example heatmap", fontsize = 8, filename = "test.pdf")
#a <- pheatmap(raw_mt, cluster_rows = True, cluster_cols = True）
#cluster_mt <- raw_mt[a$tree_row$order,a$tree_col$order]
#cluster_df = as.data.frame(cluster_mt)
#cluster_colname = colnames(cluster_df)
#col_group = ref_df[cluster_colname,]
#annotation_col = data.frame(Group = factor(col_group))
colnames(ref_col_df) = 'Group'
colnames(ref_row_df) = 'Phylum'
pheatmap(raw_mt, main="Dif genera ()", annotation_col = ref_col_df, annotation_row = ref_row_df, scale='row', cluster_col = FALSE, cutree_rows=2, filename = ofile, width = 12, height = 8)

species
library(pheatmap)
infile <- '/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data/heatmap_species_in'
ref_col <- '/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/input/heatmeta'
ref_row <- '/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data/lefse_species_group'
raw_df <-read.table(infile,sep='\t',header=TRUE,row.names = 1,encoding='utf-8')
raw_mt <- data.matrix(raw_df)
ref_col_df <- read.table(ref_col,sep='\t',header=FALSE,row.names = 1, encoding='utf-8')
ref_row_df <- read.table(ref_row,sep='\t',header=FALSE,row.names = 1, encoding='utf-8')
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/lefse/SGB/SGB_dg_adj/pic/heatmap_species_20201027.pdf'

pheatmap(raw_mt)
# #main可设置热图的标题，fontsize设置字体大小，filename可直接将热图存出，支持格式png, pdf, tiff, bmp, jpeg，并且可以通过width, height设置图片的大小；
# pheatmap(test, cellwidth = 15, cellheight = 12, main = "Example heatmap", fontsize = 8, filename = "test.pdf")
#a <- pheatmap(raw_mt, cluster_rows = True, cluster_cols = True）
#cluster_mt <- raw_mt[a$tree_row$order,a$tree_col$order]
#cluster_df = as.data.frame(cluster_mt)
#cluster_colname = colnames(cluster_df)
#col_group = ref_df[cluster_colname,]
#annotation_col = data.frame(Group = factor(col_group))
colnames(ref_col_df) = 'Group'
colnames(ref_row_df) = 'Phylum'
pheatmap(raw_mt, main="Dif species (G vs D)", annotation_col = ref_col_df, annotation_row = ref_row_df, scale='row', cluster_col = FALSE, cutree_rows=2, filename = ofile, width = 12, height = 10)


metacyc
library(pheatmap)
infile <- '/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data/heatmap_uniref50ec_metacyc_in'
ref_col <- '/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/input/heatmeta'
ref_row <- '/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data/lefse_uniref50ec_metacyc_group'
raw_df <-read.table(infile,sep='\t',header=TRUE,row.names = 1,encoding='utf-8')
raw_mt <- data.matrix(raw_df)
ref_col_df <- read.table(ref_col,sep='\t',header=FALSE,row.names = 1, encoding='utf-8')
ref_row_df <- read.table(ref_row,sep='\t',header=TRUE,row.names = 1, encoding='utf-8')
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/lefse/SGB/SGB_dg_adj/pic/heatmap_uniref50ec_metacyc_20201029.pdf'
colnames(ref_col_df) = 'Group'

pheatmap(raw_mt)
# #main可设置热图的标题，fontsize设置字体大小，filename可直接将热图存出，支持格式png, pdf, tiff, bmp, jpeg，并且可以通过width, height设置图片的大小；
# pheatmap(test, cellwidth = 15, cellheight = 12, main = "Example heatmap", fontsize = 8, filename = "test.pdf")
#a <- pheatmap(raw_mt, cluster_rows = True, cluster_cols = True）
#cluster_mt <- raw_mt[a$tree_row$order,a$tree_col$order]
#cluster_df = as.data.frame(cluster_mt)
#cluster_colname = colnames(cluster_df)
#col_group = ref_df[cluster_colname,]
#annotation_col = data.frame(Group = factor(col_group))
#colnames(ref_row_df) = 'Phylum'
pheatmap(raw_mt, main="Dif metacyc ", annotation_col = ref_col_df, annotation_row = ref_row_df, scale='row', cluster_col = FALSE, cutree_rows=2, filename = ofile, width = 18, height = 18)

kegg
library(pheatmap)
infile <- '/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data/heatmap_uniref50ec_kegg_in'
ref_col <- '/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/input/heatmeta'
ref_row <- '/work/workspace/zhurj/project/1_metadata/mouse22/SGB/GDgroup/data/lefse_kegg_group'
raw_df <-read.table(infile, sep='\t', header=TRUE, row.names = 1, encoding='utf-8')
raw_mt <- data.matrix(raw_df)
ref_col_df <- read.table(ref_col,sep='\t',header=FALSE,row.names = 1, encoding='utf-8')
ref_row_df <- read.table(ref_row,sep='\t',header=TRUE,row.names = 1, encoding='utf-8')
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/lefse/SGB/SGB_dg_adj/pic/heatmap_uniref50ec_kegg_20201029.pdf'
colnames(ref_col_df) = 'Group'

#pheatmap(raw_mt)
# #main可设置热图的标题，fontsize设置字体大小，filename可直接将热图存出，支持格式png, pdf, tiff, bmp, jpeg，并且可以通过width, height设置图片的大小；
# pheatmap(test, cellwidth = 15, cellheight = 12, main = "Example heatmap", fontsize = 8, filename = "test.pdf")
#a <- pheatmap(raw_mt, cluster_rows = True, cluster_cols = True）
#cluster_mt <- raw_mt[a$tree_row$order,a$tree_col$order]
#cluster_df = as.data.frame(cluster_mt)
#cluster_colname = colnames(cluster_df)
#col_group = ref_df[cluster_colname,]
#annotation_col = data.frame(Group = factor(col_group))
#colnames(ref_row_df) = 'Phylum'
ann_colors = list(Count=c(C1="#FFFFFF",C2="#90EE90",C3="#1E90FF",C4="#E6C3C3",C5="#F28500",C6="#FF4500"),
Group = c(HFD="#1B9E77", MNO863="#D95F02"))
pheatmap(raw_mt, main="Dif KEGG ", annotation_col = ref_col_df, annotation_row = ref_row_df, annotation_colors = ann_colors, scale='row', cluster_col = FALSE, cutree_rows=2, filename = ofile, width = 18, height = 18)

结果存在的地址：
/work/workspace/zhurj/project/1_metadata/mouse22/lefse/SGB/SGB_dg_adj/pic


### 
KEGG pathway enrichment
判断数据是否属于正态分布
test  a Seriers
shapiro.test(test)

source activate python3.6
python keggOrg.py /work/workspace/zhurj/database/kegg/20201029/ko00001_20201019.keg /work/workspace/zhurj/database/kegg/20201029

--------------------------------------------------------------------------------
### kegg enrichment --------------- 
source activate R3.6
library(clusterProfiler)
R
gene_list <- read.table("/work/workspace/zhurj/project/1_metadata/mouse22/tmp/genelist")
gene <- gene_list[,1]
background <- read.table("/work/workspace/zhurj/database/kegg/20201029/gene_path_levelC", sep="\t" ,header=FALSE)
kegg2name_data <- read.table("/work/workspace/zhurj/database/kegg/20201029/path_annotation", sep="\t" ,header=FALSE)
result <- enricher(gene,TERM2GENE=background,TERM2NAME=kegg2name_data,pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH")

# pathview: https://cloud.tencent.com/developer/article/1539928
source activate R4.0
library(pathview)
data <- read.table("/work/workspace/zhurj/project/1_metadata/mouse22/tmp/png_gene", sep="\t",header=TRUE,row.names = 1)
p <- pathview(gene.data = data, pathway.id='04726',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)

--------------------------------------------------------------------------------
source activate R3.6
R
# diff ko enriched pathway
library(clusterProfiler)
infile <- '/work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all/enrichment/ko/dif_ko_in'
ofile <- '/work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all/enrichment/ko/dif_ko_enrichedkegg'
gene_list <- read.table(infile)
gene <- gene_list[,1]
background <- read.table("/work/workspace/zhurj/database/kegg/20201029/gene_path_levelC", sep="\t" ,header=FALSE)
kegg2name_data <- read.table("/work/workspace/zhurj/database/kegg/20201029/path_annotation", sep="\t" ,header=FALSE)
result <- enricher(gene,TERM2GENE=background,TERM2NAME=kegg2name_data,pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH")
write.table(result, ofile, sep = '\t', quote = FALSE, col.names = TRUE)

# all ko enriched pathway
infile <- '/work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all/enrichment/ko/ko_all'
ofile <- '/work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all/enrichment/ko/ko_all_enrichedkegg'
gene_list <- read.table(infile)
gene <- gene_list[,1]
background <- read.table("/work/workspace/zhurj/database/kegg/20201029/gene_path_levelC", sep="\t" ,header=FALSE)
kegg2name_data <- read.table("/work/workspace/zhurj/database/kegg/20201029/path_annotation", sep="\t" ,header=FALSE)
result <- enricher(gene,TERM2GENE=background,TERM2NAME=kegg2name_data,pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH")
write.table(result, ofile, sep = '\t', quote = FALSE, col.names = TRUE)

python /work/workspace/zhurj/script/2_metapro/metapipe/koFilter.py /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all/merge_filter/filter_ko.tsv /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all/enrichment/ko

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link 'python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py {1}/humann3_uniref50_all/enrichment/ko/ko_adj {1}/input/meta_group {1}/humann3_uniref50_all/enrichment/ko/lefse_in' ::: /work/workspace/zhurj/project/1_metadata/mouse22 
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

source activate /work/workspace/zhurj/bin/miniconda3/envs/py2
parallel --link -j 1 'lefse-format_input.py {1}/lefse_in {1}/lefse_format.in -c 1 -u 2 -f r -o 1000000 && \
run_lefse.py {1}/lefse_format.in {1}/lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1 ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all/enrichment/ko 
parallel --link -j 1 'lefse-format_input.py {1}/lefse_in {1}/lefse_format.in -c 1 -u 2 -f r -o 1000000 && \
run_lefse.py {1}/lefse_format.in {1}/lefse_run_all.res -a 1 -w 1 -l 0 -y 1 ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all/enrichment/ko 
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/py2

source activate python3.6
parallel --link 'python /work/workspace/zhurj/lib/python/script/lefse_fitler.py {1}/lefse_run.res {1}/filter_lefse_run.res ' ::: /work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all/enrichment/ko
source deactivate python3.6

R4.0
# R包limma作差异基因分析 https://cloud.tencent.com/developer/article/1667505
library(limma)
#library(edgeR)
limma_gene_f <- '/work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all/enrichment/ko/limma_ko'
dif_gene_f <- '/work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all/enrichment/ko/limma_dif_ko'
infile <- '/work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all/enrichment/ko/ko_adj'
exprSet <- read.delim(infile, row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
groupf <- '/work/workspace/zhurj/project/1_metadata/mouse22/input/sgb_group'
pregroup <- read.delim(groupf, row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, header = FALSE)
group <- as.factor(pregroup[,1])
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(group)
#norm <- voom(exprSet, design, plot = TRUE)

#fit <- lmFit(norm, design, method = 'ls')
fit <- lmFit(exprSet, design, method = 'ls')
contrast <- makeContrasts('MNO863-HFD', levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

#qqt(fit2$t, df = fit2$df.prior+fit2$df.residual, pch = 16, cex = 0.2)
#abline(0,1)
#p 值校正、提取差异分析结果，详见 ?topTable
diff_gene <- topTable(fit2, number = Inf, adjust.method = 'fdr')
head(diff_gene, 10)
diff_gene_filter <- diff_gene[ with(diff_gene, abs(diff_gene$logFC) >= 1 & diff_gene$adj.P.Val < 0.05),]
write.table(diff_gene, limma_gene_f, col.names = NA, sep = '\t', quote = FALSE)
write.table(diff_gene_filter, dif_gene_f, col.names = NA, sep = '\t', quote = FALSE)

'''
#例如这里根据 |log2FC| >= 1 & FDR p-value < 0.01 定义“差异”
diff_gene[which(diff_gene$logFC >= 1 & diff_gene$adj.P.Val < 0.01),'sig'] <- 'red'
diff_gene[which(diff_gene$logFC <= -1 & diff_gene$adj.P.Val < 0.01),'sig'] <- 'blue'
diff_gene[which(abs(diff_gene$logFC) < 1 | diff_gene$adj.P.Val >= 0.01),'sig'] <- 'gray'
 
log2FoldChange <- diff_gene$logFC
FDR <- diff_gene$adj.P.Val
plot(log2FoldChange, -log10(FDR), pch = 20, col = diff_gene$sig)
abline(v = 1, lty = 2)
abline(v = -1, lty = 2)
abline(h = -log(0.01, 10), lty = 2)
'''

source activate R4.0
library(pathview)
infile = '/work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all/enrichment/ko/ko03001_in'
data <- read.table(infile, sep="\t",header=TRUE,row.names = 1)
p <- pathview(gene.data = data, pathway.id='03001',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)


library(pathview)
infile = '/work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all/enrichment/png/ko_in'
data <- read.table(infile, sep="\t",header=TRUE,row.names = 1)
p <- pathview(gene.data = data, pathway.id='00520',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00620',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='02040',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='03016',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00240',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='03011',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='02035',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='99997',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='01007',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00010',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00720',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00970',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='01011',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00230',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00250',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='02030',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00500',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00020',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00541',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00270',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00051',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00550',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00650',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='03400',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00770',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00630',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00400',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='03430',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00640',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='01502',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='01005',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00260',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00780',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00340',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='99978',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00052',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00030',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='03060',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00860',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='02048',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00710',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00670',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='99977',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='03440',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00290',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00040',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='02060',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00061',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00540',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='02024',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='99976',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='03070',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00300',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='04112',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00983',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00190',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00220',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='04122',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00680',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00730',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00521',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='03410',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00920',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00450',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='03030',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='02020',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00750',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00760',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='03012',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='03018',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00633',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='99975',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00790',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='03009',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00430',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00900',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00471',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00511',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='99985',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00330',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00561',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='03032',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='00564',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)
p <- pathview(gene.data = data, pathway.id='01501',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
diversity
USEARCH V11: https://www.drive5.com/usearch/manual/cmd_beta_div.html
input 文件需满足条件
1. input 文件 start with #OTU
2. 数据是read count数
python3.6
cd /work/workspace/zhurj/project/1_metadata/mouse22/diversity/input
import pandas as pd
import numpy as np
import os
infile = '/work/workspace/zhurj/project/1_metadata/mouse22/diversity/input/species_all_sam'
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/diversity/input/species_gb16_0.01'
df = pd.read_table(infile,header=0,index_col=0)
tmpdf = df.iloc[:,6:22]
tmpdf['count'] = (tmpdf >= 0.0001).astype(float).sum(axis=1)
resultdf_filter = tmpdf.loc[ tmpdf['count'] >= 3.0 ].iloc[:,0:16]
fdf = resultdf_filter * 10000 
final_fdf = fdf.round()
final_fdf.to_csv(ofile,header=True,index=True,sep='\t')

cd  /work/workspace/zhurj/project/1_metadata/mouse22/diversity/input
# usearch v11: https://drive5.com/usearch/manual/cmd_beta_div.html


usearch -beta_div species_gb16_0.01 -metrics unifrac,unifrac_binary
# usearch -beta_div otutable.txt -metrics jaccard,bray_curtis
# The Unifrac metric requires a tree for the OTUs, which is specified using the -tree option. The file must be in Newick format. The tree can be generated using the -cluster_agg command using the OTU FASTA file as input. If no tree is specified, the Unifrac metric will not be calcualated.
# phyloseq
R4.0
R
library(tidyverse)
library(phyloseq)
infile <- "/work/workspace/zhurj/project/1_metadata/mouse22/diversity/input/species_OUT"
df <- read.delim(infile,sep='\t',row.names = 1, header = T)
ndf <- df %>% select(c(num_range('MNA00',239:254)))
ndf_f <- ndf %>% mutate(flag = rowSums(ndf > 0.0001)) 
row.names(ndf_f) <- rownames(ndf)
nndf_f <- ndf_f %>% filter(flag >= 3) %>% select(contains('MNA'))
# ndf_f %>% summarise_all(., ~ if(is.numeric(.)) sum(.) else "add")
# data %>%
+  bind_rows(summarise_all(., ~ if (is.numeric(.)) sum(.) else "add"))
nndf_f[nndf_f < 0.0001] <- 0
ndf_ft = round(nndf_f * 10000)
# column are species, row are sample
raremax <- min(rowSums(t(ndf_ft)))
Srare <- rarefy(t(ndf_ft),raremax)
S <- specnumber(t(ndf_ft))
plot(S,Srare, xlab="Observed No. of Species", ylab="Rarefied No. of Species")
abline(0,1)
rarecurve(t(ndf_ft), step=50, cex=0.5)

#rarecurve(df,step = 200, sample = raremax, col = "blue", cex = 0.6,xlim=c(0, 20000))


beta diversity
Bray Curtis
python /work/workspace/zhurj/script/2_metapro/metapipe/distformat.py /work/workspace/zhurj/project/1_metadata/mouse22/diversity/input/bray_curtis.txt \
/work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group  /work/workspace/zhurj/project/1_metadata/mouse22/diversity/input/gb_beta_in

R4.0
library(ggplot2)
library(ggpubr)
my_comparisons <- list(c("HFD","MNO863"),c("HFD","HFD_MNO863"),c("HFD_MNO863","MNO863"))
Data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/diversity/input/gb_beta_in",header=T)
b1 <- ggplot(data=Data, aes(x=group,y=bray_curtis))+
geom_boxplot(aes(fill=group)) +
scale_fill_manual(breaks = c("HFD","MNO863","HFD_MNO863"), values=c("#00AFBB", "#FC4E07","#E7B800")) +
stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", aes(label = ..p.signif..)) + 
ggtitle("Bary Curtis Distance")+
theme(plot.title = element_text(face = "plain", vjust = 1, hjust=0.5, size = 10),
      axis.text.x = element_text(angle=0, hjust=0.5), axis.line.x=element_line(),axis.line.y=element_line(),
      axis.title.y=element_blank(), axis.text.y = element_text(angle=0, hjust=1, vjust=0.5),
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"),
      axis.text=element_text(face="plain", size=8, color = 'black'),
      panel.grid.major.y = element_blank(),,
      panel.grid.minor = element_blank(),
      panel.background=element_rect(fill='transparent', color='black',  size = 0.5),
      panel.border=element_rect(fill='transparent', color='black'),
      ) +
scale_x_discrete(name="", breaks = c("HFD", "MNO863","HFD_MNO863"),labels=c("HFD" = "HFD", "MNO863" = "MNO863", "HFD_MNO863" = "HFD_MNO863"),expand = c(0.25, 0.25)) +
scale_y_continuous(name="distance", limits=c(0.0, 0.7), breaks=seq(0.0,0.7,0.1)) +
guides(fill=FALSE)
#g1 <- ggarrange(b1,b2,b3,ncol=3,nrow=1,labels=c("A","","",""))
ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/diversity/pic/beta_div_braycurtis20201109.jpg",b1,width = 3,height = 3)

# violin plot
library(ggplot2)
library(ggpubr)
my_comparisons <- list(c("HFD","MNO863"))
Data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/diversity/input/gb_beta_in",header=T)
b1 <- ggplot(data=Data, aes(x=group,y=bray_curtis)) +
geom_violin(aes(fill=group),trim=FALSE) +
geom_boxplot(width=0.1) +
scale_x_discrete(name="", limits=c("HFD", "MNO863"),labels=c("HFD" = "HFD", "MNO863" = "MNO863"),expand = c(0.25, 0.25)) +
scale_y_continuous(name="distance", limits=c(0.0, 0.6), breaks=seq(0.0,0.6,0.1)) +
stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", aes(label = ..p.signif..)) + 
ggtitle("Bary Curtis Distance")+
theme(plot.title = element_text(face = "plain", vjust = 1, hjust=0.5, size = 10),
      axis.text.x = element_text(angle=0, hjust=0.5), axis.line.x=element_line(),axis.line.y=element_line(),
      axis.title.y=element_blank(), axis.text.y = element_text(angle=0, hjust=1, vjust=0.5),
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"),
      axis.text=element_text(face="plain", size=8, color = 'black'),
      panel.grid.major.y = element_blank(),,
      panel.grid.minor = element_blank(),
      panel.background=element_rect(fill='transparent', color='black',  size = 0.5),
      panel.border=element_rect(fill='transparent', color='black'),
      ) 
ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/diversity/pic/beta_div_braycurtisviolin_20201110.jpg",b1,width = 4,height = 4)

design = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group", header=FALSE, row.names= 1, sep="\t") 
Data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/diversity/input/bray_curtis.txt", sep="\t", header=T, row.names= 1)
idx =  rownames(design) %in% colnames(Data)
sub_design = design[idx,]
pcoa = cmdscale(Data, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points = as.data.frame(pcoa$points) # get coordinate string, format to dataframme
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig
points = cbind(points, sub_design)
p1 <- ggplot(points, aes(x=x, y=y, color=sub_design)) +
geom_point(aes(shape=sub_design),size=0.5) + 
labs(title = "PCoA - Bray-Curtis",
     x=paste("PC1: (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
     y=paste("PC2: (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="") ) +
scale_colour_manual(name="Group", values= c("#FC4E07","#00AFBB")) +
scale_y_continuous(limits = c(-0.2,0.3), breaks = c(-0.2, -0.1, 0, 0.1, 0.2, 0.3), position = "left",) + 
scale_x_continuous(limits = c(-0.3,0.5), breaks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5), position = "bottom",) +
theme(plot.title = element_text(face="plain", hjust=0.5, size = 10,vjust = .5 ),
     axis.title.x = element_text(face="plain", size=6, color = 'black',vjust = 0.5, margin=margin(0.1,0.1,0.1,0.1,"pt")),
     axis.title.y = element_text(face="plain", size=6, color = 'black',vjust = 0.5, margin=margin(0.1,0.1,0.1,0.1,"pt")),
     axis.text=element_text(face="plain", size=8, color = 'black'),
     panel.background=element_rect(fill='transparent', color='black',  size = 0.5),
     panel.border=element_rect(fill='transparent', color='black'),
     axis.ticks.length = unit(-0.15,"lines"), 
     axis.text.x = element_text(margin=margin(5,5,5,5,"pt")),
     axis.text.y = element_text(margin=margin(5,5,5,5,"pt")),
     legend.position = "right",
     ) +
geom_hline(yintercept = c(0), color = "gray", size = 0.2) +
geom_vline(xintercept = c(0), color = "gray", size = 0.2) +
stat_ellipse(linetype="dashed",size=0.2) 
#      legend.position = "none",
ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/diversity/pic/PCoA_braycurtis20201109.jpg",p1,width = 4,height = 2)


# core-pan curve
# roary: the pan genome pipeline https://sanger-pathogens.github.io/Roary/

/work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript /work/workspace/zhurj/software/proteinortho_curves-master/corePan_curve.r -p /work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/test -o /work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/corePan_curve_test

#violin plot
# http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
p<-ggplot(ToothGrowth, aes(x=dose, y=len, fill=dose)) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1)
  #stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.1 )

# http://www.maimengkong.com/kyjc/585.html
# R语言绘制花瓣图（Flower plot）示例: https://cloud.tencent.com/developer/article/1667176
# flower_plot
library(plotrix)
sample_id <- c(paste0("c",1:10), paste0("d",1:10))
otu_num <- round(runif(20,1700,2900))
core_num <- 1000 

ellipse_col <- c('#E6C3C3','#F08080','#A52A2A','#FF7F50','#4D1F00','#FFDAB9','#F4A460','#CD853F','#DAA520','#9ACD32','#73B839','#228B22','#00FFEF','#48D1CC','#008080','#00477D','#4D80E6','#87CEEB','#B57EDC','#FFB3E6')

flower_plot  <- function(sample, otu_num,core_otu,start,a,b,r,ellipse_col,circle_col){
	par(bty='n',ann=F,xaxt='n','yaxt'='n',mar=c(1,1,1,1))
	plot(c(0,10),c(0,10),type='n')
	n <- length(sample)
	deg <- 360/n
	res <- lapply(1:n,function(t){
      draw.ellipse(
      x=5 + cos((start + deg *(t-1)) * pi /180),
      y = 5 + sin((start + deg * (t-1)) * pi / 180),
      col = ellipse_col[t],
      a = a, 
      b = b, 
      angle = deg * (t-1))

      text(x=5 + 2.5 * cos((start + deg * (t-1)) * pi/180),
      y = 5 + 2.5 * sin((start + deg * (t-1)) * pi/180 ),
      otu_num[t])

      if(deg * (t-1) < 180 && deg *(t-1) > 0){
        text(
          x = 5 + 3.3 * cos((start + deg * (t-1)) * pi/180),
          y = 5 + 3.3 * sin((start + deg * (t-1)) * pi/180),
          sample[t],
          srt = deg * (t-1) - start,
          adj = 1,
          cex = 1
        )
      }else{
        text(
          x = 5 + 3.3 * cos((start + deg * (t-1)) * pi/180),
          y = 5 + 3.3 * sin((start + deg * (t-1)) * pi/180),
          sample[t],
          srt = deg * (t-1) - start,
          adj = 0,
          cex = 1
        )
      }

	})
	draw.circle(x = 5, y = 5, r = r, col = circle_col, border = NA)
	text(x = 5, y = 5, paste('Core:',core_otu))
}

png('flower.png',width = 1500, height = 1500, res = 200, units = 'px')
flower_plot(sample = sample_id, otu_num = otu_num, core_otu = core_num,
start = 90, a = 0.5, b=2, r = 1, ellipse_col = ellipse_col, circle_col = 'white')
dev.off()

# calculate core gene, and unique gene of each sample(total gene - core gene)
cd /work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out
library(dplyr)
infile <- "/work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/mergeAbundance"
flower_dat <- read.delim(infile, header = T, sep = '\t', stringsAsFactors = F, check.names = F)
df <- flower_dat %>% select(c(num_range('MNA00',239:254)))
sample_id <- colnames(df)
dfm <- data.matrix(df)
f<-function(x) sum(x != 0)
pre_otu_num <- apply(dfm,2,f)
sam_sum <- apply(dfm,1,f)
core_num <- length(which(sam_sum == 16))
otu_num <- c(t(pre_otu_num - core_num))

ellipse_col <- c('#6181BD4E','#F348004E','#64A10E4E','#9300264E','#464E044E','#049a0b4E','#4E0C664E','#D000004E','#FF6C004E','#FF00FF4E','#c7475b4E','#00F5FF4E','#BDA5004E','#A5CFED4E','#f0301c4E','#2B8BC34E')

#构建作图函数（参考自大神博客 https://www.cnblogs.com/xudongliang/p/7884667.html）
flower_plot <- function(sample, otu_num, core_otu, start, a, b, r, ellipse_col, circle_col) {
	par( bty = 'n', ann = F, xaxt = 'n', yaxt = 'n', mar = c(1,1,1,1))
	plot(c(0,10),c(0,10),type='n')
	n   <- length(sample)
	deg <- 360 / n
	res <- lapply(1:n, function(t){
		draw.ellipse(x = 5 + cos((start + deg * (t - 1)) * pi / 180),
			y = 5 + sin((start + deg * (t - 1)) * pi / 180),
			col = ellipse_col[t],
			border = ellipse_col[t],
			a = a, b = b, angle = deg * (t - 1))
		text(x = 5 + 2.5 * cos((start + deg * (t - 1)) * pi / 180),
			y = 5 + 2.5 * sin((start + deg * (t - 1)) * pi / 180),
			otu_num[t])
		
		if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
			text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
				y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
				sample[t],
				srt = deg * (t - 1) - start,
				adj = 1,
				cex = 1
				)
		} else {
			text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
				y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
				sample[t],
				srt = deg * (t - 1) + start,
				adj = 0,
				cex = 1
				)
		}			
	})
	draw.circle(x = 5, y = 5, r = r, col = circle_col, border = NA)
	text(x = 5, y = 5, paste('Core:', core_otu))
}

#调用上述函数作图，视情况修改参数
ofile <- '/work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/flower_fg16.png'
png(ofile, width = 1500, height = 1500, res = 200, units = 'px')
flower_plot(sample = sample_id, otu_num = otu_num, core_otu = core_num,
	start = 90, a = 0.5, b = 2, r = 1, ellipse_col = ellipse_col, circle_col = 'white')
dev.off()

venn 图
# http://blog.sciencenet.cn/blog-3419243-1206798.html
library(VennDiagram)
library(dplyr)
infile <- '/work/workspace/zhurj/project/1_metadata/mouse22/diversity/input/species_gb16_0.01'
pre_df <- read.delim(infile,row.names = 1, header =T, sep='\t')
df <- pre_df %>% select(c(num_range('MNA00',239:243)))
df <- t(df)
A<-names(df[1,])[df[1,]>0]
B<-names(df[2,])[df[2,]>0]
C<-names(df[3,])[df[3,]>0] 
D<-names(df[4,])[df[4,]>0]
E<-names(df[5,])[df[5,]>0]
venn.plot<-venn.diagram(x=list(MN0239=A,MN0240=B,MN0241=C,MN0242=D,MN0243=E),
filename <- 'venntest_sam5.png',
fill=c("cornflowerblue", "green", "red","darkorchid1","yellow"),
col = "transparent",alpha = 0.50,cex = 1.0,fontface = "bold",
cat.col = c("darkblue", "darkgreen","orange","darkorchid4",'pink'),
cat.cex = 1.5,cat.dist = 0.07,margin = 0.2)


source activate python3.6
cat /work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/mergeAbundance | grep "^MNA" | awk '{print $1}' > /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/in/unigene_list
split -300041 -d /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/in/unigene_list /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/in/part_ --verbose

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
python /work/workspace/zhurj/script/2_metapro/metapipe/extract_sub_data_from_fasta.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/out/cdhit.fna /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/in/part_00 /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/unigene/part_00.fasta
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
python /work/workspace/zhurj/script/2_metapro/metapipe/extract_sub_data_from_fasta.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/out/cdhit.fna /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/in/part_01 /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/unigene/part_01.fasta
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
python /work/workspace/zhurj/script/2_metapro/metapipe/extract_sub_data_from_fasta.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/out/cdhit.fna /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/in/part_02 /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/unigene/part_02.fasta
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

cd /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/p
time srun -o part_02.out -e part_02.err -N 1 -c 20 -p slurm128 -w mnclient04 bash part_02.sh &
time srun -o part_00.out -e part_00.err -N 1 -c 20 -p slurm256 -w mnclient01 bash part_00.sh &
time srun -o part_01.out -e part_01.err -N 1 -c 20 -p slurm256 -w mnclient02 bash part_01.sh &

# nr diamond db
/work/workspace/zhurj/project/1_metadata/mouse22/tmp/p
time srun -o diamondtest.out -e diamondtest.err -N 1 -c 20 -p slurm256 -w mnclient02 bash diamondtest.sh &

diamondtest.sh
```
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
diamond blastx --db /work/workspace/zhurj/database/nr/diamond/nr_microbe.dmnd --query /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/unigene/test1 --out /work/workspace/zhurj/project/1_metadata/mouse22/tmp/diamond/nr.tab --outfmt 6 --sensitive --max-target-seqs 20 --evalue 1e-5 --id 30 --block-size 20.0 --tmpdir /work/workspace/zhurj/project/1_metadata/mouse22/tmp/diamond/tmp --index-chunks 1
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

/work/workspace/zhurj/instal/megan/tools/blast2lca -i /work/workspace/zhurj/project/1_metadata/mouse22/tmp/diamond/nr.tab -f BlastTab -m BlastX -o /work/workspace/zhurj/project/1_metadata/mouse22/tmp/diamond/blast2lca -mdb /work/workspace/zhurj/database/megan/megan-map-Jul2020-2.db
```
eggnogtest.sh
```
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
diamond blastx --db /work/workspace/zhurj/database/eggnog/e5/diamond/eggnog5.dmnd --query /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/unigene/test1 --out /work/workspace/zhurj/project/1_metadata/mouse22/tmp/diamond/eggnog.tab --outfmt 6 --sensitive --max-target-seqs 20 --evalue 1e-5 --id 30 --block-size 20.0 --tmpdir /work/workspace/zhurj/project/1_metadata/mouse22/tmp/diamond/tmp --index-chunks 1
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

```
cazytest.sh
```
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
diamond blastx --db /work/workspace/zhurj/database/cazy/diamond/CAZyDB.dmnd --query /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/unigene/test1 --out /work/workspace/zhurj/project/1_metadata/mouse22/tmp/diamond/cazy.tab --outfmt 6 --sensitive --max-target-seqs 20 --evalue 1e-5 --id 30 --block-size 20.0 --tmpdir /work/workspace/zhurj/project/1_metadata/mouse22/tmp/diamond/tmp --index-chunks 1
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

```
cardtest.sh
```
source activate /work/workspace/zhurj/bin/miniconda3/envs/rgi
export MPLCONFIGDIR=/work/workspace/zhurj/project/1_metadata/mouse22/tmp/diamond/tmp
rgi main --input_sequence /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/unigene/test_all --output_file /work/workspace/zhurj/project/1_metadata/mouse22/tmp/diamond/gri_test --input_type contig --clean -a DIAMOND
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/rgi

```

cd /work/workspace/zhurj/project/1_metadata/mouse22/tmp/p
time srun -o eggnogtest.out -e eggnogtest.err -N 1 -c 20 -p slurm256 -w mnclient02 bash eggnogtest.sh & # 运行时间 0m43.892s,5条序列
time srun -o cazytest.out -e cazytest.err -N 1 -c 20 -p slurm256 -w mnclient02 bash cazytest.sh & # 50 条序列，5条序列做test时，结果文件为空
time srun -o cardtest.out -e cardtest.err -N 1 -c 20 -p slurm256 -w mnclient01 bash cardtest.sh & # 1000 条序列

diamond_part_00.sh
```
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
diamond blastx --db /work/workspace/zhurj/database/nr/diamond/nr_microbe.dmnd --query /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/unigene/part_00.fasta --out /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/part_00.tab --outfmt 6 --sensitive --max-target-seqs 20 --evalue 1e-5 --id 30 --block-size 20.0 --tmpdir /work/workspace/zhurj/project/1_metadata/mouse22/tmp/diamond/tmp --index-chunks 1

cat /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/part_00.tab | awk -F "\t" '{if($(NF-1) < 1e-10) print $0}' > /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/part_00_filter.tab

time nohup /work/workspace/zhurj/instal/megan/tools/blast2lca -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/part_00_filter.tab -f BlastTab -m BlastX -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/part_00_blast2lca.tab -mdb /work/workspace/zhurj/database/megan/megan-map-Jul2020-2.db > /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/part_00_blast2lca.log 2>&1 &

source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

```

diamond_part_01.sh
```
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
diamond blastx --db /work/workspace/zhurj/database/nr/diamond/nr_microbe.dmnd --query /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/unigene/part_01.fasta --out /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/part_01.tab --outfmt 6 --sensitive --max-target-seqs 20 --evalue 1e-5 --id 30 --block-size 20.0 --tmpdir /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/tmp --index-chunks 1

cat /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/part_01.tab | awk -F "\t" '{if($(NF-1) < 1e-10) print $0}' > /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/part_01_filter.tab

time nohup /work/workspace/zhurj/instal/megan/tools/blast2lca -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/part_01_filter.tab -f BlastTab -m BlastX -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/part_01_blast2lca.tab -mdb /work/workspace/zhurj/database/megan/megan-map-Jul2020-2.db > /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/part_01_blast2lca.log 2>&1 &

source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

```
diamond_part_02.sh
```
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
diamond blastx --db /work/workspace/zhurj/database/nr/diamond/nr_microbe.dmnd --query /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/unigene/part_02.fasta --out /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/part_02.tab --outfmt 6 --sensitive --max-target-seqs 20 --evalue 1e-5 --id 30 --block-size 20.0 --tmpdir /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/tmp --index-chunks 1

cat /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/part_02.tab | awk -F "\t" '{if($(NF-1) < 1e-10) print $0}' > /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/part_02_filter.tab

time /work/workspace/zhurj/instal/megan/tools/blast2lca -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/part_02_filter.tab -f BlastTab -m BlastX -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/part_02_blast2lca.tab -mdb /work/workspace/zhurj/database/megan/megan-map-Jul2020-2.db & 

source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

eggnog_diamond.sh
```
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
diamond blastx --db /work/workspace/zhurj/database/eggnog/e5/diamond/eggnog5.dmnd --query /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/out/cdhit.fna --out /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/eggnog/unigene_eggnog.tab --outfmt 6 --sensitive --max-target-seqs 20 --evalue 1e-5 --id 30 --block-size 20.0 --tmpdir /work/workspace/zhurj/project/1_metadata/mouse22/tmp/diamond/tmp --index-chunks 1
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

```
cazy_diamond.sh
```
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
diamond blastx --db /work/workspace/zhurj/database/cazy/diamond/CAZyDB.dmnd --query /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/out/cdhit.fna --out /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy/unigene_cazy.tab --outfmt 6 --sensitive --max-target-seqs 20 --evalue 1e-5 --id 30 --block-size 20.0 --tmpdir /work/workspace/zhurj/project/1_metadata/mouse22/tmp/diamond/tmp --index-chunks 1
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

```

CARD 抗性基因组预测
card_diamond.sh
```
source activate /work/workspace/zhurj/bin/miniconda3/envs/rgi
export MPLCONFIGDIR=/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/tmp
rgi main --input_sequence /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/out/cdhit.fna --output_file /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/gri_result --input_type contig --clean  -a DIAMOND
source deacitvate /work/workspace/zhurj/bin/miniconda3/envs/rgi

```
Uniref50_all
uniref50_diamond.sh
```
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
diamond blastx --db /work/workspace/zhurj/database/humann/201901/uniref50_all/uniref50_201901.dmnd --query /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/out/cdhit.fna --out /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/unigene_uniref50all.tab --outfmt 6 --sensitive --max-target-seqs 20 --evalue 1e-5 --id 30 --block-size 20.0 --tmpdir /work/workspace/zhurj/project/1_metadata/mouse22/tmp/diamond/tmp --index-chunks 1
source deacitvate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

```

cd /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/p
# run time: 734m54.813s
time srun -o diamond_part_00.out -e diamond_part_00.err -N 1 -c 20 -p slurm256 -w mnclient01 bash diamond_part_00.sh &
# run time: 724m54.813s
time srun -o diamond_part_01.out -e diamond_part_01.err -N 1 -c 20 -p slurm256 -w mnclient01 bash diamond_part_01.sh &
# run time :779m31.726s
time srun -o diamond_part_02.out -e diamond_part_02.err -N 1 -c 20 -p slurm256 -w mnclient02 bash diamond_part_02.sh &

# run time :15m34.990s
time srun -o cazy_diamond.out -e cazy_diamond.err -N 1 -c 20 -p slurm256 -w mnclient01 bash cazy_diamond.sh &
time srun -o eggnog_diamond.out -e eggnog_diamond.err -N 1 -c 20 -p slurm256 -w mnclient02 bash eggnog_diamond.sh &  # run time : 193m30.409s
time srun -o card_diamond.out -e card_diamond.err -N 1 -c 20 -p slurm256 -w mnclient01 bash card_diamond.sh &
time srun -o uniref50_diamond.out -e uniref50_diamond.err -N 1 -c 20 -p slurm256 -w mnclient02 bash uniref50_diamond.sh & # run time : 143m6.166s

# MDS, PCA， PCoA分析
# MDS
library(vegan)
library(ggplot2)
library(tidyverse)
iris.data <- subset(iris, select = -Species)
iris.dist <- vegdist(iris.data)
iris.mds <- monoMDS(iris.dist)
mds.df <- as.data.frame(iris.mds$points)
mds.df$Species <- iris$Species
mds.poly <- mds.df %>% group_by(Species) %>% slice(chull(MDS1,MDS2))
ggplot(mds.df, aes(MDS1,MDS2,color = Species)) +
geom_point( aes(shape = Species)) +
labs(title = "MDS - Bray-Curtis",x="MDS1",y="MDS2") +
stat_ellipse(level = 0.95, show.legend = F) + 
theme_bw() +
theme(panel.border = element_blank(),
panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), 
axis.line = element_line(color = "black"),
plot.title = element_text(face="plain", hjust=0.5, vjust = .5 ))
#aes(fill = Species) +
#geom_polygon(data = mds.poly, alpha = 0.3, linetype = "blank")

# PCoA 分析
library(vegan)
library(ggplot2)
library(tidyverse)
iris.data <- subset(iris, select = -Species)
iris.dist <- vegdist(iris.data,method="bray")
pcoa = cmdscale(iris.dist, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points = as.data.frame(pcoa$points) # get coordinate string, format to dataframme
#sub_design = iris$Species
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig
points = cbind(points, Species = iris$Species)
p1 <- ggplot(points, aes(x=x, y=y, color=Species)) +
geom_point(aes(shape=Species)) + 
labs(title = "PCoA - Bray-Curtis",
     x=paste("PC1: (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
     y=paste("PC2: (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="") ) +
stat_ellipse(level = 0.95, show.legend = F) + 
theme_bw() +
theme(panel.border = element_blank(),
panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), 
axis.line = element_line(color = "black"),
plot.title = element_text(face="plain", hjust=0.5, vjust = .5 ))

#      legend.position = "none",
ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/diversity/pic/PCoA_braycurtis20201109.jpg",p1,width = 4,height = 2)

# PCA
# https://www.plob.org/article/22240.html
library(vegan)
library(ggplot2)
library(tidyverse)
iris.data <- subset(iris, select = -Species)
iris.pca <- prcomp(iris.data, center = TRUE, scale = TRUE)
iris.pca_unscaled <- prcomp(iris.data, center = TRUE, scale = FALSE)

pca.df <- data.frame(iris.pca$x, Species = iris$Species)
#plot(iris.pca$x[,1],iris.pca$x[,2])
percentage <- round(iris.pca$sdev / sum(iris.pca$sdev) * 100,2)
percentage <- paste(colnames(pca.df)," (",as.character(percentage),"%", ")", sep="")
ggplot(pca.df,aes(x=PC1,y=PC2,color=Species)) +
geom_point(aes(shape=Species)) +
labs(title = "PCA plot",
     x=percentage[1],
     y=percentage[2]) +
stat_ellipse(level = 0.95, show.legend = F) + 
theme_bw() +
theme(panel.border = element_blank(),
panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), 
axis.line = element_line(color = "black"),
plot.title = element_text(face="plain", hjust=0.5, vjust = .5 ))

iris.anosim <- anosim(iris.data,iris$Species)
summary(iris.anosim)
plot(iris.anosim,col = c('#FFD700','#FF7F00','#EE2C2C','#D02090'),ylab="",xlab="")

# ko abundance PCA pca-20201201
library(vegan)
library(ggplot2)
library(tidyverse)
metafile <- '/work/workspace/zhurj/project/1_metadata/mouse22/input/meta_all'
infile <- '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/ko_abundance'
odir <- '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic'
mycolors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")

meta_df <- read.table(metafile,sep='\t',header=T)
cols = c('SampleId','Group')
part_meta_df = meta_df[cols]
data_df <- read.table(infile,sep='\t',header=T,row.names = 1)
dataT_df <- data.frame(t(data_df))
merge_df <- merge(dataT_df,part_meta_df, by.x = "row.names", by.y = 'SampleId', all.x = TRUE)

ko.data <- subset(merge_df, select = c(-Row.names,-Group))
ko.pca <- prcomp(ko.data, center = TRUE, scale = TRUE)
ko.pca_unscaled <- prcomp(ko.data, center = TRUE, scale = FALSE)

pca.df <- data.frame(ko.pca$x, Group = merge_df$Group)
#plot(iris.pca$x[,1],iris.pca$x[,2])
percentage <- round(ko.pca$sdev / sum(ko.pca$sdev) * 100,2)
percentage <- paste(colnames(pca.df)," (",as.character(percentage),"%", ")", sep="")
pca_pic <- ggplot(pca.df,aes(x=PC1,y=PC2,color=Group)) +
geom_point(aes(shape=Group)) +
labs(title = "PCA plot",
     x=percentage[1],
     y=percentage[2]) +
stat_ellipse(level = 0.95, show.legend = F) + 
theme_bw() +
theme(panel.border = element_blank(),
panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), 
axis.line = element_line(color = "black"),
plot.title = element_text(face="plain", hjust=0.5, vjust = .5 ))
ggsave("/work/workspace/zhurj/project/1_metadata/mouse22/diff/png/dif_phylum_20201010.jpg", p, width = 8, height =6)


used_color = mycolors[1:4]
ko.anosim <- anosim(ko.data,merge_df$Group)
summary(ko.anosim)
#pca_anosim <- plot(ko.anosim,col = c('#FFD700','#FF7F00','#EE2C2C','#D02090'),ylab="",xlab="")
pca_anosim <- plot(ko.anosim,col = used_color,ylab="",xlab="")


/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/pca_graph.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/ko_abundance -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_all -n ko

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/pca_graph.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/eggnog/level2_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_all -n eggnog_cog

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/pca_graph.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy/level2_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_all -n cazy_level2

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/pca_graph.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/aro_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_all -n card_aro

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/pca_mds_pcoa.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/ko_abundance -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_all -n ko

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/pca_mds_pcoa.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/eggnog/level2_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_all -n eggnog_cog

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/pca_mds_pcoa.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy/level2_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_all -n cazy_level2

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/pca_mds_pcoa.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/aro_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_all -n card_aro

fetch KO of each module
python /work/workspace/zhurj/lib/python/script/webfetch_KO_from_keggmoduleid.py /work/workspace/zhurj/database/kegg/20201029/test /work/workspace/zhurj/database/kegg/20201029/tmp

cd /work/workspace/zhurj/database/kegg/20201116
cat pathway_ko | grep -P "path:ko" > pathwayko_ko
cat module_ko | sort -u > unique_module_ko

data = pd.read_csv('unigene_cazy.tab',sep="\t",header=None,index_col = None)
df_data = data.loc[ (data[10] <= 0.00001) & (data[11] >= 60.0) ]


python /work/workspace/zhurj/script/2_metapro/metapipe/cazySum.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy/test /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy
python /work/workspace/zhurj/script/2_metapro/metapipe/cazySum.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy/unigene_cazy.tab /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy
GT      GT22    MNA00233_00011
CBM     CBM32   MNA00233_00042
GH      GH35    MNA00233_00042
CE      CE0     MNA00233_00055
GH      GH76    MNA00233_00056
GH      GH20    MNA00233_00058
CBM     CBM32   MNA00233_00060
GH      GH2     MNA00233_00060
GH      GH2     MNA00233_00061
GH      GH2     MNA00233_00062

data = pd.read_csv('detail_all',sep='\t',index_col = None, header = None)
cols = [0,2]
df_data = data[cols].drop_duplicates()
df_result = df_data[0].value_counts()

data = pd.read_csv('detail_all',sep='\t',index_col = None, header = None)
'''
df_group = data.groupby([0,1])[2].count() # only print once for each level1
for indexs in df_group.index:
  level1 = indexs[0]
  level2 = indexs[1]
  count - df_group[indexs].values[0]
  genelist = data[ (data[0] == level1) & (data[1] == level2) ][2].values.tolist()
  genestr = ",".join(sorted(genelist))
  print("{}\t{}\t{}".format(level1,level2,genestr))
'''
df_group = data.groupby([0,1])[2].size().reset_index(name='count')
#genelist = data[ (data[0] == 'GH') & (data[1] == 'GH0') ][2].values.tolist()
#genestr = ",".join(sorted(genelist))

for indexs in df_group.index:
  tmplist = df_group.loc[indexs].values.tolist()
  level1 = tmplist[0]
  level2 = tmplist[1]
  count = tmplist[2]
  genelist = data[ (data[0] == level1) & (data[1] == level2) ][2].values.tolist()
  genestr = ",".join(sorted(genelist))
  print("{}\t{}\t{}\t{}".format(level1,level2,count,genestr))


CARD result organization
import pandas as pd
infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/gri_result.txt'
data = pd.read_csv(infile,sep="\t",header=0,index_col = None)
df_data = data.loc[ (data['Best_Identities'] >=90.0) ]
cols = ['Contig','Best_Hit_ARO','Best_Identities','ARO','Drug Class','Resistance Mechanism']
df_part = df_data[cols]
df1 = df_part.drop(['Drug Class'],axis=1).join(df_part['Drug Class'].str.split('; ',expand=True).stack().reset_index(level = 1,drop=True).rename('Drug_Class'))
df2 = df1.drop(['Resistance Mechanism'],axis=1).join(df1['Resistance Mechanism'].str.split('; ',expand=True).stack().reset_index(level=1,drop=True).rename('Resistance_Mechanism'))
df2['Contig'] = df2['Contig'].str[0:14]

cols = ['Drug_Class','Contig']
df_drug = df2[cols].drop_duplicates()
df_druggene = df_drug.groupby(by='Drug_Class').apply(lambda x: ','.join(x['Contig']))
df_count = df_drug['Drug_Class'].value_counts()
count_df = df_count.to_frame()
count_df.columns = ['Count']
druggene_df = df_druggene.to_frame()
druggene_df.columns = ['genelist']
result_df = count_df.join(druggene_df,how='outer')
result_df.index.name = 'Drug_Class'

cols = ['Resistance_Mechanism','Contig']
df_resistance = df2[cols].drop_duplicates()
df_resgene = df_resistance.groupby(by='Resistance_Mechanism').apply(lambda x: ','.join(x['Contig']))
df_count = df_resistance['Resistance_Mechanism'].value_counts()
count_df = df_count.to_frame()
count_df.columns = ['Count']
resgene_df = df_resgene.to_frame()
resgene_df.columns = ['genelist']
result_df = count_df.join(resgene_df,how='outer')
result_df.index.name = 'Resistance_Mechanism'


# a是DataFrame格式的数据集
# a.index.name = 'date'
# a.columns.name = 'code'
# rename columns
# df.rename(columns={'two':'twotwo'},inplace=True)
# df.columns = ['newname']
# rename index
# df.rename(index={'a':'aa','b':'bb'},inplace=True)
# merge multiple dataframe
count_df.join(druggene_df,how='inner')
python /work/workspace/zhurj/script/2_metapro/metapipe/cardSum.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/gri_result.txt /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card

eggnog
# create orgid2level2
memref = '/work/workspace/zhurj/database/eggnog/e5/per_tax_level/all_members.tsv'
pre_df_mem = pd.read_csv(memref,sep='\t',header=None,index_col=None)
pre_df_mem.columns = ['taxid','level2','level2num','orgnum','level_list','org_list']
cols = ['level2','level_list']
part_df = pre_df_mem[cols]
# df_mem = part_df.drop(['level_list'],axis=1).join(part_df['level_list'].str.split(',',expand=True).stack().reset_index(level = 1,drop=True).rename('level_str')) # 内容太多，内存不够，不发运行
orgid2level2f = '/work/workspace/zhurj/database/eggnog/e5/per_tax_level/orgid2level2'
with open(orgid2level2f,'w',encoding='utf-8') as ofp:
  for indexs in part_df.index:
    level2 = part_df.loc[indexs]['level2']
    orgstr = part_df.loc[indexs]['level_list']
    orglist = orgstr.split(',')
    for orgone in orglist:
      con = "{}\t{}\n".format(orgone,level2)
      ofp.write(con)

cd /work/workspace/zhurj/database/eggnog/e5/per_tax_level/p
time srun -o slurm_orgid2level2.out -e slurm_orgid2level2.err -N 1 -c 20 -p slurm256 -w mnclient01 bash slurm_orgid2level2.sh &

```
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
import pandas as pd
annoref = '/work/workspace/zhurj/database/eggnog/e5/per_tax_level/all_annotations.tsv'
level32level2f = '/work/workspace/zhurj/database/eggnog/e5/per_tax_level/orgid2level2'
infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/eggnog/unigene_eggnog.tab'

data = pd.read_csv(infile,sep='\t',header=None,index_col=None)
df_filter = data[ (data[10] <= 1e-5) & (data[11] >= 60.0) ]
df_filter.columns = ['gene','level3','identity','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']
df_filter.drop_duplicates('gene',keep='first',inplace = True)
cols = ['gene','level3']
gene2level3_df = df_filter[cols]

df_anno = pd.read_csv(annoref,sep='\t',header=None,index_col=None)
df_anno.columns = ['taxid','level2','level1','note_']
cols = ['level2','level1']
anno_df = df_anno[cols].drop_duplicates()

level3to2_df = pd.read_csv(level32level2f,sep='\t',header=None,index_col=None)
level3to2_df.columns = ['level3','level2']

df1 = pd.merge(gene2level3_df,level3to2_df,how='left',on='level3')
all_df = pd.merge(df1,anno_df,how='left',on='level2')

# eggnog level2 abundance eggnog-20201201
import pandas as pd
detailf = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/eggnog/detail_all'
samabunf = '/work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/mergeAbundance'
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/eggnog/level2_abun'
levels_df = pd.read_csv(detailf,sep='\t',header=0)
data_df = pd.read_csv(samabunf,sep='\t',header=0)
level2_df = levels_df[['level2','gene']].dropna(how='any')
cog_level2_df = level2_df[level2_df['level2'].str.contains("COG")]
data_df.rename(columns={'Sample':'gene'},inplace=True)
merge_df = pd.merge(cog_level2_df,data_df,on='gene',how='inner')
samcols = [x for x in merge_df.columns if 'MNA' in x ]
select_cols = ['level2'] + samcols
part_merge_df = merge_df[select_cols]
group_eggnog_level2_df = part_merge_df.groupby('level2').sum()
group_eggnog_level2_df.to_csv(ofile,sep='\t',header=True,index=True)


# cazy level2 abundance
import pandas as pd
detailf = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy/detail_all'
samabunf = '/work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/mergeAbundance'
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy/level2_abun'
levels_df = pd.read_csv(detailf,sep='\t',header=None)
levels_df.columns = ['level1','level2','gene']
data_df = pd.read_csv(samabunf,sep='\t',header=0)
level2_df = levels_df[['level2','gene']].dropna(how='any')
data_df.rename(columns={'Sample':'gene'},inplace=True)
merge_df = pd.merge(level2_df,data_df,on='gene',how='inner')
samcols = [x for x in merge_df.columns if 'MNA' in x ]
select_cols = ['level2'] + samcols
part_merge_df = merge_df[select_cols]
group_cazy_level2_df = part_merge_df.groupby('level2').sum()
group_cazy_level2_df.to_csv(ofile,sep='\t',header=True,index=True)

# card argo abundance
import pandas as pd
detailf = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/detail_all'
samabunf = '/work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/mergeAbundance'
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/aro_abun'

python /work/workspace/zhurj/script/2_metapro/metapipe/pca_in.py -l /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/card_detail -i /work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/mergeAbundance -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/aro_abun --listkey ARO Contig

source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

```

cd /work/workspace/zhurj/script/2_metapro/metapipe
time srun -o test.out -e test.err -N 1 -c 20 -p slurm256 -w mnclient01 bash test.sh &

python /work/workspace/zhurj/script/2_metapro/metapipe/eggnogSum.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/eggnog/test /work/workspace/zhurj/database/eggnog/e5/per_tax_level/all_annotations.tsv /work/workspace/zhurj/database/eggnog/e5/per_tax_level/orgid2level2 /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/eggnog

python /work/workspace/zhurj/script/2_metapro/metapipe/eggnogSum.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/eggnog/unigene_eggnog.tab /work/workspace/zhurj/database/eggnog/e5/per_tax_level/all_annotations.tsv /work/workspace/zhurj/database/eggnog/e5/per_tax_level/orgid2level2 /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/eggnog


kegg reference: /work/workspace/zhurj/database/kegg/20201116
all_ko
all_module
all_pathway
module_ko
pathway_ko
pathwayko_ko
readme
unique_module_ko

/work/workspace/zhurj/database/humann/201901/misc
map_ko_name.txt.gz
map_ko_uniref50.txt.gz

/work/workspace/zhurj/database/kegg/20201029
gene_annotation
gene_path_levelA
gene_path_levelB
gene_path_levelC
path_annotation
path_levels
uniref2ko
module_ko

# read uniref to ko
import pandas as pd
import gzip
import sys
import os

uniref2kof = '/work/workspace/zhurj/database/humann/201901/misc/map_ko_uniref50.txt.gz'
uniref2koformatf = '/work/workspace/zhurj/database/kegg/20201029/uniref2ko'

contents = []
koid = ''
uniref50id = ''
with gzip.open(uniref2kof,'rt',encoding='utf-8') as fp:
  contents = fp.readlines()

with open(uniref2koformatf,'w',encoding='utf-8') as ofp:
  for line in contents:
    if not line or line.startswith('#'):
      continue
    tmpcontents = line.strip().split("\t")
    koid = tmpcontents[0]
    for uniref50id in tmpcontents[1:]:
      con = "{}\t{}\n".format(uniref50id,koid)
      ofp.write(con)


import pandas as pd
import gzip
import sys
import os

def unirefname(x):
  pos = x.find('|')
  return x[0:pos]

odir = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all'
infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/unigene_uniref50all.tab'
uniref2koformatf = '/work/workspace/zhurj/database/kegg/20201029/uniref2ko'
samabunf = '/work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/mergeAbundance'

levelAf = '/work/workspace/zhurj/database/kegg/20201029/gene_path_levelA'
levelBf = '/work/workspace/zhurj/database/kegg/20201029/gene_path_levelB'
levelCf = '/work/workspace/zhurj/database/kegg/20201029/gene_path_levelC'
level_all = '/work/workspace/zhurj/database/kegg/20201029/gene_path_all_level'
modulef = '/work/workspace/zhurj/database/kegg/20201029/module_ko'
pathannof = '/work/workspace/zhurj/database/kegg/20201029/path_annotation' 

omoduleabunf = os.path.join(odir,'kegg_module_abundance')
olevelCabunf = os.path.join(odir,'kegg_levelC_abundance')
olevelBabunf = os.path.join(odir,'kegg_levelB_abundance')
olevelAabunf = os.path.join(odir,'kegg_levelA_abundance')
okoabunf = os.path.join(odir,'ko_abundance')
gene2kof = os.path.join(odir,'gene2kof')
opathgenecountf = os.path.join(odir,'kegg_path_genecount')
olevelAcountf = os.path.join(odir,'kegg_levelA_genecount')
olevelBcountf = os.path.join(odir,'kegg_levelB_genecount')
olevelCcountf = os.path.join(odir,'kegg_levelC_genecount')
con = ''

data = pd.read_csv(infile,sep='\t',header=None,index_col=None)
data.columns = ['gene','uniref50','identity','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']
df_filter = data[ (data['evalue'] <= 1e-5) & (data['bitscore'] >= 60.0) ]
df_filter.drop_duplicates('gene',keep='first',inplace = True)
df_filter['uniref50'] = df_filter['uniref50'].apply(lambda x: unirefname(x))
cols = ['gene','uniref50']
data_df = df_filter[cols].drop_duplicates()


data = pd.read_csv(uniref2koformatf,sep='\t',header=None,index_col=None)
data.columns = ['uniref50','koid']
data.drop_duplicates()


df1 = pd.merge(data_df,data,how='left',on='uniref50')
#df1.dropna(axis=0,how)
df1_filter = df1.dropna(axis=0,subset=['koid'])
df1_filter.to_csv(gene2kof,sep='\t',header=None,index=None)

data = pd.read_csv(samabunf,sep='\t',header=0,index_col=None)
abun_df = data.rename(columns={'Sample':'gene'})

gka_df = pd.merge(df1_filter,abun_df,how='left',on='gene')
gka_filter = gka_df.dropna(axis=0,how='any')
secols = [x for x in gka_filter.columns if 'MNA' in x]
secols.insert(0,'koid')
ka_df = gka_filter[secols]
groupko_df = ka_df.groupby('koid').sum()
groupko_df.to_csv(okoabunf,sep='\t',header=True,index=True)

# pathway summary
# levela 
data = pd.read_csv(levelAf,sep='\t',header=None,index_col=None)
data.columns = ['path','koid']
data.drop_duplicates()
levelA_df = pd.merge(data,groupko_df,how='left',on='koid')
levelA_filter = levelA_df.dropna(axis=0,how='any')
secols[0] = 'path'
levelA_path = levelA_filter[secols]
groupA_df = levelA_path.groupby('path').sum()
groupA_df.to_csv(olevelAabunf,sep='\t',header=True,index=True)

# levelB
data = pd.read_csv(levelBf,sep='\t',header=None,index_col=None)
data.columns = ['path','koid']
data.drop_duplicates()
levelB_df = pd.merge(data,groupko_df,how='left',on='koid')
levelB_filter = levelB_df.dropna(axis=0,how='any')
secols[0] = 'path'
levelB_path = levelB_filter[secols]
groupB_df = levelB_path.groupby('path').sum()
groupB_df.to_csv(olevelBabunf,sep='\t',header=True,index=True)

# levelC
data = pd.read_csv(levelCf,sep='\t',header=None,index_col=None)
data.columns = ['path','koid']
data.drop_duplicates()
levelC_df = pd.merge(data,groupko_df,how='left',on='koid')
levelC_filter = levelC_df.dropna(axis=0,how='any')
secols[0] = 'path'
levelC_path = levelC_filter[secols]
groupC_df = levelC_path.groupby('path').sum()
groupC_df.to_csv(olevelCabunf,sep='\t',header=True,index=True)

# module
data = pd.read_csv(modulef,sep='\t',header=None,index_col=None)
data.columns = ['module','koid']
data.drop_duplicates()
module_df = pd.merge(data,groupko_df,how='left',on='koid')
module_filter = module_df.dropna(axis=0,how='any')
secols[0] = 'module'
module_ko = module_filter[secols]
module_result_df = module_ko.groupby('module').sum()
module_result_df.to_csv(omoduleabunf,sep='\t',header=True,index=True)

sam_cols =  [x for x in gka_filter.columns if 'MNA' in x]
module_result_df[sam_cols].apply(lambda x : sum(x),axis=0)
#groupA_df[sam_cols].apply(lambda x : sum(x),axis=0)
#groupB_df[sam_cols].apply(lambda x : sum(x),axis=0)
#groupC_df[sam_cols].apply(lambda x : sum(x),axis=0)

# 统计各pathway包含的基因数目
# gka_filter
# levelA
data = pd.read_csv(level_all,sep='\t',header=None,index_col=None)
data.columns = ['path','koid']
data_level_filter = data.drop_duplicates()
path_sam_list = ['path'] + sam_cols

path_list = list(data_level_filter['path'].unique())
pathone = path_list[0]
path_df = data_level_filter[ data_level_filter['path'] == pathone ]
tmp_pkg_df = pd.merge(gka_filter,path_df,how='left',on='koid')
tmp_pkg_filter = tmp_pkg_df.dropna(axis=0,subset=['path'])
pkg_one_df = tmp_pkg_filter[sam_cols]
pkg_sum_series = (pkg_one_df > 0.0).astype(int).sum(axis=0)
pkg_sum_one_df = pkg_sum_series.to_frame()
pkg_sum_one_df.columns = [pathone]
result_df = pkg_sum_one_df
for pathone in path_list[1:]:
  path_df = data_level_filter[ data_level_filter['path'] == pathone ]
  tmp_pkg_df = pd.merge(gka_filter,path_df,how='left',on='koid')
  tmp_pkg_filter = tmp_pkg_df.dropna(axis=0,subset=['path'])
  pkg_one_df = tmp_pkg_filter[sam_cols]
  pkg_sum_series = (pkg_one_df > 0.0).astype(int).sum(axis=0)
  pkg_sum_one_df = pkg_sum_series.to_frame()
  pkg_sum_one_df.columns = [pathone]
  result_df = pd.concat([result_df,pkg_sum_one_df],axis=1)

final_result_df = result_df.loc[:,~(result_df == 0).all(axis=0)]
fT_result = final_result_df.T
fT_result['path'] = fT_result.index

path_anno = pd.read_csv(pathannof,sep='\t',header=None,index_col=None)
path_anno.columns = ['path','anno']
result_ann_df = pd.merge(fT_result,path_anno,on='path',how='left')
new_cols = ['path','anno'] + sam_cols
result_ann_df = result_ann_df[new_cols]
result_ann_df.to_csv(opathgenecountf,sep='\t',header=True,index=None)

python /work/workspace/zhurj/script/2_metapro/metapipe/keggSum.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/unigene_uniref50all.tab /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test
python /work/workspace/zhurj/script/2_metapro/metapipe/keggUnigeneSum.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/unigene_uniref50all.tab /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/unigene

做bar图
# 查看配置信息，mpl.rcParams
#plt.savefig("filename.png")
#plt.savefig('squares_plot.png', bbox_inches='tight')
#plt.show()

查了一下，原来pyplot有一个interactive开关：
#ioff() # Turn interactive plotting off
#ion() # Turn interactive plotting on
source activate ipython
# 进入命令行
ipython --pylab
%matplotlib inline
import numpy as np
import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt

x = np.random.randn(100)
fig = plt.figure(figsize=(8,6))
ax1 = fig.add_subplot(2,2,1)
ax1.hist(x,bins=100)
ax1.set_title('default')

ax2 = fig.add_subplot(2,2,2)
ax2.hist(x,bins=100,density=True)
ax2.set_title('normed')

ax3 = fig.add_subplot(2,2,3)
ax3.hist(x,bins=100,density=True,histtype='step')
ax3.set_title('histtype')

ax4 = fig.add_subplot(2,2,4)
rsult = ax4.hist(x,bins=100,density=True,rwidth=0.6,label='aaa',color='r')
ax4.set_title('rwidth&color')

%pylab

ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/test.png'
plt.savefig(ofile)
# plt.show()

df = pd.DataFrame({'lab': ['A', 'B', 'C'], 'val': [10, 30, 20]})
ax = df.plot.barh(x='lab', y='val')
plt.savefig(ofile)

infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelA_genecount'
path_gene_count_df = pd.read_csv(infile,sep='\t',header=0,index_col=None)

cols = ['anno','MNA00233']
part_df = path_gene_count_df[cols]
ax = path_gene_count_df.plot.barh(x='anno',figsize=(15,10))
samcols =  [x for x in path_gene_count_df.columns if 'MNA' in x]
part_df = path_gene_count_df[samcols]
part_df.index = path_gene_count_df['anno']
plt.figure()
part_df.plot(figsize=(15,10))
plt.xticks(rotation=50)
plt.savefig(ofile, dpi=300, bbox_inches="tight")

# plt.xlabel("x轴自定义描述")
# plt.ylable('y轴自定义描述')
# plt.xlim(tuple) # 设置x轴范围 plt.xlim((-5,5))
# plt.ylim((tuple)) # 设置y轴的范围 plt.ylim((-5,5))
# plt.xticks(arr) # 设置x坐标轴刻度显示设置 plt.xticks([-2.0, -1.0, 0, 1, 2, 4, 8])
# plt.yticks(arr) # 设置y坐标轴刻度显示设置 plt.yticks(np.linspace(-3, 3, 10))

'''
# 移动坐标轴

figure = plt.figure(num=100)

# x,y
x = np.linspace(-4, 4, 50)
y = x ** 2 - 4

# 获取到坐标轴
ax = plt.gca()

# 隐藏右边、上边的spine
ax.spines["right"].set_color("none")
ax.spines["top"].set_color("none")

# 移动两个spine到0,0，达到把坐标轴移动的目的
ax.spines["bottom"].set_position(("data", 0))
ax.spines["left"].set_position(("data", 0))

xx_label = r"y = x ^ 2 - 4"
x_label = r"y = x"
plt.title("here is title")

# 绘图
plt.plot(x, y, color="#ff0000")

# 显示图例
plt.legend()

plt.plot(x, x)

# 显示网格
plt.grid(True)

# 显示图例
plt.legend(labels=[xx_label, x_label])

plt.show()
'''

'''
plt.bar(x_arr, y_arr, facecolor="#999999", edgecolor="white")
x_arr表示直方图x数据组成的数组；
y_arr表示直方图y数据组成的数组；
facecolor="#999999"表示直方图的颜色；
edgecolor="white"表示直方图的边框颜色；
'''
'''
plt.text(-3, 40, "function: y = x * x", size = 15,\
         family = "fantasy", color = "r", style = "italic", weight = "light",\
         bbox = dict(facecolor = "r", alpha = 0.2))
# 第一个参数是x轴坐标
# 第二个参数是y轴坐标
# 第三个参数是要显式的内容
# alpha 设置字体的透明度
# family 设置字体
# size 设置字体的大小
# style 设置字体的风格
# wight 字体的粗细
# bbox 给字体添加框，alpha 设置框体的透明度， facecolor 设置框体的颜色
'''

'''
# list ,array 相互转换

import pandas as pd
import numpy as np

array=numpy.array(list)
list=array.tolist()

### dict
# dict to dataframe
df = pd.DataFrame(dict)
# dict to series
series = pd.Series(dict) # series.shape

### list
# list to series
series = pd.Series(list,index=list('abcde')) # list 为（5,3）
# list to dataframe
df = pd.DataFrame(list,index=list('abcde'),columns=['c1','c2','c3'])
# [i + j for i in 'ABC' for j in '123'] , result: ['A1', 'A2', 'A3', 'B1', 'B2', 'B3', 'C1', 'C2', 'C3']
# list to array
ndarray = np.array(list)

## array
# array to dataframe
# list
list = [[2000, 'Ohino', 1.5],
        [2001, 'Ohino', 1.7],
        [2002, 'Ohino', 3.6],
        [2001, 'Nevada', 2.4],
        [2002, 'Nevada', 2.9]]  # type(data) 为 list

ndarray = np.array(list)
df = pd.DataFrame(ndarray,index = list('abcde'),columns=['C' + j for j in '123'])


# series, dataframe to array
array = np.array(df)


# dataframe to dict
df1 = df.to_dict(orient='dict')
df2 = df.to_dict(orient='list')
df3 = df.to_dict(orient='series')
df4 = df.to_dict(orient='records')
'''


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
plt.cla()
pathlevelf = '/work/workspace/zhurj/database/kegg/20201029/path_levels'
infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/unigene/kegg_levelB_genecount'
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/test.png'
level_df = pd.read_csv(pathlevelf,sep='\t',index_col=None,header=None)
data = pd.read_csv(infile,sep='\t',index_col=None,header=0)
level_df.columns = ['level3','level2','level1']

value = data['MNAall']
Y = np.arange(0,value.size,1)
height = value.to_list()
width = np.array([0.8] * value.size)

name = data['anno'].to_list()
name.reverse()

#plt.barh(Y,height,width, align ="edge", color ="g",edgecolor="k",linewidth=2,
#        tick_label=["A","B","C","D","E","F","G","H"],xerr =xerr,ecolor="b",capsize=6)
plt.barh(Y,height,width, align ="center", color ="g",edgecolor="w",
        tick_label=name,capsize=6)
plt.title("levelB gene count",fontproperties="SimHei",fontsize = 18)
plt.xlabel("gene count")
plt.ylabel("levelB pathways")

#array = np.array(value)
#plt.barh(np.arange(array.size), array, label='test')
plt.legend()
plt.yticks(size = 3)
plt.xticks(size = 3)
plt.savefig(ofile, dpi=300, bbox_inches="tight")



# dcols = ['MNA00' + str(j) for j in range(239,247,1)]
# gcols = ['MNA00' + str(j) for j in range(247,255,1)]
# dgavg_df = part_level21_df.join([d_df['davg'],g_df['gavg']])
# dgavg_sort_df = dgavg_df.sort_values(by="level1",ascending=True,inplace=False, axis=0)
# list(dgavg_df['level1'].unique())

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
plt.cla()
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/test.png'
men_means, men_std = (20, 35, 30, 35, 27), (2, 3, 4, 1, 2)
women_means, women_std = (25, 32, 34, 20, 25), (3, 5, 2, 3, 3)

ind = np.arange(len(men_means))  # the x locations for the groups
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.barh(ind - width/2, men_means, width, xerr=men_std, label='Men')
rects2 = ax.barh(ind + width/2, women_means, width, xerr=women_std, label='Women')


colorlist = ['red','red','green','blue','blue']

redstr='green'
# Add some text for labels, title and custom x-axis tick labels, etc.
colorstr = 'red'
ax.set_xlabel('Scores')
ax.set_title('Scores by group and gender')
ax.set_yticks(ind)
ax.set_yticklabels(('G1', 'G2', 'G3', 'G4', 'G5'))
ax.get_yticklabels()[2].set_color(redstr)
for i in range(0,len(colorlist)):
  colorstr = colorlist[i]
  ax.get_yticklabels()[i].set_color(colorstr)

ax.legend()

def autolabel(rects, xpos='middle'):
  va = {'middle':'middle','top':'top','bottom':'bottom'}
  offset = {'middle':0,'top':1,'bottom':-1}
  
  for rect in rects:
    width = rect.get_width()
    ax.annotate('{}'.format(width),
    xy=(width,rect.get_y() + rect.get_height() / 2),
    xytext=(offset[xpos]*3, 3),  # use 3 points offset
    textcoords="offset points",  # in both directions
    va=va[xpos], ha='left')

autolabel(rects1, "middle")
autolabel(rects2, "middle")

fig.tight_layout()
plt.savefig(ofile, dpi=300, bbox_inches="tight")


###### levelB gene count
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
plt.cla()
pathlevelf = '/work/workspace/zhurj/database/kegg/20201029/path_levels'
infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/unigene/kegg_levelB_genecount'
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/test.png'
level_df = pd.read_csv(pathlevelf,sep='\t',index_col=None,header=None)
data = pd.read_csv(infile,sep='\t',index_col=None,header=0)
level_df.columns = ['level3','level2','level1']
level21_df = level_df[['level2','level1']].drop_duplicates()
level21_df.rename(columns={'level2':'path'},inplace=True)
mycolor_m = ["#0048BA","#B0BF1A","#7CB9E8","#C0E8D5","#B284BE","#72A0C1","#EDEAE0","#F0F8FF","#C46210","#EFDECD",
"#E52B50","#9F2B68","#F19CBB","#AB274F","#D3212D","#3B7A57","#FFBF00","#FF7E00","#9966CC","#A4C639","#CD9575",
"#665D1E","#915C83","#841B2D","#FAEBD7","#008000","#8DB600","#FBCEB1","#00FFFF","#7FFFD4","#D0FF14","#4B5320",
"#8F9779","#E9D66B","#B2BEB5","#87A96B","#FF9966","#A52A2A","#FDEE00","#568203"]
mycolor_17 = ["maroon","green","navy","teal","olive","indigo","coral","orange","bronze","ocher","sepia","emerald","cerise","dimgray","salmon","lavender","pink"]
mycolor_10 = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",]

dcols = ['MNA00' + str(j) for j in range(239,247,1)]
gcols = ['MNA00' + str(j) for j in range(247,255,1)]
data['dmean'] = data[dcols].mean(axis=1)
data['dstd'] = data[dcols].std(axis=1)
data['gmean'] = data[gcols].mean(axis=1)
data['gstd'] = data[gcols].std(axis=1)

secols = ['path','anno','dmean','dstd','gmean','gstd']
sec_df = data[secols]
merge_df = pd.merge(sec_df,level21_df,on='path',how='left')
merge_sort_df = merge_df.sort_values(by="level1",ascending=True,inplace=False, axis=0)
level1_unique_list = list(merge_sort_df['level1'].unique())
level1_unique_list.sort()
sec_cols = mycolor[0:len(level1_unique_list)]
level1_cols_df = pd.DataFrame({'level1':level1_unique_list,'color':sec_cols})
path_mean_std_color_df = pd.merge(merge_df,level1_cols_df,on='level1',how='left')
path_mean_std_color_sort_df = path_mean_std_color_df.sort_values(by="level1",ascending=True,inplace=False, axis=0)

d_mean =  list(path_mean_std_color_sort_df['dmean'])
d_std = list(path_mean_std_color_sort_df['dstd'])
g_mean =  list(path_mean_std_color_sort_df['gmean'])
g_std = list(path_mean_std_color_sort_df['gstd'])
colorlist = list(path_mean_std_color_sort_df['color'])

ind = np.arange(len(d_mean))  # the x locations for the groups
width = 0.35  # the width of the bars

figw = 6.4
figh = 4.8
fig, ax = plt.subplots(figsize=(figw,figh))
plt.xticks(fontsize=6)
plt.yticks(fontsize=4)
#rects1 = ax.barh(ind - width/2, d_mean, width, xerr=d_std, label='HFD')
#rects2 = ax.barh(ind + width/2, g_mean, width, xerr=g_std, label='MNO863')
rects1 = ax.barh(ind - width/2, d_mean, width, label='HFD')
rects2 = ax.barh(ind + width/2, g_mean, width, label='MNO863')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_xlabel('Gene count',size=8)
ax.set_ylabel('Kegg pathway',size=8)
ax.set_title('KEGG levelB (MNO863 && HFD)',size=10)
ax.set_yticks(ind)
ax.set_yticklabels(path_mean_std_color_df['anno'])
for i in range(0,len(colorlist)):
  colorstr = colorlist[i]
  ax.get_yticklabels()[i].set_color(colorstr)

plt.legend(loc='best',fontsize='6')
nbins = len(ind)
ax.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))

'''
def autolabel(rects, xpos='middle'):
  va = {'top':'top','bottom':'bottom'}
  offset = {'middle':0,'top':1,'bottom':-1}
  
  for rect in rects:
    width = rect.get_width()
    ax.annotate('{}'.format(width),
    size = 3,
    xy=(width,rect.get_y() + rect.get_height() / 2),
    xytext=(offset[xpos]*3, 3),  # use 3 points offset
    textcoords="offset points",  # in both directions
    va=va[xpos], ha='left')

autolabel(rects1, "top")
autolabel(rects2, "top")
'''

fig.tight_layout()
plt.savefig(ofile, dpi=500, bbox_inches="tight")



python /work/workspace/zhurj/script/2_metapro/metapipe/kegg_genecount_bargraph_dvsg.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/unigene/kegg_levelB_genecount /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/keggLevelB.png level2
python /work/workspace/zhurj/script/2_metapro/metapipe/kegg_genecount_bargraph_dvsg.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/unigene/kegg_levelC_genecount /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/keggLevelC.png level3
python /work/workspace/zhurj/script/2_metapro/metapipe/kegg_genecount_bargraph_dvsg.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/unigene/kegg_levelA_genecount /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/keggLevelA.png level1

# cazy gene count
import pandas as pd
import numpy as np
import sys
import os

defined_level = 'level1'
level1annof = '/work/workspace/zhurj/database/cazy/20201111/cazyAnno'
samabunf = '/work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/mergeAbundance'
incazydetailf = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy/detail_all'
data = pd.read_csv(samabunf,sep='\t',header=0,index_col=None)
abun_df = data.rename(columns={'Sample':'gene'})
cazy_df = pd.read_csv(incazydetailf,sep='\t',header=None, index_col=None)
cazy_df.columns = ['level1','level2','gene']
cazy_abun_df = pd.merge(abun_df,cazy_df,on='gene',how='left')
level_anno_df = pd.read_csv(level1annof,sep='\t',header=None,index_col=None)
level_anno_df.columns = ['level1','anno']
# filter genes without cazy annotation
cazy_abun_fiter_df = cazy_abun_df.dropna(axis=0,subset=['level1'])
sam_cols =  [x for x in cazy_abun_fiter_df.columns if 'MNA' in x]

# 统计各cazy level1, level2 包含的基因数目
# level1
cazy_sam_list = [defined_level] + sam_cols

path_list = list(cazy_abun_fiter_df[defined_level].unique())
pathone = path_list[0]
path_df = cazy_abun_fiter_df[ cazy_abun_fiter_df[defined_level] == pathone ]
pkg_one_df = path_df[sam_cols]
pkg_sum_series = (pkg_one_df > 0.0).astype(int).sum(axis=0)
pkg_sum_one_df = pkg_sum_series.to_frame()
pkg_sum_one_df.columns = [pathone]
result_df = pkg_sum_one_df
for pathone in path_list[1:]:
  path_df = cazy_abun_fiter_df[ cazy_abun_fiter_df[defined_level] == pathone ]
  pkg_one_df = path_df[sam_cols]
  pkg_sum_series = (pkg_one_df > 0.0).astype(int).sum(axis=0)
  pkg_sum_one_df = pkg_sum_series.to_frame()
  pkg_sum_one_df.columns = [pathone]
  result_df = pd.concat([result_df,pkg_sum_one_df],axis=1)

final_result_df = result_df.loc[:,~(result_df == 0).all(axis=0)]
fT_result = final_result_df.T
fT_result[defined_level] = fT_result.index

result_ann_df = pd.merge(fT_result,level_anno_df,on=defined_level,how='left')
new_cols = [defined_level,'anno'] + sam_cols
result_ann_df = result_ann_df[new_cols]
result_ann_df.to_csv(opathgenecountf,sep='\t',header=True,index=None)


python /work/workspace/zhurj/script/2_metapro/metapipe/cazy_gene_count_sum.py /work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/mergeAbundance /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy/detail_all /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy/cazy_genecount
python /work/workspace/zhurj/script/2_metapro/metapipe/cazy_genecount_bargraph_dvsg.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy/cazy_genecount /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy/cazy_genecount_bar.png
python /work/workspace/zhurj/script/2_metapro/metapipe/cazy_gene_abun_sum.py /work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/mergeAbundance /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy/detail_all /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy/cazy_abun


# CARD Drug and Resistance mechanism summary
import pandas as pd
import numpy as np
import sys
import os

defined_level = 'Resistance_Mechanism'

samabunf = '/work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/mergeAbundance'
infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/card_detail'
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/card_resistance_genecount'
data = pd.read_csv(samabunf,sep='\t',header=0,index_col=None)
abun_df = data.rename(columns={'Sample':'gene'})
card_df = pd.read_csv(infile,sep='\t',header=0,index_col=None)
card_df.rename(columns={'Contig':'gene'}, inplace=True)
part_card_df = card_df[['gene',defined_level]].drop_duplicates()
merge_df = pd.merge(abun_df,part_card_df,on='gene',how='left')
card_merge_df = merge_df.dropna(axis=0,subset=[defined_level])

sam_cols =  [x for x in card_merge_df.columns if 'MNA' in x]
card_sam_list = [defined_level] + sam_cols

path_list = list(card_merge_df[defined_level].unique())
pathone = path_list[0]
path_df = card_merge_df[ card_merge_df[defined_level] == pathone ]
pkg_one_df = path_df[sam_cols]
pkg_sum_series = (pkg_one_df > 0.0).astype(int).sum(axis=0)
pkg_sum_one_df = pkg_sum_series.to_frame()
pkg_sum_one_df.columns = [pathone]
result_df = pkg_sum_one_df
for pathone in path_list[1:]:
  path_df = card_merge_df[ card_merge_df[defined_level] == pathone ]
  pkg_one_df = path_df[sam_cols]
  pkg_sum_series = (pkg_one_df > 0.0).astype(int).sum(axis=0)
  pkg_sum_one_df = pkg_sum_series.to_frame()
  pkg_sum_one_df.columns = [pathone]
  result_df = pd.concat([result_df,pkg_sum_one_df],axis=1)

final_result_df = result_df.loc[:,~(result_df == 0).all(axis=0)]
fT_result = final_result_df.T
fT_result[defined_level] = fT_result.index

new_cols = [defined_level] + sam_cols
result_ann_df = fT_result[new_cols]
result_ann_df.to_csv(ofile,sep='\t',header=True,index=None)

python /work/workspace/zhurj/script/2_metapro/metapipe/card_gene_count_sum.py /work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/mergeAbundance /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/card_detail /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/card_resistance_genecount Resistance_Mechanism
python /work/workspace/zhurj/script/2_metapro/metapipe/card_gene_count_sum.py /work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/mergeAbundance /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/card_detail /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/card_Drug_genecount Drug_Class
python /work/workspace/zhurj/script/2_metapro/metapipe/card_genecount_bargraph_dvsg.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/card_resistance_genecount /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/card_resistance_genecount.png Resistance_Mechanism
python /work/workspace/zhurj/script/2_metapro/metapipe/card_genecount_bargraph_dvsg.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/card_Drug_genecount /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/card_Drug_genecount.png Drug_Class
python /work/workspace/zhurj/script/2_metapro/metapipe/card_gene_abun_sum.py /work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/mergeAbundance /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/card_detail /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/card_resistance_abun Resistance_Mechanism
python /work/workspace/zhurj/script/2_metapro/metapipe/card_gene_abun_sum.py /work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/mergeAbundance /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/card_detail /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/card_Drug_abun Drug_Class

# eggnog genecount summary
import pandas as pd
import numpy as np
import sys
import os

defined_level = 'level1'
levelannof = '/work/workspace/zhurj/database/eggnog/e5/per_tax_level/level1_anno'
samabunf = '/work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/mergeAbundance'
infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/eggnog/detail_all'
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/eggnog/eggnog_genecount'
data = pd.read_csv(samabunf,sep='\t',header=0,index_col=None)
abun_df = data.rename(columns={'Sample':'gene'})
eggnog_df = pd.read_csv(infile,sep='\t',header=0,index_col=None)
level_anno_df = pd.read_csv(levelannof,sep='\t',header=None,index_col=None)
level_anno_df.columns = [defined_level,'anno']
part_eggnog_df = eggnog_df[['gene',defined_level]].drop_duplicates()
merge_df = pd.merge(abun_df,part_eggnog_df,on='gene',how='left')
eggnog_merge_df = merge_df.dropna(axis=0,subset=[defined_level])

sam_cols =  [x for x in eggnog_merge_df.columns if 'MNA' in x]
eggnog_sam_list = [defined_level] + sam_cols

path_list = list(level_anno_df[defined_level].unique())
pathone = path_list[0]
path_df = eggnog_merge_df[ eggnog_merge_df[defined_level] == pathone ]
pkg_one_df = path_df[sam_cols]
pkg_sum_series = (pkg_one_df > 0.0).astype(int).sum(axis=0)
pkg_sum_one_df = pkg_sum_series.to_frame()
pkg_sum_one_df.columns = [pathone]
result_df = pkg_sum_one_df
for pathone in path_list[1:]:
  path_df = eggnog_merge_df[ eggnog_merge_df[defined_level] == pathone ]
  pkg_one_df = path_df[sam_cols]
  pkg_sum_series = (pkg_one_df > 0.0).astype(int).sum(axis=0)
  pkg_sum_one_df = pkg_sum_series.to_frame()
  pkg_sum_one_df.columns = [pathone]
  result_df = pd.concat([result_df,pkg_sum_one_df],axis=1)

#final_result_df = result_df.loc[:,~(result_df == 0).all(axis=0)]
fT_result = result_df.T
fT_result[defined_level] = fT_result.index
merge_result_df = pd.merge(fT_result,level_anno_df,on=defined_level,how='left')

new_cols = [defined_level,'anno'] + sam_cols
result_ann_df = merge_result_df[new_cols]
result_ann_df.to_csv(ofile,sep='\t',header=True,index=None)


python /work/workspace/zhurj/script/2_metapro/metapipe/eggnog_gene_count_sum.py /work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/mergeAbundance /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/eggnog/detail_all /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/eggnog/eggnog_genecount
python /work/workspace/zhurj/script/2_metapro/metapipe/eggnog_genecount_bargraph_dvsg.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/eggnog/eggnog_genecount /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/eggnog/eggnog_genecount.png
python /work/workspace/zhurj/script/2_metapro/metapipe/eggnog_gene_abun_sum.py /work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/mergeAbundance /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/eggnog/detail_all /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/eggnog/eggnog_abun


# kegg abundance bar graph
import pandas as pd
import matplotlib.pyplot as plt

inannof = '/work/workspace/zhurj/database/kegg/20201029/path_annotation'
infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelA_abundance'
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/kegg_levelA_abundance.png'
metaf = '/work/workspace/zhurj/project/1_metadata/mouse22/input/meta_all'

meta_df = pd.read_csv(metaf,sep='\t',header=0,index_col=None)
data = pd.read_csv(infile,sep='\t',header=0,index_col=None)
anno_df = pd.read_csv(inannof,sep='\t',header=None,index_col=None)
anno_df.columns = ['path','anno']
merge_df = pd.merge(data,anno_df,on='path',how='left')
samcols = [x for x in data.columns if 'MNA' in x]
data_df = merge_df[samcols]
data_df.index = merge_df['anno']
data_df['avg'] = data_df.sum(axis=1)
data_df.sort_values(by='avg',ascending=False,inplace=True)
pre_df = data_df[samcols]
pre_df.index = data_df.index
df = pre_df.T
xnames = meta_df['Expeid'].tolist()

plt.cla()
fig,ax = plt.subplots(figsize=(7.0,5.0))
df.plot(ax=ax,kind='bar',stacked=True)
ax.legend(bbox_to_anchor=(1.05,1.0),fontsize='6')
ax.set_title('kegg level 1',fontsize=12)
ax.set_xlabel('samples', fontsize=10)
ax.set_ylabel('Relative Abundance',fontsize=10)
ax.set_xticklabels(meta_df['Expeid'])
ax.xaxis.set_tick_params(rotation=90,labelsize=6,colors='black')
ax.yaxis.set_tick_params(labelsize=6,colors='black')
plt.savefig(ofile, dpi=500, bbox_inches="tight")

python /work/workspace/zhurj/script/2_metapro/metapipe/kegg_level1_abun_bar.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelA_abundance /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/kegg_levelA_abundance.png

## card drug and resistance mechanism
import sys
import pandas as pd
import matplotlib.pyplot as plt

defined_level = 'Drug_Class'
infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/card_Drug_abun'
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/card_Drug_abun.png'
metaf = '/work/workspace/zhurj/project/1_metadata/mouse22/input/meta_all'

meta_df = pd.read_csv(metaf,sep='\t',header=0,index_col=None)
data = pd.read_csv(infile,sep='\t',header=0,index_col=None)
samcols = [x for x in data.columns if 'MNA' in x]
data_df = data[samcols]
data_df.index = data[defined_level]
data_df['avg'] = data_df.sum(axis=1)
data_df.sort_values(by='avg',ascending=False,inplace=True)
pre_df = data_df[samcols]
pre_df.index = data_df.index
df = pre_df.T
xnames = meta_df['Expeid'].tolist()

plt.cla()
fig,ax = plt.subplots(figsize=(7.0,5.0))
df.plot(ax=ax,kind='bar',stacked=True)
ax.legend(bbox_to_anchor=(1.05,1.0),fontsize='6')
ax.set_title(defined_level,fontsize=12)
ax.set_xlabel('samples', fontsize=10)
ax.set_ylabel('Relative Abundance',fontsize=10)
ax.set_xticklabels(meta_df['Expeid'])
ax.xaxis.set_tick_params(rotation=90,labelsize=6,colors='black')
ax.yaxis.set_tick_params(labelsize=6,colors='black')
plt.savefig(ofile, dpi=500, bbox_inches="tight")

python /work/workspace/zhurj/script/2_metapro/metapipe/card_abun_bar.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/card_Drug_abun /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/card_Drug_abun.png Drug_Class
python /work/workspace/zhurj/script/2_metapro/metapipe/card_abun_bar.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/card_resistance_abun /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/card_resistance_abun.png Resistance_Mechanism

# eggnog and cazy
import sys
import pandas as pd
import matplotlib.pyplot as plt

infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy/cazy_abun'
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/cazy_abun.png'
metaf = '/work/workspace/zhurj/project/1_metadata/mouse22/input/meta_all'

meta_df = pd.read_csv(metaf,sep='\t',header=0,index_col=None)
data = pd.read_csv(infile,sep='\t',header=0,index_col=None)
samcols = [x for x in data.columns if 'MNA' in x]
data_df = data[samcols]
data_df.index = data['level1'].str.cat(data['anno'],sep=' : ')
data_df['avg'] = data_df.sum(axis=1)
data_df.sort_values(by='avg',ascending=False,inplace=True)
pre_df = data_df[samcols]
pre_df.index = data_df.index
df = pre_df.T
xnames = meta_df['Expeid'].tolist()

plt.cla()
fig,ax = plt.subplots(figsize=(7.0,5.0))
df.plot(ax=ax,kind='bar',stacked=True)
ax.legend(bbox_to_anchor=(1.05,1.0),fontsize='6')
ax.set_title(defined_level,fontsize=12)
ax.set_xlabel('samples', fontsize=10)
ax.set_ylabel('Relative Abundance',fontsize=10)
ax.set_xticklabels(meta_df['Expeid'])
ax.xaxis.set_tick_params(rotation=90,labelsize=6,colors='black')
ax.yaxis.set_tick_params(labelsize=6,colors='black')
plt.savefig(ofile, dpi=500, bbox_inches="tight")

python /work/workspace/zhurj/script/2_metapro/metapipe/eggnog_cazy_abun_bar.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy/cazy_abun /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/cazy_abun.png Cazy_level1

python /work/workspace/zhurj/script/2_metapro/metapipe/eggnog_cazy_abun_bar.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/eggnog/eggnog_abun /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/eggnog_abun.png Eggnog_level1


core-pan gene curve graph
/work/workspace/liangzj/script/gene_analysis/example_of_cdhitclstr2profile2curve.sh
python cdhitclstr2profile01.py -i Pd_162.ffn -c Pd_162.ffn.clstr -m MNID_list -o genes_profile
python /work/workspace/liangzj/script/Pd_analysis/profile01_2_curve.py -i genes_profile.profile -l 10 -o curve_data.txt
# /work/workspace/liangzj/script/gene_analysis/genes_profile.profile

import pandas as pd
import numpy as np
infile = '/work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/mergeAbundance'
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/corepan/geneprofile.profile'
data = pd.read_csv(infile,sep='\t',header=0,index_col=0)
data[data > 0] = 1
result_df = data.astype(int)
result_df.to_csv(ofile,sep='\t',header=True, index=True)

python /work/workspace/liangzj/script/Pd_analysis/profile01_2_curve.py -i geneprofile.profile -l 10 -o curve_data.txt


python /work/workspace/zhurj/script/2_metapro/metapipe/corePan_profile.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/corepan/newtest -l 10 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/corepan/curve_data.txt

python /work/workspace/zhurj/script/2_metapro/metapipe/corePan_profile.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/corepan/geneprofile.profile -l 10 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/corepan/geneprofile_curve_in

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/core_pan_curve.R
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/corepan/core_pan_curve.R

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/core_pan.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/corepan/test_in -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/corepan/geneprofile_curve_in.png


# R语言绘制带聚类树的堆叠柱形图 20201202
import pandas as pd
import sys
infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy/cazy_abun'
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy/new_cazy_abun'

data = pd.read_csv(infile,sep='\t',header=0)
col_names = data.columns.values.tolist()
key1 = col_names[0]
key2 = col_names[1]
data.index =  data[key1].str.cat(data[key2],sep=' : ')
select_cols = col_names[2:]
o_data = data[select_cols]
o_data.to_csv(ofile,sep='\t',header=True,index=True)

python /work/workspace/zhurj/script/2_metapro/metapipe/reformat_abun.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy/cazy_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy/new_cazy_abun 
python /work/workspace/zhurj/script/2_metapro/metapipe/reformat_abun.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/eggnog/eggnog_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/eggnog/new_eggnog_abun

# bar cluster group file creation
import pandas as pd
import sys
infile = '/work/workspace/zhurj/project/1_metadata/mouse22/input/meta_all'
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/input/meta_bar_cluster'
keylist = ['SampleId','Group']

data = pd.read_csv(infile,sep='\t',header=0)
key1 = keylist[0]
key2 = keylist[1]
part_data = data[keylist]
unique_key2_list = list(data[key2].unique())
unique_num = len(unique_key2_list) + 1
num_list = [x for x in range(1,unique_num,1)]
unique_key2_df = pd.DataFrame({key2:unique_key2_list,'Group_num':num_list})
merge_df = pd.merge(part_data,unique_key2_df,on=key2,how='left')
merge_df.to_csv(ofile,sep='\t',header=True,index=False)

python /work/workspace/zhurj/script/2_metapro/metapipe/bar_cluster_meta_creation.py -i /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_all -o /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_bar_cluster --keylist SampleId Group


R
library(vegan)
library(RColorBrewer)
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#mycolors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
library(ggsci)

infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy/new_cazy_abun'
metaf = '/work/workspace/zhurj/project/1_metadata/mouse22/input/meta_bar_cluster'
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/bar_cluster_cazy.png'
oname = 'Cazy level1'
otitle = paste0(oname,' relative abundance')

mycolors <- unique(c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                      pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
                      pal_locuszoom("default")(7),pal_igv("default")(51),
                      pal_uchicago("default")(9),pal_startrek("uniform")(7),
                      pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
                      pal_simpsons("springfield")(16),pal_gsea("default")(12)))


#mycolors = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")


predata <- read.delim(infile,row.names = 1, sep = '\t', header = TRUE, check.names = FALSE) 
barx = -0.1 * (max(colSums(predata)))
predata$sum = rowSums(predata)
tmpdata <- predata[order(-predata$sum),]
data <- subset(tmpdata, select = -c(sum))

data.dist <- vegdist(t(data),method='bray')
tree <- hclust(data.dist,method='average')

group <- read.delim(metaf,sep='\t',row.names = 1,header=TRUE,check.names = FALSE,  stringsAsFactors = FALSE)
grp <- group[2]
unique_num <- nrow(unique(group[2]))
group_col <- mycolors[1:unique_num]
names(group_col) <- as.vector(as.matrix(unique(group[2])))
group_name <- as.vector(as.matrix(unique(group[1])))

#样本分组标签
png(ofile,width = 1200, height = 500, units = "px")
layout(t(c(1, 1, 2, 2, 2, 3)))
par(mar = c(5, 2, 0, 0))

plot(0, type = 'n', hang=-1, xaxt = 'n', yaxt = 'n', frame.plot = FALSE, xlab = '', ylab = '', xlim = c(-max(tree$height), 0), ylim = c(0, length(tree$order)))
legend('topleft', legend = group_name, pch = 15, col = group_col, bty = 'n', cex = 1)

#聚类树绘制，按分组给分支上色
treeline <- function(pos1, pos2, height, col1, col2) {
    meanpos = (pos1[1] + pos2[1]) / 2
    segments(y0 = pos1[1] - 0.4, x0 = -pos1[2], y1 = pos1[1] - 0.4, x1 = -height,  col = col1,lwd = 2)
    segments(y0 = pos1[1] - 0.4, x0 = -height,  y1 = meanpos - 0.4, x1 = -height,  col = col1,lwd = 2)
    segments(y0 = meanpos - 0.4, x0 = -height,  y1 = pos2[1] - 0.4, x1 = -height,  col = col2,lwd = 2)
    segments(y0 = pos2[1] - 0.4, x0 = -height,  y1 = pos2[1] - 0.4, x1 = -pos2[2], col = col2,lwd = 2)
}
 
meanpos = matrix(rep(0, 2 * length(tree$order)), ncol = 2)
meancol = rep(0, length(tree$order))
for (step in 1:nrow(tree$merge)) {
    if(tree$merge[step, 1] < 0){
        pos1 <- c(which(tree$order == -tree$merge[step, 1]), 0)
        col1 <- group_col[as.character(grp[tree$labels[-tree$merge[step, 1]],1])]
    } else {
        pos1 <- meanpos[tree$merge[step, 1], ]
        col1 <- meancol[tree$merge[step, 1]]
    }
    if (tree$merge[step, 2] < 0) {
        pos2 <- c(which(tree$order == -tree$merge[step, 2]), 0)
        col2 <- group_col[as.character(grp[tree$labels[-tree$merge[step, 2]],1])]
    } else {
        pos2 <- meanpos[tree$merge[step, 2], ]
        col2 <- meancol[tree$merge[step, 2]]
    }
    height <- tree$height[step]
    treeline(pos1, pos2, height, col1, col2)
    meanpos[step, ] <- c((pos1[1] + pos2[1]) / 2, height)
    if (col1 == col2) meancol[step] <- col1 else meancol[step] <- 'grey'
}

##堆叠柱形图
#样本顺序调整为和聚类树中的顺序一致
data <- data[ ,tree$order]
rownum = nrow(data)
 
#物种颜色设置
item_color <- mycolors[1:rownum]
names(item_color) <- rownames(data)
 
#堆叠柱形图
par(mar = c(5, 6, 0, 0))

bar <- barplot(as.matrix(data), col = item_color, space = 0.4, width = 0.7, cex.axis = 1, horiz = TRUE, cex.lab = 1.2, xlab = otitle, yaxt = 'n', las = 1, ylim = c(0, ncol(data)), family = 'mono')
text(x = barx, y = bar, labels = colnames(data), col = group_col[group[tree$order, 2]], xpd = TRUE)
#mtext('cazy level1 relative abundance', side = 3, line = 1, cex = 1)

#柱形图图例
par(mar = c(5, 1, 0, 0))
plot(0, type = 'n', xaxt = 'n', yaxt = 'n', bty = 'n', xlab = '', ylab = '')
legend('topleft', pch = 15, col = item_color, legend = names(item_color), bty = 'n', cex = 1)

dev.off()

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/bar_cluster.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy/new_cazy_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/bar_cluster_cazy.png -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_bar_cluster -n 'Cazy level1'

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/bar_cluster.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/card_Drug_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/bar_cluster_card_drug.png -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_bar_cluster -n 'Card drug'

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/bar_cluster.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/card_resistance_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/bar_cluster_card_resistance.png -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_bar_cluster -n 'Card resistance'

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/bar_cluster.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/eggnog/new_eggnog_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/bar_cluster_eggnog_level1.png -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_bar_cluster -n 'Eggnog level1'

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/bar_cluster.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/new_kegg_levelA_abundance -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/bar_cluster_kegg_levelA.png -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_bar_cluster -n 'Kegg levelA'

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/bar_cluster.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/new_kegg_levelB_abundance -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/bar_cluster_kegg_levelB.png -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_bar_cluster -n 'Kegg levelB' --plot_height 1000

# kegg pathway + annotation, create file for bar cluster plot 20201203
import pandas as pd
import os
import sys

infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelA_abundance'
metaf = '/work/workspace/zhurj/database/kegg/20201029/path_annotation'
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/new_kegg_levelA_abundance'

data = pd.read_csv(infile,sep='\t',header=0)
ref_df = pd.read_csv(metaf,sep='\t',header=None)
ref_df.columns = ['path','anno']
merge_df = pd.merge(data,ref_df,on='path',how='left')

sam_cols = [x for x in data.columns if 'MN' in x]
selected_cols = ['anno'] + sam_cols
result_df = merge_df[selected_cols]
result_df.to_csv(ofile,sep='\t',header=True,index=False)

python /work/workspace/zhurj/script/2_metapro/metapipe/kegg_abun_add_anno.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelA_abundance -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/new_kegg_levelA_abundance

python /work/workspace/zhurj/script/2_metapro/metapipe/kegg_abun_add_anno.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelB_abundance -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/new_kegg_levelB_abundance

# lefse 差异分析 20201203
# Lefse in
source activate /work/workspace/zhurj/bin/miniconda3/envs/py2
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelA_abundance /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelA_abundance_lefsein
lefse-format_input.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelA_abundance_lefsein /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelA_abundance_lefseformatin -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelA_abundance_lefseformatin /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelA_abundance_lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse_fitler.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelA_abundance_lefse_run.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelA_abundance_lefse_run_filter.res

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelB_abundance /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelB_abundance_lefsein
lefse-format_input.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelB_abundance_lefsein /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelB_abundance_lefseformatin -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelB_abundance_lefseformatin /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelB_abundance_lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse_fitler.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelB_abundance_lefse_run.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelB_abundance_lefse_run_filter.res

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelC_abundance /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelC_abundance_lefsein
lefse-format_input.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelC_abundance_lefsein /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelC_abundance_lefseformatin -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelC_abundance_lefseformatin /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelC_abundance_lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1
lefse-plot_res.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelC_abundance_lefse_run.res  /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelC_abundance_lefse_run.png --feature_font_size 8 --width 10 --format png
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse_fitler.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelC_abundance_lefse_run.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelC_abundance_lefse_run_filter.res

source deactivate /work/workspace/zhurj/bin/miniconda3/envs/py2


# lefse 差异分析 HFD vs MNO863
# 提取 HFD and MNO863
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/extract_cols.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelA_abundance -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelA_dg_abundance -r MNA00233 MNA00234 MNA00235 MNA00236 MNA00237 MNA00238
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/extract_cols.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelB_abundance -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelB_dg_abundance -r MNA00233 MNA00234 MNA00235 MNA00236 MNA00237 MNA00238
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/extract_cols.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelC_abundance -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelC_dg_abundance -r MNA00233 MNA00234 MNA00235 MNA00236 MNA00237 MNA00238

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/extract_cols.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/card_Drug_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/card_Drug_dg_abun -r MNA00233 MNA00234 MNA00235 MNA00236 MNA00237 MNA00238
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/extract_cols.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/card/card_resistance_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/card_resistance_dg_abun -r MNA00233 MNA00234 MNA00235 MNA00236 MNA00237 MNA00238


/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/extract_cols.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy/cazy_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_dg_abun -r anno MNA00233 MNA00234 MNA00235 MNA00236 MNA00237 MNA00238
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/extract_cols.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/cazy/level2_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_level2_dg_abun -r MNA00233 MNA00234 MNA00235 MNA00236 MNA00237 MNA00238


/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/extract_cols.py -i/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/eggnog/eggnog_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/eggnog_dg_abun -r anno MNA00233 MNA00234 MNA00235 MNA00236 MNA00237 MNA00238
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/extract_cols.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/eggnog/level2_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/eggnog_level2_dg_abun -r MNA00233 MNA00234 MNA00235 MNA00236 MNA00237 MNA00238

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/extract_cols.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/ko_abundance -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/ko_dg_abundance -r MNA00233 MNA00234 MNA00235 MNA00236 MNA00237 MNA00238


source activate py2
# kegg level1
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelA_dg_abundance /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelA_dg_abundance_lefsein
lefse-format_input.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelA_dg_abundance_lefsein /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelA_dg_abundance_lefseformatin -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelA_dg_abundance_lefseformatin /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelA_dg_abundance_lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse_fitler.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelA_dg_abundance_lefse_run.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelA_dg_abundance_lefse_run_filter.res
# 差异信息的heatmap
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/toheadmapin.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelA_dg_abundance /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelA_dg_abundance_lefse_run_filter.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelA_dg_preheatmap
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/kegg_abun_add_pathanno.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelA_dg_preheatmap -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelA_dg_heatmapin

# kegg level2
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelB_dg_abundance /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelB_dg_abundance_lefsein
lefse-format_input.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelB_dg_abundance_lefsein /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelB_dg_abundance_lefseformatin -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelB_dg_abundance_lefseformatin /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelB_dg_abundance_lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse_fitler.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelB_dg_abundance_lefse_run.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelB_dg_abundance_lefse_run_filter.res
# 差异信息的heatmap
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/toheadmapin.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelB_dg_abundance /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelB_dg_abundance_lefse_run_filter.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelB_dg_preheatmap
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/kegg_abun_add_pathanno.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelB_dg_preheatmap -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelB_dg_heatmapin

# kegg level3
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelC_dg_abundance /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelC_dg_abundance_lefsein
lefse-format_input.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelC_dg_abundance_lefsein /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelC_dg_abundance_lefseformatin -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelC_dg_abundance_lefseformatin /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelC_dg_abundance_lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse_fitler.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelC_dg_abundance_lefse_run.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelC_dg_abundance_lefse_run_filter.res
# 差异信息的heatmap
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/toheadmapin.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelC_dg_abundance /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelC_dg_abundance_lefse_run_filter.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelC_dg_preheatmap
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/kegg_abun_add_pathanno.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelC_dg_preheatmap -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelC_dg_heatmapin

# CARD DRUG
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/card_Drug_dg_abun /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/card_Drug_dg_abun_lefsein
lefse-format_input.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/card_Drug_dg_abun_lefsein /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/card_Drug_dg_abun_lefseformatin -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/card_Drug_dg_abun_lefseformatin /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/card_Drug_dg_abun_lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse_fitler.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/card_Drug_dg_abun_lefse_run.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/card_Drug_dg_abun_lefse_run_filter.res
# 无差异card drug

# card resistance
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/card_resistance_dg_abun /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/card_resistance_dg_abun_lefsein
lefse-format_input.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/card_resistance_dg_abun_lefsein /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/card_resistance_dg_abun_lefseformatin -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/card_resistance_dg_abun_lefseformatin /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/card_resistance_dg_abun_lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse_fitler.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/card_resistance_dg_abun_lefse_run.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/card_resistance_dg_abun_lefse_run_filter.res
# 无差异card resistance

# cazy_dg_abun
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_dg_abun /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_dg_abun_lefsein
lefse-format_input.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_dg_abun_lefsein /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_dg_abun_lefseformatin -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_dg_abun_lefseformatin /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_dg_abun_lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse_fitler.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_dg_abun_lefse_run.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_dg_abun_lefse_run_filter.res
# 差异信息的heatmap
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/toheadmapin.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_dg_abun /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_dg_abun_lefse_run_filter.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_dg_preheatmap
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/cazy_abun_add_anno.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_dg_preheatmap -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_dg_heatmapin

# cazy level 2
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_level2_dg_abun /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_level2_dg_abun_lefsein
lefse-format_input.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_level2_dg_abun_lefsein /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_level2_dg_abun_lefseformatin -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_level2_dg_abun_lefseformatin /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_level2_dg_abun_lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse_fitler.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_level2_dg_abun_lefse_run.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_level2_dg_abun_lefse_run_filter.res
# 差异信息的heatmap
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/toheadmapin.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_level2_dg_abun /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_level2_dg_abun_lefse_run_filter.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_level2_dg_heatmapin

# eggnog level 1
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/eggnog_dg_abun /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/eggnog_dg_abun_lefsein
lefse-format_input.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/eggnog_dg_abun_lefsein /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/eggnog_dg_abun_lefseformatin -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/eggnog_dg_abun_lefseformatin /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/eggnog_dg_abun_lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse_fitler.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/eggnog_dg_abun_lefse_run.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/eggnog_dg_abun_lefse_run_filter.res
# 差异信息的heatmap
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/toheadmapin.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/eggnog_dg_abun /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/eggnog_dg_abun_lefse_run_filter.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/eggnog_level1_dg_preheatmap
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/item_abun_add_pathanno.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/eggnog_level1_dg_preheatmap -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/eggnog_level1_dg_heatmapin -m /work/workspace/zhurj/database/eggnog/e5/per_tax_level/level1_anno

# ko


source deactivate py2


library(pheatmap)
library(ggsci)

mycolors <- unique(c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                      pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
                      pal_locuszoom("default")(7),pal_igv("default")(51),
                      pal_uchicago("default")(9),pal_startrek("uniform")(7),
                      pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
                      pal_simpsons("springfield")(16),pal_gsea("default")(12)))

infile <- '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelA_dg_heatmapin'
ref_col <- '/work/workspace/zhurj/project/1_metadata/mouse22/input/dg_group'
intitle <- 'Dif KEGG level1'
figwidth = 10
figheight = 5

raw_df <-read.table(infile, sep='\t', header=TRUE, row.names = 1, encoding='utf-8')
raw_mt <- data.matrix(raw_df)
ref_col_df <- read.table(ref_col,sep='\t',header=FALSE,row.names = 1, encoding='utf-8')
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/dg_group/kegg_levelA_dg_heatmapin.pdf'
colnames(ref_col_df) = 'Group'
group_v = as.vector(unique(ref_col_df$Group))
select_col = mycolors[1:length(group_v)]
names(select_col) = group_v

#pheatmap(raw_mt)
# #main可设置热图的标题，fontsize设置字体大小，filename可直接将热图存出，支持格式png, pdf, tiff, bmp, jpeg，并且可以通过width, height设置图片的大小；
ann_colors = list(Group=select_col)
pheatmap(raw_mt, main=intitle, annotation_col = ref_col_df, angle_col = 45, annotation_colors = ann_colors, scale='row', cluster_col = FALSE, cutree_rows=2, filename = ofile, width = figwidth, height = figheight)


/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/dif_pheatmap.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelA_dg_heatmapin -r /work/workspace/zhurj/project/1_metadata/mouse22/input/dg_group -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/dg_group/kegg_levelA_dg_heatmapin.pdf -n 'KEGG levelA' --plot_width 15 --plot_height 2
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/dif_pheatmap.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelB_dg_heatmapin -r /work/workspace/zhurj/project/1_metadata/mouse22/input/dg_group -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/dg_group/kegg_levelB_dg_heatmapin.pdf -n 'KEGG levelB' --plot_width 15 --plot_height 6
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/dif_pheatmap.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/kegg_levelC_dg_heatmapin -r /work/workspace/zhurj/project/1_metadata/mouse22/input/dg_group -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/dg_group/kegg_levelC_dg_heatmapin.pdf -n 'KEGG levelC' --plot_width 15 --plot_height 12

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/dif_pheatmap.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_dg_heatmapin -r /work/workspace/zhurj/project/1_metadata/mouse22/input/dg_group -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/dg_group/cazy_dg_heatmapin.pdf -n 'Cazy level1' --plot_width 15 --plot_height 2

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/dif_pheatmap.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/cazy_level2_dg_heatmapin -r /work/workspace/zhurj/project/1_metadata/mouse22/input/dg_group -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/dg_group/cazy_level2_dg_heatmapin.pdf -n 'Cazy level2' --plot_width 15 --plot_height 15

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/dif_pheatmap.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/eggnog_level1_dg_heatmapin -r /work/workspace/zhurj/project/1_metadata/mouse22/input/dg_group -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/dg_group/eggnog_level1_dg_heatmapin.pdf -n 'Eggnog level1' --plot_width 15 --plot_height 6


# ko list
import pandas as pd
import sys
import os

infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/ko_dg_abun_lefse_run_filter.res'
odir = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg'
name = 'dif_ko'
filename = name + '_all'
ofile = os.path.join(odir, filename)

data = pd.read_csv(infile,sep='\t',header=None,index_col=None)
data.columns = ['item','v1','group','foldchange','pvalue']
data['item'].drop_duplicates().to_csv(ofile,sep='\t',header=False,index=False)

unique_group = list(data['group'].unique())
for groupone  in unique_group:
  filename = "{}_{}".format(name,groupone)
  ofile = os.path.join(odir,filename)
  filter_df = data[data['group'] == groupone]
  filter_df['item'].drop_duplicates().to_csv(ofile,sep='\t', header=False, index=False)

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_to_ko_list.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/ko_dg_abun_lefse_run_filter.res -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg -n dif_ko



# enrichment kegg pathway - 20201204
source activate R3.6
# diff ko enriched pathway
library(clusterProfiler)
infile <- '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/dif_ko_all'
ofile <- '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/dif_ko_all_enrichedkegg'
path_kof <- '/work/workspace/zhurj/database/kegg/20201029/gene_path_levelC'
path_annf <- '/work/workspace/zhurj/database/kegg/20201029/path_annotation'

gene_list <- read.table(infile)
gene <- gene_list[,1]
background <- read.table(path_kof, sep="\t" ,header=FALSE)
kegg2name_data <- read.table(path_annf, sep="\t" ,header=FALSE)
result <- enricher(gene,TERM2GENE=background,TERM2NAME=kegg2name_data,pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH")
write.table(result, ofile, sep = '\t', quote = FALSE, col.names = TRUE)

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/from_ko_to_enriched_pathway.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/dif_ko_all -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/dif_ko_all_enrichedkegg

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/from_ko_to_enriched_pathway.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/dif_ko_MNO863 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/dif_ko_MNO863_enrichedkegg

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/from_ko_to_enriched_pathway.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/dif_ko_HFD -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/dif_ko_HFD_enrichedkegg


/work/workspace/liangzj/project/multi_omics/MNO0863mice/metagenomics_metabolome/run.sh
/work/workspace/liangzj/project/multi_omics/MNO0863mice/metagenomics_metabolome/kegg_pathway

source activate R4.0
library(pathview)
infile = '/work/workspace/zhurj/project/1_metadata/mouse22/humann3_uniref50_all/enrichment/ko/ko03001_in'
data <- read.table(infile, sep="\t",header=TRUE,row.names = 1)
p <- pathview(gene.data = data, pathway.id='03001',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T)


import pandas as pd
import argparse
import os
import sys
import math

infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/dif_ko_all' 
abunf = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/ko_dg_abun'
metaf = '/work/workspace/zhurj/project/1_metadata/mouse22/input/dg_group'
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/dif_ko_foldchange'


# fold change

meta_df = pd.read_csv(metaf,sep='\t',header=None,index_col=None)
meta_df.columns = ['sample','group']

ko_list_df = pd.read_csv(infile,sep='\t',header=None,index_col=None,names=['koid'])

data = pd.read_csv(abunf,sep='\t',header=0,index_col = 0)
not_zero_min = data[data > 0].min().min()
multiply_index = math.pow(10,math.ceil(math.log10(1/not_zero_min)) + 2)
data = data * multiply_index + 1
data['count'] = (data > 1).astype(int).sum(axis=1)
data_filter_df = data[data['count'] > 2]
data_t = data_filter_df.T
item_cols = data_t.columns.values.tolist()
data_t['sample'] = data_t.index
merge_df = pd.merge(data_t,meta_df,on='sample',how='inner')
select_cols = ['group'] + item_cols
part_merge_df = merge_df[select_cols]
foldchange_df = part_merge_df.groupby('group').mean().T
key_list = foldchange_df.columns.values.tolist()
foldchange_df['Fc'] = foldchange_df[key_list[1]] / foldchange_df[key_list[0]]
foldchange_df['logFC'] = foldchange_df['Fc'].apply(lambda x: math.log2(x))
foldchange_df['FLAG'] = 0
foldchange_df['koid'] = foldchange_df.index
part_foldchange_df = foldchange_df[['koid','logFC','FLAG']]
fc_merge_df = pd.merge(foldchange_df,ko_list_df,on='koid',how='inner')
result_cols = ['koid','logFC','FLAG']
result_df = fc_merge_df[result_cols]

result_df.to_csv(ofile,sep='\t',header=True,index=False)


/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/ko_foldchange.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/dif_ko_all --metaf /work/workspace/zhurj/project/1_metadata/mouse22/input/dg_group -r /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/ko_dg_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/dif_ko_foldchange


R4.0

library(pathview)
library(dplyr)

path_creation <- function(x,f_dir,f_data) {
  pathid <- substr(x,3,nchar(x))
  pathview(gene.data = f_data, pathway.id = pathid, species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T, kegg.dir = f_dir)
}

enrichpathf <- '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/dif_ko_MNO863_enrichedkegg'
infile <- '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/dif_ko_foldchange'
odir <- '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/pathway/MNO863'
enrichpathf <- '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/test_path'

data <- read.table(infile, sep="\t",header=TRUE,row.names = 1)
path_df <- read.table(enrichpathf,sep='\t',header=TRUE,row.names = 1)
part_path_df <- filter(path_df,grepl('PATH:ko',Description))

lapply(part_path_df$ID, path_creation, f_dir = odir, f_data = data)



#pathview(gene.data = data, pathway.id='00970',species = "ko", gene.idtype = "kegg", keys.align = "y", kegg.native = T, multi.state = T, same.layer = T, kegg.dir = odir)

cd /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/pathway/MNO863 && /work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/pathview.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/dif_ko_foldchange -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/pathway/MNO863 -p /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/test_path
cd /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/pathway/MNO863 && /work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/pathview.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/dif_ko_foldchange -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/pathway/MNO863 -p /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/dif_ko_MNO863_enrichedkegg

cd /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/pathway/HFD && /work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/pathview.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/dif_ko_foldchange -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/pic/pathway/HFD -p /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/lefse_dg/dif_ko_HFD_enrichedkegg



# ko reference
/work/database/kegg/index/2019.Mar27/ko.compress.faa

# number of gene detected in each sample
import pandas as pd
import re
import os


def gene_taxo_define(infile_,ofile_):
  con = ''
  with open(ofile_,'w', encoding='utf-8') as fo, open(infile_,'r',encoding='utf-8') as fp:
    con = "Gene\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n"
    fo.write(con)
    line = fp.readline()
    while line:
      tmplist = []
      #line = re.sub(r"\s+","", line, flags=re.UNICODE)
      tmplist = line.strip().split(';')
      line = fp.readline()
      geneid = tmplist[0]
      domain = tmplist[2]
      domain_score = int(tmplist[3].strip())
      phylum = tmplist[4]
      phylum_score = int(tmplist[5].strip())
      class = tmplist[6]
      class_score = int(tmplist[7].strip())
      order = tmplist[8]
      order_score = int(tmplist[9].strip())
      family = tmplist[10]
      family_score = int(tmplist[11].strip())
      genus = tmplist[12]
      genus_score = int(tmplist[13].strip())
      species = tmplist[14]
      species_score = int(tmplist[15].strip())
      if domain_score < 50:
        con = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(geneid,'','','','','','','')
        fo.write(con)
        continue
      if phylum_score < 50:
        con = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(geneid,domain,'','','','','','')
        fo.write(con)
        continue
      if class_score < 50:
        con = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(geneid,domain,phylum,'','','','','')
        fo.write(con)
        continue
      if order_score < 50:
        con = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(geneid,domain,phylum,class,'','','','')
        fo.write(con)
        continue
      if family_score < 50:
        con = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(geneid,domain,phylum,class,order,'','','')
        fo.write(con)
        continue
      if genus_score < 50:
        con = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(geneid,domain,phylum,class,order,family,'','')
        fo.write(con)
        continue
      if species_score < 50:
        con = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(geneid,domain,phylum,class,order,family,genus,'')
        fo.write(con)
      else:
        con = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(geneid,domain,phylum,class,order,family,genus,species)
        fo.write(con)

def levels_abun(data_df_,sam_cols_,type_):
  level_cols = [type_] + sam_cols_
  level_df = data_df_[level_cols].dropna(how = 'any')
  level_abun_df = data_df_.groupby(type_).sum()
  name_ = type_ + '_abun'
  ofile = os.path.join(odir,name_)
  level_abun_df.to_csv(ofile,sep='\t',header=True,index=True)

def levels_gene_count(data_df_,sam_cols_,type_):
  level_cols = [type_] + sam_cols_
  level_df = data_df_[level_cols].dropna(how = 'any')
  level_gene_count_df = data_df_.groupby(type_).agg(lambda x: x.gt(0).sum())
  name_ = type_ + '_gene_count'
  ofile = os.path.join(odir,name_)
  level_gene_count_df.to_csv(ofile,sep='\t',header=True,index=True)


#sentence = re.sub(r"\s+","", sentence, flags=re.UNICODE)

abunf = '/work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/mergeAbundance'
in_taxof = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/all_blast2lca.tab'
odir = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond'
genetaxof = os.path.join(odir,'gene_taxo')
gene_taxo_abunf = os.path.join(odir,'gene_taxo_abun')
ofile = ''
#gene_taxo_define(in_taxof,genetaxof)

abun_df = pd.read_csv(abunf,sep='\t',header=0,index_col=None)
genetaxo_df = pd.read_csv(genetaxof,sep='\t',header=0,index_col=None)
merge_df = pd.merge(genetaxo_df,abun_df,left_on='Gene',right_on='Sample', how='inner')
select_cols = genetaxo_df.columns.values.tolist() + [x for x in merge_df.columns if 'MN' in x]
part_merge_df = merge_df[select_cols]
part_merge_df.to_csv(gene_taxo_abunf,sep='\t',header=True,index=False)

sam_cols = [x for x in merge_df.columns if 'MN' in x]
levels_abun(merge_df,sam_cols,'Species')
levels_abun(merge_df,sam_cols,'Genus')
levels_abun(merge_df,sam_cols,'Family')
levels_abun(merge_df,sam_cols,'Order')
levels_abun(merge_df,sam_cols,'Class')
levels_abun(merge_df,sam_cols,'Phylum')
levels_abun(merge_df,sam_cols,'Domain')

levels_gene_count(merge_df,sam_cols,'Species')
levels_gene_count(merge_df,sam_cols,'Genus')
levels_gene_count(merge_df,sam_cols,'Family')
levels_gene_count(merge_df,sam_cols,'Order')
levels_gene_count(merge_df,sam_cols,'Class')
levels_gene_count(merge_df,sam_cols,'Phylum')
levels_gene_count(merge_df,sam_cols,'Domain')



/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/taxo_from_gene.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/all_blast2lca.tab -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond --abunf /work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/mergeAbundance --cutoff 50


287185 + 288806 + 285728 = 861719les 

krona in creation - 20201209
import pandas as pd
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import os
import subprocess


infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/gene_taxo_abun'
data = pd.read_csv(infile,sep='\t',header=0,index_col=0)
odir = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/krona/'
"""
tmpfile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/krona/input/test'
krona_in = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/krona/input/MNA00233_in'
"""
tmpfile = os.path.join(odir,'tmpfile')


kronapro = '/work/workspace/zhurj/bin/miniconda3/envs/metaAna/bin/ktImportText'
con = ''

sam_cols = [x for x in data.columns if 'MN' in x]
name_cols = data.columns.values.tolist()
level_cols = [x for x in name_cols if x not in sam_cols]

for sam in sam_cols:
  name = sam + '_in'
  krona_in = os.path.join(odir,name)
  name = sam + '.krona.html'
  krona_of = os.path.join(odir,name)
  exam_cols = [sam] + level_cols
  sam_df = data[exam_cols]
  sam_filter_df = sam_df[ sam_df[sam] > 0 ]
  sam_filter_df.dropna(subset=['Domain'],inplace=True)
  sam_filter_df.fillna('',inplace=True)
  sam_group_df = sam_filter_df.groupby(level_cols).sum()
  result_df = sam_group_df[sam_group_df[sam] >= 0.0001]
  result_df.to_csv(tmpfile,sep='\t',header=True,index=True)
  
  re_data_df = pd.read_csv(tmpfile,sep='\t',header=0,index_col=None)
  final_df = re_data_df[exam_cols]
  final_df.fillna('',inplace=True)
  final_df.to_csv(krona_in,sep='\t',header=False,index=False)
  
  con = "{} {} -o {} -n {} ".format(kronapro,krona_in,krona_of,'root')
  subprocess.call(con,shell=True)


/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/krona_taxo_in.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/gene_taxo_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/krona


#/work/workspace/zhurj/bin/miniconda3/envs/metaAna/bin/ktImportText /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/krona/input/MNA00233_in -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/krona/output/MNA00233_krona.html -n root


taxo bar graph 

taxo bar graph input -- sample

import pandas as pd
import numpy as np
import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

# normlization to sum = 1 for each column
infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Phylum_abun'
metaf = '/work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group'
odir = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar'
top_n = 10
name = 'phylum' 
ofile = filename = ''

df = pd.read_csv(infile,sep='\t',header=0, index_col=0)
col_names = df.columns.values.tolist()

# sort according to sum of each row
norm_df = df.apply(lambda x: x/x.sum(), axis=0)
norm_df['sum'] = norm_df.apply(lambda x: x.sum(),axis=1)
sort_df = norm_df.sort_values(by='sum', ascending=False)
part_df = sort_df[col_names]
topn_df = part_df.head(top_n)
result_df = topn_df.copy()
result_df.loc['others'] = result_df.apply(lambda x: 1- x.sum(), axis=0)
filename =  "{}{}{}{}".format('barin_', name, '_top', top_n)
ofile = os.path.join(odir,filename)
result_df.to_csv(ofile,sep='\t',header=True,index=True)

# barin_phylum_top10

# group in
# variable taxonomy value group
meta_df = pd.read_csv(metaf,sep='\t',header=None, index_col=None)
meta_df.columns = ['variable','group']
pre_df = pd.DataFrame(result_df.stack())
pre_df.columns = ['value']
pre_df['tmp'] = pre_df.index
pre_df['variable'] = pre_df['tmp'].apply(lambda x: x[1])
pre_df['Taxonomy'] = pre_df['tmp'].apply(lambda x: x[0])
select_cols = ['variable','Taxonomy','value','group']
merge_df = pd.merge(pre_df,meta_df,on='variable',how='left')
result_df = merge_df[select_cols]
filename = "{}{}{}".format('barin_',name,'_group')
ofile = os.path.join(odir,filename)
result_df.to_csv(ofile,sep='\t',header=True,index=False)
taxo_cols = list(result_df['Taxonomy'].drop_duplicates())
#--metaf /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group 

# group taxonomy abundance 
pre_group_col = ['Taxonomy','value','group']
group_part_df = result_df[pre_group_col]
group_groupby_df = group_part_df.groupby(['Taxonomy','group']).mean()
group_groupby_df['tmp'] = group_groupby_df.index
group_groupby_df['group'] = group_groupby_df['tmp'].apply(lambda x: x[1])
group_groupby_df['Taxonomy'] = group_groupby_df['tmp'].apply(lambda x: x[0])
group_result_df = group_groupby_df[pre_group_col]
rownum = group_result_df.shape[0]
group_result_df.index = [x for x in range(0,rownum)]
group_result_df['Taxonomy'] = group_result_df['Taxonomy'].astype('category')
group_result_df['Taxonomy'].cat.reorder_categories(taxo_cols,inplace=True)
group_result_df.sort_values('Taxonomy',inplace=True)
filename = "{}{}{}".format('bargroup_',name,'_in')
ofile = os.path.join(odir,filename)
group_result_df.to_csv(ofile,sep='\t',header=True,index=False)
# 

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/bar_taxo_sample_in.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Phylum_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar --name phylum --metaf /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/bar_taxo_sample_in.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Class_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar --name class --metaf /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/bar_taxo_sample_in.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Order_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar --name order --metaf /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/bar_taxo_sample_in.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Family_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar --name family --metaf /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/bar_taxo_sample_in.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Genus_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar --name genus --metaf /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/bar_taxo_sample_in.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Species_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar --name species --metaf /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group

R
library(ggsci)
library(reshape2)
library(ggplot2)


mycolors <- unique(c(pal_gsea("default")(12),pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                      pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
                      pal_locuszoom("default")(7),pal_igv("default")(51),
                      pal_uchicago("default")(9),pal_startrek("uniform")(7),
                      pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
                      pal_simpsons("springfield")(16)))


#mycolors = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")


#metaf = '/work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group'
infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/barin_phylum_group'
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/bar_graph_phylum.png'

data <- read.delim(infile, stringsAsFactors = FALSE, check.names = FALSE, header = TRUE, sep = "\t")
taxo_num = length(unique(data$Taxonomy))
group_num = length(unique(data$group))
taxo_col = mycolors[1:taxo_num]
data$Taxonomy = factor(data$Taxonomy,levels=rev(unique(data$Taxonomy)))
p <- ggplot(data, aes(variable, 100 * value, fill = Taxonomy)) +
geom_col(position = 'stack', width = 0.6) +
labs(x = 'Sample Name', y = 'Relative Abundance(%)') +
theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
theme(legend.text = element_text(size = 11)) +
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
theme(legend.title = element_blank()) +
scale_fill_manual(values =  rev(taxo_col)) +
facet_wrap(~group, scales = 'free_x', ncol = group_num) +
theme(strip.text = element_text(size = 12))

ggsave(ofile, p, width = 10, height = 6)

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/bar_taxo_sample.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/barin_phylum_group -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/bar_graph_phylum.png
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/bar_taxo_sample.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/barin_class_group -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/bar_graph_class.png
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/bar_taxo_sample.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/barin_order_group -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/bar_graph_order.png
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/bar_taxo_sample.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/barin_family_group -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/bar_graph_family.png
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/bar_taxo_sample.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/barin_genus_group -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/bar_graph_genus.png
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/bar_taxo_sample.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/barin_species_group -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/bar_graph_species.png


R group
library(ggsci)
library(reshape2)
library(ggplot2)
library(dplyr)


mycolors <- unique(c(pal_gsea("default")(12),pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                      pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
                      pal_locuszoom("default")(7),pal_igv("default")(51),
                      pal_uchicago("default")(9),pal_startrek("uniform")(7),
                      pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
                      pal_simpsons("springfield")(16)))


#mycolors = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")


#metaf = '/work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group'
infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/bargroup_phylum_in'
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/bar_group_phylum.png'

data <- read.delim(infile, stringsAsFactors = FALSE, check.names = FALSE, header = TRUE, sep = "\t")
sort_col <- c('Taxonomy','value')
sort_df <- data[sort_col]
sort_group_df <- groupby(sort_df,'Taxonomy').sum()
taxo_num = length(unique(data$Taxonomy))
group_num = length(unique(data$group))
taxo_col = mycolors[1:taxo_num]
data$Taxonomy = factor(data$Taxonomy,levels=rev(unique(data$Taxonomy)))
p <- ggplot(data, aes(variable, 100 * value, fill = Taxonomy)) +
geom_col(position = 'stack', width = 0.6) +
labs(x = 'Sample Name', y = 'Relative Abundance(%)') +
theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
theme(legend.text = element_text(size = 11)) +
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
theme(legend.title = element_blank()) +
scale_fill_manual(values =  rev(taxo_col)) +
facet_wrap(~group, scales = 'free_x', ncol = group_num) +
theme(strip.text = element_text(size = 12))

ggsave(ofile, p, width = 10, height = 6)

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/bar_taxo_group.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/bargroup_phylum_in -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/bar_group_phylum.png --plot_height 6 --plot_width 8
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/bar_taxo_group.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/bargroup_class_in -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/bar_group_class.png --plot_height 6 --plot_width 8
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/bar_taxo_group.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/bargroup_order_in -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/bar_group_order.png --plot_height 6 --plot_width 8
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/bar_taxo_group.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/bargroup_family_in -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/bar_group_family.png --plot_height 6 --plot_width 8
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/bar_taxo_group.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/bargroup_genus_in -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/bar_group_genus.png --plot_height 6 --plot_width 8
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/bar_taxo_group.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/bargroup_species_in -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_bar/bar_group_species.png --plot_height 6 --plot_width 8


top35, gene abundance 
heatmap in - 20201211
import pandas as pd
import numpy as np
import argparse
import os
import sys
import math

infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Phylum_abun' 
odir = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap'
taxonomy = 'phylum'
topn = 35
filename = "{}{}{}".format(taxonomy,'_geneabun_heatmapin_top',topn)
ofile = os.path.join(odir,filename)

# normalize the relative abundance
#  log10(data * match.pow(10,magnitude(1/min(all data > 0) + 2) + 1)
data = pd.read_csv(infile,sep='\t',header=0,index_col = 0)
not_zero_min = data[data > 0].min().min()
multiply_index = math.pow(10,math.ceil(math.log10(1/not_zero_min)) + 2)
data = data * multiply_index + 1
norm_df = data.apply(np.log10)
sam_cols = norm_df.columns.values.tolist()

# topn
tmp_df = norm_df
tmp_df['sum'] = tmp_df.apply(lambda x: x.sum(),axis=1)
sort_df = tmp_df.sort_values('sum',ascending=False)
part_sort_df = sort_df.head(topn)
result_df = part_sort_df[sam_cols]
result_df.to_csv(ofile,sep='\t',header=True,index=True)

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/heatmapin_abun_topn.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Phylum_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap --name phylum -n 35
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/heatmapin_abun_topn.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Class_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap --name class -n 35
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/heatmapin_abun_topn.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Order_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap --name order -n 35
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/heatmapin_abun_topn.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Family_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap --name family -n 35
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/heatmapin_abun_topn.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Genus_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap --name genus -n 35
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/heatmapin_abun_topn.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Species_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap --name species -n 35


top35, gene count
heatmap in - 20201211
import pandas as pd
import numpy as np
import argparse
import os
import sys
import math

infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Phylum_gene_count' 
odir = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap'
taxonomy = 'phylum'
topn = 35
filename = "{}{}{}".format(taxonomy,'_genecount_heatmapin_top',topn)
ofile = os.path.join(odir,filename)

# normalize the relative abundance
#  log10(data + 1)
data = pd.read_csv(infile,sep='\t',header=0,index_col = 0)
data = data + 1
norm_df = data.apply(np.log10)
sam_cols = norm_df.columns.values.tolist()

# topn
tmp_df = norm_df
tmp_df['sum'] = tmp_df.apply(lambda x: x.sum(),axis=1)
sort_df = tmp_df.sort_values('sum',ascending=False)
part_sort_df = sort_df.head(topn)
result_df = part_sort_df[sam_cols]
result_df.to_csv(ofile,sep='\t',header=True,index=True)


/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/heatmapin_genecount_topn.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Phylum_gene_count -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap --name phylum -n 35
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/heatmapin_genecount_topn.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Class_gene_count -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap --name class -n 35
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/heatmapin_genecount_topn.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Order_gene_count -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap --name order -n 35
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/heatmapin_genecount_topn.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Family_gene_count -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap --name family -n 35
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/heatmapin_genecount_topn.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Genus_gene_count -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap --name genus -n 35
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/heatmapin_genecount_topn.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Species_gene_count -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap --name species -n 35



gene abun, gene count pheatmap
library(pheatmap)
library(ggsci)

mycolors <- unique(c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                      pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
                      pal_locuszoom("default")(7),pal_igv("default")(51),
                      pal_uchicago("default")(9),pal_startrek("uniform")(7),
                      pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
                      pal_simpsons("springfield")(16),pal_gsea("default")(12)))

infile <- '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap/phylum_geneabun_heatmapin_top35'
reff <- '/work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group'
odir <- '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap'
infilename <- basename(infile)
name_suffix <- unlist(strsplit(infilename,'.',fixed=TRUE))
prename <- name_suffix[1]
scalerow_ofile <- paste0(odir,'/',prename,'_scalerow.pdf')
noscalerow_ofile <- paste0(odir,'/',prename,'.pdf')
intitle <- 'phylum abundance'
figwidth = 10
figheight = 5

raw_df <-read.table(infile, sep='\t', header=TRUE, row.names = 1, encoding='utf-8')
raw_mt <- data.matrix(raw_df)
ref_col_df <- read.table(reff,sep='\t',header=FALSE,row.names = 1, encoding='utf-8')
colnames(ref_col_df) = 'Group'
group_v = as.vector(unique(ref_col_df$Group))
group_n = length(unique(ref_col_df$Group))
select_col = mycolors[1:length(group_v)]
names(select_col) = group_v

#pheatmap(raw_mt)
# #main可设置热图的标题，fontsize设置字体大小，filename可直接将热图存出，支持格式png, pdf, tiff, bmp, jpeg，并且可以通过width, height设置图片的大小；
ann_colors = list(Group=select_col)
pheatmap(raw_mt, main=intitle, annotation_col = ref_col_df, angle_col = 45, annotation_colors = ann_colors, scale='row', cluster_col = FALSE, cutree_rows=group_n, filename = scalerow_ofile, width = figwidth, height = figheight)
pheatmap(raw_mt, main=intitle, annotation_col = ref_col_df, angle_col = 45, annotation_colors = ann_colors, cluster_col = FALSE, cutree_rows=group_n, filename = noscalerow_ofile, width = figwidth, height = figheight)

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/heatmap_taxo_topn.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap/phylum_geneabun_heatmapin_top35 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap -n 'phylum abundance top35' -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/heatmap_taxo_topn.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap/class_geneabun_heatmapin_top35 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap -n 'class abundance top35' -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/heatmap_taxo_topn.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap/order_geneabun_heatmapin_top35 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap -n 'order abundance top35' -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/heatmap_taxo_topn.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap/family_geneabun_heatmapin_top35 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap -n 'family abundance top35' -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/heatmap_taxo_topn.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap/genus_geneabun_heatmapin_top35 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap -n 'genus abundance top35' -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/heatmap_taxo_topn.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap/species_geneabun_heatmapin_top35 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap -n 'species abundance top35' -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/heatmap_taxo_topn.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap/phylum_genecount_heatmapin_top35 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap -n 'phylum gene count top35' -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/heatmap_taxo_topn.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap/class_genecount_heatmapin_top35 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap -n 'class gene count top35' -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/heatmap_taxo_topn.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap/order_genecount_heatmapin_top35 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap -n 'order gene count top35' -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/heatmap_taxo_topn.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap/family_genecount_heatmapin_top35 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap -n 'family gene count top35' -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/heatmap_taxo_topn.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap/genus_genecount_heatmapin_top35 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap -n 'genus gene count top35' -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/heatmap_taxo_topn.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap/species_genecount_heatmapin_top35 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_genecountabun_heatmap -n 'species gene count top35' -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group


taxo pca_mds_pcoa
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/pca_mds_pcoa.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Phylum_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_pca_mds_pcoa -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_all -n phylum

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/pca_mds_pcoa.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Class_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_pca_mds_pcoa -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_all -n class

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/pca_mds_pcoa.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Order_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_pca_mds_pcoa -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_all -n order

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/pca_mds_pcoa.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Family_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_pca_mds_pcoa -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_all -n family

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/pca_mds_pcoa.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Species_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_pca_mds_pcoa -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_all -n species

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/pca_mds_pcoa.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Genus_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_pca_mds_pcoa -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_all -n genus

2020.12.17
filtered taxo
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/pca_mds_pcoa.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/phylum/phylum_lefse_pre -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_pca_mds_pcoa_dg/phylum -r /work/workspace/zhurj/project/1_metadata/mouse22/input/sgb_group_title -n phylum

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/pca_mds_pcoa.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/class/class_lefse_pre -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_pca_mds_pcoa_dg/class -r /work/workspace/zhurj/project/1_metadata/mouse22/input/sgb_group_title -n class

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/pca_mds_pcoa.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/order/order_lefse_pre -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_pca_mds_pcoa_dg/order -r /work/workspace/zhurj/project/1_metadata/mouse22/input/sgb_group_title -n order

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/pca_mds_pcoa.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/family/family_lefse_pre -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_pca_mds_pcoa_dg/family -r /work/workspace/zhurj/project/1_metadata/mouse22/input/sgb_group_title -n family

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/pca_mds_pcoa.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/genus/genus_lefse_pre -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_pca_mds_pcoa_dg/genus -r /work/workspace/zhurj/project/1_metadata/mouse22/input/sgb_group_title -n genus

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/pca_mds_pcoa.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_lefse_pre -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_pca_mds_pcoa_dg/species  -r /work/workspace/zhurj/project/1_metadata/mouse22/input/sgb_group_title -n species


/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/bar_cluster.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_group_bar/barin_phylum_top10 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_bar_cluster/bar_cluster_phylum.png -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_bar_cluster -n 'Phylum'

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/bar_cluster.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_group_bar/barin_class_top10 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_bar_cluster/bar_cluster_class.png -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_bar_cluster -n 'Class'

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/bar_cluster.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_group_bar/barin_order_top10 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_bar_cluster/bar_cluster_order.png -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_bar_cluster -n 'Order'

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/bar_cluster.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_group_bar/barin_family_top10 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_bar_cluster/bar_cluster_family.png -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_bar_cluster -n 'Family'

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/bar_cluster.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_group_bar/barin_genus_top10 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_bar_cluster/bar_cluster_genus.png -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_bar_cluster -n 'Genus'

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/bar_cluster.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/sample_group_bar/barin_species_top10 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_bar_cluster/bar_cluster_species.png -r /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_bar_cluster -n 'Species'



# taxo lefse differential analysis
import pandas as pd
import numpy as np
import os


infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/gene_taxo_abun'
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse/taxo_lefse_pre'
filter_n = 3

data = pd.read_csv(infile,sep='\t',header=0,index_col = 0)
data.fillna('',inplace=True)
data['taxo'] = data['Domain'].str.cat([data['Phylum'],data['Class'],data['Order'],data['Family'],data['Genus'],data['Species']], sep='|')
sam_col = [x for x in data.columns if 'MN' in x]
select_col = ['taxo'] + sam_col
part_df = data[select_col]
group_df = part_df.groupby('taxo').sum()
group_df['count'] = (group_df[sam_col] >= 0.0001 ).astype(int).sum(axis=1)

# select taxo with >= 3 sample with abundance
group_filter_df = group_df[group_df['count'] >= filter_n]
result_df = group_filter_df[sam_col]
result_df['taxos'] = result_df.index
result_df['taxos'] = result_df['taxos'].replace('\|+$','',regex=True)
result_col = ['taxos'] + sam_col
final_df = result_df[result_col]
final_df.to_csv(ofile,sep='\t',header=True,index=False)


source activate py2
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse/taxo_lefse_pre /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse/taxo_lefse_lefsein
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/lefse-format_input.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse/taxo_lefse_lefsein /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse/taxo_lefse_lefseformatin -c 1 -u 2 -f r -o 1000000
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/run_lefse.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse/taxo_lefse_lefseformatin /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse/taxo_lefse_lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1
# lefse-plot_res.py, xuyao matplotlib=1.5.0 , must notice it
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/lefse-plot_res.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse/taxo_lefse_lefse_run.res  /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse/taxo_lefse_lefse_run.pdf --feature_font_size 8 --width 15 --format pdf
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/lefse-plot_cladogram.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse/taxo_lefse_lefse_run.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse/taxo_lefse_cladogram.pdf --format pdf --dpi 200

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse_fitler.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelC_abundance_lefse_run.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/test/kegg_levelC_abundance_lefse_run_filter.res
source deactivate py2

source activate py2
taxo comparation between HFD vs MNO863
1. copy one files from /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefs to /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg
taxo_lefse_pre: merge taxo, abundance in each sample
2. remove 6 samples from /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/taxo_lefse_pre
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/extract_cols.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/taxo_lefse_pre -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/taxo_lefse_pre_dg -r MNA00233 MNA00234 MNA00235 MNA00236 MNA00237 MNA00238
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/taxo_lefse_pre_dg /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/taxo_lefse_dg_lefsein
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/lefse-format_input.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/taxo_lefse_dg_lefsein /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/taxo_lefse_dg_lefseformatin -c 1 -u 2 -f r -o 1000000
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/run_lefse.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/taxo_lefse_dg_lefseformatin /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/taxo_lefse_dg_run.res -a 0.05 -w 0.05 -l 2 -y 1
# lefse-plot_res.py, xuyao matplotlib=1.5.0 , must notice it
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/lefse-plot_res.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/taxo_lefse_dg_run.res  /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/taxo_lefse_dg_run.pdf --feature_font_size 8 --width 10 --format pdf
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/lefse-plot_cladogram.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/taxo_lefse_dg_run.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/taxo_lefse_dg_cladogram.pdf --format pdf --dpi 200

source deactivate py2

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/extract_cols.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Species_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/Species_dg_abun -r MNA00233 MNA00234 MNA00235 MNA00236 MNA00237 MNA00238

import pandas as pd
import os
import sys

"""
description_comment = "create input of lefse, criteria, 1: <= cutoff as 0, 2: "
"""
sam_cutoff = 3
cutoff = 0.0001
infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/Species_dg_abun'
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_lefse_pre'
data = pd.read_csv(infile,sep='\t',header=0,index_col=0)
data[data < 0.0001] = 0
sam_cols = data.columns.values.tolist()
data['count'] = (data > 0).astype(int).sum(axis=1)
filter_df = data[ data['count'] >= sam_cutoff]
result_df = filter_df[sam_cols]
result_df['taxo'] = result_df.index
result_df['taxo'] = result_df['taxo'].str.replace('(','')
result_df['taxo'] = result_df['taxo'].str.replace(')','')
result_df['taxo'] = result_df['taxo'].str.replace(']','')
result_df['taxo'] = result_df['taxo'].str.replace('[','')
result_df['taxo'] = result_df['taxo'].str.replace(':','_')
result_df['taxo'] = result_df['taxo'].str.replace('sp.','sp')
result_df['taxo'] = result_df['taxo'].str.replace('-','_')
result_df['taxo'] = result_df['taxo'].str.replace(' ','_')
select_cols = ['taxo'] + sam_cols
final_df = result_df[select_cols]
final_df.to_csv(ofile,sep='\t',header=True,index=False)

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/heatmap_taxo_formation_abun2count.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/Species_dg_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_lefse_pre

source activate py2
#taxo comparation between HFD vs MNO863, species strain level > 0.0001, >=3 sample with abundance >= 0.0001
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_lefse_pre /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_lefse_dg_lefsein
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/lefse-format_input.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_lefse_dg_lefsein /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_lefse_dg_lefseformatin -c 1 -u 2 -f r -o 1000000
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/run_lefse.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_lefse_dg_lefseformatin /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_lefse_dg_run.res -a 0.05 -w 0.05 -l 2 -y 1
# lefse-plot_res.py, xuyao matplotlib=1.5.0 , must notice it
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/lefse-plot_res.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_lefse_dg_run.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_lefse_dg_run.pdf --feature_font_size 8 --width 10 --format pdf

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse_fitler.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_lefse_dg_run.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_lefse_dg_run_filter.res
python /work/workspace/zhurj/script/2_metapro/metapipe/toheadmapin.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_lefse_pre /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_lefse_dg_run_filter.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_heatmap_prein

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/heatmapin_abun_topn.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_heatmap_prein -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species --name species -n 150

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/heatmap_taxo_topn.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_geneabun_heatmapin_top150 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species -n 'dif species abundance' -r /work/workspace/zhurj/project/1_metadata/mouse22/input/sgb_group --plot_width 12  --plot_height 18
# heatmap add cols
create columns (up and down)

import pandas as pd
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_lefse_dg_run_filter.res'
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/status'
item_index_in = 1
group_index_in = 3
status_up = 'MNO863'
status_down = 'HFD'
item_index = item_index_in - 1
group_index = group_index_in -1

data = pd.read_csv(infile,sep='\t',header=None,index_col=None)
select_col = [item_index,group_index]
result_df = data[select_col]
result_df[group_index] = result_df[group_index].str.replace(status_up,'UP')
result_df[group_index] = result_df[group_index].str.replace(status_down,'DOWN')
result_df.to_csv(ofile,sep='\t',header=False,index=False)

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/extract_item_updown_2group.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_lefse_dg_run_filter.res -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/status --status_up MNO863 --status_down HFD --item_index 1 --group_index 3


/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/heatmap_taxo_topn_addrowanno.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_geneabun_heatmapin_top150 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species -n 'dif species abundance' -r /work/workspace/zhurj/project/1_metadata/mouse22/input/sgb_group --plot_width 12  --plot_height 18 --rowf /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/status



source activate py2
# genus
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/extract_cols.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Genus_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/genus/Genus_dg_abun -r MNA00233 MNA00234 MNA00235 MNA00236 MNA00237 MNA00238
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/heatmap_taxo_formation_abun2count.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/genus/Genus_dg_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/genus/genus_lefse_pre
#taxo comparation between HFD vs MNO863, species strain level > 0.0001, >=3 sample with abundance >= 0.0001
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/genus/genus_lefse_pre /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/genus/genus_lefsein
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/lefse-format_input.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/genus/genus_lefsein /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/genus/genus_lefseformatin -c 1 -u 2 -f r -o 1000000
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/run_lefse.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/genus/genus_lefseformatin /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/genus/genus_lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1
# lefse-plot_res.py, xuyao matplotlib=1.5.0 , must notice it
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/lefse-plot_res.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/genus/genus_lefse_run.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/genus/genus_lefse_run.pdf --feature_font_size 8 --width 10 --format pdf
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse_fitler.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/genus/genus_lefse_run.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/genus/genus_lefse_run_filter.res
python /work/workspace/zhurj/script/2_metapro/metapipe/toheadmapin.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/genus/genus_lefse_pre /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/genus/genus_lefse_run_filter.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/genus/genus_heatmap_prein
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/heatmapin_abun_topn.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/genus/genus_heatmap_prein -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/genus --name genus -n 150
# heatmap add cols
create columns (up and down)
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/extract_item_updown_2group.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/genus/genus_lefse_run_filter.res -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/genus/status --status_up MNO863 --status_down HFD --item_index 1 --group_index 3
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/heatmap_taxo_topn_addrowanno.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/genus/genus_geneabun_heatmapin_top150 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/genus -n 'dif genus abundance' -r /work/workspace/zhurj/project/1_metadata/mouse22/input/sgb_group --plot_width 12  --plot_height 15 --rowf /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/genus/status

# family
source activate py2
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/extract_cols.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Family_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/family/Family_dg_abun -r MNA00233 MNA00234 MNA00235 MNA00236 MNA00237 MNA00238
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/heatmap_taxo_formation_abun2count.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/family/Family_dg_abun -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/family/family_lefse_pre
#taxo comparation between HFD vs MNO863, species strain level > 0.0001, >=3 sample with abundance >= 0.0001
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/family/family_lefse_pre /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/family/family_lefsein
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/lefse-format_input.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/family/family_lefsein /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/family/family_lefseformatin -c 1 -u 2 -f r -o 1000000
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/run_lefse.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/family/family_lefseformatin /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/family/family_lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1
# lefse-plot_res.py, xuyao matplotlib=1.5.0 , must notice it
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/lefse-plot_res.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/family/family_lefse_run.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/family/family_lefse_run.pdf --feature_font_size 8 --width 10 --format pdf
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse_fitler.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/family/family_lefse_run.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/family/family_lefse_run_filter.res
python /work/workspace/zhurj/script/2_metapro/metapipe/toheadmapin.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/family/family_lefse_pre /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/family/family_lefse_run_filter.res /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/family/family_heatmap_prein
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/heatmapin_abun_topn.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/family/family_heatmap_prein -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/family --name family -n 150
# heatmap add cols
create columns (up and down)
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/extract_item_updown_2group.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/family/family_lefse_run_filter.res -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/family/status --status_up MNO863 --status_down HFD --item_index 1 --group_index 3
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/heatmap_taxo_topn_addrowanno.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/family/family_geneabun_heatmapin_top150 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/family -n 'dif family abundance' -r /work/workspace/zhurj/project/1_metadata/mouse22/input/sgb_group --plot_width 12  --plot_height 8 --rowf /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/family/status
source deactivate py2

# species
source activate py2
parallel --link -j 1 '{1} {2}/extract_cols.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Species_abun -o {3}/{4}/{4}_dg_abun -r MNA00233 MNA00234 MNA00235 MNA00236 MNA00237 MNA00238 && \
{1} {2}/heatmap_taxo_formation_abun2count.py -i {3}/{4}/{4}_dg_abun -o {3}/{4}/{4}_lefse_pre && \
{1} {2}/lefse_in.py {3}/{4}/{4}_lefse_pre /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group {3}/{4}/{4}_lefsein && \
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/lefse-format_input.py {3}/{4}/{4}_lefsein {3}/{4}/{4}_lefseformatin -c 1 -u 2 -f r -o 1000000 && \
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/run_lefse.py {3}/{4}/{4}_lefseformatin {3}/{4}/{4}_lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1 && \
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/lefse-plot_res.py {3}/{4}/{4}_lefse_run.res {3}/{4}/{4}_lefse_run.pdf --feature_font_size 8 --width 10 --format pdf  && \
{1} /work/workspace/zhurj/lib/python/script/lefse_fitler.py {3}/{4}/{4}_lefse_run.res {3}/{4}/{4}_lefse_run_filter.res && \
{1} {2}/toheadmapin.py {3}/{4}/{4}_lefse_pre {3}/{4}/{4}_lefse_run_filter.res {3}/{4}/{4}_heatmap_prein && \
{1} {2}/heatmapin_abun_topn.py -i {3}/{4}/{4}_heatmap_prein -o {3}/{4}  --name {4} -n 150 && \
{1} {2}/extract_item_updown_2group.py -i {3}/{4}/{4}_lefse_run_filter.res -o {3}/{4}/status --status_up MNO863 --status_down HFD --item_index 1 --group_index 3 && \
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript {2}/heatmap_taxo_topn_addrowanno.r -i {3}/{4}/{4}_geneabun_heatmapin_top150 -o {3}/{4} -n dif_{5}_abundance -r /work/workspace/zhurj/project/1_metadata/mouse22/input/sgb_group --plot_width 12  --plot_height {5} --rowf {3}/{4}/status  \
' ::: /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python ::: /work/workspace/zhurj/script/2_metapro/metapipe ::: /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg ::: species ::: 18


# genus
source activate py2
parallel --link -j 1 '{1} {2}/extract_cols.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Genus_abun -o {3}/{4}/{4}_dg_abun -r MNA00233 MNA00234 MNA00235 MNA00236 MNA00237 MNA00238 && \
{1} {2}/heatmap_taxo_formation_abun2count.py -i {3}/{4}/{4}_dg_abun -o {3}/{4}/{4}_lefse_pre && \
{1} {2}/lefse_in.py {3}/{4}/{4}_lefse_pre /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group {3}/{4}/{4}_lefsein && \
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/lefse-format_input.py {3}/{4}/{4}_lefsein {3}/{4}/{4}_lefseformatin -c 1 -u 2 -f r -o 1000000 && \
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/run_lefse.py {3}/{4}/{4}_lefseformatin {3}/{4}/{4}_lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1 && \
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/lefse-plot_res.py {3}/{4}/{4}_lefse_run.res {3}/{4}/{4}_lefse_run.pdf --feature_font_size 8 --width 10 --format pdf  && \
{1} /work/workspace/zhurj/lib/python/script/lefse_fitler.py {3}/{4}/{4}_lefse_run.res {3}/{4}/{4}_lefse_run_filter.res && \
{1} {2}/toheadmapin.py {3}/{4}/{4}_lefse_pre {3}/{4}/{4}_lefse_run_filter.res {3}/{4}/{4}_heatmap_prein && \
{1} {2}/heatmapin_abun_topn.py -i {3}/{4}/{4}_heatmap_prein -o {3}/{4}  --name {4} -n 150 && \
{1} {2}/extract_item_updown_2group.py -i {3}/{4}/{4}_lefse_run_filter.res -o {3}/{4}/status --status_up MNO863 --status_down HFD --item_index 1 --group_index 3 && \
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript {2}/heatmap_taxo_topn_addrowanno.r -i {3}/{4}/{4}_geneabun_heatmapin_top150 -o {3}/{4} -n dif_{5}_abundance -r /work/workspace/zhurj/project/1_metadata/mouse22/input/sgb_group --plot_width 12  --plot_height {5} --rowf {3}/{4}/status  \
' ::: /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python ::: /work/workspace/zhurj/script/2_metapro/metapipe ::: /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg ::: genus ::: 15

# family
source activate py2
parallel --link -j 1 '{1} {2}/extract_cols.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Family_abun -o {3}/{4}/{4}_dg_abun -r MNA00233 MNA00234 MNA00235 MNA00236 MNA00237 MNA00238 && \
{1} {2}/heatmap_taxo_formation_abun2count.py -i {3}/{4}/{4}_dg_abun -o {3}/{4}/{4}_lefse_pre && \
{1} {2}/lefse_in.py {3}/{4}/{4}_lefse_pre /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group {3}/{4}/{4}_lefsein && \
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/lefse-format_input.py {3}/{4}/{4}_lefsein {3}/{4}/{4}_lefseformatin -c 1 -u 2 -f r -o 1000000 && \
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/run_lefse.py {3}/{4}/{4}_lefseformatin {3}/{4}/{4}_lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1 && \
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/lefse-plot_res.py {3}/{4}/{4}_lefse_run.res {3}/{4}/{4}_lefse_run.pdf --feature_font_size 8 --width 10 --format pdf  && \
{1} /work/workspace/zhurj/lib/python/script/lefse_fitler.py {3}/{4}/{4}_lefse_run.res {3}/{4}/{4}_lefse_run_filter.res && \
{1} {2}/toheadmapin.py {3}/{4}/{4}_lefse_pre {3}/{4}/{4}_lefse_run_filter.res {3}/{4}/{4}_heatmap_prein && \
{1} {2}/heatmapin_abun_topn.py -i {3}/{4}/{4}_heatmap_prein -o {3}/{4}  --name {4} -n 150 && \
{1} {2}/extract_item_updown_2group.py -i {3}/{4}/{4}_lefse_run_filter.res -o {3}/{4}/status --status_up MNO863 --status_down HFD --item_index 1 --group_index 3 && \
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript {2}/heatmap_taxo_topn_addrowanno.r -i {3}/{4}/{4}_geneabun_heatmapin_top150 -o {3}/{4} -n dif_{5}_abundance -r /work/workspace/zhurj/project/1_metadata/mouse22/input/sgb_group --plot_width 12  --plot_height {5} --rowf {3}/{4}/status  \
' ::: /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python ::: /work/workspace/zhurj/script/2_metapro/metapipe ::: /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg ::: family ::: 8

# order
source activate py2
parallel --link -j 1 '{1} {2}/extract_cols.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Order_abun -o {3}/{4}/{4}_dg_abun -r MNA00233 MNA00234 MNA00235 MNA00236 MNA00237 MNA00238 && \
{1} {2}/heatmap_taxo_formation_abun2count.py -i {3}/{4}/{4}_dg_abun -o {3}/{4}/{4}_lefse_pre && \
{1} {2}/lefse_in.py {3}/{4}/{4}_lefse_pre /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group {3}/{4}/{4}_lefsein && \
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/lefse-format_input.py {3}/{4}/{4}_lefsein {3}/{4}/{4}_lefseformatin -c 1 -u 2 -f r -o 1000000 && \
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/run_lefse.py {3}/{4}/{4}_lefseformatin {3}/{4}/{4}_lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1 && \
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/lefse-plot_res.py {3}/{4}/{4}_lefse_run.res {3}/{4}/{4}_lefse_run.pdf --feature_font_size 8 --width 10 --format pdf  && \
{1} /work/workspace/zhurj/lib/python/script/lefse_fitler.py {3}/{4}/{4}_lefse_run.res {3}/{4}/{4}_lefse_run_filter.res && \
{1} {2}/toheadmapin.py {3}/{4}/{4}_lefse_pre {3}/{4}/{4}_lefse_run_filter.res {3}/{4}/{4}_heatmap_prein && \
{1} {2}/heatmapin_abun_topn.py -i {3}/{4}/{4}_heatmap_prein -o {3}/{4}  --name {4} -n 150 && \
{1} {2}/extract_item_updown_2group.py -i {3}/{4}/{4}_lefse_run_filter.res -o {3}/{4}/status --status_up MNO863 --status_down HFD --item_index 1 --group_index 3 && \
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript {2}/heatmap_taxo_topn_addrowanno.r -i {3}/{4}/{4}_geneabun_heatmapin_top150 -o {3}/{4} -n dif_{5}_abundance -r /work/workspace/zhurj/project/1_metadata/mouse22/input/sgb_group --plot_width 12  --plot_height {5} --rowf {3}/{4}/status  \
' ::: /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python ::: /work/workspace/zhurj/script/2_metapro/metapipe ::: /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg ::: order ::: 5 

# class
parallel --link -j 1 '{1} {2}/extract_cols.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Class_abun -o {3}/{4}/{4}_dg_abun -r MNA00233 MNA00234 MNA00235 MNA00236 MNA00237 MNA00238 && \
{1} {2}/heatmap_taxo_formation_abun2count.py -i {3}/{4}/{4}_dg_abun -o {3}/{4}/{4}_lefse_pre && \
{1} {2}/lefse_in.py {3}/{4}/{4}_lefse_pre /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group {3}/{4}/{4}_lefsein && \
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/lefse-format_input.py {3}/{4}/{4}_lefsein {3}/{4}/{4}_lefseformatin -c 1 -u 2 -f r -o 1000000 && \
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/run_lefse.py {3}/{4}/{4}_lefseformatin {3}/{4}/{4}_lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1 && \
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/lefse-plot_res.py {3}/{4}/{4}_lefse_run.res {3}/{4}/{4}_lefse_run.pdf --feature_font_size 8 --width 10 --format pdf  && \
{1} /work/workspace/zhurj/lib/python/script/lefse_fitler.py {3}/{4}/{4}_lefse_run.res {3}/{4}/{4}_lefse_run_filter.res && \
{1} {2}/toheadmapin.py {3}/{4}/{4}_lefse_pre {3}/{4}/{4}_lefse_run_filter.res {3}/{4}/{4}_heatmap_prein && \
{1} {2}/heatmapin_abun_topn.py -i {3}/{4}/{4}_heatmap_prein -o {3}/{4}  --name {4} -n 150 && \
{1} {2}/extract_item_updown_2group.py -i {3}/{4}/{4}_lefse_run_filter.res -o {3}/{4}/status --status_up MNO863 --status_down HFD --item_index 1 --group_index 3 && \
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript {2}/heatmap_taxo_topn_addrowanno.r -i {3}/{4}/{4}_geneabun_heatmapin_top150 -o {3}/{4} -n dif_{5}_abundance -r /work/workspace/zhurj/project/1_metadata/mouse22/input/sgb_group --plot_width 12  --plot_height {5} --rowf {3}/{4}/status  \
' ::: /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python ::: /work/workspace/zhurj/script/2_metapro/metapipe ::: /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg ::: class ::: 5 

# phylum
parallel --link -j 1 '{1} {2}/extract_cols.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Phylum_abun -o {3}/{4}/{4}_dg_abun -r MNA00233 MNA00234 MNA00235 MNA00236 MNA00237 MNA00238 && \
{1} {2}/heatmap_taxo_formation_abun2count.py -i {3}/{4}/{4}_dg_abun -o {3}/{4}/{4}_lefse_pre && \
{1} {2}/lefse_in.py {3}/{4}/{4}_lefse_pre /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group {3}/{4}/{4}_lefsein && \
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/lefse-format_input.py {3}/{4}/{4}_lefsein {3}/{4}/{4}_lefseformatin -c 1 -u 2 -f r -o 1000000 && \
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/run_lefse.py {3}/{4}/{4}_lefseformatin {3}/{4}/{4}_lefse_run.res -a 0.05 -w 0.05 -l 2 -y 1 && \
/work/workspace/zhurj/bin/miniconda3/envs/py2/bin/lefse-plot_res.py {3}/{4}/{4}_lefse_run.res {3}/{4}/{4}_lefse_run.pdf --feature_font_size 8 --width 10 --format pdf  && \
{1} /work/workspace/zhurj/lib/python/script/lefse_fitler.py {3}/{4}/{4}_lefse_run.res {3}/{4}/{4}_lefse_run_filter.res && \
{1} {2}/toheadmapin.py {3}/{4}/{4}_lefse_pre {3}/{4}/{4}_lefse_run_filter.res {3}/{4}/{4}_heatmap_prein && \
{1} {2}/heatmapin_abun_topn.py -i {3}/{4}/{4}_heatmap_prein -o {3}/{4}  --name {4} -n 150 && \
{1} {2}/extract_item_updown_2group.py -i {3}/{4}/{4}_lefse_run_filter.res -o {3}/{4}/status --status_up MNO863 --status_down HFD --item_index 1 --group_index 3 && \
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript {2}/heatmap_taxo_topn_addrowanno.r -i {3}/{4}/{4}_geneabun_heatmapin_top150 -o {3}/{4} -n dif_{5}_abundance -r /work/workspace/zhurj/project/1_metadata/mouse22/input/sgb_group --plot_width 12  --plot_height {5} --rowf {3}/{4}/status  \
' ::: /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python ::: /work/workspace/zhurj/script/2_metapro/metapipe ::: /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg ::: phylum ::: 3 

source deactivate py2

KO statistic
import pandas as pd

infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/uniref50_all/unigene/ko_abundance'
data = pd.read_csv(infile,sep='\t',header=0,index_col = 0)
(data > 0).astype(int).sum()

# taxo statistics
infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Phylum_abun'
data = pd.read_csv(infile,sep='\t',header=0,index_col = 0)
(data > 0).astype(int).sum()

infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Class_abun'
data = pd.read_csv(infile,sep='\t',header=0,index_col = 0)
(data > 0).astype(int).sum()

infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Order_abun'
data = pd.read_csv(infile,sep='\t',header=0,index_col = 0)
(data > 0).astype(int).sum()

infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Family_abun'
data = pd.read_csv(infile,sep='\t',header=0,index_col = 0)
(data > 0).astype(int).sum()

infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Genus_abun'
data = pd.read_csv(infile,sep='\t',header=0,index_col = 0)
(data > 0).astype(int).sum()

infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Species_abun'
data = pd.read_csv(infile,sep='\t',header=0,index_col = 0)
(data > 0).astype(int).sum()

infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Species_abun'
data = pd.read_csv(infile,sep='\t',header=0,index_col = 0)
(data >= 0.0001).astype(int).sum()

2020.12.15
classification of high-quality assembled genomes
source activate python3.6
export GTDBTK_DATA_PATH=/work/workspace/zhurj/reference/GTDB/release89
gtdbtk classify_wf --batchfile /work/workspace/zhurj/project/1_metadata/mouse22/metabat/input/check_completeness90_contamination5 --out_dir /work/workspace/zhurj/project/1_metadata/mouse22/metabat/gtdbtk_hq_genome --cpus 20

muscle -in /work/workspace/zhurj/project/1_metadata/mouse22/metabat/gtdbtk_hq_genome/gtdbtk.bac120.user_msa.fasta -out /work/workspace/zhurj/project/1_metadata/mouse22/metabat/phylo/mouse22_genome249.fa -maxiters 16  

from Bio import SeqIO
records = SeqIO.parse("/work/workspace/zhurj/project/1_metadata/mouse22/metabat/phylo/mouse22_genome249.fa", "fasta")
count = SeqIO.write(records, "/work/workspace/zhurj/project/1_metadata/mouse22/metabat/phylo/mouse22_genome249.phylip", "phylip")
print("Converted %i records" % count)

nohup raxmlHPC-PTHREADS-SSE3 -s /work/workspace/zhurj/project/1_metadata/mouse22/metabat/phylo/mouse22_genome249.phylip -n tree -m PROTGAMMAJTT -T 30 -N 1000 -p 20170808 -f a -x 20170808 >  genome249.log &


MNA00245, MN00252 prokka 预测




cd /work/workspace/zhurj/project/1_metadata/mouse22/com/merge/spades/MNA00245 && /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/seqtk seq -L 500 scaffolds.fasta > scaffolds_500.fasta
cd /work/workspace/zhurj/project/1_metadata/mouse22/com/merge/spades/MNA00252 && /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/seqtk seq -L 500 scaffolds.fasta > scaffolds_500.fasta

source activate /work/workspace/zhurj/bin/miniconda3/envs/prokka
parallel -j 1 --link 'prokka --outdir {1}/prokka/{2} -cpus 20 --prefix {2} --force --metagenome --compliant --locustag {2} {1}/spades/{2}/scaffolds_500.fasta' ::: /work/workspace/zhurj/project/1_metadata/mouse22/com/merge ::: MNA00245
parallel -j 1 --link 'prokka --outdir {1}/prokka/{2} -cpus 20 --prefix {2} --force --metagenome --compliant --locustag {2} {1}/spades/{2}/scaffolds_500.fasta' ::: /work/workspace/zhurj/project/1_metadata/mouse22/com/merge ::: MNA00252
echo `date +"%Y-%m-%d %H:%M:%S"`
echo "Process8: gene prediction finished "
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/prokka

cd /work/workspace/zhurj/project/1_metadata/p/p
time srun -o prokka245_252.out -e prokka245_252.err -N 1 -c 20 -p slurm128 -w mnclient04 bash prokka245_252.sh &

cat /work/workspace/zhurj/project/1_metadata/mouse22/metabat/input/check_completeness90_contamination5 | awk -F "\t" '{print "ln -s "$1" /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/phylophlan/input/"$2".fna"}' > /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/phylophlan/create_phylophlanin.sh

source activate /work/workspace/zhurj/bin/miniconda3/envs/phylophlan
#phylophlan_write_config_file -d a -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/phylophlan/config.cfg --db_aa diamond --map_dna diamond --map_aa diamond --msa mafft --trim trimal --tree1 raxml --verbose 
phylophlan -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/phylophlan/input -d phylophlan --diversity medium --accurate -f /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/phylophlan/config.cfg -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/phylophlan/output --nproc 20 -t a --verbose 

source deactivate /work/workspace/zhurj/bin/miniconda3/envs/phylophlan
scontrol show job JobId
time srun -o phylophlan_249.out -e phylophlan_249.err -N 1 -c 18 -p slurm256 -w mnclient02 bash phylophlan_249.sh &

assembled family bar
library(ggplot2)
library(dplyr)
library(forcats)

infile <- '/work/workspace/zhurj/project/1_metadata/mouse22/metabat/phylo/bar_family_in'
ofile <- '/work/workspace/zhurj/project/1_metadata/mouse22/metabat/phylo/family_bar.pdf'
title <- 'Family of high-quality assembled genome 248'
data <- read.table(infile,sep='\t',header=TRUE,comment.char="!")

used_color = data$color
names(used_color) = data$family
fb <- ggplot(data,aes(x=family,y=value)) + 
geom_bar(stat="identity",fill=used_color) +
ggtitle(title) + 
theme(plot.title = element_text(hjust = 0.5)) +
coord_flip()
ggsave(ofile,fb,width = 10,height = 6)


/work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/item_value_color_to_bargraph.r -i /work/workspace/zhurj/project/1_metadata/mouse22/metabat/phylo/bar_family_in -o /work/workspace/zhurj/project/1_metadata/mouse22/metabat/phylo/family_bar.pdf --plot_title "Family-19 of high-quality assembled genome 248"
/work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/item_value_color_to_bargraph.r -i /work/workspace/zhurj/project/1_metadata/mouse22/metabat/phylo/bar_genus_in -o /work/workspace/zhurj/project/1_metadata/mouse22/metabat/phylo/genus_bar.pdf --plot_title "Genus-46 of high-quality assembled genome 248"

beta diversity - assembled result
gene count

beta diversity

# abundance to count
import pandas as pd
import os
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_lefse_pre'
ofile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_count_100000'
multi_time = 100000

data = pd.read_csv(infile,sep='\t',header=0,index_col=0)
new_data = data * multi_time
new_data = new_data.astype(int)
new_data.to_csv(ofile,sep='\t',header=True,index=True)



Bray Curtis
library(vegan)

infile <- '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_lefse_pre'
ofile <-  '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_braycurtis_dist'
data <- read.table(infile,sep = '\t',header=T,row.names = 1)
data.dist <- vegdist(t(data),method='bray',diag=TRUE,upper=TRUE)
write.table(as.matrix(data.dist), file = ofile,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE, fileEncoding = "utf-8")


R4.0
library(ggplot2)
library(ggpubr)
infile <- '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_beta_diversity_in'
odir <- '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species'
name <- 'species_beta'
title = "Bray_Curtis_Distance"
ofile <- paste0(odir,'/',name,'_diversity.jpg')
Data <- read.table(infile,header=T)
unique_group <- as.vector(unique(Data$group))
my_comparisons <- list(c(unique_group[1],unique_group[3]),c(unique_group[1],unique_group[2]),c(unique_group[2],unique_group[3]))
b1 <- ggplot(data=Data, aes(x=group,y=bray_curtis))+
geom_boxplot(aes(fill=group)) +
scale_fill_manual(breaks = unique_group, values=c("#00AFBB", "#FC4E07","#E7B800")) +
stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", aes(label = ..p.signif..)) + 
ggtitle(title)+
theme(plot.title = element_text(face = "plain", vjust = 1, hjust=0.5, size = 10),
      axis.text.x = element_text(angle=0, hjust=0.5), axis.line.x=element_line(),axis.line.y=element_line(),
      axis.title.y=element_blank(), axis.text.y = element_text(angle=0, hjust=1, vjust=0.5),
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"),
      axis.text=element_text(face="plain", size=8, color = 'black'),
      panel.grid.major.y = element_blank(),,
      panel.grid.minor = element_blank(),
      panel.background=element_rect(fill='transparent', color='black',  size = 0.5),
      panel.border=element_rect(fill='transparent', color='black'),
      ) +
scale_x_discrete(name="", breaks = unique_group, labels= unique_group) +
scale_y_continuous(name="distance", limits=c(0.0, 0.7), breaks=seq(0.0,0.7,0.1))  + 
 guides(fill=FALSE)
ggsave(ofile,b1,width = 3,height = 3)

ofile <- paste0(odir,'/',name,'_diversity_violin.jpg')
my_comparisons <- list(c(unique_group[1],unique_group[3]))
Data = read.table(infile,header=T)
b1 <- ggplot(data=Data, aes(x=group,y=bray_curtis)) +
geom_violin(aes(fill=group),trim=FALSE) +
geom_boxplot(width=0.1) +
scale_x_discrete(name="", limits=c(unique_group[1],unique_group[3]),labels=c(unique_group[1],unique_group[3]),expand = c(0.25, 0.25)) +
scale_y_continuous(name="distance", limits=c(0.0, 0.6), breaks=seq(0.0,0.6,0.1)) +
stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", aes(label = ..p.signif..)) + 
ggtitle(title)+
theme(plot.title = element_text(face = "plain", vjust = 1, hjust=0.5, size = 10),
      axis.text.x = element_text(angle=0, hjust=0.5), axis.line.x=element_line(),axis.line.y=element_line(),
      axis.title.y=element_blank(), axis.text.y = element_text(angle=0, hjust=1, vjust=0.5),
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"),
      axis.text=element_text(face="plain", size=8, color = 'black'),
      panel.grid.major.y = element_blank(),,
      panel.grid.minor = element_blank(),
      panel.background=element_rect(fill='transparent', color='black',  size = 0.5),
      panel.border=element_rect(fill='transparent', color='black'),
      ) 
ggsave(ofile,b1,width = 3,height = 3)


/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/vedist_to_dist.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_count_100000 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_braycurtis_dist -m bray
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/distformat.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_braycurtis_dist /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group  /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_beta_diversity_in
/work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/beta_diversity_creation.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_beta_diversity_in -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species -n species_beta -t Bray_Curtis_Distance --plot_width 4 --plot_height 4

library(phyloseq)
library(ape)
# ref: https://blog.csdn.net/woodcorpse/article/details/106554382

infile <- '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_count_100000'
odir <- '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species'
name <- 'species'
data <- read.table(infile,sep='\t', header=TRUE, row.name = 1)
OTU <- otu_table(data, taxa_are_rows = TRUE)
physeq <- phyloseq(OTU)
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
physeq1 = phyloseq(OTU, random_tree)
ofile <- paste0(odir,'/',name,'_unweighted_unifrac_dist')
dist <- UniFrac(physeq1, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE)
mat <- as.matrix(dist)
write.table(mat, file = ofile,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE, fileEncoding = "utf-8")

ofile <- paste0(odir,'/',name,'_weighted_unifrac_dist')
dist <- UniFrac(physeq1, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
mat <- as.matrix(dist)
write.table(mat, file = ofile,sep='\t',quote=FALSE,row.names = TRUE,col.names = TRUE, fileEncoding = "utf-8")




/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/from_abun_to_count.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_lefse_pre -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_count_100000 --multi_index 100000

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/vedist_to_dist.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_count_100000 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_braycurtis_dist -m bray
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/distformat.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_braycurtis_dist /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group  /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_beta_diversity_in
/work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/beta_diversity_creation.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_beta_diversity_in -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species -n species_beta -t Bray_Curtis_Distance --plot_width 4 --plot_height 4

/work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/from_count_to_unifrac_dist.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_count_100000 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species -n species
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/distformat.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_unweighted_unifrac_dist /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group  /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_unifrac_beta_diversity_in
/work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/beta_diversity_creation.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_unifrac_beta_diversity_in -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species -n species_unifrac_beta -t Unifrac_Distance --plot_width 4 --plot_height 4 --y_max 0.4

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/distformat.py /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_weighted_unifrac_dist /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group  /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_weighted_unifrac_beta_diversity_in
/work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/beta_diversity_creation.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_weighted_unifrac_beta_diversity_in -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species -n species_weighted_unifrac_beta -t Weighted_Unifrac_Distance --plot_width 4 --plot_height 4 --y_max 0.4



# genus
parallel --link -j 1 '{1} {2}/from_abun_to_count.py -i {3}/{4}/{4}_lefse_pre -o {3}/{4}/{4}_count_100000 --multi_index 100000 && \
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript {2}/vedist_to_dist.r -i {3}/{4}/{4}_count_100000 -o {3}/{4}/{4}_braycurtis_dist -m bray && \
{1} {2}/distformat.py {3}/{4}/{4}_braycurtis_dist {5} {3}/{4}/{4}_beta_diversity_in && \
/work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript {2}/beta_diversity_creation.r -i {3}/{4}/{4}_beta_diversity_in -o {3}/{4} -n {4}_beta -t Bray_Curtis_Distance --plot_width {6} --plot_height {6} --y_max 0.5 && \

/work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript {2}/from_count_to_unifrac_dist.r -i {3}/{4}/{4}_count_100000 -o {3}/{4} -n {4} && \
{1} {2}/distformat.py {3}/{4}/{4}_unweighted_unifrac_dist {5} {3}/{4}/{4}_unifrac_beta_diversity_in && \
/work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript {2}/beta_diversity_creation.r -i {3}/{4}/{4}_unifrac_beta_diversity_in -o {3}/{4} -n {4}_unifrac_beta -t Unifrac_Distanc --plot_width {6} --plot_height {6} --y_max 0.2 && \
{1} {2}/distformat.py {3}/{4}/{4}_weighted_unifrac_dist {5} {3}/{4}/{4}_weighted_unifrac_beta_diversity_in && \
/work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript {2}/beta_diversity_creation.r -i {3}/{4}/{4}_weighted_unifrac_beta_diversity_in -o {3}/{4} -n {4}_weighted_unifrac_beta -t Weighted_Unifrac_Distanc --plot_width {6} --plot_height {6} --y_max 0.3 \
' ::: /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python ::: /work/workspace/zhurj/script/2_metapro/metapipe ::: /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg ::: genus ::: /work/workspace/zhurj/project/1_metadata/mouse22/input/meta_group ::: 4


import pandas as pd
cutoff = 0.0001
infile = '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/diamond/Genus_abun'
data = pd.read_csv(infile,sep='\t',header=0,index_col=0)
data['max'] = data.apply(lambda x: x.max(),axis=1)
data[data['max'] >= cutoff]



co-abundance
species
#/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_count_100000
#'fastspar':'/work/workspace/liangzj/program/Miniconda3/envs/local/bin/fastspar',
#'fastspar_bootstrap':'/work/workspace/liangzj/program/Miniconda3/envs/local/bin/fastspar_bootstrap',
#'fastspar_pvalues':'/work/workspace/liangzj/program/Miniconda3/envs/local/bin/fastspar_pvalues'

mkdir bootstrap_counts
mkdir bootstrap_correlation
/work/workspace/liangzj/program/Miniconda3/envs/local/bin/fastspar -c species_count_100000 -r cor_sparcc.fastspar.txt -a cov_sparcc.fastspar.txt -t 20 
/work/workspace/liangzj/program/Miniconda3/envs/local/bin/fastspar_bootstrap --otu_table species_count_100000 --n 100 --prefix bootstrap_counts/bootstrap -t 20
parallel /work/workspace/liangzj/program/Miniconda3/envs/local/bin/fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 50 -t 20 ::: bootstrap_counts/*
/work/workspace/liangzj/program/Miniconda3/envs/local/bin/fastspar_pvalues --otu_table species_count_100000 --correlation cor_sparcc.fastspar.txt --prefix bootstrap_correlation/cor_bootstrap --permutations 100 --outfile pvalues.tsv -t 20


parallel -j 1 'cp {1}/taxo_lefse_dg/{2}/{2}_count_100000 {1}/coabundance/{2}' ::: /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo ::: genus
R3.6
library(igraph)

in_cor_file <- '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/coabundance/species/cor_sparcc.fastspar.txt'
in_pvalue_file <- '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/coabundance/species/pvalues.tsv'
odir <- '/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/coabundance/species'
cor_cutoff <- 0.6
pvalue_cutoff <- 0.05

cor_sparcc <- read.delim(in_cor_file, row.names = 1, sep = '\t', check.names = FALSE)
#伪 p 值矩阵
pvals <- read.delim(in_pvalue_file, row.names = 1, sep = '\t', check.names = FALSE)
#保留 |相关性|≥0.8 且 p<0.01的值
cor_sparcc[abs(cor_sparcc) <= cor_cutoff] <- 0
pvals[pvals >= pvalue_cutoff] <- 0
pvals[pvals < pvalue_cutoff & pvals > 0] <- 1
 
#筛选后的邻接矩阵
adj <- as.matrix(cor_sparcc) * as.matrix(pvals)
diag(adj) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
ofile <- paste0(odir,'/network.adj.txt')
write.table(data.frame(adj, check.names = FALSE), ofile, col.names = NA, sep = '\t', quote = FALSE)
##网络格式转换
#输入数据，邻接矩阵
tmp_infile <- ofile
network_adj <- read.delim(tmp_infile, row.names = 1, sep = '\t', check.names = FALSE)
#邻接矩阵 -> igraph 的邻接列表，获得含权的无向网络
g <- graph_from_adjacency_matrix(as.matrix(network_adj), mode = 'undirected', weighted = TRUE, diag = FALSE)
#g    #igraph 的邻接列表
#这种转换模式下，默认的边权重代表了 sparcc 计算的相关性（存在负值）
#由于边权重通常为正值，因此最好取个绝对值，相关性重新复制一列作为记录
E(g)$sparcc <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)
edge <- data.frame(as_edgelist(g))
 
edge_list <- data.frame(
    source = edge[[1]],
    target = edge[[2]],
    weight = E(g)$weight,
    sparcc = E(g)$sparcc
)
ofile <- paste0(odir,'/network.edge_list.txt')
write.table(edge_list, ofile, sep = '\t', row.names = FALSE, quote = FALSE) 
#节点属性列表，对应边列表，记录节点属性，例如
node_list <- data.frame(
    nodes_id = V(g)$name,    #节点名称
    degree = degree(g)    #节点度
)
ofile <- paste0(odir,'/network.node_list.txt')
write.table(node_list, ofile, sep = '\t', row.names = FALSE, quote = FALSE)



/usr/bin/bash  /work/workspace/zhurj/script/2_metapro/metapipe/fastspar_to_cor_pvalue /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/coabundance/species/species_count_100000 /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/coabundance/species
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/cor_pvalue_to_edge_node_list.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/coabundance/species/cor_sparcc.fastspar.txt -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/coabundance/species --pfile /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/coabundance/species/pvalues.tsv --p_cutoff 0.05 --cor_cutoff 0.6

/usr/bin/bash  /work/workspace/zhurj/script/2_metapro/metapipe/fastspar_to_cor_pvalue /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/coabundance/genus/genus_count_100000 /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/coabundance/genus
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/cor_pvalue_to_edge_node_list.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/coabundance/genus/cor_sparcc.fastspar.txt -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/coabundance/genus --pfile /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/coabundance/genus/pvalues.tsv --p_cutoff 0.05 --cor_cutoff 0.6

parallel -j 1 '/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/from_abun_to_count.py -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/{1}/{1}_lefse_pre -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/coabundance/{1}/{1}_count_100000 --multi_index 100000' ::: phylum class order family
# {1}_count_100000 taxo 改为 "#OTU ID"

parallel -j 1 '/usr/bin/bash  /work/workspace/zhurj/script/2_metapro/metapipe/fastspar_to_cor_pvalue /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/coabundance/{1}/{1}_count_100000 /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/coabundance/{1} && \
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/cor_pvalue_to_edge_node_list.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/coabundance/{1}/cor_sparcc.fastspar.txt -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/coabundance/{1} --pfile /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/coabundance/{1}/pvalues.tsv --p_cutoff 0.05 --cor_cutoff 0.6 \
' ::: phylum class order family


pyfastx stat /work/rawdata/fastq/metagenome/MNA/MNA002/MNA00233/E1/L1/S1/MNA00233.1.fq.gz

find /work/workspace/zhurj/project/1_metadata/mouse22/clean | sort -u | grep "gz$"  > /work/workspace/zhurj/project/1_metadata/mouse22/input/clean_r1r2_1c.txt
cat /work/rawdata/fastq/metagenome/MNA/MNA002/MNA00233/E1/L1/S1/MNA00233.1.fq.gz | head -n 100 > /work/workspace/zhurj/project/1_metadata/mouse22/tmp/test1.fq
cat /work/rawdata/fastq/metagenome/MNA/MNA002/MNA00233/E1/L1/S1/MNA00233.1.fq.gz | head -n 1000 > /work/workspace/zhurj/project/1_metadata/mouse22/tmp/test2.fq
cat /work/rawdata/fastq/metagenome/MNA/MNA002/MNA00233/E1/L1/S1/MNA00233.1.fq.gz | head -n 10000 > /work/workspace/zhurj/project/1_metadata/mouse22/tmp/test3.fq

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/pyfastx stat /work/workspace/zhurj/project/1_metadata/mouse22/tmp/test.fq

import linecache
import os
import sys
import subprocess
from argparse import ArgumentParser,ArgumentDefaultsHelpFormmater

#infile = '/work/workspace/zhurj/project/1_metadata/mouse22/input/clean_r1r2_1c.txt'
infile = '/work/workspace/zhurj/project/1_metadata/mouse22/tmp/inlist'
odir ='/work/workspace/zhurj/project/1_metadata/mouse22/tmp'
content_line_num = 2
title_line_num = 1
name = 'output'
gc_calculate_pro = '/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/pyfastx'
command = ''

ofilename = "{}.txt".format(name)
tmpfile = os.path.join(odir,'tmpfile')
ofile = os.path.join(odir,ofilename)

files = []
with open(infile,'r',encoding='utf-8') as fp:
  files = fp.readlines()

with open(ofile,'w',encoding='utf-8') as ofp:
  file_path = files[0].strip()
  command = "{} stat {} > {}".format(gc_calculate_pro,file_path,tmpfile)
  subprocess.call(command,shell=True)
  ofp.write(linecache.getline(tmpfile, title_line_num))
  ofp.write(linecache.getline(tmpfile, content_line_num))

  for file_path in files[1:]:
    file_path = file_path.strip()
    command = "{} stat {} > {}".format(gc_calculate_pro,file_path,tmpfile)
    subprocess.call(command,shell=True)
    ofp.write(linecache.getline(tmpfile, content_line_num))

/work/workspace/zhurj/project/1_metadata/mouse22/tmp

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/fastq_stat_gc_length.py -o /work/workspace/zhurj/project/1_metadata/mouse22/qc/GCstat -i /work/workspace/zhurj/project/1_metadata/mouse22/input/clean_r1r2_1c.txt -n sam22GC
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

cd /work/workspace/zhurj/project/1_metadata/p/p
time srun -o qcGCstat.out -e qcGCstat.err -N 1 -c 20 -p slurm128 -w mnclient03 bash qcGCstat.sh &


library(corrplot)
library("PerformanceAnalytics")
pdf(file = "/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/corr_MNO863.pdf",width = 400, height = 350)
Data = read.table("/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species//corr_MNO863_bgin",header=T,row.names=1)
res <- cor(Data)
corrplot(res, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = res, col = col, symm = TRUE)
dev.off()

chart.Correlation(Data, histogram=TRUE, method="spearman", pch=19)

chart.Correlation(Data, histogram=TRUE, method="pearson", pch=19)
dev.off()
#
# https://stackoverflow.com/questions/32220217/correlation-matrix-with-significance-testing-in-r
library(scales)
library(ggplot2)
infile <- '/work/workspace/zhurj/temp/corr_in'
longformData = read.table(infile,header=T)
longformData$stars <- cut(longformData$pValue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), 
                           label=c("***", "**", "*", ""))  # Create column of significance labels
muted_red = muted("red")
muted_blue = muted("blue")
ggplot(longformData, aes(X,Var2))+
  geom_tile(data=longformData, aes(fill=CorrValue), color="white")+
  geom_text(aes(label=stars), color="black", size=5,vjust=-1.5)+
  geom_text(aes(fill = longformData$value, label = round(longformData$CorrValue, 2)))+
  scale_colour_gradient(low = 'blue', high = 'red', limit=c(-1,1), name="Correlation\n(Pearson)")+
  theme(axis.text.x = element_text(size=12,  colour='black'),
        axis.text.y=element_text(colour='black'),
        panel.background=element_rect(colour="black", fill=NA))+
  coord_equal()
================================================================================================================================================================================
Data filter:
$fastpro -i $samplefile{$gr_sample}{'r1'} -o $gr_oread1 -I $samplefile{$gr_sample}{'r2'} -O $gr_oread2 --poly_g_min_len 10 --poly_x_min_len 10 -q 20 -u 40 -n 5 -l 50 -w $thread

python3.6
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/fastqc -o /work/workspace/zhurj/project/1_metadata/mouse22/qc/clean/fastqc -t 10 --quiet /work/workspace/zhurj/project/1_metadata/mouse22/clean/MNA00233.1.fq.gz /work/workspace/zhurj/project/1_metadata/mouse22/clean/MNA00233.2.fq.gz

/work/workspace/liangyj/bin/conda_env/RNASeq/bin/multiqc -o /work/workspace/zhurj/project/1_metadata/mouse22/qc/clean/multiqc -n test /work/workspace/zhurj/project/1_metadata/mouse22/qc/clean/fastqc -f
/work/workspace/liangyj/bin/conda_env/RNASeq/bin/multiqc -o /work/workspace/zhurj/project/1_metadata/mouse22/qc/clean/multiqc -n test -f /work/workspace/zhurj/project/1_metadata/mouse22/qc/clean/fastqc

# qiime2 taxa -> normalization
import pandas as pd
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import os
import sys

infile = '/work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa/genus_tab'
odir = '/work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm'
sam_num = 52
abun_cutoff = 0.005
remove_cols_from_last = 2
samnum_cutoff = int(sam_num/10)
name = "{}_sam{}_abun{}".format('genus',samnum_cutoff,abun_cutoff)
ofile = os.path.join(odir,name)

pre_data = pd.read_csv(infile,sep='\t',header=0,index_col=0)
#*genus_columns,_,_ = pre_data.columns.tolist()
tmp_columns_list = pre_data.columns.tolist()
genus_last_index = -1 * remove_cols_from_last
genus_columns = tmp_columns_list[0:genus_last_index]
mid_data = pre_data[genus_columns]
norm_df = mid_data.div(mid_data.sum(axis=1),axis=0)
norm_df.loc['mean'] = norm_df.mean(axis=0)
norm_df.loc['count'] = (norm_df.iloc[:-1] >= abun_cutoff).astype(int).sum()
# norm_df.loc[:,norm_df.loc['count'] >= samnum_cutoff]
# df.loc[(df['Salary_in_1000']>=100) & (df['Age']< 60) & (df['FT_Team'].str.startswith('S')),['Name','FT_Team']]
filter_norm_df = norm_df.loc[:,(norm_df.loc['count'] >= samnum_cutoff) | (norm_df.loc['mean'] >= abun_cutoff)]
filter_rows = filter_norm_df.index.tolist()
select_rows = filter_rows[0:-2]
result_df = filter_norm_df.loc[select_rows,:]
final_result_df = result_df.T

final_result_df.to_csv(ofile,sep='\t',header=True,index=True)


/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/qiime2_levels_to_filtered_frequency.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa/genus_tab -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm --sample_num 52

# taxa name modification
1. g__xxx 保留 xxx
2. g__, f__ 删除
3. g__，f__xxx 改为 f__xxx_uncl.
4. __(last) 删除
5. g__[xxx] 删除


import os
import sys
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter



infile = "/work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm/genus_sam5_abun0.005"
odir = "/work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm"
title_row = 1
infile_name = os.path.basename(infile)
ofile_name = "name_filter_{}".format(infile_name)
ofile = os.path.join(odir,ofile_name)

title = con = ''
uselines = []

ofp = open(ofile,'w')

with open(infile,'r',encoding='utf-8') as fp:
  lines = fp.readlines()
  uselines = lines;
  if title_row:
    ofp.write(lines[0])
    uselines = lines[1:]
  
  for line in uselines:
    # didn't remove the line break
    pre_levels,*tmp_items = line.split("\t")
    tmp_levels = pre_levels.split(";")
    last_level_mark,last_level_name = tmp_levels[-1].split('__')
    one_before_last_level_mark,one_before_last_level_name = tmp_levels[-2].split('__')
    if not last_level_mark:
      continue
    else:
      if not last_level_name:
        if not one_before_last_level_name:
          continue
        elif one_before_last_level_name.startswith('['):
          continue
        else:
          con = "{}_{}_uncl.\t{}".format(one_before_last_level_mark,one_before_last_level_name,"\t".join(tmp_items))
          ofp.write(con)
      elif last_level_name.startswith('['):
        continue
      else:
        con = "{}\t{}".format(last_level_name,"\t".join(tmp_items))
        ofp.write(con)


ofp.close()


# 转录组数据处理
cd /work/workspace/zhurj/project/1_metadata/mouse22/transcription/input
sed -i '1d' /work/workspace/zhurj/project/1_metadata/mouse22/transcription/input/hfd_mnc863_colon_counts.txt

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values

library(dplyr)
library(Hmisc)

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}

columnf <- '/work/workspace/zhurj/project/1_metadata/mouse22/transcription/corr/input/col_in'
rowf <- '/work/workspace/zhurj/project/1_metadata/mouse22/transcription/corr/input/row_in'
infile <- '/work/workspace/zhurj/project/1_metadata/mouse22/transcription/corr/input/data_in'

data <- read.table(infile,header = TRUE, sep = "\t",row.names= 1)
row_df <- read.table(rowf,sep="\t")
column_df <- read.table(columnf,sep="\t")

#res2 <- rcorr(as.matrix(data),type="spearman")
res2 <- rcorr(as.matrix(data),type="pearson")
format_df <- flattenCorrMatrix(res2$r, res2$P)
format_df$stars <- cut(format_df$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
plot.data <- subset(format_df, row %in% row_df$V1 & column %in% column_df$V1)
#plot.data <- filter(plot.data,stars %in% c('*','***','**'))

# 0.05
plot.data <- filter(plot.data,stars %in% c('***','**','*'))
ofile <- '/work/workspace/zhurj/project/1_metadata/mouse22/transcription/corr/output/corr_meta2_pearson_0.05'
write.table(plot.data,ofile,sep="\t",col.names = TRUE, row.names = FALSE, quote = FALSE)

# 0.01
plot.data <- filter(plot.data,stars %in% c('***','**'))
ofile <- '/work/workspace/zhurj/project/1_metadata/mouse22/transcription/corr/output/corr_meta2_pearson_0.01'
write.table(plot.data,ofile,sep="\t",col.names = TRUE, row.names = FALSE, quote = FALSE)

# 0.001
plot.data <- filter(plot.data,stars %in% c('***'))
ofile <- '/work/workspace/zhurj/project/1_metadata/mouse22/transcription/corr/output/corr_meta2_pearson_0.001'
write.table(plot.data,ofile,sep="\t",col.names = TRUE, row.names = FALSE, quote = FALSE)





/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/33_corr_heatmap/corr_heatmap_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/sam52_taxa_CBG_clinic_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/sam52_taxa_CBG_clinic.pdf -r /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/sam52_taxa_clinic_row -c /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/sam52_taxa_clinic_column --method spearman --plot_width 10 --plot_height 4 --text_size 12


library(dplyr)
kof <- '/work/workspace/zhurj/database/humann/current/kegg/ko_list'
ko_id <- '/work/workspace/zhurj/project/1_metadata/mouse22/transcription/corr/input/ko_id'
select_of <- '/work/workspace/zhurj/project/1_metadata/mouse22/transcription/corr/output/select_Kgene'

ko_ref <- read.table(kof,header = TRUE, sep = "\t",quote = "")
koid_df <- read.table(ko_id,sep="\t")
ko_subset_ref <- ko_ref[c("knum","definition")]
ko_select_df <- filter(ko_subset_ref, knum %in% koid_df$V1)
write.table(ko_select_df,select_of,sep="\t",col.names = TRUE, row.names = FALSE, quote = FALSE)


##tran_liver
library(dplyr)
library(Hmisc)

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}

columnf <- '/work/workspace/zhurj/project/1_metadata/mouse22/transcription/corr/input/col_in'
rowf <- '/work/workspace/zhurj/project/1_metadata/mouse22/transcription/corr/input/tran_liver_row_in'
infile <- '/work/workspace/zhurj/project/1_metadata/mouse22/transcription/corr/input/tran_liver_in'

data <- read.table(infile,header = TRUE, sep = "\t",row.names= 1)
row_df <- read.table(rowf,sep="\t")
column_df <- read.table(columnf,sep="\t")

#res2 <- rcorr(as.matrix(data),type="spearman")
res2 <- rcorr(as.matrix(data),type="pearson")
format_df <- flattenCorrMatrix(res2$r, res2$P)
format_df$stars <- cut(format_df$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
plot.data <- subset(format_df, row %in% row_df$V1 & column %in% column_df$V1)
#plot.data <- filter(plot.data,stars %in% c('*','***','**'))

# 0.05
plot.data <- filter(plot.data,stars %in% c('***','**','*'))
ofile <- '/work/workspace/zhurj/project/1_metadata/mouse22/transcription/corr/output/corr_tran_liver_pearson_0.05'
write.table(plot.data,ofile,sep="\t",col.names = TRUE, row.names = FALSE, quote = FALSE)

# 0.01
plot.data <- filter(plot.data,stars %in% c('***','**'))
ofile <- '/work/workspace/zhurj/project/1_metadata/mouse22/transcription/corr/output/corr_tran_liver_pearson_0.01'
write.table(plot.data,ofile,sep="\t",col.names = TRUE, row.names = FALSE, quote = FALSE)

# 0.001
plot.data <- filter(plot.data,stars %in% c('***'))
ofile <- '/work/workspace/zhurj/project/1_metadata/mouse22/transcription/corr/output/corr_tran_liver_pearson_0.001'
write.table(plot.data,ofile,sep="\t",col.names = TRUE, row.names = FALSE, quote = FALSE)


##tran_testicle
library(dplyr)
library(Hmisc)

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}

columnf <- '/work/workspace/zhurj/project/1_metadata/mouse22/transcription/corr/input/col_in'
rowf <- '/work/workspace/zhurj/project/1_metadata/mouse22/transcription/corr/input/tran_testicle_row_in'
infile <- '/work/workspace/zhurj/project/1_metadata/mouse22/transcription/corr/input/tran_testicle_in'

data <- read.table(infile,header = TRUE, sep = "\t",row.names= 1)
row_df <- read.table(rowf,sep="\t")
column_df <- read.table(columnf,sep="\t")

#res2 <- rcorr(as.matrix(data),type="spearman")
res2 <- rcorr(as.matrix(data),type="pearson")
format_df <- flattenCorrMatrix(res2$r, res2$P)
format_df$stars <- cut(format_df$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
plot.data <- subset(format_df, row %in% row_df$V1 & column %in% column_df$V1)
#plot.data <- filter(plot.data,stars %in% c('*','***','**'))

# 0.05
plot.data <- filter(plot.data,stars %in% c('***','**','*'))
ofile <- '/work/workspace/zhurj/project/1_metadata/mouse22/transcription/corr/output/corr_tran_testicle_pearson_0.05'
write.table(plot.data,ofile,sep="\t",col.names = TRUE, row.names = FALSE, quote = FALSE)

# 0.01
plot.data <- filter(plot.data,stars %in% c('***','**'))
ofile <- '/work/workspace/zhurj/project/1_metadata/mouse22/transcription/corr/output/corr_tran_testicle_pearson_0.01'
write.table(plot.data,ofile,sep="\t",col.names = TRUE, row.names = FALSE, quote = FALSE)

# 0.001
plot.data <- filter(plot.data,stars %in% c('***'))
ofile <- '/work/workspace/zhurj/project/1_metadata/mouse22/transcription/corr/output/corr_tran_testicle_pearson_0.001'
write.table(plot.data,ofile,sep="\t",col.names = TRUE, row.names = FALSE, quote = FALSE)

