# tools-readme
Record readme of often used software
2021.05.21
中山肿瘤43个宏基因组测序样品数据分析



1. input create
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
perl /work/workspace/zhurj/script/5_wgspro/1_filter/mnid2rawdata.pl -i /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid -t 0 -r metagenome -o /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input -n sam43

2. raw data QC (real    88m23.977s)
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
perl /work/workspace/zhurj/script/2_metapro/metapipe/dataqc.pl -i /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/sam53_r1r2_2c.txt -o /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/QC --thread 36 --basecutoff 2.7  --cutoffq20 85 --cutoffq30 80
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

cd /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/p
time srun -o qc_sam43.out -e qc_sam43.err -N 1 -c 20 -p slurm256 -w mnclient01 bash qc_sam43.sh &

3. data filter (real    140m11.605s)
# automatically generate folder: clean
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
perl /work/workspace/zhurj/script/5_wgspro/1_filter/datafilter.pl -i /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/sam43_r1r2_2c.txt -o /work/workspace/zhurj/project/1_metadata/zzuc43_20210521
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

cd /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/p
time srun -o datafilter_sam43.out -e datafilter_sam43.err -N 1 -c 20 -p slurm256 -w mnclient02 bash datafilter_sam43.sh &

4. dehost
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
bwa index -a bwtsw /work/workspace/zhurj/reference/human/GRCh38/hg38.fa.gz
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

cd /work/workspace/zhurj/reference/mouse/GRCm38
time srun -o index.out -e index.err -N 1 -c 20 -p slurm128 -w bash index.sh &

# create folder
/usr/bin/bash /work/workspace/zhurj/script/1_tools/41_linux_tools/createfolder_dir_filename.sh /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/dehost /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid
/usr/bin/bash /work/workspace/zhurj/script/1_tools/41_linux_tools/createfolder_dir_filename.sh /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/kneaddata /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid
/usr/bin/bash /work/workspace/zhurj/script/1_tools/41_linux_tools/createfolder_dir_filename.sh /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/bwarmhost /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid
/usr/bin/bash /work/workspace/zhurj/script/1_tools/41_linux_tools/createfolder_dir_filename.sh /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/metaphlan3 /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid
/usr/bin/bash /work/workspace/zhurj/script/1_tools/41_linux_tools/createfolder_dir_filename.sh /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/bowtie2rmhost /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid
/usr/bin/bash /work/workspace/zhurj/script/1_tools/41_linux_tools/createfolder_dir_filename.sh /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/bwa /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid
/usr/bin/bash /work/workspace/zhurj/script/1_tools/41_linux_tools/createfolder_dir_filename.sh /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/metabat /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid
/usr/bin/bash /work/workspace/zhurj/script/1_tools/41_linux_tools/createfolder_dir_filename.sh /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/prokka /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid
parallel -j 1 'mkdir -p {1}/prokka/{2} ' ::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521 ::: faa ffn
parallel -j 1 'mkdir -p {1}/metaphlan3/{2} ' ::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521 ::: all phylum genus species merge
/usr/bin/bash /work/workspace/zhurj/script/1_tools/41_linux_tools/createfolder_dir_filename.sh /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/humann3 /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid
parallel -j 1 'mkdir -p {1}/humann3/{2} ' ::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521 ::: merge


# remove host DNA with bowtie2 http://www.metagenomics.wiki/tools/short-read/remove-host-sequences
# samtools flag search https://broadinstitute.github.io/picard/explain-flags.html
# bwarmhost (real    752m31.592s)
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel --link -j 2 'bwa mem {4} {1}/{2}.clean.1.fq.gz {1}/{2}.clean.2.fq.gz -t 20 -o {3}/{2}/{2}.sam && \
samtools view -@ 20 -bS {3}/{2}/{2}.sam -o {3}/{2}/{2}.bam && \
samtools view -b -f 12 -F 256  -@ 20 {3}/{2}/{2}.bam -o {3}/{2}/{2}_PEunmapped.bam && \
samtools sort -n -@ 20 {3}/{2}/{2}_PEunmapped.bam -o {3}/{2}/{2}_PEunmapped_sorted.bam && \
bedtools bamtofastq -i {3}/{2}/{2}_PEunmapped_sorted.bam -fq {3}/{2}/{2}_r1.fastq -fq2 {3}/{2}/{2}_r2.fastq && \
rm {3}/{2}/*.sam {3}/{2}/*.bam -f ' \
::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/clean :::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid \
::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/bwarmhost ::: /work/workspace/zhurj/reference/human/GRCh38/hg38.fa.gz
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

bwarmdehost ()
cd /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/p
time srun -o bwarmhost.out -e bwarmhost.err -N 2 -c 20 -p slurm128 -w mnclient03 bash bwarmhost.sh &
-rw-r--r-- 1 zhurj bioinfor   64397286 May 24 22:21 MNC00262_PEunmapped_sorted.bam.tmp.0020.bam
-rw-r--r-- 1 zhurj bioinfor   66010757 May 24 22:21 MNC00262_PEunmapped_sorted.bam.tmp.0021.bam
-rw-r--r-- 1 zhurj bioinfor   65331773 May 24 22:21 MNC00262_PEunmapped_sorted.bam.tmp.0022.bam
-rw-r--r-- 1 zhurj bioinfor   64170304 May 24 22:21 MNC00262_PEunmapped_sorted.bam.tmp.0023.bam
-rw-r--r-- 1 zhurj bioinfor   64113932 May 24 22:21 MNC00262_PEunmapped_sorted.bam.tmp.0024.bam
-rw-r--r-- 1 zhurj bioinfor   64337455 May 24 22:21 MNC00262_PEunmapped_sorted.bam.tmp.0025.bam
-rw-r--r-- 1 zhurj bioinfor   66066876 May 24 22:21 MNC00262_PEunmapped_sorted.bam.tmp.0026.bam
-rw-r--r-- 1 zhurj bioinfor   64486740 May 24 22:21 MNC00262_PEunmapped_sorted.bam.tmp.0027.bam
-rw-r--r-- 1 zhurj bioinfor   64190202 May 24 22:21 MNC00262_PEunmapped_sorted.bam.tmp.0028.bam
-rw-r--r-- 1 zhurj bioinfor   64016283 May 24 22:21 MNC00262_PEunmapped_sorted.bam.tmp.0029.bam
-rw-r--r-- 1 zhurj bioinfor   66134108 May 24 22:21 MNC00262_PEunmapped_sorted.bam.tmp.0030.bam
-rw-r--r-- 1 zhurj bioinfor   65273952 May 24 22:21 MNC00262_PEunmapped_sorted.bam.tmp.0031.bam
-rw-r--r-- 1 zhurj bioinfor   63877785 May 24 22:21 MNC00262_PEunmapped_sorted.bam.tmp.0032.bam
-rw-r--r-- 1 zhurj bioinfor   63891616 May 24 22:21 MNC00262_PEunmapped_sorted.bam.tmp.0033.bam
-rw-r--r-- 1 zhurj bioinfor   65592611 May 24 22:21 MNC00262_PEunmapped_sorted.bam.tmp.0034.bam
-rw-r--r-- 1 zhurj bioinfor   66388810 May 24 22:21 MNC00262_PEunmapped_sorted.bam.tmp.0035.bam
-rw-r--r-- 1 zhurj bioinfor   64061058 May 24 22:21 MNC00262_PEunmapped_sorted.bam.tmp.0036.bam
-rw-r--r-- 1 zhurj bioinfor   63627361 May 24 22:21 MNC00262_PEunmapped_sorted.bam.tmp.0037.bam
-rw-r--r-- 1 zhurj bioinfor   63858588 May 24 22:21 MNC00262_PEunmapped_sorted.bam.tmp.0038.bam
-rw-r--r-- 1 zhurj bioinfor   71061156 May 24 22:21 MNC00262_PEunmapped_sorted.bam.tmp.0039.bam
-rw-r--r-- 1 zhurj bioinfor 7579116567 May 24 22:21 MNC00262_r1.fastq
rm MNC00262_PEunmapped_sorted* -f
-rw-r--r-- 1 zhurj bioinfor   50132456 May 24 22:23 MNC00263_PEunmapped_sorted.bam.tmp.0020.bam
-rw-r--r-- 1 zhurj bioinfor   50099699 May 24 22:23 MNC00263_PEunmapped_sorted.bam.tmp.0021.bam
-rw-r--r-- 1 zhurj bioinfor   50209593 May 24 22:23 MNC00263_PEunmapped_sorted.bam.tmp.0022.bam
-rw-r--r-- 1 zhurj bioinfor   50916561 May 24 22:23 MNC00263_PEunmapped_sorted.bam.tmp.0023.bam
-rw-r--r-- 1 zhurj bioinfor   51464407 May 24 22:23 MNC00263_PEunmapped_sorted.bam.tmp.0024.bam
-rw-r--r-- 1 zhurj bioinfor   50271327 May 24 22:23 MNC00263_PEunmapped_sorted.bam.tmp.0025.bam
-rw-r--r-- 1 zhurj bioinfor   50141963 May 24 22:23 MNC00263_PEunmapped_sorted.bam.tmp.0026.bam
-rw-r--r-- 1 zhurj bioinfor   50095869 May 24 22:23 MNC00263_PEunmapped_sorted.bam.tmp.0027.bam
-rw-r--r-- 1 zhurj bioinfor   50814335 May 24 22:23 MNC00263_PEunmapped_sorted.bam.tmp.0028.bam
-rw-r--r-- 1 zhurj bioinfor   52338861 May 24 22:23 MNC00263_PEunmapped_sorted.bam.tmp.0029.bam
-rw-r--r-- 1 zhurj bioinfor   50050310 May 24 22:23 MNC00263_PEunmapped_sorted.bam.tmp.0030.bam
-rw-r--r-- 1 zhurj bioinfor   49902732 May 24 22:23 MNC00263_PEunmapped_sorted.bam.tmp.0031.bam
-rw-r--r-- 1 zhurj bioinfor   49910212 May 24 22:23 MNC00263_PEunmapped_sorted.bam.tmp.0032.bam
-rw-r--r-- 1 zhurj bioinfor   51221454 May 24 22:23 MNC00263_PEunmapped_sorted.bam.tmp.0033.bam
-rw-r--r-- 1 zhurj bioinfor   52462236 May 24 22:23 MNC00263_PEunmapped_sorted.bam.tmp.0034.bam
-rw-r--r-- 1 zhurj bioinfor   50136573 May 24 22:23 MNC00263_PEunmapped_sorted.bam.tmp.0035.bam
-rw-r--r-- 1 zhurj bioinfor   49808776 May 24 22:23 MNC00263_PEunmapped_sorted.bam.tmp.0036.bam
-rw-r--r-- 1 zhurj bioinfor   49683944 May 24 22:23 MNC00263_PEunmapped_sorted.bam.tmp.0037.bam
-rw-r--r-- 1 zhurj bioinfor   50842854 May 24 22:23 MNC00263_PEunmapped_sorted.bam.tmp.0038.bam
-rw-r--r-- 1 zhurj bioinfor   55710529 May 24 22:23 MNC00263_PEunmapped_sorted.bam.tmp.0039.bam
-rw-r--r-- 1 zhurj bioinfor 6967386125 May 24 22:23 MNC00263_r1.fastq
-rw-r--r-- 1 zhurj bioinfor 6967366365 May 24 22:23 MNC00263_r2.fastq
rm MNC00263_PEunmapped_sorted.bam* -f

source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel -j 2 ' kneaddata  -i {1}/{2}.clean.1.fq.gz -i {1}/{2}.clean.2.fq.gz --output /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/kneaddata/{2} -db /work/workspace/zhurj/reference/human/kneaddata/20200921  --bypass-trim --run-trf -t 20 -p 20 --output-prefix  {2} --reorder  --bowtie2-options "--very-sensitive --dovetail" --remove-intermediate-output'  ::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/clean :::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid
parallel -j 1 ' ln -s /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/kneaddata/{1}/{1}.log /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/kneaddata/log/{1}.log ' :::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid
parallel -j 2 ' cat {1}/{2}/{2}.repeats.removed.1.fastq {1}/{2}/{2}.repeats.removed.2.fastq > {1}/{2}/{2}.fastq ' ::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/kneaddata :::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid
kneaddata_read_count_table --input /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/kneaddata/log --output /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/kneaddata/merge/kneaddata_read_counts.txt
echo `date +"%Y-%m-%d %H:%M:%S"`
echo "Process-4: kneaddata remove host DNA finished"
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

# kneaddata (real    1824m24.649s)
cd /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/p
time srun -o kneaddata.out -e kneaddata.err -N 1 -c 20 -p slurm256 -w mnclient01 bash kneaddata.sh &

# QC clean, bwarmhost
find /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/clean | grep "clean.1" | sort  > /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/clean_read1
find /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/clean | grep "clean.2" | sort  > /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/clean_read2
paste /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/clean_read1 /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/clean_read2 > /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/clean_qc_r1r2.txt
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
perl /work/workspace/zhurj/script/2_metapro/metapipe/dataqc.pl -i /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/clean_qc_r1r2.txt -o /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/QC --thread 36 --basecutoff 2.7  --cutoffq20 85 --cutoffq30 80
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

# clean_qc (real    127m14.262s)
cd /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/p
time srun -o clean_qc.out -e clean_qc.err -N 1 -c 20 -p slurm256 -w mnclient01 bash clean_qc.sh &


find /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/bwarmhost | grep "_r1" > /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/bwarmhost_read1
find /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/bwarmhost | grep "_r2" > /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/bwarmhost_read2
paste /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/bwarmhost_read1 /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/bwarmhost_read2 > /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/bwarmhost_qc_r1r2.txt

# bwarmhost_QC(real    37m24.216s)
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
perl /work/workspace/zhurj/script/2_metapro/metapipe/dataqc.pl -i /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/bwarmhost_qc_r1r2.txt -o /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/QC/bwarmhost --thread 36 --basecutoff 2.7  --cutoffq20 85 --cutoffq30 80
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
bwarmhost_qc(real    33m37.440s)
cd /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/p
time srun -o bwarmhost_qc.out -e bwarmhost_qc.err -N 1 -c 20 -p slurm256 -w mnclient02 bash bwarmhost_qc.sh &

find /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/bowtie2rmhost | grep _R1 | sort -u > /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/bowtie2rmhost_read1
find /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/bowtie2rmhost | grep _R2 | sort -u > /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/bowtie2rmhost_read2
paste /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/bowtie2rmhost_read1 /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/bowtie2rmhost_read2 > /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/bowtie2rmhost_qc_r1r2.txt
perl /work/workspace/zhurj/script/2_metapro/metapipe/dataqc.pl -i /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/bowtie2rmhost_qc_r1r2.txt -o /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/QC/bowtie2rmhost --thread 36 --basecutoff 2.7  --cutoffq20 85 --cutoffq30 80


bowtiermhost
Remove host sequences ref: https://www.metagenomics.wiki/tools/short-read/remove-host-sequences
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/bowtie2 -p 20 -x /work/workspace/zhurj/reference/human/kneaddata/20200921/hg37dec_v0.1 -1 /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/clean/MNC00224.clean.1.fq.gz -2 /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/clean/MNC00224.clean.2.fq.gz --un-conc-gz /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/bowtie2rmhost/MNC00224/MNC00224_host_removed > /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/bowtie2rmhost/MNC00224/MNC00224_mapped_and_unmapped.sam

```
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel -j 1 'bowtie2 -p 20 -x {3} -1 {1}/clean/{2}.clean.1.fq.gz -2 {1}/clean/{2}.clean.2.fq.gz \
--un-conc-gz {1}/bowtie2rmhost/{2}/{2}_host_free > {1}/bowtie2rmhost/{2}/{2}.sam && \
rm {1}/bowtie2rmhost/{2}/{2}.sam -f && \
mv {1}/bowtie2rmhost/{2}/{2}_host_free.1 {1}/bowtie2rmhost/{2}/{2}_host_free_R1.fastq.gz && \
mv {1}/bowtie2rmhost/{2}/{2}_host_free.2 {1}/bowtie2rmhost/{2}/{2}_host_free_R2.fastq.gz ' \
::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521 :::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid \
::: /work/workspace/zhurj/reference/human/kneaddata/20200921/hg37dec_v0.1
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

```
bowtie2rmhost ()
cd /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/p
time srun -o bowtie2rmhost.out -e bowtie2rmhost.err -N 1 -c 20 -p slurm128 -w mnclient03 bash bowtie2rmhost.sh &

# metaphlan3 taxonomy 2021.05.24
# metaphlan3 (real    274m17.789s)
source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel -j 2 --link 'metaphlan {1}/kneaddata/{2}/{2}.repeats.removed.1.fastq,{1}/kneaddata/{2}/{2}.repeats.removed.2.fastq --bowtie2out {1}/metaphlan3/{2}/metagenome.bowtie2.bz2 --input_type fastq --nproc 12 --sample_id_key {2} -o {1}/metaphlan3/{2}/all.txt && \
metaphlan {1}/metaphlan3/{2}/metagenome.bowtie2.bz2 --input_type bowtie2out --tax_lev "p" --nproc 12 --sample_id_key {2} -o {1}/metaphlan3/{2}/phylum.txt && \
metaphlan {1}/metaphlan3/{2}/metagenome.bowtie2.bz2 --input_type bowtie2out --tax_lev "g" --nproc 12 --sample_id_key {2} -o {1}/metaphlan3/{2}/genus.txt && \
metaphlan {1}/metaphlan3/{2}/metagenome.bowtie2.bz2 --input_type bowtie2out --tax_lev "s" --nproc 12 --sample_id_key {2} -o {1}/metaphlan3/{2}/species.txt \
' ::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521 :::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid
parallel -j 1 'ln -s {1}/metaphlan3/{2}/all.txt {1}/metaphlan3/all/{2}.txt && \
ln -s {1}/metaphlan3/{2}/phylum.txt {1}/metaphlan3/phylum/{2}.txt && \
ln -s {1}/metaphlan3/{2}/genus.txt {1}/metaphlan3/genus/{2}.txt && \
ln -s {1}/metaphlan3/{2}/species.txt {1}/metaphlan3/species/{2}.txt \
' ::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521 :::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid
parallel -j 1 --link ' merge_metaphlan_tables.py  {1}/{2}/*.txt > {1}/merge/{2}_merge.txt ' ::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/metaphlan3 ::: all phylum genus species
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

cd /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/p
time srun -o metaphlan3.out -e metaphlan3.err -N 1 -c 20 -p slurm256 -w mnclient02 bash metaphlan3.sh &

# humann
# https://forum.biobakery.org/t/high-proportion-of-unmapped-uniref90-reads-and-very-few-kos-after-regroup/785
# for human tool use uniref90 (/work/workspace/zhurj/database/humann/201901/uniref)

source activate /work/workspace/zhurj/bin/miniconda3/envs/humann
parallel -j 1 'humann --threads 20 --input {1}/kneaddata/{2}/{2}.fastq --output {1}/humann3/{2} --output-basename {2}_metacyc --protein-database {3} && \
rm {1}/humann3/{2}/{2}_metacyc_humann_temp -rf && \
humann_regroup_table -i {1}/humann3/{2}/{2}_metacyc_genefamilies.tsv -g uniref90_ko --output  {1}/humann3/{2}/{2}_uniref90_ko.tsv && \
humann_regroup_table -i {1}/humann3/{2}/{2}_metacyc_genefamilies.tsv -g uniref90_go --output  {1}/humann3/{2}/{2}_uniref90_go.tsv && \
humann_regroup_table -i {1}/humann3/{2}/{2}_metacyc_genefamilies.tsv -g uniref90_level4ec --output  {1}/humann3/{2}/{2}_uniref90_level4ec.tsv && \
humann_regroup_table -i {1}/humann3/{2}/{2}_metacyc_genefamilies.tsv -g uniref90_pfam --output  {1}/humann3/{2}/{2}_uniref90_pfam.tsv && \
humann_regroup_table -i {1}/humann3/{2}/{2}_metacyc_genefamilies.tsv -g uniref90_eggnog --output  {1}/humann3/{2}/{2}_uniref90_eggnog.tsv \
' ::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521 :::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid ::: /work/workspace/zhurj/database/humann/201901/uniref
parallel -j 1 'humann -i {1}/humann3/{2}/{2}_uniref90_ko.tsv --pathways-database {3} --output {1}/humann3/{2} --threads 20 --output-basename {2}_kegg' \
 ::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521 :::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid ::: /work/workspace/zhurj/database/humann/current/kegg/keggc
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/humann

# uniref 需要diamond 0.9.24 需要将 0.8.38 改为 0.9.24
# /work/workspace/zhurj/bin/diamond --version => diamond version 0.9.24
# /work/workspace/zhurj/software/diamondv0.8.38/diamond --version => diamond version 0.8.38
cd /work/workspace/zhurj/database/humann/201901/uniref/ 
/work/workspace/zhurj/bin/diamond getseq --db uniref90_201901.dmnd > uniref90_201901.faa
/work/workspace/zhurj/software/diamondv0.8.38/diamond getseq --db uniref90_201901.dmnd > uniref90_201901.faa
/work/workspace/zhurj/bin/diamond makedb --in uniref90_201901.faa --db uniref90_201901

humann（real    3353m28.632s）
cd /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/p
time srun -o humann3.out -e humann3.err -N 1 -c 20 -p slurm256 -w mnclient01 bash humann3.sh &

# QC, merge result of read1, and read2
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/3_qcstat/qcstat_r1r2_merge.py -i /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/QC/cleanqc.txt -o /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/QC/merge_cleanqc.txt --titleline --namepos -5
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/3_qcstat/qcstat_r1r2_merge.py -i /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/QC/rawdata.txt -o /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/QC/merge_rawdata.txt --titleline --namepos -4
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/3_qcstat/qcstat_r1r2_merge.py -i /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/QC/bwarmhost.txt -o /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/QC/merge_bwarmhost.txt --titleline --namepos -3
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/3_qcstat/qcstat_r1r2_merge.py -i /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/QC/bwarmhost/output.txt -o /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/QC/bwarmhost/merge_bwarmhost.txt --titleline --namepos -3

find /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/bwarmhost | grep "_r1" > /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/bwarmhost_read1
find /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/bwarmhost | grep "_r2" > /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/bwarmhost_read2
paste /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/bwarmhost_read1 /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/bwarmhost_read2 > /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/bwarmhost_qc_r1r2.txt

# bwarmhost_QC(real    37m24.216s)
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
perl /work/workspace/zhurj/script/2_metapro/metapipe/dataqc.pl -i /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/bwarmhost_qc_r1r2.txt -o /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/QC/bwarmhost --thread 36 --basecutoff 2.7  --cutoffq20 85 --cutoffq30 80
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

## 组装物种鉴定

# create folder
/usr/bin/bash /work/workspace/zhurj/script/1_tools/41_linux_tools/createfolder_dir_filename.sh /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/spades /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid
mkdir -p /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/spades/tmp


1. contig assemble
spades.sh
删除中间文件只保留每个样本的 contigs.fasta scaffolds.fasta

```
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
mkdir -p /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/spades/tmp
parallel -j 1 ' \
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/moonis/bin/genoassemble/soft/SPAdes-3.13.1-Linux/bin/spades.py --meta -1 {1}/bowtie2rmhost/{2}/{2}_host_free_R1.fastq.gz \
-2 {1}/bowtie2rmhost/{2}/{2}_host_free_R2.fastq.gz -t 20 -o {1}/spades/{2}  && \
rm {1}/spades/tmp/* -rf && \
mv {1}/spades/{2}/contigs.fasta  {1}/spades/tmp && \
mv {1}/spades/{2}/scaffolds.fasta  {1}/spades/tmp && \
rm {1}/spades/{2}/* -rf && \
mv {1}/spades/tmp/* {1}/spades/{2} ' \
::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521 :::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

test.sh (单个样品 spades的测试 real    156m12.625s)
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
mkdir -p /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/spades/tmp
parallel -j 1 ' \
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/moonis/bin/genoassemble/soft/SPAdes-3.13.1-Linux/bin/spades.py --meta -1 {1}/bowtie2rmhost/{2}/{2}_host_free_R1.fastq.gz \
-2 {1}/bowtie2rmhost/{2}/{2}_host_free_R2.fastq.gz -t 20 -o {1}/spades/{2}  && \
rm {1}/spades/tmp/* -rf && \
mv {1}/spades/{2}/contigs.fasta  {1}/spades/tmp && \
mv {1}/spades/{2}/scaffolds.fasta  {1}/spades/tmp && \
rm {1}/spades/{2}/* -rf && \
mv {1}/spades/tmp/* {1}/spades/{2} ' \
::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521 ::: MNC00223
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

```
cd /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/p
time srun -o spades_kneaddata.out -e spades_kneaddata.err -N 1 -c 20 -p slurm256 -w mnclient02 bash spades_kneaddata.sh &

# 单个样品spades组装测试
# 33个样品运行时间（61h14m）
time srun -o test.out -e test.err -N 1 -c 20 -p slurm256 -w mnclient02 bash test.sh &
time srun -o spades_bowtie2.out -e spades_bowtie2.err -N 1 -c 20 -p slurm256 -w mnclient02 bash spades_bowtie2.sh &
# spades_bowtie2.sh 未运行275-285
# 从新spades运行275-285， 20210604 (real    833m37.664s) 
# metaSPAdes stuck on read error correction: https://github.com/ablab/spades/issues/152
# spades add “--meta --only-assembler” to solve the problem "stucking on read error correction"
time srun -o spades_bowtie2_275to285.out -e spades_bowtie2_275to285.err -N 1 -c 20 -p slurm256 -w mnclient02 bash spades_bowtie2_275to285.sh &


# metabat 基因组组组装
metabat (real    589m25.422s)

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel -j 1 'bwa index {1}/spades/{2}/contigs.fasta && \
bwa mem -t {3} {1}/spades/{2}/contigs.fasta {1}/bowtie2rmhost/{2}/{2}_host_free_R1.fastq.gz {1}/bowtie2rmhost/{2}/{2}_host_free_R2.fastq.gz > {1}/bwa/{2}/{2}.sam && \
samtools view -@ {3} -S -b {1}/bwa/{2}/{2}.sam -o {1}/bwa/{2}/{2}.bam && \
samtools sort -@ {3} {1}/bwa/{2}/{2}.bam -o {1}/bwa/{2}/{2}.sorted.bam && \
rm {1}/bwa/{2}/{2}.sam {1}/bwa/{2}/{2}.bam -f  && \
jgi_summarize_bam_contig_depths --outputDepth {1}/metabat/{2}/depth.txt  {1}/bwa/{2}/{2}.sorted.bam && \
metabat2 -i {1}/spades/{2}/contigs.fasta -a {1}/metabat/{2}/depth.txt -o {1}/metabat/{2}/{2} -m 2000 ' \
::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521 :::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid ::: 20
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

time srun -o metabat.out -e metabat.err -N 1 -c 20 -p slurm256 -w mnclient01 bash metabat.sh &

gene prediction 如何加速****************
gene prediction
prokka-cdhit 
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel -j 1 'seqtk seq -L 500 {1}/spades/{2}/scaffolds.fasta > {1}/spades/{2}/scaffolds_500.fasta ' \
::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521 :::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid
source deactivate  /work/workspace/zhurj/bin/miniconda3/envs/python3.6
source activate /work/workspace/zhurj/bin/miniconda3/envs/prokka
parallel -j 1 'prokka --outdir {1}/prokka/{2} -cpus 20 --prefix {2} --force --metagenome --compliant --locustag {2} {1}/spades/{2}/scaffolds_500.fasta && \
ln -s {1}/prokka/{2}/{2}.ffn {1}/prokka/ffn ' \
::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521 :::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/prokka

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
parallel -j 1 'mkdir -p {1}/cdhit/{2} ' ::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521 ::: in out
cat /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/prokka/ffn/*.ffn > /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/cdhit/in/pangene.ffn
parallel -j 1 'cd-hit-est -i {1}/cdhit/in/pangene.ffn -o {1}/cdhit/out/cdhit.fna -c 0.95 -G 0 -aS 0.9 -g 1 -d 0 -T 20 '  ::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

time srun -o prokka-cdhit.out -e prokka-cdhit.err -N 1 -c 20 -p slurm128 -w mnclient03 bash prokka-cdhit.sh &

# prokka运行速度太慢，增加参数 --fast 加速prokka，并将样品拆分为3个tasks提交

/usr/bin/bash /work/workspace/zhurj/script/1_tools/41_linux_tools/createfolder_dir_filename.sh /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/prokka261_272 /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid261_272
# prokka261_272.sh
source activate /work/workspace/zhurj/bin/miniconda3/envs/prokka
parallel -j 1 'prokka --outdir {1}/prokka261_272/{2} -cpus 20 --prefix {2} --force --metagenome --fast --compliant --locustag {2} {1}/spades/{2}/scaffolds_500.fasta && \
ln -s {1}/prokka261_272/{2}/{2}.ffn {1}/prokka/ffn ' \
::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521 :::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid261_272
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/prokka

/usr/bin/bash /work/workspace/zhurj/script/1_tools/41_linux_tools/createfolder_dir_filename.sh /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/prokka273_285 /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid275_285
# prokka273_285.sh (real    4075m11.983s)
source activate /work/workspace/zhurj/bin/miniconda3/envs/prokka
parallel -j 1 'prokka --outdir {1}/prokka273_285/{2} -cpus 20 --prefix {2} --force --metagenome --fast --compliant --locustag {2} {1}/spades/{2}/scaffolds_500.fasta && \
ln -s {1}/prokka273_285/{2}/{2}.ffn {1}/prokka/ffn ' \
::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521 :::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid275_285
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/prokka

/work/workspace/zhurj/project/1_metadata/zzuc43_20210521/p
# prokka261_272.sh(real    5154m33.265s)
time srun -o prokka261_272.out -e prokka261_272.err -N 1 -c 20 -p slurm256 -w mnclient01 bash prokka261_272.sh &
time srun -o prokka273_285.out -e prokka273_285.err -N 1 -c 20 -p slurm256 -w mnclient02 bash prokka273_285.sh &

prokka280_285.sh(real    4442m27.002s)
parallel -j 1 'seqtk seq -L 500 {1}/spades/{2}/scaffolds.fasta > {1}/spades/{2}/scaffolds_500.fasta ' \
::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521 ::: MNC00280 MNC00281 MNC00283 MNC00284 MNC00285
time srun -o prokka280_285.out -e prokka280_285.err -N 1 -c 20 -p slurm256 -w mnclient01 bash prokka280_285.sh &

parallel -j 1 'mkdir -p {1}/{2}' ::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/prokka230_260 :::: /work/workspace/zhurj/project/1_metadata/zzuc43_20210521/input/mnid230_260
time srun -o prokka230_260.out -e prokka230_260.err -N 1 -c 20 -p slurm256 -w mnclient01 bash prokka230_260.sh &
time srun -o prokka273_274.out -e prokka273_274.err -N 1 -c 20 -p slurm128 -w mnclient03 bash prokka273_274.sh &
time srun -o prokka223_253.out -e prokka223_253.err -N 1 -c 20 -p slurm256 -w mnclient01 bash prokka223_253.sh &
time srun -o prokka254_260.out -e prokka254_260.err -N 1 -c 20 -p slurm256 -w mnclient02 bash prokka254_260.sh &
