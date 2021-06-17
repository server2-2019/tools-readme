# tools-readme
Record readme of often used software
20200914
nanopore 介绍： http://tiramisutes.github.io/2019/03/12/Oxford-Nanopore.html
Nanopore -long sequencing - data analysis

#conda remove -n nanopore --all

conda create -n nanopore python=3.6
/work/workspace/zhurj/bin/miniconda3/envs/nanopore

source activate py2
# https://anaconda.org/bioconda/poretools
conda install -c bioconda poretools
 v0.6.1a1

source activate nanopore
# https://anaconda.org/bioconda/porechop
conda install -c bioconda porechop
v0.2.4

source activate nanopore
conda install -c bioconda filtlong
# https://anaconda.org/bioconda/filtlong
v0.2.0

source activate nanopore
conda install -c bioconda flye
# https://anaconda.org/bioconda/flye
v2.8.1

source activate nanopore
conda install -c bioconda canu
# https://anaconda.org/bioconda/canu
v2.0

source activate nanopore
conda install -c bioconda unicycler
# https://anaconda.org/bioconda/unicycler
v0.4.8

source activate nanopore
conda install -c bioconda nanoplot
# https://anaconda.org/bioconda/nanoplot
v1.20.0

source activate nanopore
conda install -c bioconda nanoqc nanostat nanofilt

source activate nanopore
conda install -c bioconda minimap2
# https://anaconda.org/bioconda/minimap2
v2.17


source activate nanopore
# https://anaconda.org/bioconda/racon
# https://github.com/lbcb-sci/racon
conda install -c bioconda racon
v1.4.13


source activate nanopore
# https://anaconda.org/bioconda/medaka
# https://github.com/nanoporetech/medaka
conda install -c bioconda medaka
v1.0.3

nanopore测序技术分析专栏
# https://zhuanlan.zhihu.com/p/102392867
nanopore测序技术专题（十四）：nanopore测序质量怎么样
# https://zhuanlan.zhihu.com/p/100217613
nanopore测序技术专题（十五）：利用NanoPlot进行数据质控
https://zhuanlan.zhihu.com/p/100218081
conda update -n nanopore NanoPlot

cd 
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/pigz SRR8651554_1.fastq
mkdir -p /work/workspace/zhurj/project/13_nanopore/kejing20200914/MN214516/nanoQC
source activate nanopore
nanoQC  -o /work/workspace/zhurj/project/13_nanopore/kejing20200914/MN214516/nanoQC/raw /work/rawdata/run/guangzhou/novogene/2020/09/20200908/run00062/rawdata/X101SC20050850-Z01-J002/Result/01.Cleandata/MN214516/barcode21.fastq.gz

#nanoQC  -l 1000 -o /work/workspace/zhurj/project/13_nanopore/kejing20200914/MN214516/nanoQC/min1000 /work/rawdata/run/guangzhou/novogene/2020/09/20200908/run00062/rawdata/X101SC20050850-Z01-J002/Result/01.Cleandata/MN214516/barcode21.fastq.gz
#NanoStat --fastq /work/rawdata/run/guangzhou/novogene/2020/09/20200908/run00062/rawdata/X101SC20050850-Z01-J002/Result/01.Cleandata/MN214516/barcode21.fastq.gz -o /work/workspace/zhurj/project/13_nanopore/kejing20200914/MN214516/nanoStat/all --tsv -t 36  -n MN214516_rawdata

NanoPlot --fastq /work/rawdata/run/guangzhou/novogene/2020/09/20200908/run00062/rawdata/X101SC20050850-Z01-J002/Result/01.Cleandata/MN214516/barcode21.fastq.gz  -t 36 --plots hex dot -o /work/workspace/zhurj/project/13_nanopore/kejing20200914/MN214516/nanoPlot/raw

# nanopore测序技术专题（十六）：利用NanoFilt对数据进行过滤
# https://mp.weixin.qq.com/s?__biz=MzI2MjA1MDQxMg==&mid=2649709950&idx=1&sn=3fd58d3b8c24e036d65d90b00ed0a829&chksm=f24afa7dc53d736b1ddf6bb223d8a1e5d64c21d9f31f3a12c058b4916a0c5b985959aecc6e30&scene=21#wechat_redirect
gunzip -c /work/rawdata/run/guangzhou/novogene/2020/09/20200908/run00062/rawdata/X101SC20050850-Z01-J002/Result/01.Cleandata/MN214516/barcode21.fastq.gz | NanoFilt -q 10 -l 1000 --headcrop 50 --tailcrop  50| gzip > /work/workspace/zhurj/project/13_nanopore/kejing20200914/MN214516/clean/MN214516.NanoFilt.fastq.gz
# clean data QC
nanoQC  -o /work/workspace/zhurj/project/13_nanopore/kejing20200914/MN214516/nanoQC/clean /work/workspace/zhurj/project/13_nanopore/kejing20200914/MN214516/clean/MN214516.NanoFilt.fastq.gz
NanoPlot --fastq /work/workspace/zhurj/project/13_nanopore/kejing20200914/MN214516/clean/MN214516.NanoFilt.fastq.gz  -t 36 --plots hex dot -o /work/workspace/zhurj/project/13_nanopore/kejing20200914/MN214516/nanoPlot/clean

# 文献笔记三十三：结合二代三代测序数据组装叶绿体基因组
# https://cloud.tencent.com/developer/article/1593316

# 结合二代三代数据
#unicycler -1 short_reads_1.fastq -2 short_reads_2.fastq -l long_reads_high_depth.fastq -o output_dir -t 16

flye --nano-raw /work/workspace/zhurj/project/13_nanopore/kejing20200914/MN214516/clean/MN214516.NanoFilt.fastq.gz -g 6m -o /work/workspace/zhurj/project/13_nanopore/kejing20200914/MN214516/flye/assemble -t 36 

# 三种长片段测序方法对比
#https://www.mgitech.cn/news/255/

source activate python3.6
quast -t 16 -o /work/workspace/zhurj/project/13_nanopore/kejing20200914/MN214516/flye/quast /work/workspace/zhurj/project/13_nanopore/kejing20200914/MN214516/flye/assemble/assembly.fasta
checkm lineage_wf /work/workspace/zhurj/project/13_nanopore/kejing20200914/MN214516/flye/assemble /work/workspace/zhurj/project/13_nanopore/kejing20200914/MN214516/flye/checkm -x fasta -t 36 --tab_table -f /work/workspace/zhurj/project/13_nanopore/kejing20200914/MN214516/flye/checkm/MN214516_checkm.txt

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/pigz -p 36 MN210742.fast5
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/pigz -p 36 MN210742.fastq
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/pigz -p 36 barcode19.fast5
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/pigz -p 36 barcode19.fastq
cd /work/rawdata/run/guangzhou/novogene/2020/09/20200908/run00062/rawdata/X101SC20050850-Z01-J002/Result/01.Cleandata/MN27565
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/pigz -p 36 barcode20.fast5
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/pigz -p 36 barcode20.fastq
cd /work/rawdata/run/guangzhou/novogene/2020/09/20200908/run00062/rawdata/X101SC20050850-Z01-J002/Result/01.Cleandata/MN235798
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/pigz -p 36 barcode16.fast5
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/pigz -p 36 barcode16.fastq
cd /work/rawdata/run/guangzhou/novogene/2020/09/20200908/run00062/rawdata/X101SC20050850-Z01-J002/Result/01.Cleandata/MN230415
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/pigz -p 36 barcode18.fast5
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/pigz -p 36 barcode18.fastq
cd /work/rawdata/run/guangzhou/novogene/2020/09/20200908/run00062/rawdata/X101SC20050850-Z01-J002/Result/01.Cleandata/MN230352
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/pigz -p 36 barcode17.fast5
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/pigz -p 36 barcode17.fastq
cd /work/rawdata/run/guangzhou/novogene/2020/09/20200908/run00062/rawdata/X101SC20050850-Z01-J002/Result/01.Cleandata/MN214516
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/pigz -p 36 barcode21.fast5

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/pigz -p 36 /work/rawdata/run/guangzhou/novogene/2020/09/20200908/run00062/rawdata/X101SC20050850-Z01-J003/Result/01.Cleandata/MN216370/MN216370.fast5
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/pigz -p 36 /work/rawdata/run/guangzhou/novogene/2020/09/20200908/run00062/rawdata/X101SC20050850-Z01-J003/Result/01.Cleandata/MN221662/MN221662.fast5
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/pigz -p 36 /work/rawdata/run/guangzhou/novogene/2020/09/20200908/run00062/rawdata/X101SC20050850-Z01-J003/Result/01.Cleandata/MN214635/MN214635.fast5
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/pigz -p 36 /work/rawdata/run/guangzhou/novogene/2020/09/20200908/run00062/rawdata/X101SC20050850-Z01-J003/Result/01.Cleandata/MN216370/MN216370.fastq
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/pigz -p 36 /work/rawdata/run/guangzhou/novogene/2020/09/20200908/run00062/rawdata/X101SC20050850-Z01-J003/Result/01.Cleandata/MN221662/MN221662.fastq
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/pigz -p 36 /work/rawdata/run/guangzhou/novogene/2020/09/20200908/run00062/rawdata/X101SC20050850-Z01-J003/Result/01.Cleandata/MN214635/MN214635.fastq

/work/workspace/zhurj/lib/sh/nanoporeV1 /work/workspace/zhurj/project/13_nanopore/kejing20200914/input/test 36

srun -o MN230352.out -e MN230352.err -N 1 -c 20 -p slurm256 bash /work/workspace/zhurj/project/13_nanopore/kejing20200914/MN230352/p/process.sh &

/work/workspace/zhurj/lib/sh/nanoporeV1 /work/workspace/zhurj/project/13_nanopore/kejing20200914/input/sam11 36

srun -o sam10.out -e sam10.err -N 1 -c 20 -p slurm256 bash /work/workspace/zhurj/project/13_nanopore/kejing20200914/p/nanopore10.sh &

source activate python3.6
fastANI -q /work/workspace/zhurj/project/13_nanopore/kejing20200914/Q10/MN214516/flye/assemble/assembly.fasta -r /work/database/ncbi/genome/all/GCA/GCA_900/GCA_900103/GCA_900103445/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/13_nanopore/kejing20200914/Q10/MN214516/fastANI/MN214516_GCA_900103445

gtdbtk classify_wf --batchfile /work/workspace/zhurj/project/13_nanopore/kejing20200914/input/gtdbtk_test --out_dir /work/workspace/zhurj/project/13_nanopore/kejing20200914/GTDBTK/test --cpus 36

/work/workspace/zhurj/project/13_nanopore/kejing20200914/Q10/MN214516/flye/assemble/assembly.fasta	MN214516
/work/workspace/zhurj/project/13_nanopore/kejing20200914/MN230352/flye/assemble/assembly.fasta	MN230352

/work/workspace/zhurj/project/13_nanopore/kejing20200914/MN230352/flye/assemble/assembly.fasta	MN230352
/work/workspace/zhurj/project/13_nanopore/kejing20200914/MN214715/flye/assemble/assembly.fasta	MN214715
/work/workspace/zhurj/project/13_nanopore/kejing20200914/MN230415/flye/assemble/assembly.fasta	MN230415
/work/workspace/zhurj/project/13_nanopore/kejing20200914/MN29329/flye/assemble/assembly.fasta	MN29329
/work/workspace/zhurj/project/13_nanopore/kejing20200914/MN235798/flye/assemble/assembly.fasta	MN235798
/work/workspace/zhurj/project/13_nanopore/kejing20200914/MN27565/flye/assemble/assembly.fasta	MN27565
/work/workspace/zhurj/project/13_nanopore/kejing20200914/MN216370/flye/assemble/assembly.fasta	MN216370
/work/workspace/zhurj/project/13_nanopore/kejing20200914/MN221662/flye/assemble/assembly.fasta	MN221662
/work/workspace/zhurj/project/13_nanopore/kejing20200914/MN210742/flye/assemble/assembly.fasta	MN210742
/work/workspace/zhurj/project/13_nanopore/kejing20200914/MN214635/flye/assemble/assembly.fasta	MN214635
/work/workspace/zhurj/project/13_nanopore/kejing20200914/MN214516/flye/assemble/assembly.fasta	MN214516

gtdbtk classify_wf --batchfile /work/workspace/zhurj/project/13_nanopore/kejing20200914/input/gtdbtk_sam11 --out_dir /work/workspace/zhurj/project/13_nanopore/kejing20200914/GTDBTK/sam11_20200916 --cpus 36

srun -o gtdbtksam11.out -e gtdbtksam11.err -N 1 -c 20 -p slurm256 bash /work/workspace/zhurj/project/13_nanopore/kejing20200914/p/gtdbtk_sam11.sh &


MN230352
MN214715
MN230415
MN29329
MN235798
MN27565
MN216370
MN221662
MN210742
MN214635
MN214516
2020.09.15
nanopore 数据质控，及基因组组装当前可以完成，需要对结果进行整理


find `pwd` | grep NanoStats.txt | grep -v Q10 | awk -F "/" '{print $0"\t"$(NF-3)"\t"$(NF-1)}' > /work/workspace/zhurj/project/13_nanopore/kejing20200914/input/in_mergenanoPlot

python /work/workspace/zhurj/lib/python/script/mergenanoPlot.py /work/workspace/zhurj/project/13_nanopore/kejing20200914/input/test.txt

python /work/workspace/zhurj/lib/python/script/mergenanoPlot.py /work/workspace/zhurj/project/13_nanopore/kejing20200914/input/in_mergenanoPlot /work/workspace/zhurj/project/13_nanopore/kejing20200914/merge/nanoPlot/sam11_nanoplot

python /work/workspace/zhurj/lib/python/script/mergequast.py /work/workspace/zhurj/project/13_nanopore/kejing20200914/input/test /work/workspace/zhurj/project/13_nanopore/kejing20200914/merge/quast/test

python /work/workspace/zhurj/lib/python/script/mergequast.py /work/workspace/zhurj/project/13_nanopore/kejing20200914/input/in_mergequast /work/workspace/zhurj/project/13_nanopore/kejing20200914/merge/quast/quast_sam11

cd /work/workspace/zhurj/project/13_nanopore/kejing20200914
find `pwd` | grep check.txt | grep -v Q10 | awk -F "/" '{print $0"\t"$(NF-3)}' > /work/workspace/zhurj/project/13_nanopore/kejing20200914/input/in_mergecheckm

python /work/workspace/zhurj/lib/python/script/mergecheckm.py /work/workspace/zhurj/project/13_nanopore/kejing20200914/input/in_mergecheckm /work/workspace/zhurj/project/13_nanopore/kejing20200914/merge/checkm/checkm_sam11


2020.12.18
40 nanopore sequencing data analysis

source activate /work/workspace/zhurj/bin/miniconda3/envs/nanopore
mkdir -p /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN211200/nanoQC/raw
mkdir -p /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN211200/nanoQC/clean
mkdir -p /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN211200/clean
mkdir -p /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN211200/flye/assemble
mkdir -p /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN211200/flye/quast
mkdir -p /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN211200/flye/checkm

/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/nanoQC -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN211200/nanoQC/raw /work/rawdata/run/guangzhou/novogene/2020/10/20201022/run00072/rawdata/MN211200/MN211200.fastq
/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/NanoPlot --fastq /work/rawdata/run/guangzhou/novogene/2020/10/20201022/run00072/rawdata/MN211200/MN211200.fastq -t 36 --plots hex dot -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN211200/nanoQC/raw -p MN211200
/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/NanoFilt -q 8 -l 1000 --headcrop 50 --tailcrop 50 /work/rawdata/run/guangzhou/novogene/2020/10/20201022/run00072/rawdata/MN211200/MN211200.fastq | /work/workspace/zhurj/bin/pigz > /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN211200/clean/nanoFilt.fastq.gz
/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/nanoQC -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN211200/nanoQC/clean /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN211200/clean/nanoFilt.fastq.gz
/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/NanoPlot --fastq /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN211200/clean/nanoFilt.fastq.gz -t 36 --plots hex dot -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN211200/nanoQC/clean -p MN211200
/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/flye --nano-raw /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN211200/clean/nanoFilt.fastq.gz -g 6m -t 36 -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN211200/flye/assemble

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/quast -t 36 -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN211200/flye/quast /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN211200/flye/assemble/assembly.fasta
/work/workspace/zhurj/bin/miniconda3/bin/checkm lineage_wf /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN211200/flye/assemble /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN211200/flye/checkm -x fasta -t 36 --tab_table -f /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN211200/flye/checkm/check.txt

/work/rawdata/run/guangzhou/novogene/2020/10/20201022/run00072/rawdata/MN231199/barcode21.fastq
mkdir -p /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231199/nanoQC/raw
mkdir -p /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231199/nanoQC/clean
mkdir -p /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231199/clean
mkdir -p /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231199/flye/assemble
mkdir -p /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231199/flye/quast
mkdir -p /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231199/flye/checkm

/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/nanoQC -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231199/nanoQC/raw /work/rawdata/run/guangzhou/novogene/2020/10/20201022/run00072/rawdata/MN231199/barcode21.fastq
/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/NanoPlot --fastq /work/rawdata/run/guangzhou/novogene/2020/10/20201022/run00072/rawdata/MN231199/barcode21.fastq -t 36 --plots hex dot -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231199/nanoQC/raw -p MN231199
/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/NanoFilt -q 8 -l 1000 --headcrop 50 --tailcrop 50 /work/rawdata/run/guangzhou/novogene/2020/10/20201022/run00072/rawdata/MN231199/barcode21.fastq | /work/workspace/zhurj/bin/pigz > /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231199/clean/nanoFilt.fastq.gz
/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/nanoQC -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231199/nanoQC/clean /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231199/clean/nanoFilt.fastq.gz
/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/NanoPlot --fastq /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231199/clean/nanoFilt.fastq.gz -t 36 --plots hex dot -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231199/nanoQC/clean -p MN231199
/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/flye --nano-raw /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231199/clean/nanoFilt.fastq.gz -g 6m -t 36 -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231199/flye/assemble

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/quast -t 36 -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231199/flye/quast /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231199/flye/assemble/assembly.fasta
/work/workspace/zhurj/bin/miniconda3/bin/checkm lineage_wf /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231199/flye/assemble /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231199/flye/checkm -x fasta -t 36 --tab_table -f /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231199/flye/checkm/check.txt

source deactivate /work/workspace/zhurj/bin/miniconda3/envs/nanopore


find /work/rawdata/run/guangzhou/novogene/2020/10/20201022/run00072/rawdata | grep -P "fastq$" | grep -v "MN211200" | grep -v "MN231199" | awk -F "/" '{print $0"\t"$(NF-1)}' > /work/workspace/zhurj/project/13_nanopore/kejing20201218/input/input38

from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
import os
import sys

def createfolder(odir,second,third):
  if third:
    folder = os.path.join(odir,second,third)
  else:
    folder = os.path.join(odir,second)
  com = "mkdir -p {}\n".format(folder)
  ofp.write(com)

infile = '/work/workspace/zhurj/project/13_nanopore/kejing20201218/input/test2'
odir = '/work/workspace/zhurj/project/13_nanopore/kejing20201218'
obashf = '/work/workspace/zhurj/project/13_nanopore/kejing20201218/p/sam38.sh'
thread = 36
flye_g = '6m'
nanofilt_q = 8
nanofilt_l = 1000
nanofilt_headcrop = 50
nanofilt_tailcrop = 50

tools = {
  'nanoQC':'/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/nanoQC',
  'NanoPlot':'/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/NanoPlot',
  'NanoFilt':'/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/NanoFilt',
  'flye':'/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/flye',
  'quast':'/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/quast',
  'checkm':'/work/workspace/zhurj/bin/miniconda3/bin/checkm',
  'pigz':'/work/workspace/zhurj/bin/pigz'
}

with open(obashf,'w',encoding='utf-8') as ofp, open(infile,'r',encoding='utf-8') as fp:
  command = "source activate /work/workspace/zhurj/bin/miniconda3/envs/nanopore\n"
  ofp.write(command)
  lines = fp.readlines()
  for line in lines:
    if not line or line.startswith('#'):
      continue
    [ifile,sample] = line.strip().split('\t')
    # mkdir folder
    
    createfolder(odir,'nanoQC','raw')
    createfolder(odir,'nanoQC','clean')
    createfolder(odir,'clean','')
    createfolder(odir,'flye','assemble')
    createfolder(odir,'flye','quast')
    createfolder(odir,'flye','checkm')
    
    # raw data nanoQC
    rawqcfolder = os.path.join(odir,'nanoQC','raw')
    command = "{} {} {}\n".format(tools['nanoQC'],rawqcfolder,ifile)
    ofp.write(command)

    # raw nanoplot
    command = "{program} --fastq {in} -t {thread} --plots hex dot -o {out} -p {pre}\n".format(program=tools['NanoPlot'],in=ifile,thread=thread,out=rawqcfolder,pre=sample)
    ofp.write(command)

    # data filter
    clean_fq = os.path.join(odir,'clean','nanoFilt.fastq.gz')
    command = "{program} -q {qual} -l {length} --headcrop {headcut} --tailcrop {tailcut} {in_} | {pigz} > {out}\n".format(program=tools['NanoFilt'],qual=nanofilt_q,length=nanofilt_l,headcut=nanofilt_headcrop,tailcut=nanofilt_tailcrop,in_=ifile,pigz=tools['NanoFilt'],out=clean_fq)
    ofp.write(command)

    # clean data nanoQC
    cleanqcfolder = os.path.join(odir,'nanoQC','clean')
    command = "{} {} {}\n".format(tools['nanoQC'],cleanqcfolder,clean_fq)
    ofp.write(command)

    # cleand data nanoplot
    command = "{program} --fastq {in_} -t {thread} --plots hex dot -o {out} -p {pre}\n".format(program=tools['NanoPlot'],in_=clean_fq,thread=thread,out=cleanqcfolder,pre=sample)
    ofp.write(command)

    # genome assemble
    assemblefolder = os.path.join(odir,'flye','assemble')
    command = "{program} --nano-raw {in_} -g {flye_g} -t {thread} -o {out}\n".format(program=tools['flye'],in_=clean_fq,flye_g=flye_g,thread=thread,out=assemblefolder)
    ofp.write(command)

    # assembled qualification - quast
    out = os.path.join(odir,'flye','quast')
    in_ = os.path.join(assemblefolder,'assembly.fasta')
    command = "{program} -t {thread} -o {out} {in}\n".format(program=tools['quast'],thread=thread,out=out,in_=in_)
    ofp.write(command)

    # assembled qualification - checkm
    out = os.path.join(odir,'flye','checkm')
    ocheckmf = os.path.join(odir,'flye','checkm','check.txt')
    command = "{program} lineage_wf {assemblef} {out} -x fasta -t {thread} --tab_table -f {ifile_}\n".format(program=tools['checkm'],assemblef=assemblefolder,out=out,thread=thread,ifile_=ocheckmf)
    ofp.write(command)
    ofp.write("\n")
  
  command = "source deactivate /work/workspace/zhurj/bin/miniconda3/envs/nanopore\n"  
  ofp.write(command)

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/6_nanoporepro/nanopore_only_flye_assemble.py -i /work/workspace/zhurj/project/13_nanopore/kejing20201218/input/test -o /work/workspace/zhurj/project/13_nanopore/kejing20201218 --obashf /work/workspace/zhurj/project/13_nanopore/kejing20201218/p/runbash.sh

find /work/rawdata/run/guangzhou/novogene/2020/10/20201022/run00072/rawdata | grep -P "fastq$" | grep -v "MN211200" | grep -v "MN231199" | grep -v "MN212517"| awk -F "/" '{print $0"\t"$(NF-1)}' > /work/workspace/zhurj/project/13_nanopore/kejing20201218/input/input37
cat input37 | tail -n 27 > input27

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/6_nanoporepro/nanopore_only_flye_assemble.py -i /work/workspace/zhurj/project/13_nanopore/kejing20201218/input/input37 -o /work/workspace/zhurj/project/13_nanopore/kejing20201218 --obashf /work/workspace/zhurj/project/13_nanopore/kejing20201218/p/sam37.sh -t 36
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/6_nanoporepro/nanopore_only_flye_assemble.py -i /work/workspace/zhurj/project/13_nanopore/kejing20201218/input/input27 -o /work/workspace/zhurj/project/13_nanopore/kejing20201218 --obashf /work/workspace/zhurj/project/13_nanopore/kejing20201218/p/sam27.sh -t 36
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/6_nanoporepro/nanopore_only_flye_assemble.py -i /work/workspace/zhurj/project/13_nanopore/kejing20201218/input/input_11_1 -o /work/workspace/zhurj/project/13_nanopore/kejing20201218 --obashf /work/workspace/zhurj/project/13_nanopore/kejing20201218/p/sam11_1.sh -t 36
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/6_nanoporepro/nanopore_only_flye_assemble.py -i /work/workspace/zhurj/project/13_nanopore/kejing20201218/input/input_11_2 -o /work/workspace/zhurj/project/13_nanopore/kejing20201218 --obashf /work/workspace/zhurj/project/13_nanopore/kejing20201218/p/sam11_2.sh -t 36
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/6_nanoporepro/nanopore_only_flye_assemble.py -i /work/workspace/zhurj/project/13_nanopore/kejing20201218/input/input_5 -o /work/workspace/zhurj/project/13_nanopore/kejing20201218 --obashf /work/workspace/zhurj/project/13_nanopore/kejing20201218/p/sam5_20201224.sh -t 36
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/6_nanoporepro/nanopore_only_flye_assemble.py -i /work/workspace/zhurj/project/13_nanopore/kejing20201218/input/input_8 -o /work/workspace/zhurj/project/13_nanopore/kejing20201218 --obashf /work/workspace/zhurj/project/13_nanopore/kejing20201218/p/sam8_20201224.sh -t 36

cd /work/workspace/zhurj/project/13_nanopore/kejing20201218/p
time srun -o sam37.out -e sam37.err -N 1 -c 20 -p slurm256 -w mnclient02 bash sam37.sh &
time srun -o sam27.out -e sam27.err -N 1 -c 20 -p slurm256 -w mnclient02 bash sam27.sh &
time srun -o sam11_1.out -e sam11_1.err -N 1 -c 20 -p slurm256 -w mnclient01 bash sam11_1.sh &
time srun -o sam11_2.out -e sam11_2.err -N 1 -c 20 -p slurm256 -w mnclient02 bash sam11_2.sh &
time srun -o sam5_20201224.out -e sam5_20201224.err -N 1 -c 20 -p slurm256 -w mnclient01 bash sam5_20201224.sh &
time srun -o sam8_20201224.out -e sam8_20201224.err -N 1 -c 20 -p slurm256 -w mnclient02 bash sam8_20201224.sh &
time srun -o MN233887_29385.out -e MN233887_29385.err -N 1 -c 20 -p slurm256 -w mnclient01 bash MN233887_29385.sh &
time srun -o MN231218.out -e MN231218.err -N 1 -c 20 -p slurm256 -w mnclient02 bash MN231218.sh &

MN230323 need input something when running flye of MN230323, so run this ample manually
/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/flye --nano-raw /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN230323/clean/nanoFilt.fastq.gz -g 6m -t 10 -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN230323/flye/assemble
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/quast -t 36 -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN230323/flye/quast /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN230323/flye/assemble/assembly.fasta
/work/workspace/zhurj/bin/miniconda3/bin/checkm lineage_wf /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN230323/flye/assemble /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN230323/flye/checkm -x fasta -t 36 --tab_table -f /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN230323/flye/checkm/check.txt


source activate /work/workspace/zhurj/bin/miniconda3/envs/nanopore
echo "MN29385"
mkdir -p /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN29385/nanoQC/raw
mkdir -p /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN29385/nanoQC/clean
mkdir -p /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN29385/clean
mkdir -p /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN29385/flye/assemble
mkdir -p /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN29385/flye/quast
mkdir -p /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN29385/flye/checkm
/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/nanoQC -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN29385/nanoQC/raw /work/rawdata/run/guangzhou/novogene/2020/10/20201022/run00072/rawdata/MN29385/MN29385.fastq
/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/NanoPlot --fastq /work/rawdata/run/guangzhou/novogene/2020/10/20201022/run00072/rawdata/MN29385/MN29385.fastq -t 36 --plots hex dot -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN29385/nanoQC/raw -p MN29385
/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/NanoFilt -q 8 -l 1000 --headcrop 50 --tailcrop 50 /work/rawdata/run/guangzhou/novogene/2020/10/20201022/run00072/rawdata/MN29385/MN29385.fastq | /work/workspace/zhurj/bin/pigz > /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN29385/clean/nanoFilt.fastq.gz
/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/nanoQC -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN29385/nanoQC/clean /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN29385/clean/nanoFilt.fastq.gz
/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/NanoPlot --fastq /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN29385/clean/nanoFilt.fastq.gz -t 36 --plots hex dot -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN29385/nanoQC/clean -p MN29385
/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/flye --nano-raw /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN29385/clean/nanoFilt.fastq.gz -g 6m -t 36 -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN29385/flye/assemble
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/quast -t 36 -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN29385/flye/quast /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN29385/flye/assemble/assembly.fasta
/work/workspace/zhurj/bin/miniconda3/bin/checkm lineage_wf /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN29385/flye/assemble /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN29385/flye/checkm -x fasta -t 36 --tab_table -f /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN29385/flye/checkm/check.txt

echo "MN233887"
/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/flye --nano-raw /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN233887/clean/nanoFilt.fastq.gz -g 6m -t 36 -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN233887/flye/assemble
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/quast -t 36 -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN233887/flye/quast /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN233887/flye/assemble/assembly.fasta
/work/workspace/zhurj/bin/miniconda3/bin/checkm lineage_wf /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN233887/flye/assemble /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN233887/flye/checkm -x fasta -t 36 --tab_table -f /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN233887/flye/checkm/check.txt
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/nanopore


source activate /work/workspace/zhurj/bin/miniconda3/envs/nanopore
echo "MN231218"
/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/flye --nano-raw /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231218/clean/nanoFilt.fastq.gz -g 6m -t 36 -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231218/flye/assemble
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/quast -t 36 -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231218/flye/quast /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231218/flye/assemble/assembly.fasta
/work/workspace/zhurj/bin/miniconda3/bin/checkm lineage_wf /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231218/flye/assemble /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231218/flye/checkm -x fasta -t 36 --tab_table -f /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN231218/flye/checkm/check.txt
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/nanopore


source activate /work/workspace/zhurj/bin/miniconda3/envs/nanopore
echo "MN241760"
/work/workspace/zhurj/bin/miniconda3/envs/nanopore/bin/flye --nano-raw /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN241760/clean/nanoFilt.fastq.gz -g 6m -t 36 -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN241760/flye/assemble
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/quast -t 36 -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN241760/flye/quast /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN241760/flye/assemble/assembly.fasta
/work/workspace/zhurj/bin/miniconda3/bin/checkm lineage_wf /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN241760/flye/assemble /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN241760/flye/checkm -x fasta -t 36 --tab_table -f /work/workspace/zhurj/project/13_nanopore/kejing20201218/MN241760/flye/checkm/check.txt
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/nanopore
cd /work/workspace/zhurj/project/13_nanopore/kejing20201218/p
srun -o MN241760.out -e MN241760.err -N 1 -c 20 -p slurm256 -w mnclient01 bash MN241760.sh &

nanoQC result -- merged

find /work/workspace/zhurj/project/13_nanopore/kejing20201218 | grep "NanoStats" | grep "clean" | awk -F "/" '{print $0"\t"$(NF-3)}' > /work/workspace/zhurj/project/13_nanopore/kejing20201218/input/clean_nanostats40
cat /work/workspace/zhurj/project/13_nanopore/kejing20201218/input/clean_nanostats40 | head -n 2 > /work/workspace/zhurj/project/13_nanopore/kejing20201218/input/clean_nanotest2
NanoStats.txt
clean

import os
import sys
import re
import pandas as pd

def read_content_extract(str_):
  temp = []
  temp = line.split(':')
  return temp[1].strip()

def qual_content_extract(str_):
  temp = []
  temp = re.findall(r' [(](.*)[)] ',str1)
  return temp[0]

con = ''

infile = '/work/workspace/zhurj/project/13_nanopore/kejing20201218/input/clean_nanotest2'
ofile = '/work/workspace/zhurj/project/13_nanopore/kejing20201218/tmp/output'

# get files names
data = pd.read_csv(infile,sep='\t',names=['value','id'],index_col=None)
name2files_dict = data.set_index('id')['value'].to_dict()

ofp = open(ofile,'w')
con = "Sample\tMean read length\tMean read quality\tMedian read length\tMedian read quality\tNumber of reads\tRead length N50\tTotal bases\t>Q5 percentage of read\t>Q7 percentage of read\tQ10 percentage of read\tQ12 percentage of read\tQ15 percentage of read\n"
ofp.write(con)

for name,fileone in name2files_dict.items():
  contemp = []
  contemp.append(name)
  mean_read_length = mean_read_quality = median_read_length = median_read_quality = number_of_reads = read_length_N50 = total_bases = q5 = q7 = q10 = q12 = q15 = ''
  with open(fileone,'r',encoding='utf-8') as fp:
    for line in fp.readlines():
      if line.startswith('Mean read length'):
        mean_read_length = read_content_extract(line)
        print(mean_read_length)
      if line.startswith('Mean read quality'):
        mean_read_quality = read_content_extract(line)
      if line.startswith('Median read length'):
        median_read_length = read_content_extract(line)
      if line.startswith('Median read quality'):
        median_read_quality = read_content_extract(line)
      if line.startswith('Number of reads'):
        number_of_reads = read_content_extract(line)
      if line.startswith('Read length N50'):
        read_length_N50 = read_content_extract(line)
      if line.startswith('Total bases'):
        total_bases = read_content_extract(line)
      if line.startswith('>Q5'):
        q5 = qual_content_extract(line)
      if line.startswith('>Q7'):
        q7 = qual_content_extract(line)
      if line.startswith('>Q10'):
        q10 = qual_content_extract(line)
      if line.startswith('>Q12'):
        q12 = qual_content_extract(line)
      if line.startswith('>Q15'):
        q15 = qual_content_extract(line)
  contemp.append(mean_read_length)
  contemp.append(mean_read_quality) 
  contemp.append(median_read_length)     
  contemp.append(median_read_quality) 
  contemp.append(number_of_reads) 
  contemp.append(read_length_N50) 
  contemp.append(total_bases) 
  contemp.append(q5) 
  contemp.append(q7) 
  contemp.append(q10) 
  contemp.append(q12) 
  contemp.append(q15) 
  con = "\t".join(contemp)
  con = con + "\n"
  ofp.write(con)

ofp.close()

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/6_nanoporepro/merge_nanoQC_output.py -i /work/workspace/zhurj/project/13_nanopore/kejing20201218/input/clean_nanostats40 -o /work/workspace/zhurj/project/13_nanopore/kejing20201218/output/nanoQC_merge


