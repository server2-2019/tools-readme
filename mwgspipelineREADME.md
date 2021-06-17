'''
# tools-readme
Record readme of often used software
#===============================================================================
#
#      PROJECT: Microbial Whole genome sequencing analysis 
#
#  DESCRIPTION: pipeline development and test
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: --
#        NOTES: ---
#       AUTHOR: Ruijuan Zhu (Zhu RJ), ruijuan6@gmail.com
# ORGANIZATION: MOONBIOTECH 
#      VERSION: 1.0
#      CREATED: 2019年 02月 21日 星期四 14:13:31 CST
#     REVISION: ---
#===============================================================================

pipeline 1: /mnt/d/home/yuanjie/script/isolate-strain2analysis/fastq-clean-bwa-flash.soapdenovo.pl
source deactivate MWGS --all
conda remove -n MWGS --all
conda create -n MWGS # /mnt/d/home/ruijuan/.conda/envs/MWGS
source activate MWGS
# circos Tutorial： http://circos.ca/documentation/tutorials/
# bioconda source: https://anaconda.org/bioconda/circos
conda install --name MWGS -c bioconda circos
circos -version
circos -module
# To see a list of paths that are searched, use the -debug_group flag.
circos -debug_group io
# can specify the conf file by "circos -conf myimage.conf"
# add samtools
conda install -c bioconda bwa
conda install -c bioconda samtools
conda uninstall samtools
conda install -c conda-forge r-gplots
conda install -c r r-rcolorbrewer
conda install -c bioconda seqkit
conda install -c bioconda fastqc
conda install -c bioconda multiqc
conda install -c bioconda gatk4 
conda install -c bioconda realphy
conda install -c bioconda bowtie2
conda install -c bioconda raxml
conda install -c bioconda phyml 
# pyani
# ref: https://github.com/widdowquinn/pyani && https://anaconda.org/bioconda/pyani
# install pip3
sudo yum -y install python36-pip # install pip3 on server
sudo pip3 install pyani
# pyani require the blast in bashrc



ANI: OrthoANI https://www.ezbiocloud.net/tools/orthoaniu

ANI命令行软件
ref: http://www.zhengyue90.com/?p=127

/mnt/d/home/ruijuan/workflows/software/pyani-master/bin/average_nucleotide_identity.py -i /work/workspace/ruijuan/project/2_mwgs/run00020-25/result/ksnpdata -o /work/workspace/ruijuan/project/2_mwgs/run00020-25/result/anitmp/tmp -m ANIb -g

ref: https://ngs-data-for-pathogen-analysis.readthedocs.io/zh_CN/latest/chapter_03/snp.html
BOWTIE2 /mnt/d/home/ruijuan/.conda/envs/MWGS/bin/bowtie2
BOWTIE2BUILDER /mnt/d/home/ruijuan/.conda/envs/MWGS/bin/bowtie2-build
/mnt/d/home/ruijuan/workflows/software/tree-puzzle-5.3.rc16-linux
./configure --prefix /mnt/d/home/ruijuan/.conda/envs/MWGS && make install
TREEPUZZLE /mnt/d/home/ruijuan/.conda/envs/MWGS/bin/puzzle
RAXML /mnt/d/home/ruijuan/.conda/envs/MWGS/bin/raxmlHPC
Rscript /mnt/d/home/ruijuan/.conda/envs/MWGS/bin/Rscript
MaxPars /mnt/d/home/ruijuan/workflows/software/phylip-3.697/exe/dnapars
PhyML /mnt/d/home/ruijuan/workflows/software/PhyML-3.1/PhyML-3.1_linux64

/mnt/d/home/ruijuan/.conda/envs/MWGS/bin/realphy
java -Xmx18g -jar /mnt/d/home/ruijuan/.conda/envs/MWGS/bin/RealPhy_v112.jar [Sequence folder] [Output folder] [Options] -ref 

download phylip-3.697.tar.gz from http://evolution.gs.washington.edu/phylip/getme-new1.html
tar zxvf phylip-3.697.tar.gz
make -f Makefile.unx install
sequencing data download /mnt/d/home/ruijuan/.conda/envs/MWGS/bin/bowtie2
① download data of 49 samples to "/mnt/d/home/x/data/W/MoonNGS/raw/20190225Novogene49samples"
file structure
cd /mnt/d/home/x/data/W/MoonNGS/raw/20190225Novogene49samples
find ./ -type f -name *.txt > file.txt
#create md5check.sh based on file.txt
s10_BDMS190003326-1a
  |_MD5_s10_BDMS190003326-1a.txt
  |_s10_BDMS190003326-1a_1.fq.gz
  |_s10_BDMS190003326-1a_2.fq.gz
# md5sum check
cd /mnt/d/home/x/data/W/MoonNGS/raw/20190225Novogene49samples
bash md5check.sh
RESULT: TRUE

② download data of 18 samples to "/mnt/d/home/x/data/W/MoonNGS/raw/20190225Novogene18samples"
cd /mnt/d/home/x/data/W/MoonNGS/raw/20190225Novogene18samples
find ./ -type f -name *.txt > file.txt
# create md5check.sh based on file 
1.rawdata
   |_1019251_BDMS190002985-1a
       |_MD5_1019251_BDMS190002985-1a.txt
       |_1019251_BDMS190002985-1a_1.fq.gz
       |_1019251_BDMS190002985-1a_2.fq.gz
# md5check
cd /mnt/d/home/x/data/W/MoonNGS/raw/20190225Novogene18samples 
bash md5check.sh
RESULT: TRUE



reference download
ascp -T -k 1 -i /home/p/ascp/connect/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/genomes/GENOME_REPORTS/prokaryotes.txt  /work/workspace/ruijuan/script/2_MWGS/ref
perl 1_refdownload.pl -i /work/workspace/ruijuan/script/2_MWGS/input/seqbacteriainfo.txt -f /work/workspace/ruijuan/script/2_MWGS/ref/prokaryotes.txt -o /work/workspace/ruijuan/script/2_MWGS/output

ascp -T -k 1 -i /home/p/ascp/connect/etc/asperaweb_id_dsa.openssh anonftp@ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/154/205/GCF_000154205.1_ASM15420v1/GCF_000154205.1_ASM15420v1_genomic.fna.gz  /work/workspace/ruijuan/project/5_tmp

format sequencing data (run00021)
creat sampleinfo file
account x
cd /mnt/d/home/x/data/W/MoonNGS/raw/20190225Novogene49samples
find ./ -type f | grep ".gz" | awk -F "/" '{print $NF}' | awk -F "_[12]{1}" ' {print $1}' | uniq > sample.txt
# sample.txt include the sample name
# test command
perl 3_rawdataorganize.pl -i /work/rawdata/test/rawdata/run/novogene/20190225/run00020 -o test -r input/run00021.txt
# enter root account
cd /work
perl /work/workspace/ruijuan/script/3_rawdataOrg/3_rawdataorganize.pl -i /work/rawdata/run/guangzhou/novogene/20190225/run00021 -o /work/rawdata -r /work/workspace/ruijuan/script/3_rawdataOrg/input/run00021.txt

format sequencing data (run00020)
account root
mkdir /work/rawdata/run/guangzhou/novogene/20190225/run00020
cp /mnt/d/home/x/data/W/MoonNGS/raw/20190225Novogene18samples/1.rawdata /work/rawdata/run/guangzhou/novogene/20190225/run00020 -rf
create sampleinfo file
find /mnt/d/home/x/data/W/MoonNGS/raw/20190225Novogene18samples -type f | grep ".gz" | awk -F "/" '{ print $NF }' | awk -F "_[12]{1}[_\.]{1}" ' {print $1} '| uniq
create /work/workspace/ruijuan/script/3_rawdataOrg/input/run00020.txt
perl /work/workspace/ruijuan/script/3_rawdataOrg/3_rawdataorganize.pl -i /work/rawdata/run/guangzhou/novogene/20190225/run00020 -o /work/rawdata -r /work/workspace/ruijuan/script/3_rawdataOrg/input/run00020.txt


referebce download
info in prokaryotes.txt
# species : 64247 (cat prokaryotes.txt | awk -F "\t" ' {print $1}  ' | sort | uniq  | wc -l )
  taxid: 64247 (cat prokaryotes.txt | awk -F "\t" ' {print $2}  ' | sort | uniq  | wc -l)
  ftp: 190080 (cat prokaryotes.txt | awk -F "\t" ' {print $21}  ' | sort | uniq  | wc -l)
  assembly_id: 190080
  status: Complete Genome, Complete, Chromosome, Scaffold, Contig (cat prokaryotes.txt | awk -F "\t" ' {print $16}  ' | sort | uniq  | wc -l)

  perl 1_refdownload_going.pl

  $finalsaminfo{$fsf_spenm} = "$fsf_spenm\tContig\t".$selectedsampleinfo{$fsf_spenm}{'Contig'}{'size'}."\t".$selectedsampleinfo{$fsf_spenm}{'Contig'}{'md'}."\t".$selectedsampleinfo{$fsf_spenm}{'Contig'}{'ftp'};

  perl 1_refdownload_going.pl -i input/speinfo.txt -r ref/prokaryotes.txt -o test # test
  account: root
  perl /work/workspace/ruijuan/script/3_rawdataOrg/1_refdownload_going.pl -i /work/workspace/ruijuan/script/3_rawdataOrg/input/speinfo.txt -r /work/workspace/ruijuan/script/3_rawdataOrg/ref/prokaryotes.txt -o /work/workspace/ruijuan/script/3_rawdataOrg/test # test
  perl /work/workspace/ruijuan/script/3_rawdataOrg/1_refdownload_going.pl -i /work/workspace/ruijuan/script/3_rawdataOrg/input/speinfo.txt -r /work/workspace/ruijuan/script/3_rawdataOrg/ref/prokaryotes.txt -o /mnt/d/work/database/ncbi/genome/all
  ## directory for the reference: /mnt/d/work/database/ncbi/genome/all


2019.03.01
cluParabacteroides_distasonis 
1. sequencing data organization
mkdir /work/workspace/ruijuan/project/2_mwgs/2019.Mar1.Parabacteroides_distasonis
find /work/workspace/yajun/project/2019.Feb27.run00021.genome_assembly/1.Cleandata/ -name "*.gz" |awk -F "[/_]" '{print $9"\t"$0}' >1.Cleandata.fastqlist

2. reference creation
mkdir /work/workspace/ruijuan/project/2_mwgs/2019.Mar1.Parabacteroides_distasonis/1.Reference
cd /work/workspace/ruijuan/project/2_mwgs/2019.Mar1.Parabacteroides_distasonis/
cp /mnt/d/work/database/ncbi/genome/all/GCA/GCA_000/GCA_000012/GCA_000012845/genomic.fna.gz 1.Reference
mv 1.Reference/genomic.fna.gz 1.Reference/GCA_000012845.fna.gz
gunzip 1.Reference/GCA_000012845.fna.gz
source activate MWGS
bwa index 1.Reference/GCA_000012845.fna
samtools faidx 1.Reference/GCA_000012845.fna

3.call SNPs
cd /work/workspace/ruijuan/project/2_mwgs/2019.Mar1.Parabacteroides_distasonis
perl /work/workspace/yuanjie/script/isolate-strain2analysis/fastqlist-clean-bwa2snp.pl -i 1.Cleandata.fastqlist -d ./ -p 2019.Mar1.Parabacteroides_distasonis -r 1.Reference/GCA_000012845.fna -c
sh work.sh 1>/dev/null 2>work.sh.err

4. 生成VCF.LIST 得到SNP table
find 6.Analysis/2019.Mar1.Parabacteroides_distasonis/ -name "*.vcf" |awk -F "[/.]" '{print $(NF-1)"\t"$0}' >6.Analysis.vcflist
perl /work/workspace/yuanjie/script/isolate-strain2analysis/vcflist-filter2snptable.pl --vcflist 6.Analysis.vcflist --refset PD --result 6.Analysis/2019.Mar1.Parabacteroides_distasonis/2019.Mar1.Parabacteroides_distasonis.PD

5. SNP 过滤
perl /work/workspace/yuanjie/script/isolate-strain2analysis/filter-vcf-snp.pl 6.Analysis/2019.Mar1.Parabacteroides_distasonis/2019.Mar1.Parabacteroides_distasonis.PD /work/workspace/ruijuan/project/2_mwgs/2019.Mar1.Parabacteroides_distasonis/1.Reference/GCA_000012845.fna 2019.Mar1.Parabacteroides_distasonis.PD.clean

6. 绘图
perl /work/workspace/yuanjie/script/isolate-strain2analysis/snp_table_matrix.pl 2019.Mar1.Parabacteroides_distasonis.PD.clean 2019.Mar1.Parabacteroides_distasonis.PD.clean.matrix

perl /work/workspace/yuanjie/script/usefulKit/general-heatmap.pl 2019.Mar1.Parabacteroides_distasonis.PD.clean.matrix 2019.Mar1.Parabacteroides_distasonis.PD.clean.matrix

R -f 2019.Mar1.Parabacteroides_distasonis.PD.clean.matrix.R >/dev/null

library(RColorBrewer)
library(gplots)

MEGA X
cd  /work/workspace/ruijuan/project/2_mwgs/2019.Mar1.Parabacteroides_distasonis
perl /work/workspace/yuanjie/script/isolate-strain2analysis/snp_table_fas.pl 2019.Mar1.Parabacteroides_distasonis.PD.clean 2019.Mar1.Parabacteroides_distasonis.PD.clean.fas

#-------------------------------------------------------------------------------------------------------------------
Akkermansia muciniphila ATCC BAA-835
1. reference creation
reference creation
mkdir -p /work/workspace/ruijuan/project/2_mwgs/2019.Mar1.Akkermansia_muciniphila/1.Reference
cd /work/workspace/ruijuan/project/2_mwgs/2019.Mar1.Akkermansia_muciniphila/
cp /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000020/GCA_000020225.1/genomic.fna.gz 1.Reference
mv 1.Reference/genomic.fna.gz 1.Reference/GCA_000020225.fna.gz
gunzip 1.Reference/GCA_000020225.fna.gz
source activate MWGS
bwa index 1.Reference/GCA_000020225.fna
samtools faidx 1.Reference/GCA_000020225.fna

2. sequencing data organization
find /work/workspace/yajun/project/2019.Feb28.run00020.Akkermansia_muciniphila/1.Cleandata -name "*.gz" |awk -F "[/_]" '{print $9"\t"$0}' >1.Cleandata.fastqlist

3.call SNPs
cd /work/workspace/ruijuan/project/2_mwgs/2019.Mar1.Akkermansia_muciniphila
perl /work/workspace/yuanjie/script/isolate-strain2analysis/fastqlist-clean-bwa2snp.pl -i 1.Cleandata.fastqlist -d ./ -p 2019.Mar1.Akkermansia_muciniphila -r 1.Reference/GCA_000020225.fna -c
bash work.sh 1>/dev/null 2>work.sh.err

4. 生成VCF.LIST 得到SNP table
find 6.Analysis/2019.Mar1.Akkermansia_muciniphila/ -name "*.vcf" |awk -F "[/.]" '{print $(NF-1)"\t"$0}' >6.Analysis.vcflist
perl /work/workspace/yuanjie/script/isolate-strain2analysis/vcflist-filter2snptable.pl --vcflist 6.Analysis.vcflist --refset AKK --result 6.Analysis/2019.Mar1.Akkermansia_muciniphila/2019.Mar1.Akkermansia_muciniphila.AKK

5. SNP 过滤
perl /work/workspace/yuanjie/script/isolate-strain2analysis/filter-vcf-snp.pl 6.Analysis/2019.Mar1.Akkermansia_muciniphila/2019.Mar1.Akkermansia_muciniphila.AKK /work/workspace/ruijuan/project/2_mwgs/2019.Mar1.Akkermansia_muciniphila/1.Reference/GCA_000020225.fna 2019.Mar1.Akkermansia_muciniphila.AKK.clean

6. 绘图
perl /work/workspace/yuanjie/script/isolate-strain2analysis/snp_table_matrix.pl 2019.Mar1.Akkermansia_muciniphila.AKK.clean 2019.Mar1.Akkermansia_muciniphila.AKK.clean.matrix

perl /work/workspace/yuanjie/script/usefulKit/general-heatmap.pl 2019.Mar1.Akkermansia_muciniphila.AKK.clean.matrix 2019.Mar1.Akkermansia_muciniphila.AKK.clean.matrix

R -f 2019.Mar1.Akkermansia_muciniphila.AKK.clean.matrix.R >/dev/null

AKKsnpCluster.sh
perl /work/workspace/yuanjie/script/isolate-strain2analysis/fastqlist-clean-bwa2snp.pl -i 1.Cleandata.fastqlist -d ./ -p 2019.Mar1.Akkermansia_muciniphila -r 1.Reference/GCA_000020225.fna -c
bash work.sh 1>/dev/null 2>work.sh.err
find 6.Analysis/2019.Mar1.Akkermansia_muciniphila/ -name "*.vcf" |awk -F "[/.]" '{print $(NF-1)"\t"$0}' >6.Analysis.vcflist
perl /work/workspace/yuanjie/script/isolate-strain2analysis/vcflist-filter2snptable.pl --vcflist 6.Analysis.vcflist --refset AKK --result 6.Analysis/2019.Mar1.Akkermansia_muciniphila/2019.Mar1.Akkermansia_muciniphila.AKK
perl /work/workspace/yuanjie/script/isolate-strain2analysis/filter-vcf-snp.pl 6.Analysis/2019.Mar1.Akkermansia_muciniphila/2019.Mar1.Akkermansia_muciniphila.AKK /work/workspace/ruijuan/project/2_mwgs/2019.Mar1.Akkermansia_muciniphila/1.Reference/GCA_000020225.fna 2019.Mar1.Akkermansia_muciniphila.AKK.clean
perl /work/workspace/yuanjie/script/isolate-strain2analysis/snp_table_matrix.pl 2019.Mar1.Akkermansia_muciniphila.AKK.clean 2019.Mar1.Akkermansia_muciniphila.AKK.clean.matrix

perl /work/workspace/yuanjie/script/usefulKit/general-heatmap.pl 2019.Mar1.Akkermansia_muciniphila.AKK.clean.matrix 2019.Mar1.Akkermansia_muciniphila.AKK.clean.matrix

R -f 2019.Mar1.Akkermansia_muciniphila.AKK.clean.matrix.R >/dev/null

mega x TREE
cd /work/workspace/ruijuan/project/2_mwgs/2019.Mar1.Akkermansia_muciniphila
perl /work/workspace/yuanjie/script/isolate-strain2analysis/snp_table_fas.pl 2019.Mar1.Akkermansia_muciniphila.AKK.clean 2019.Mar1.Akkermansia_muciniphila.AKK.clean.fas


2019.03.06
无参基因组组装

 python GenomeAssemblingJupyterPipeline_rj190306.py /work/rawdata/run/guangzhou/2019/02/20190225/run00021/fastq/MNH00499.1.fq.gz /work/rawdata/run/guangzhou/2019/02/20190225/run00021/fastq/MNH00499.2.fq.gz  /mnt/d/linux/W/MoonNGS/result/MWGS/run00021/MNH00499 MNH00499 

 MNH04496
mkdir -p /mnt/d/linux/W/MoonNGS/result/MWGS/run00021/MNH04496 
nohup python GenomeAssemblingJupyterPipeline_rj190306.py /work/rawdata/run/guangzhou/2019/02/20190225/run00021/fastq/MNH04496.1.fq.gz /work/rawdata/run/guangzhou/2019/02/20190225/run00021/fastq/MNH04496.2.fq.gz  /mnt/d/linux/W/MoonNGS/result/MWGS/run00021/MNH04496 MNH04496  > /mnt/d/linux/W/MoonNGS/result/MWGS/run00021/MNH04496/run.log 2>&1 &



root
mkdir -p /work/rawdata/test/run/novogene/20190225/run00020/fastq
cp /mnt/d/home/x/data/W/MoonNGS/raw/20190225Novogene49samples/s1_BDMS190003317-1a/*.gz /work/rawdata/test/run/novogene/20190225/run00020/fastq

mkdir -p /work/rawdata/test/rawdata/fastq/current/MNH/MNH004/MNH00499/E1/L1/S1
ln -s  /work/rawdata/test/rawdata/run/novogene/20190225/run00020/fastq/s1_BDMS190003317-1a_* /work/rawdata/test/rawdata/fastq/current/MNH/MNH004/MNH00499/E1/L1/S1 && cd /work/rawdata/test/rawdata/fastq/current/MNH/MNH004/MNH00499/E1/L1/S1 && rename s1_BDMS190003317-1a_ MNH00499. *.gz

mkdir -p /work/database/test/all/
mkdir -p /work/rawdata/run/novogene/20190225/run00021
cp /mnt/d/home/x/data/W/MoonNGS/raw/20190225Novogene49samples/* /work/rawdata/run/novogene/20190225/run00021 -rf

find /work/rawdata/run/novogene/20190225/run00021 -type d | awk -F "/" '{print $NF}'

perl 1_refdownload_going.pl -i input/speinfo.txt -r ref/prokaryotes.txt -o test

=------------------------------------------------------------------------------------------------------------------
序列拼装软件的选择
if (micro diversity is not a major issue&& the primary research goal is to bin && reconstruct representative bacterial genomes from a given environment){

  metaSPAdes should clearly be the assembler of choice. # This assembler yields the best contig size statistics  while capturing a high degree of community diversity, even at high complexity and low read coverage;

}elsif(mico diversity is however an issue || the degree of
  captured diversity is far more important than contig
  lengths){

  then IDBA-UD or Megahit should be preferred. #  The sensitivity of these assemblers, both for diversity as  well as micro diversity, makes them optimal choices when trying to discover novel species in complex habitats. Whenever computational resources become limiting, 
  Megahit becomes the most attractive option, due to its good compromise between contig size statistics, captured diversity and required memory.
}

 However, the bias of Megahit towards relatively low coverage genomes may provide a disadvantage for very large datasets, leading to a suboptimal assembly of high abundant community member genomes.
 In such cases, Megahit may provide better results when assembling subsets of the sequencing data in a “divide and conquer” approach.


 20190306
 download 55 sequencing data to /work/rawdata/run/guangzhou/novogene/20190306/run00022
 find ./ -type f -name *.txt | awk -F "[/]" '{print "cd "$2" && md5sum -c "$3" && cd .."}' > md5checklist.sh
 bash md5checklist.sh
 MNH12268-1
MNH02743-1
MNH04191-1
MNH09210-1
MNH07080-1
MNH02752-1
MNH02683-1
MNH08990-1
MNH20130-1
MNH09022-1
MNH03588-1
MNH20105-1
MNH09691-1
MNH13391-1
MNH01988-1
MNH03582-1
MNH07037-1
MNH09206-1
MNH12866-1
MNH09207-1
MNH19666-1
MNH06245-1
MNH20128-1
MNH03633-1
MNH19831-1
MNH04174-1
MNH03859-1
MNH12274-1
MNH14763-1
MNH03718-1
MNH19820-1
MNH12921-1
MNH07060-1
MNH20160-1
MNH13780-1
MNH12934-1
MNH13313-1
MNH13979-1
MNH12881-1
MNH12853-1
MNH13392-1
MNH03654-1
MNH12868-1
MNH08149-1
MNH04531-1
MNH20069-1
MNH13400-1
MNH02271-1
MNH07066-1
MNH03922-1
MNH12498-1
MNH12477-1
MNH06395-1
MNH13272-1
MNH08509-1

root account
perl /work/workspace/ruijuan/script/3_rawdataOrg/3_rawdataorganize.pl -i /work/rawdata/run/guangzhou/novogene/20190306/run00022/1.rawdata -o /work/rawdata -r /work/workspace/ruijuan/script/3_rawdataOrg/input/run00022.txt


reference download
cd /work/workspace/ruijuan/script/3_rawdataOrg
perl 1_refdownload_done.pl -i /work/workspace/ruijuan/script/3_rawdataOrg/input/speinfo_run00022.txt \
            -r /work/workspace/ruijuan/script/3_rawdataOrg/ref/prokaryotes.txt -o /mnt/d/work/database/ncbi/genome/all

/home/p/prokka/prokka-master/bin/prokka --kingdom Bacteria --cpus 20  MNH04496.genome.fa --outdir .  --prefix MNH04496 --force --metagenome --locustag MNH04496
tbl2asn -V b -a r10k -l paired-ends -M n -N 1 -y 'Annotated using prokka 1.14-dev from https://github.com/tseemann/prokka' -Z ./MNH04496.err -i ./MNH04496.fsa 2> /dev/null


# pipeline
1. quality control
   1) READ, BASE, Q20, Q30,
      seqkit stats -a /work/rawdata/run/guangzhou/2019/02/20190225/run00021/fastq/MNH06648.1.fq.gz -j 18
      ## cd  
      ## perl qcstat.pl -i qcstat/input/infile.txt -o qcstat/output --basecutoff 0.8 --cutoffq20 90 --cutoffq30 90
      ## vi qcstat/output/qcstat.txt
   2) FASTQC
      fastqc -o result/test -t 18 --quiet /work/rawdata/run/guangzhou/2019/02/20190225/run00021/fastq/MNH06648.1.fq.gz /work/rawdata/run/guangzhou/2019/02/20190225/run00021/fastq/MNH06648.2.fq.gz
   3) multiqc
      multiqc -o result/multiqc -n run00020 result/test

      ## perl /work/workspace/ruijuan/script/2_MWGS/mfastqc.pl -i /work/workspace/ruijuan/script/2_MWGS/qcstat/input/infile.txt -o /work/workspace/ruijuan/script/2_MWGS/mfqcresult -n test

 2. data filter (adapter, low quality reads, short size reads) 
    fastp -i /work/rawdata/run/guangzhou/2019/02/20190225/run00021/fastq/MNH06648.1.fq.gz -o result/clean/MNH06648.clean.1.fq.gz -I /work/rawdata/run/guangzhou/2019/02/20190225/run00021/fastq/MNH06648.2.fq.gz -O result/clean/MNH06648.clean.2.fq.gz --poly_g_min_len 10 --poly_x_min_len 10 -q 15 -u 40 -n 5 -l 50 -w 16
    ## perl datafilter.pl -i /work/workspace/ruijuan/script/2_MWGS/qcstat/input/infile.txt -o /work/workspace/ruijuan/script/2_MWGS/datafilter


3. alignment - bwa (1254.123 sec/MNH06648)
   VARIANTS CALL WITH GATK: (https://approachedinthelimit.wordpress.com/2015/10/09/variant-calling-with-gatk/)
   ref: 
   cd /work/workspace/ruijuan/script/2_MWGS/ref
   cp /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000012/GCA_000012845.1/genomic.fna.gz GCA_000012845.fa.gz && gunzip GCA_000012845.fa.gz
   bwa index GCA_000012845.fa
   samtools faidx GCA_000012845.fa
   gatk CreateSequenceDictionary -R /work/workspace/ruijuan/script/2_MWGS/ref/GCA_000012845.fa
   
   bwa mem -t 18  -R '@RG\tID:MNH06648\tPL:ILLUMINA\tLB:MNH06648\tSM:PA01' ref/GCA_000012845.fa /work/workspace/ruijuan/script/2_MWGS/result/clean/MNH06648.clean.1.fq.gz /work/workspace/ruijuan/script/2_MWGS/result/clean/MNH06648.clean.2.fq.gz -o test/result/variant/MNH06648.sam

   gatk SortSam -I result/variant/MNH06648.sam -O result/variant/MNH06648.sorted.bam --SORT_ORDER coordinate

   gatk MarkDuplicates -I result/variant/MNH06648.sorted.bam -O result/variant/MNH06648.dedup.bam --REMOVE_DUPLICATES true --METRICS_FILE metrics.txt 

   gatk BuildBamIndex -I result/variant/MNH06648.dedup.bam

   gatk  HaplotypeCaller -R ref/GCA_000012845.fa -I result/variant/MNH06648.dedup.bam -O result/variant/MNH06648.raw.vcf -stand-call-conf 30.0 -mbq 10 --QUIET true -ploidy 1 --dont-use-soft-clipped-bases true -ERC GVCF

   gatk   HaplotypeCallerSpark -R ref/GCA_000012845.fa -I result/variant/MNH06648.dedup.bam -O result/variant/MNH06648.raw.vcf --spark-master local[20] -stand-call-conf 30.0 -mbq 10 --QUIET true -ploidy 1 --dont-use-soft-clipped-bases true -ERC GVCF
  
   OPTIONS:
   gatk FilterVcf --MIN_DP 10 -I result/variant/MNH06648.raw.vcf -O result/variant/MNH06648.filter.vcf
   gatk SelectVariants -R ref/GCA_000012845.fa -V result/variant/MNH06648.raw.vcf -select-type SNP -O  result/variant/MNH06648.SNP.vcf

   gatk SelectVariants -R ref/GCA_000012845.fa -V result/variant/MNH06648.raw.vcf -select-type INDEL -O  result/variant/MNH06648.INDEL.vcf


   realphy (ref: https://ngs-data-for-pathogen-analysis.readthedocs.io/zh_CN/latest/chapter_03/snp.html)
   java -Xmx18g -jar /mnt/d/home/ruijuan/.conda/envs/MWGS/bin/RealPhy_v112.jar /work/workspace/ruijuan/script/2_MWGS/data /work/workspace/ruijuan/script/2_MWGS/output -ref GCA_000012845 -treeBuilder 4 -d -config /work/workspace/ruijuan/script/2_MWGS/input/config.txt

   java -Xmx18g -jar /mnt/d/home/ruijuan/.conda/envs/MWGS/bin/RealPhy_v112.jar /work/workspace/ruijuan/script/2_MWGS/data /work/workspace/ruijuan/script/2_MWGS/realphytest -ref GCA_000012845 -treeBuilder 4 -d -config /work/workspace/ruijuan/script/2_MWGS/input/config.txt

   
   RAxML 来构建基于 ML 的进化树
   cd /work/workspace/ruijuan/script/2_MWGS/realphytest/GCA_000012845/PolySeqOut_NoGenes
   /mnt/d/home/ruijuan/.conda/envs/MWGS/bin/raxmlHPC -f a -x 12345 -p 12345 -# 100 -m GTRCAT -s polymorphisms_move.phy -n raxml -T 20

   java -Xmx20g  -jar /mnt/d/home/ruijuan/workflows/software/FigTree_v1.4.4/lib/figtree.jar -graphic PDF RAxML_bestTree.raxml RAxML_bestTree.pdf
   

for tmpfile in `ls | grep gz` 
do 
file=${tmpfile%@} 
newfile=`echo $file | sed 's/\.R/\_R/g'`
mv $file $newfile
done 

for tmpfile in `ls | grep gz` 
do 
file=${tmpfile%@} 
newfile=`echo $file | sed 's/fq/fastq/g'`
mv $file $newfile
done 


kSNP3 test

/mnt/d/home/ruijuan/workflows/software/kSNP3.1_Linux_package/kSNP3/kSNP3 -k 21 -in /work/workspace/ruijuan/script/2_MWGS/kSNPtest/ksnpinput/input.txt -outdir /work/workspace/ruijuan/script/2_MWGS/kSNPtest/output -core -CPU 10

/mnt/d/home/ruijuan/workflows/software/kSNP3.1_Linux_package/kSNP3/MakeFasta /work/workspace/ruijuan/script/2_MWGS/kSNPtest/ksnpinput/ksnpinput.txt /work/workspace/ruijuan/script/2_MWGS/kSNPtest/ksnpinput/in_kchooser.fa

/mnt/d/home/ruijuan/workflows/software/kSNP3.1_Source/Kchooser /work/workspace/ruijuan/script/2_MWGS/kSNPtest/ksnpinput/in_kchooser.fa
q
/mnt/d/home/ruijuan/workflows/software/kSNP3.1_Linux_package/kSNP3/jellyfish

/mnt/d/home/ruijuan/workflows/software/kSNP3.1_Linux_package/kSNP3/kSNP3 -k 21 -in /work/workspace/ruijuan/script/2_MWGS/kSNPtest/ksnpinput/file.txt -outdir /work/workspace/ruijuan/script/2_MWGS/kSNPtest/output1 -core -CPU 10

java -Xmx10g  -jar /mnt/d/home/ruijuan/workflows/software/FigTree_v1.4.4/lib/figtree.jar -graphic PDF /work/workspace/ruijuan/script/2_MWGS/kSNPtest/output1/tree_AlleleCounts.parsimony.tre /work/workspace/ruijuan/script/2_MWGS/kSNPtest/output1/tree_AlleleCounts.parsimony.pdf

java -Xmx10g  -jar /mnt/d/home/ruijuan/workflows/software/FigTree_v1.4.4/lib/figtree.jar -graphic PDF /work/workspace/ruijuan/script/2_MWGS/kSNPtest/output1/tree_tipAlleleCounts.parsimony.tre /work/workspace/ruijuan/script/2_MWGS/kSNPtest/output1/tree_tipAlleleCounts.parsimony.pdf

java -Xmx10g  -jar /mnt/d/home/ruijuan/workflows/software/FigTree_v1.4.4/lib/figtree.jar -graphic PDF /work/workspace/ruijuan/script/2_MWGS/kSNPtest/output1/tree.parsimony.tre /work/workspace/ruijuan/script/2_MWGS/kSNPtest/output1/tree.parsimony.pdf

java -Xmx10g  -jar /mnt/d/home/ruijuan/workflows/software/FigTree_v1.4.4/lib/figtree.jar -graphic PDF /work/workspace/ruijuan/script/2_MWGS/kSNPtest/output1/tree.core.tre /work/workspace/ruijuan/script/2_MWGS/kSNPtest/output1/tree.core.pdf

java -Xmx10g  -jar /mnt/d/home/ruijuan/workflows/software/FigTree_v1.4.4/lib/figtree.jar -graphic PDF /work/workspace/ruijuan/script/2_MWGS/kSNPtest/output1/tree_AlleleCounts.core.tre /work/workspace/ruijuan/script/2_MWGS/kSNPtest/output1/tree_AlleleCounts.core.pdf

java -Xmx10g  -jar /mnt/d/home/ruijuan/workflows/software/FigTree_v1.4.4/lib/figtree.jar -graphic PDF /work/workspace/ruijuan/script/2_MWGS/kSNPtest/output1/tree_tipAlleleCounts.core.tre /work/workspace/ruijuan/script/2_MWGS/kSNPtest/output1/tree_tipAlleleCounts.core.pdf

 perl genoPhyana.pl -i /work/workspace/ruijuan/script/2_MWGS/kSNPtest/groupTest/input/group.txt -o /work/workspace/ruijuan/project/2_mwgs/run00020-25/result -d /work/workspace/ruijuan/project/2_mwgs/run00020-25/result/finalAssemble 

  perl genoPhyana.pl -i /work/workspace/ruijuan/script/2_MWGS/kSNPtest/groupTest/input/group.txt -o /work/workspace/ruijuan/project/2_mwgs/run00020-25/result -d /work/workspace/ruijuan/project/2_mwgs/run00020-25/result/finalAssemble --ani

 java -Xmx10g  -jar /mnt/d/home/ruijuan/workflows/software/FigTree_v1.4.4/lib/figtree.jar -graphic PDF /work/workspace/ruijuan/project/2_mwgs/run00020-25/result/ksnptmp/Akkermansia_muciniphila_ATCC_BAA-835/tree.core.tre /work/workspace/ruijuan/project/2_mwgs/run00020-25/result/ksnptmp/Akkermansia_muciniphila_ATCC_BAA-835/tree.core.pdf


# test ref /mnt/d/work/workspace/yajun/project/2019.Feb28.run00020.Akkermansia_muciniphila/2.Assembly/

Figtree

java -Xmx10g  -jar /mnt/d/home/ruijuan/workflows/software/FigTree_v1.4.4/lib/figtree.jar -graphic PDF /work/workspace/ruijuan/script/2_MWGS/kSNPtest/output/tree_tipAlleleCounts.core.tre /work/workspace/ruijuan/script/2_MWGS/kSNPtest/output/tree_tipAlleleCounts.core.pdf
java -Xmx10g  -jar /mnt/d/home/ruijuan/workflows/software/FigTree_v1.4.4/lib/figtree.jar -graphic PDF /work/workspace/ruijuan/script/2_MWGS/kSNPtest/output1/tree_tipAlleleCounts.parsimony.tre /work/workspace/ruijuan/script/2_MWGS/kSNPtest/output1/tree_tipAlleleCounts.parsimony.pdf

java -Xmx10g  -jar /mnt/d/home/ruijuan/workflows/software/FigTree_v1.4.4/lib/figtree.jar -graphic PDF /work/workspace/ruijuan/script/2_MWGS/realphytest/GCA_000012845/PolySeqOut_NoGenes/polymorphisms_move.phy /work/workspace/ruijuan/script/2_MWGS/realphytest/GCA_000012845/PolySeqOut_NoGenes/polymorphisms_move.pdf

java -Xmx10g  -jar /mnt/d/home/ruijuan/workflows/software/FigTree_v1.4.4/lib/figtree.jar -graphic PDF /work/workspace/ruijuan/project/2_mwgs_kSNP3/2019.Mar18.Akkermansia_muciniphila/output/tree_tipAlleleCounts.core.tre /work/workspace/ruijuan/project/2_mwgs_kSNP3/2019.Mar18.Akkermansia_muciniphila/output/tree_tipAlleleCounts.core.pdf

AKK Amuc_1100 mutation check (ref: https://www.ncbi.nlm.nih.gov/protein/187425647)
cd /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000020/GCA_000020225.1/
zcat cds_from_genomic.fna.gz| grep Amuc_1100
result: >lcl|CP001071.1_cds_ACD04926.1_1079 [locus_tag=Amuc_1100] [protein=hypothetical protein] [protein_id=ACD04926.1] [location=complement(1314411..1315364)] [gbkey=CDS]
cd /work/workspace/ruijuan/project/2_mwgs/2019.Mar1.Akkermansia_muciniphila/

perl /work/workspace/yuanjie/script/alignKit/quick-blast.pl -q 2.Assembly/MNH17503.fas -s 5.Reference/GCA_000020225.Amuc_1100.fas -o 6.Analysis/2019.Mar18.Amuc_1100/MNH17503
perl /work/workspace/yuanjie/script/alignKit/quick-blast.pl -q 2.Assembly/MNH17505.fas -s 5.Reference/GCA_000020225.Amuc_1100.fas -o 6.Analysis/2019.Mar18.Amuc_1100/MNH17505
perl /work/workspace/yuanjie/script/alignKit/quick-blast.pl -q 2.Assembly/MNH19250.fas -s 5.Reference/GCA_000020225.Amuc_1100.fas -o 6.Analysis/2019.Mar18.Amuc_1100/MNH19250
perl /work/workspace/yuanjie/script/alignKit/quick-blast.pl -q 2.Assembly/MNH19631.fas -s 5.Reference/GCA_000020225.Amuc_1100.fas -o 6.Analysis/2019.Mar18.Amuc_1100/MNH19631
perl /work/workspace/yuanjie/script/alignKit/quick-blast.pl -q 2.Assembly/MNH19632.fas -s 5.Reference/GCA_000020225.Amuc_1100.fas -o 6.Analysis/2019.Mar18.Amuc_1100/MNH19632
perl /work/workspace/yuanjie/script/alignKit/quick-blast.pl -q 2.Assembly/MNH19714.fas -s 5.Reference/GCA_000020225.Amuc_1100.fas -o 6.Analysis/2019.Mar18.Amuc_1100/MNH19714

 AKK kSNP3
for tmpfile in `ls | grep fas$` 
do 
file=${tmpfile%@} 
newfile=`echo $file | sed 's/fas$/fasta/g'`
mv $file $newfile
done 

/mnt/d/home/ruijuan/workflows/software/kSNP3.1_Linux_package/kSNP3/kSNP3 -k 19 -in /work/workspace/ruijuan/project/2_mwgs_kSNP3/2019.Mar18.Akkermansia_muciniphila/input/input.txt -outdir /work/workspace/ruijuan/project/2_mwgs_kSNP3/2019.Mar18.Akkermansia_muciniphila/output -core -NJ -CPU 20 -ML

java -Xmx10g  -jar /mnt/d/home/ruijuan/workflows/software/FigTree_v1.4.4/lib/figtree.jar -graphic PDF /work/workspace/ruijuan/project/2_mwgs_kSNP3/2019.Mar18.Akkermansia_muciniphila/output/tree_tipAlleleCounts.core.tre /work/workspace/ruijuan/project/2_mwgs_kSNP3/2019.Mar18.Akkermansia_muciniphila/output/tree_tipAlleleCounts.core.pdf

java -Xmx10g  -jar /mnt/d/home/ruijuan/workflows/software/FigTree_v1.4.4/lib/figtree.jar -graphic PDF /work/workspace/ruijuan/project/2_mwgs_kSNP3/2019.Mar18.Akkermansia_muciniphila/output/tree_tipAlleleCounts.ML.tre /work/workspace/ruijuan/project/2_mwgs_kSNP3/2019.Mar18.Akkermansia_muciniphila/output/tree_tipAlleleCounts.ML.pdf

java -Xmx10g  -jar /mnt/d/home/ruijuan/workflows/software/FigTree_v1.4.4/lib/figtree.jar -graphic PDF /work/workspace/ruijuan/project/2_mwgs_kSNP3/2019.Mar18.Akkermansia_muciniphila/output/tree_tipAlleleCounts.NJ.tre /work/workspace/ruijuan/project/2_mwgs_kSNP3/2019.Mar18.Akkermansia_muciniphila/output/tree_tipAlleleCounts.NJ.pdf

java -Xmx10g  -jar /mnt/d/home/ruijuan/workflows/software/FigTree_v1.4.4/lib/figtree.jar -graphic PDF /work/workspace/ruijuan/project/2_mwgs_kSNP3/2019.Mar18.Akkermansia_muciniphila/output/tree_tipAlleleCounts.parsimony.tre /work/workspace/ruijuan/project/2_mwgs_kSNP3/2019.Mar18.Akkermansia_muciniphila/output/tree_tipAlleleCounts.parsimony.pdf

# docker
su root

docker run hello-world
docker image ls
docker container ls --all
(ref:https://blog.csdn.net/u010164190/article/details/80065555)
docker pull ubuntu
docker run -ti ubuntu bash # run ubuntu command
ctrl+p, then ctrl+q # exit ubuntu 
docker run -it -v /work/workspace/ruijuan/script/5_docker/ubuntu:/home ubuntu bash
docker ps -a
docker attach containerID
cat /etc/issue # check ubuntu version
cp /etc/apt/sources.list /etc/apt/sources.list.backup
apt-get update
apt-get upgrade

apt-get instal openssh-server openssh-client # install ssh

# assemble 
/home/p/anaconda/anaconda3_5.2.0/bin/SOAPdenovo-127mer all -s 
/work/workspace/ruijuan/script/2_MWGS/datafilter/clean/MNH00499.clean.1.fq.gz
/work/workspace/ruijuan/script/2_MWGS/datafilter/clean/MNH00499.clean.2.fq.gz

md5check
perl md5checkpro.pl -i /work/rawdata/run/guangzhou/novogene/20190321/run00023 -o ./

GCE: predict genome size
cd /work/workspace/ruijuan/script/2_MWGS/GCE
/mnt/d/home/ruijuan/workflows/software/GCE/gce-1.0.0/kmerfreq/kmer_freq_hash/kmer_freq_hash -k 21 -l read.list -t 20 -i 50000000 -o 0 -p test &> kmer_freq.log

/mnt/d/home/ruijuan/workflows/software/GCE/gce-1.0.0/gce \
  -f test.freq.stat -m 1 -D 8 -b 1 > test.table 2> test.log

  # K > 19 时，设置 -b 1；

run00023
sudo perl /work/workspace/ruijuan/script/3_rawdataOrg/3_rawdataorganize.pl -i /work/rawdata/run/guangzhou/novogene/20190321/run00023 -o /work/rawdata -r /work/workspace/ruijuan/script/3_rawdataOrg//input/run00023.txt
sudo perl /work/workspace/ruijuan/script/3_rawdataOrg/3_rawdataorganize.pl -i /work/rawdata/run/guangzhou/novogene/20190322/run00024 -o /work/rawdata -r /work/workspace/ruijuan/script/3_rawdataOrg//input/run00024.txt
sudo perl /work/workspace/ruijuan/script/3_rawdataOrg/3_rawdataorganize.pl -i /work/rawdata/run/guangzhou/novogene/20190322/run00025 -o /work/rawdata -r /work/workspace/ruijuan/script/3_rawdataOrg//input/run00025.txt
#samples	run		dir
4samples	run00023	/work/rawdata/run/guangzhou/2019/03/20190321/run00023/fastq
93samples	run00024	/work/rawdata/run/guangzhou/2019/03/20190322/run00024/fastq
32samples	run00025	/work/rawdata/run/guangzhou/2019/03/20190322/run00025/fastq

run00023
# input -o dir
mkdir -p /work/workspace/ruijuan/project/2_mwgs/run00023/md5checkresult
mkdir -p /work/workspace/ruijuan/project/2_mwgs/run00023/qcstat
mkdir -p /work/workspace/ruijuan/project/2_mwgs/run00023/fastqc
mkdir -p /work/workspace/ruijuan/project/2_mwgs/run00023/multiqc
mkdir -p /work/workspace/ruijuan/project/2_mwgs/run00023/clean


perl /work/workspace/ruijuan/script/2_MWGS/md5checkpro.pl -i /work/rawdata/run/guangzhou/novogene/20190321/run00023 -o /work/workspace/ruijuan/project/2_mwgs/run00023/md5checkresult # result in  md5checkresult.txt
notes: "/work/workspace/ruijuan/project/2_mwgs/run00023/md5checkresult" should exist
perl /work/workspace/ruijuan/script/2_MWGS/qcstat.pl -i /work/workspace/ruijuan/project/2_mwgs/run00023/input/sampledir.txt -o /work/workspace/ruijuan/project/2_mwgs/run00023/qcstat --basecutoff 0.8 --cutoffq20 90 --cutoffq30 90
perl /work/workspace/ruijuan/script/2_MWGS/datafilter.pl -i /work/workspace/ruijuan/project/2_mwgs/run00023/input/sampledir.txt -o /work/workspace/ruijuan/project/2_mwgs/run00023/clean







qcstat.pl -i qcstat/input/test.txt -o qcstat/output --readcutoff 2000 --basecutoff 0.00001 --cutoffq20 90 --cutoffq30 90
sudo perl 1_refdownload_done.pl -i /work/workspace/ruijuan/script/3_rawdataOrg/input/speinfo_all20190325.txt \
            -r /work/workspace/ruijuan/script/3_rawdataOrg/ref/prokaryotes.txt -o /mnt/d/work/database/ncbi/genome/all


# spades
mkdir -p /work/workspace/ruijuan/project/2_mwgs/run00023/assemble

/home/p/anaconda/anaconda3_5.2.0/bin/spades.py -1 /work/workspace/ruijuan/project/2_mwgs/run00023/clean/MNH03591.clean.1.fq.gz -2 /work/workspace/ruijuan/project/2_mwgs/run00023/clean/MNH03591.clean.2.fq.gz -t 20 --careful --cov-cutoff 5 -o /work/workspace/ruijuan/project/2_mwgs/run00023/assemble

GCE: predict genome size
cd /work/workspace/ruijuan/project/2_mwgs/run00023/gce
create read.list # content as the following two lines
/work/workspace/ruijuan/project/2_mwgs/run00023/clean/MNH03591.clean.1.fq.gz
/work/workspace/ruijuan/project/2_mwgs/run00023/clean/MNH03591.clean.2.fq.gz
 
/mnt/d/home/ruijuan/workflows/software/GCE/gce-1.0.0/kmerfreq/kmer_freq_hash/kmer_freq_hash -k 21 -l gce_read.list -t 20 -i 50000000 -o 0 -p MNH03591 &> kmer_freq.log

/mnt/d/home/ruijuan/workflows/software/GCE/gce-1.0.0/gce \
  -f MNH03591.freq.stat -m 1 -D 8 -b 1 > MNH03591.table 2> test.log
raw_peak        now_node        low_kmer        now_kmer        cvg     genome_size     a[1]    b[1]
206     5163040 20576065        1074704270      204.365 5.28453e+06     1       1.00394  

# soapdenovo

kmergenie
mkdir -p /work/workspace/ruijuan/project/2_mwgs/run00023/kmergenie
/home/p/KmerGenie/kmergenie-1.7048/kmergenie /work/workspace/ruijuan/project/2_mwgs/run0-o 0023/gce/gce_read.list -o /work/workspace/ruijuan/project/2_mwgs/run00023/kmergenie/kmergenie -l 21 -k 121 -s 10 -t 20 

mkdir -p /work/workspace/ruijuan/project/2_mwgs/run00023/soapdenovo-assemble
/home/p/anaconda/anaconda3_5.2.0/bin/SOAPdenovo-127mer -s /work/workspace/ruijuan/project/2_mwgs/run00023/input/soapdenovo_config.txt -u -p 20 -d 1 -D 1 -K 71 -R -o MNH03591
/home/p/anaconda/anaconda3_5.2.0/bin/SOAPdenovo-127mer all -s /work/workspace/ruijuan/project/2_mwgs/run00023/input/soapdenovo_config.txt -u -p 20 -d 1 -D 1 -K 97 -R -o MNH03591

quast
mkdir -p /work/workspace/ruijuan/project/2_mwgs/run00023/quast
/home/p/quast/quast-4.6.3/quast.py

quast.py ~/Seqs/SPAdesout_7942/contigs.fasta -o ~/Seqs/SPAdesout_7942/quast_out
quast.py -o compare_spa_velvet ./SPAdesout_7942_new/contigs.fasta ./velvet_out/contigs.fa
/home/p/quast/quast-4.6.3/quast.py -t 20 -o /work/workspace/ruijuan/project/2_mwgs/run00023/quast /work/workspace/ruijuan/script/2_MWGS/soapdenovo-yj/MNH03591.fas /work/workspace/ruijuan/project/2_mwgs/run00023/assemble/scaffolds.fasta /work/workspace/ruijuan/project/2_mwgs/run00023/soapdenovo-assemble/MNH03591.scafSeq
open report

/home/p/quast/quast-4.6.3/quast.py -t 20 /work/workspace/ruijuan/script/2_MWGS/kSNPtest/ksnpref/P5.fasta

perl /work/workspace/ruijuan/script/git/gitMWGS/genoassemble.pl -i /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00023/clean -o /work/workspace/ruijuan/project/2_mwgs/run00023/assemble_tmp -t 20

perl /work/workspace/ruijuan/script/git/gitMWGS/genoassemble.pl -i /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00023/clean/test -o /work/workspace/ruijuan/project/2_mwgs/run00023/assemble_tmp -t 20

perl genoassemble.pl -i /work/workspace/ruijuan/project/2_mwgs/run00023/clean/test/ -o /work/workspace/ruijuan/project/2_mwgs/run00023/assemble_tmp/ -t 20 -l 21 --method 2 

/home/p/quast/quast-4.6.3/quast.py /work/workspace/ruijuan/project/2_mwgs/run00023/assemble_tmp/finalAssemble/MNH03591.fasta /work/workspace/ruijuan/project/2_mwgs/run00023/assemble_tmp/finalAssemble/MNH03591.scafSeq -o /work/workspace/ruijuan/project/2_mwgs/run00023/assemble_tmp/quastResult/MNH03591


ANI 
/mnt/d/home/yajun/.conda/envs/yajun-env/bin/average_nucleotide_identity.py

run00020 - run00025 assemble
dir='/work/workspace/ruijuan/project/2_mwgs/'
for i in {20..25}
do
	file=$dir'run000'$i
	if [[ ! -d "$file" ]]; then
		mkdir $file
	else

	fi
done

find /work/rawdata/run -type l | grep gz > /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00020-25/input/run20-25file.txt

 perl /work/workspace/ruijuan/script/git/gitMWGS/datafilter.pl -i /work/workspace/ruijuan/project/2_mwgs/run00020-25/input/testfile.txt -o /work/workspace/ruijuan/project/2_mwgs/run00020-25/test

 nohup perl /work/workspace/ruijuan/script/git/gitMWGS/datafilter.pl -i /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00020-25/input/run20-25file.txt -o /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00020-25/result -t 16 > /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00020-25/log/datafilter-run20-25.log &
 nohup perl /work/workspace/ruijuan/script/git/gitMWGS/genoassemble.pl -i /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00020-25/result/clean -o /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00020-25/result -t 20 > /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00020-25/log/assemble-run20-25.log &

 cd /work/rawdata/run/guangzhou/novogene/20190402/run00026/1.rawdata
 ls | awk -F "[-_]" 'BEGIN {print "# run_acc\tseq_acc\tmn_acc\teid\tlid\tsid\tdate"}{print "run00026\t"$0"\t"$1"\t\t\t\t20190402"}' |sed 's/\///' > /work/workspace/ruijuan/script/3_rawdataOrg/input/run00026.txt

 cd /work/rawdata/run/guangzhou/novogene/20190402/run00027/1.rawdata
 ls | awk -F "[-_]" 'BEGIN {print "# run_acc\tseq_acc\tmn_acc\teid\tlid\tsid\tdate"}{print "run00027\t"$0"\t"$1"\t\t\t\t20190402"}' |sed 's/\///' > /work/workspace/ruijuan/script/3_rawdataOrg/input/run00027.txt

 nohup perl /work/workspace/ruijuan/script/3_rawdataOrg/3_rawdataorganize.pl -i /work/rawdata/run/guangzhou/novogene/20190402/run00026 -o /work/rawdata -r /work/workspace/ruijuan/script/3_rawdataOrg/input/run00026.txt > /work/workspace/ruijuan/script/3_rawdataOrg/log/run00026.txt &

 nohup perl /work/workspace/ruijuan/script/3_rawdataOrg/3_rawdataorganize.pl -i /work/rawdata/run/guangzhou/novogene/20190402/run00027 -o /work/rawdata -r /work/workspace/ruijuan/script/3_rawdataOrg/input/run00027.txt > /work/workspace/ruijuan/script/3_rawdataOrg/log/run00027.txt &

 nohup perl /work/workspace/ruijuan/script/3_rawdataOrg/3_rawdataorganize.pl -i /work/rawdata/run/guangzhou/novogene/20190225/run00020 -o /work/rawdata -r /work/workspace/ruijuan/script/3_rawdataOrg/input/run00020_mod.txt > /work/workspace/ruijuan/script/3_rawdataOrg/log/run00020_mod.txt &

 nohup perl /work/workspace/ruijuan/script/3_rawdataOrg/3_rawdataorganize.pl -i /work/rawdata/run/guangzhou/novogene/20190322/run00024 -o /work/rawdata -r /work/workspace/ruijuan/script/3_rawdataOrg/input/run00024_mod.txt > /work/workspace/ruijuan/script/3_rawdataOrg/log/run00024_mod.txt &

 perl /work/workspace/ruijuan/script/3_rawdataOrg/1_refdownload_done.pl -i /work/workspace/ruijuan/script/3_rawdataOrg/input/speinfo_20190308.txt -r /work/workspace/ruijuan/script/3_rawdataOrg/ref/prokaryotes.txt -o /mnt/d/work/database/ncbi/genome/all > log/log_ref_20190308.txt

 find /work/rawdata/run/guangzhou/2019/04/20190402 -type l | grep gz > /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00026-27/input/run26-27file.txt
 nohup perl /work/workspace/ruijuan/script/git/gitMWGS/datafilter.pl -i /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00026-27/input/run26-27file.txt -o /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00026-27/result -t 16 > /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00026-27/log/datafilter-run26-27.log &

 nohup perl /work/workspace/ruijuan/script/git/gitMWGS/genoassemble.pl -i /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00026-27/result/clean -o /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00026-27/result -t 20 > /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00026-27/log/assemble-run26-27.log &

 MNH19252, MNH19632 genome assemble
 find /work/rawdata/run/guangzhou -type l | grep -P "MNH19252|MNH19632" > /work/workspace/ruijuan/project/2_mwgs/run00021_4twosample/input/run21-24twosamples.txt
 nohup perl /work/workspace/ruijuan/script/git/gitMWGS/datafilter.pl -i /work/workspace/ruijuan/project/2_mwgs/run00021_4twosample/input/run21-24twosamples.txt -o /work/workspace/ruijuan/project/2_mwgs/run00021_4twosample/result -t 16 > /work/workspace/ruijuan/project/2_mwgs/run00021_4twosample/log/datafilter-run21-24twosamples.log &

 nohup perl /work/workspace/ruijuan/script/git/gitMWGS/genoassemble.pl -i /work/workspace/ruijuan/project/2_mwgs/run00021_4twosample/result/clean -o /work/workspace/ruijuan/project/2_mwgs/run00021_4twosample/result -t 20 > /work/workspace/ruijuan/project/2_mwgs/run00021_4twosample/log/assemble-run21_24twosamples.log &

cd /work/rawdata/run/guangzhou/novogene/20190415/run00028
 ls | awk -F "[-_]" 'BEGIN {print "# run_acc\tseq_acc\tmn_acc\teid\tlid\tsid\tdate"}{print "run00028\t"$0"\t"$1"\t\t\t\t20190415"}' |sed 's/\///' > /work/workspace/ruijuan/script/3_rawdataOrg/input/run00028.txt

 nohup perl /work/workspace/ruijuan/script/3_rawdataOrg/3_rawdataorganize.pl -i /work/rawdata/run/guangzhou/novogene/20190415/run00028 -o /work/rawdata -r /work/workspace/ruijuan/script/3_rawdataOrg/input/run00028.txt > /work/workspace/ruijuan/script/3_rawdataOrg/log/run00028.txt &

 mkdir -p /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00028/input
 mkdir -p /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00028/log
 find /work/rawdata/run/guangzhou/2019/04/20190415 -type l | grep gz > /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00028/input/run28file.txt
 nohup perl /work/workspace/ruijuan/script/git/gitMWGS/datafilter.pl -i /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00028/input/run28file.txt -o /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00028/result -t 16 > /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00028/log/datafilter-run28.log &

 nohup perl /work/workspace/ruijuan/script/git/gitMWGS/genoassemble.pl -i /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00028/result/clean -o /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00028/result -t 20 > /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00028/log/assemble-run28.log &

 mkdir /work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble
 mkdir /work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/in_genome
 mkdir /work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/result
 ln -s /work/workspace/ruijuan/project/2_mwgs/run00028/result/finalAssemble/*.fasta /work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/in_genome
 ln -s /work/workspace/ruijuan/project/2_mwgs/run00020-25/result/finalAssemble/*.fasta /work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/in_genome
 ln -s /work/workspace/ruijuan/project/2_mwgs/run00021_4twosample/result/finalAssemble/*.fasta /work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/in_genome
 ln -s /work/workspace/ruijuan/project/2_mwgs/run00026-27/result/finalAssemble/*.fasta /work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/in_genome

 perl assembleOrg.pl -i /work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/filter_genome -r /work/workspace/ruijuan/script/3_rawdataOrg/input/run00020.txt -o /work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/test

 sudo perl /work/workspace/ruijuan/script/git/gitMWGS/assembleOrg.pl -i /work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/filter_genome -r /work/workspace/ruijuan/script/3_rawdataOrg/input/run00020.txt -o /work/assembly
 sudo perl /work/workspace/ruijuan/script/git/gitMWGS/assembleOrg.pl -i /work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/filter_genome -r /work/workspace/ruijuan/script/3_rawdataOrg/input/run00021.txt -o /work/assembly
 sudo perl /work/workspace/ruijuan/script/git/gitMWGS/assembleOrg.pl -i /work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/filter_genome -r /work/workspace/ruijuan/script/3_rawdataOrg/input/run00022.txt -o /work/assembly
 sudo perl /work/workspace/ruijuan/script/git/gitMWGS/assembleOrg.pl -i /work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/filter_genome -r /work/workspace/ruijuan/script/3_rawdataOrg/input/run00023.txt -o /work/assembly
 sudo perl /work/workspace/ruijuan/script/git/gitMWGS/assembleOrg.pl -i /work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/filter_genome -r /work/workspace/ruijuan/script/3_rawdataOrg/input/run00024.txt -o /work/assembly
 sudo perl /work/workspace/ruijuan/script/git/gitMWGS/assembleOrg.pl -i /work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/filter_genome -r /work/workspace/ruijuan/script/3_rawdataOrg/input/run00025.txt -o /work/assembly
 sudo perl /work/workspace/ruijuan/script/git/gitMWGS/assembleOrg.pl -i /work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/filter_genome -r /work/workspace/ruijuan/script/3_rawdataOrg/input/run00026.txt -o /work/assembly
 sudo perl /work/workspace/ruijuan/script/git/gitMWGS/assembleOrg.pl -i /work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/filter_genome -r /work/workspace/ruijuan/script/3_rawdataOrg/input/run00027.txt -o /work/assembly
 sudo perl /work/workspace/ruijuan/script/git/gitMWGS/assembleOrg.pl -i /work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/filter_genome -r /work/workspace/ruijuan/script/3_rawdataOrg/input/run00028.txt -o /work/assembly

 nohup perl /work/workspace/ruijuan/script/git/gitMWGS/genoPhyana.pl -i /mnt/d/work/workspace/ruijuan/script/3_rawdataOrg/input/groupInfo349sample.txt -o /work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/result -d /work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/in_genome --ani > /work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/log/run20-28phylogeneticAni.log &

 nohup perl /work/workspace/ruijuan/script/git/gitMWGS/genoPhyana_nodel.pl -i /mnt/d/work/workspace/ruijuan/script/3_rawdataOrg/input/groupInfo349sample.txt -o /work/workspace/ruijuan/project/2_mwgs/run00020-28phyana/result -d /work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/in_genome --ani > /work/workspace/ruijuan/project/2_mwgs/run00020-28phyana/log/run20-28phylogeneticAni.log &

 2019.May29-18:50
 run00030, run00031
 原始数据存储
 su root
 perl /work/workspace/ruijuan/script/3_rawdataOrg/3_rawdataorganize.pl -i /work/rawdata/run/guangzhou/novogene/20190529/run00030 -o /work/rawdata -r /mnt/d/work/workspace/ruijuan/script/3_rawdataOrg/input/run00030.txt
 perl /work/workspace/ruijuan/script/3_rawdataOrg/3_rawdataorganize.pl -i /work/rawdata/run/guangzhou/novogene/20190529/run00031 -o /work/rawdata -r /mnt/d/work/workspace/ruijuan/script/3_rawdataOrg/input/run00031.txt

 find /work/rawdata/run/guangzhou/2019/05/20190529 -type l | grep gz | sort -n > /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00030-31/input/run30-31file.txt

 bash 运行 datafilter 和 genoassemble
 cd /work/workspace/ruijuan/project/2_mwgs/run00030-31/program
 filterAssemble.sh 
 nohup bash /work/workspace/ruijuan/project/2_mwgs/run00030-31/program/filterAssemble.sh > /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00030-31/log/filterAssebmle30-31.txt &
perl /work/workspace/ruijuan/script/git/gitMWGS/datafilter.pl -i /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00030-31/input/run30-31file.txt -o /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00030-31/result -t 16 > /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00030-31/log/datafilter-run30-31.log 

perl /work/workspace/ruijuan/script/git/gitMWGS/genoassemble.pl -i /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00030-31/result/clean -o /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00030-31/result -t 20 > /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00028/log/assemble-run30-31.log 

2019.Jun02-14:38
find /work/workspace/ruijuan/project/2_mwgs/run00030-31/result/finalAssemble -type f | grep fasta | awk -F "[/.]" '{print "ln -s "$0" /work/workspace/ruijuan/project/2_mwgs/phyana/genome/ingenome_run20-31/"$(NF-1)".fasta"}' > /work/workspace/ruijuan/project/2_mwgs/p/run30-31genomeSoftlink.sh
bash /work/workspace/ruijuan/project/2_mwgs/p/run30-31genomeSoftlink.sh

nohup perl /work/workspace/ruijuan/script/git/gitMWGS/genoPhyana_nodel.pl -i /mnt/d/work/workspace/ruijuan/script/3_rawdataOrg/input/groupInfo396sample.txt -o /work/workspace/ruijuan/project/2_mwgs/phyana/run20-31 -d /work/workspace/ruijuan/project/2_mwgs/phyana/genome/ingenome_run20-31 --ani > /work/workspace/ruijuan/project/2_mwgs/phyana/log/run20-31phylogeneticAni.log &
 
$vresultdir = "$resultdir/fasta/$dbversion";
 95 $runresultdir = "$resultdir/run/$city";
 96 $gresultdir = "$resultdir/fasta/genome";
 97 $currentdir = "$resultdir/current";

test:
perl /work/workspace/ruijuan/script/git/assembleOrg/assembleOrg.pl -i /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00030-31/result/finalAssemble -o /work/workspace/ruijuan/project/2_mwgs/test/assembleTest -r /work/workspace/ruijuan/script/git/assembleOrg/input/run00030.txt

run00020-run00031
find run/ | grep fas | wc -l
su root
perl /work/workspace/ruijuan/script/git/assembleOrg/assembleOrg.pl -i /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00030-31/result/finalAssemble -o /work/assembly -r /work/workspace/ruijuan/script/git/assembleOrg/input/run00030.txt
perl /work/workspace/ruijuan/script/git/assembleOrg/assembleOrg.pl -i /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00030-31/result/finalAssemble -o /work/assembly -r /work/workspace/ruijuan/script/git/assembleOrg/input/run00031.txt
perl /work/workspace/ruijuan/script/git/assembleOrg/assembleOrg.pl -i /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/filter_genome -o /work/assembly -r /work/workspace/ruijuan/script/git/assembleOrg/input/run00020.txt
perl /work/workspace/ruijuan/script/git/assembleOrg/assembleOrg.pl -i /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/filter_genome -o /work/assembly -r /work/workspace/ruijuan/script/git/assembleOrg/input/run00021.txt
perl /work/workspace/ruijuan/script/git/assembleOrg/assembleOrg.pl -i /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/filter_genome -o /work/assembly -r /work/workspace/ruijuan/script/git/assembleOrg/input/run00022.txt
perl /work/workspace/ruijuan/script/git/assembleOrg/assembleOrg.pl -i /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/filter_genome -o /work/assembly -r /work/workspace/ruijuan/script/git/assembleOrg/input/run00023.txt
perl /work/workspace/ruijuan/script/git/assembleOrg/assembleOrg.pl -i /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/filter_genome -o /work/assembly -r /work/workspace/ruijuan/script/git/assembleOrg/input/run00024.txt
perl /work/workspace/ruijuan/script/git/assembleOrg/assembleOrg.pl -i /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/filter_genome -o /work/assembly -r /work/workspace/ruijuan/script/git/assembleOrg/input/run00025.txt
perl /work/workspace/ruijuan/script/git/assembleOrg/assembleOrg.pl -i /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/filter_genome -o /work/assembly -r /work/workspace/ruijuan/script/git/assembleOrg/input/run00026.txt
perl /work/workspace/ruijuan/script/git/assembleOrg/assembleOrg.pl -i /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/filter_genome -o /work/assembly -r /work/workspace/ruijuan/script/git/assembleOrg/input/run00027.txt
perl /work/workspace/ruijuan/script/git/assembleOrg/assembleOrg.pl -i /mnt/d/work/workspace/ruijuan/project/2_mwgs/run00020-28Assemble/filter_genome -o /work/assembly -r /work/workspace/ruijuan/script/git/assembleOrg/input/run00028.txt

## prokka annotation
/home/p/prokka/prokka-master/bin/prokka --kingdom Bacteria --cpus 20  --outdir /work/workspace/ruijuan/project/2_mwgs/test/annotation/result  --prefix MN211086 --force --metagenome --locustag MN211086 /work/workspace/ruijuan/project/2_mwgs/test/assembleTest/run/guangzhou/novogene/2019/05/20190529/run00030/MN211086.fas 
tbl2asn -V b -a r10k -l paired-ends -M n -N 1 -y 'Annotated using prokka 1.14-dev from https://github.com/tseemann/prokka' -Z ./MNH04496.err -i ./MNH04496.fsa 2> /dev/null

/home/p/quast/quast-4.6.3/quast.py
/home/p/blast/bin/blastn
/home/p/blast/bin/makeblastdb
/home/p/quast/quast-4.6.3/quast_libs/MUMmer/nucmer
python {QUAST} -t {thread} -o quast ./prokka/{bacteriaStrain}.genome.fa ./CISA/abyss.fa.p.fa ./CISA/velvet.fa.p.fa ./CISA/spades.fa.p.fa ./CISA/IDBA-UD.fa.p.fa ./CISA/soapdenovo.fa.p.fa

2019.Jun11-16:27
农业组样品进化分析
nohup perl /work/workspace/ruijuan/script/git/gitMWGS/genoPhyana_nodel.pl -i /work/workspace/ruijuan/project/2_mwgs/run00030-31/input/methy_algricultureGroup_20190611.txt -o /work/workspace/ruijuan/project/2_mwgs/phyana/run20-31/groupMN37sample20190611 -d /work/workspace/ruijuan/project/2_mwgs/phyana/genome/ingenome_run20-31 --ani > /work/workspace/ruijuan/project/2_mwgs/phyana/log/groupmn37sample20190611.log &

run00038
mkdir -p /work/workspace/ruijuan/project/2_mwgs/run00038/tracking
mkdir -p /work/workspace/ruijuan/project/2_mwgs/run00038/input
mkdir -p /work/workspace/ruijuan/project/2_mwgs/run00038/result
sudo mkdir -p /work/rawdata/run/guangzhou/ novogene/2019/07/20190729/run00038/
数据下载：  sudo bash /mnt/d/home/x/software/ossdownload.sh -k LTAI6hHDkkxiP5U4 -s q1ABZDVUr2czJhEI2gc2zEII2cVCBA -c oss://novo-data-nj/customer-znZW7VCU/ -d /work/rawdata/run/guangzhou/novogene/2019/07/20190729/run00038/rawdata --location huangzhou
原始数据存储：  sudo perl /work/workspace/ruijuan/script/git/rawdataOrg/rawdataOrg.pl -i /work/rawdata/run/guangzhou/novogene/2019/07/20190729/run00038/rawdata -o /work/rawdata -r /work/rawdata/run/guangzhou/ novogene/2019/07/20190729/run00038/run00038.txt -p genome
数据过滤所需文件creation： 
find /work/rawdata/run/guangzhou/2019/07/20190729/run00038/fastq | grep gz$ | sort -u | awk -F "." '{if($(NF-2) == 1){tmp=$0;getline;print tmp"\t"$0}}' > /work/workspace/ruijuan/project/2_mwgs/run00038/input/infile.txt
数据过滤： perl /work/workspace/ruijuan/script/git/gitMWGS/datafilter.pl -i /work/workspace/ruijuan/project/2_mwgs/run00038/input/infile.txt -o /work/workspace/ruijuan/project/2_mwgs/run00038/result
序列拼装：
序列拼装： nohup perl /work/workspace/ruijuan/script/git/gitMWGS/genoassemble.pl -i /work/workspace/ruijuan/project/2_mwgs/run00038/result/clean -o /work/workspace/ruijuan/project/2_mwgs/run00038/result/assemble -t 20 > /work/workspace/ruijuan/project/2_mwgs/run00038/tracking/run00038_assemble_20190730.txt 2>&1 &
拼装结果存储：

# 2019.Sep20-11:26
# install resfinder
# ref: https://bitbucket.org/genomicepidemiology/resfinder/src/master/
mkdir -p /work/workspace/zhurj/software/vfanno
cd /work/workspace/zhurj/software/vfanno
git clone https://git@bitbucket.org/genomicepidemiology/resfinder.git
git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git

cd cd /work/workspace/zhurj/download
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.9.0+-x64-linux.tar.gz


# install pip3 without root
# https://unix.stackexchange.com/questions/445906/how-to-get-pip3-without-sudo-privileges
source activate base
wget https://bootstrap.pypa.io/get-pip.py
python3 get-pip.py --user
pip3 --version

python3 resfinder.py -i /work/workspace/zhurj/reference/strainest/AKK/ref/fna/MNH20651.fas -o ./test -p /work/workspace/zhurj/software/vfanno/resfinder_db -mp /work/workspace/zhurj/bin/miniconda3/bin/blastn -d aminoglycoside -t 0.80 -l 0.60
python3 resfinder.py -i /work/workspace/zhurj/reference/strainest/AKK/ref/fna/MNH17505.fas -o ./test -p /work/workspace/zhurj/software/vfanno/resfinder_db \
-mp /work/workspace/zhurj/bin/miniconda3/bin/blastn -d aminoglycoside -t 0.80 -l 0.60
python3 resfinder.py -i /work/workspace/zhurj/reference/strainest/AKK/ref/fna/MNH19250.fas -o ./test -p /work/workspace/zhurj/software/vfanno/resfinder_db \
-mp /work/workspace/zhurj/bin/miniconda3/bin/blastn -d aminoglycoside -t 0.80 -l 0.60
python3 resfinder.py -i /work/workspace/zhurj/reference/strainest/AKK/ref/fna/MNH19252.fas -o ./test -p /work/workspace/zhurj/software/vfanno/resfinder_db \
-mp /work/workspace/zhurj/bin/miniconda3/bin/blastn -d aminoglycoside -t 0.80 -l 0.60
python3 resfinder.py -i /work/workspace/zhurj/reference/strainest/AKK/ref/fna/MNH19631.fas -o ./test -p /work/workspace/zhurj/software/vfanno/resfinder_db \
-mp /work/workspace/zhurj/bin/miniconda3/bin/blastn -d aminoglycoside -t 0.80 -l 0.60
python3 resfinder.py -i /work/workspace/zhurj/reference/strainest/AKK/ref/fna/MNH19632.fas -o ./test -p /work/workspace/zhurj/software/vfanno/resfinder_db \
-mp /work/workspace/zhurj/bin/miniconda3/bin/blastn -d aminoglycoside -t 0.80 -l 0.60
python3 resfinder.py -i /work/workspace/zhurj/reference/strainest/AKK/ref/fna/MNH19714.fas -o ./test -p /work/workspace/zhurj/software/vfanno/resfinder_db \
-mp /work/workspace/zhurj/bin/miniconda3/bin/blastn -d aminoglycoside -t 0.80 -l 0.60
python3 resfinder.py -i /work/workspace/zhurj/reference/strainest/AKK/ref/fna/MNH17503.fas -o ./test -p /work/workspace/zhurj/software/vfanno/resfinder_db \
-mp /work/workspace/zhurj/bin/miniconda3/bin/blastn -d aminoglycoside -t 0.80 -l 0.60




python3 resfinder.py -i test.fsa -o ./test -p /work/workspace/zhurj/software/vfanno/resfinder_db \
-mp /work/workspace/zhurj/bin/miniconda3/bin/blastn -d aminoglycoside -t 0.90 -l 0.60

python3 resfinder.py -i /work/workspace/zhurj/reference/strainest/PD/ref/fna/MNH09897.fas -o ./test -p /work/workspace/zhurj/software/vfanno/resfinder_db -mp /work/workspace/zhurj/bin/miniconda3/bin/blastn -d aminoglycoside -t 0.90 -l 0.60
python3 resfinder.py -i /work/workspace/zhurj/reference/strainest/PD/ref/fna/MNH07338.fas -o ./test -p /work/workspace/zhurj/software/vfanno/resfinder_db -mp /work/workspace/zhurj/bin/miniconda3/bin/blastn -d aminoglycoside -t 0.90 -l 0.60


/work/workspace/zhurj/bin/miniconda3/envs/tcsh/bin/tcsh /work/workspace/zhurj/software/kSNP3.1_Linux_package/kSNP3
java -jar /work/workspace/zhurj/software/FigTree_v1.4.4/lib/figtree.jar 


/work/workspace/zhurj/bin/miniconda3/bin/mash sketch -o /work/workspace/zhurj/project/2_swgs/test/mash/output/akk -l /work/workspace/zhurj/project/2_swgs/test/mash/input/input.txt

/work/workspace/zhurj/bin/miniconda3/bin/mash dist -t /work/workspace/zhurj/project/2_swgs/test/mash/output/akk.msh /work/workspace/zhurj/project/2_swgs/test/mash/output/akk.msh > /work/workspace/zhurj/project/2_swgs/test/mash/output/mash.dist

Clustering: https://uc-r.github.io/hc_clustering
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms

df <- read.table("/work/workspace/zhurj/project/2_swgs/test/mash/output/mash_adj.dist",sep = "\t",header=T,row.names = 1)
dm <- as.matrix(df)
ds <- as.dist(dm, diag = FALSE, upper = FALSE)
hc <- hclust(ds, method = "complete")
dev.new()
png("test.png",width=800,height=400)
png("test.png",width=4000,height=1000)
plot(hc,hang=-1) 
dev.off()

# new server

bwa mem -t 16 /work/workspace/zhurj/reference/NCBI/genome/GCF_000154205/bwa/GCF_000154205.fna.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH31610.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH31610.clean.2.fq.gz -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH31610/MNH31610.sam
samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH31610/MNH31610.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH31610/MNH31610.sam
samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH31610/MNH31610.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/stat/MNH31610.txt

samtools idxstats --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH31610/MNH31610.sorted.sam
OUTPUT of samtools idxstats
ref.name ref.length mapped.reads unmapped.reads


mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH31671
bwa mem -t 16 /work/workspace/zhurj/reference/NCBI/genome/GCF_000154205/bwa/GCF_000154205.fna.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH31671.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH31671.clean.2.fq.gz -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH31671/MNH31671.sam
samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH31671/MNH31671.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH31671/MNH31671.sam
samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH31671/MNH31671.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/stat/MNH31671.txt

mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH33080
bwa mem -t 16 /work/workspace/zhurj/reference/NCBI/genome/GCF_000154205/bwa/GCF_000154205.fna.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33080.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33080.clean.2.fq.gz -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH33080/MNH33080.sam
samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH33080/MNH33080.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH33080/MNH33080.sam
samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH33080/MNH33080.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/stat/MNH33080.txt

mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH33156
bwa mem -t 16 /work/workspace/zhurj/reference/NCBI/genome/GCF_000154205/bwa/GCF_000154205.fna.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33156.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33156.clean.2.fq.gz -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH33156/MNH33156.sam
samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH33156/MNH33156.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH33156/MNH33156.sam
samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH33156/MNH33156.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/stat/MNH33156.txt

mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH33403
bwa mem -t 16 /work/workspace/zhurj/reference/NCBI/genome/GCF_000154205/bwa/GCF_000154205.fna.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33403.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33403.clean.2.fq.gz -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH33403/MNH33403.sam
samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH33403/MNH33403.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH33403/MNH33403.sam
samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH33403/MNH33403.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/stat/MNH33403.txt

mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH36188
bwa mem -t 16 /work/workspace/zhurj/reference/NCBI/genome/GCF_000154205/bwa/GCF_000154205.fna.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36188.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36188.clean.2.fq.gz -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH36188/MNH36188.sam
samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH36188/MNH36188.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH36188/MNH36188.sam
samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH36188/MNH36188.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/stat/MNH36188.txt

mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH36362
bwa mem -t 16 /work/workspace/zhurj/reference/NCBI/genome/GCF_000154205/bwa/GCF_000154205.fna.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36362.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36362.clean.2.fq.gz -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH36362/MNH36362.sam
samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH36362/MNH36362.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH36362/MNH36362.sam
samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/sam/MNH36362/MNH36362.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/stat/MNH36362.txt

perl /work/workspace/zhurj/script/git/tools/getfile/getfile.pl -i /work/workspace/zhurj/project/5_tools/getfile/input/input.txt -n test -o /work/workspace/zhurj/project/5_tools/getfile/output 
perl /work/workspace/zhurj/script/git/tools/getfile/getfile.pl -i /work/workspace/zhurj/project/5_tools/getfile/input/input.txt -n test -o /work/workspace/zhurj/project/5_tools/getfile/output --quiet

# /work/classify/current/MNH/MNH048/MNH04835/MNH04835.profile
# /work/classify/current/MNH/MNH104/MNH10420/MNH10420.profile
# /work/classify/current/MNH/MNH180/MNH18018/MNH18018.profile
# /work/classify/current/MNH/MNH187/MNH18780/MNH18780.profile
# /work/classify/current/MNH/MNH023/MNH02324/MNH02324.profile
# /work/classify/current/MNH/MNH027/MNH02737/MNH02737.profile

checkm ssu_finder /work/assembly/current/MNH/MNH047/MNH04744/genomic.fna /work/assembly/current/MNH/MNH047/MNH04744 ./extract16S -x fna -t 20
checkm lineage_wf /work/assembly/current/MNH/MNH047/MNH04744 /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH04744 -x fna -t 20

quast under python3.6
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/quast/MNH06145 /work/assembly/current/MNH/MNH061/MNH06145/genomic.fna




checkm ssu_finder /work/assembly/current/MNH/MNH061/MNH06145/genomic.fna /work/assembly/current/MNH/MNH061/MNH06145 /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH06145 -x fna -t 20 
checkm ssu_finder /work/assembly/current/MNH/MNH047/MNH04744/genomic.fna /work/assembly/current/MNH/MNH047/MNH04744 /work/workspace/zhurj/project/2_swgs/taxcheck/extract16S/MNH04744 -x fna -t 20 -c 500
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/quast/MNH06145 /work/assembly/current/MNH/MNH061/MNH06145/genomic.fna

MNH26050
GCA_000020225： 2.76M, 55.3%
quast: 4.83M, 50.69%
AKK: ANI 98.86%, 59.7%(836     1400)
lactobacillus ruminis: ANI 96.26% 29.86% (418     1400)
checkm contamination: 94.46
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/quast/MNH26050
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/quast/MNH26050 /work/assembly/current/MNH/MNH260/MNH26050/genomic.fna

mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH26050
checkm ssu_finder /work/assembly/current/MNH/MNH260/MNH26050/genomic.fna /work/assembly/current/MNH/MNH260/MNH26050 /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH26050 -x fna -t 20 
checkm lineage_wf /work/assembly/current/MNH/MNH260/MNH26050 /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH26050 -x fna -t 20

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH26050/ssu_filter.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH26050/MNH26050.out -max_target_seqs 5 -num_threads 16 -outfmt 7 

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "BCVU01000117|CP001071|PJKF01000002|PJKB01000002"
BCVU01000117    Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;Lactobacillus ruminis;NBRC 102161(T)
CP001071        Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;Akkermansia muciniphila;ATCC BAA-835(T)
PJKB01000002    Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;PJKB_s;GP22
PJKF01000002    Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;PJKF_s;GP15


fastANI -q /work/assembly/current/MNH/MNH260/MNH26050/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000020/GCA_000020225/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH26050_AKK
fastANI -q /work/assembly/current/MNH/MNH260/MNH26050/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000224/GCA_000224985/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH26050_lactobacillus_ruminis



MNH27359
GCA_000020225： 2.76M, 55.3%
quast: 5.17M, 50.82%
AKK: ANI 98.49%, 66.83%(828/1239)
ruminococcus gnavus: ANI 96.09% 20.82% (258     1239)
checkm contamination: 74.31


mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/quast/MNH27359
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/quast/MNH27359 /work/assembly/current/MNH/MNH273/MNH27359/genomic.fna

mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH27359
checkm ssu_finder /work/assembly/current/MNH/MNH273/MNH27359/genomic.fna /work/assembly/current/MNH/MNH273/MNH27359 /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH27359 -x fna -t 20 
checkm lineage_wf /work/assembly/current/MNH/MNH273/MNH27359 /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH27359 -x fna -t 20

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH27359/ssu_filter.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH27359/MNH27359.out -max_target_seqs 5 -num_threads 16 -outfmt 7 

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "AAYG02000025|CP001071|PJKF01000002|PJKB01000002"
AAYG02000025    Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Mediterraneibacter;Ruminococcus gnavus;ATCC 29149(T)
CP001071        Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;Akkermansia muciniphila;ATCC BAA-835(T)
PJKB01000002    Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;PJKB_s;GP22
PJKF01000002    Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;PJKF_s;GP15

fastANI -q /work/assembly/current/MNH/MNH273/MNH27359/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000020/GCA_000020225/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH27359_AKK
fastANI -q /work/assembly/current/MNH/MNH273/MNH27359/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000526/GCA_000526735/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH27359_ruminococcus_gnavus

MNH27614 
GCA_000020225： 2.76M, 55.3%
quast: 6.987M, GC 50.86%
AKK: ANI 97.47% 65.3%
ruminococcus gnavus: ANI 96.62% 22.49% (280     1245)
checkm: contamination 120.13


split -l 100 -e --additional-suffix _20200110 -d sam293_20200110.sh
srun -o x00.out -e x00.err -N 1 -c 20 -p slurm256 bash x00_20200110 &
srun -o x01.out -e x01.err -N 1 -c 20 -p slurm256 bash x01_20200110 &
srun -o x02.out -e x02.err -N 1 -c 20 -p slurm128 bash x02_20200110 &

quast -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/quast/MNH27614 /work/assembly/current/MNH/MNH276/MNH27614/genomic.fna
checkm ssu_finder /work/assembly/current/MNH/MNH276/MNH27614/genomic.fna /work/assembly/current/MNH/MNH276/MNH27614 /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH27614 -x fna -t 20 
checkm lineage_wf /work/assembly/current/MNH/MNH276/MNH27614 /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH27614 -x fna -t 20 --tab_table -f /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH27614/MNH27614_checkm.txt
mkdir -p /work/workspace/zhurj/project/2_swgs/run00050/checkM/MNH32731 && checkm ssu_finder /work/assembly/current/MNH/MNH327/MNH32731/genomic.fna /work/assembly/current/MNH/MNH327/MNH32731 /work/workspace/zhurj/project/2_swgs/run00050/checkM/MNH32731 -x fna -t 20 && checkm lineage_wf /work/assembly/current/MNH/MNH327/MNH32731  /work/workspace/zhurj/project/2_swgs/run00050/checkM/MNH32731 -x fna -t 20 --tab_table -f /work/workspace/zhurj/project/2_swgs/run00050/checkM/MNH32731/MNH32731_checkm.txt




fastANI -q /work/assembly/current/MNH/MNH276/MNH27614/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000020/GCA_000020225/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH27614_AKK
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH27614/ssu_filter.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH27614/MNH27614.out -max_target_seqs 5 -num_threads 16 -outfmt 7 
fastANI -q /work/assembly/current/MNH/MNH276/MNH27614/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000010/GCA_000010425/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH27614_bifi_adolescentis
fastANI -q /work/assembly/current/MNH/MNH276/MNH27614/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000526/GCA_000526735/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH27614_ruminococcus_gnavus

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "AP009256|KF990498|CP001071|PJKF01000002|PJKB01000002|AAYG02000025|CP022823|AJJI01000018|CP010523|HG933296|CBZR010000040"


/work/database/ezbio/16S/current/Ezbio_16S_seqs
/work/database/ncbi/16S/current/16SMicrobial
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH06145/ssu.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH06145/MNH06145.out -max_target_seqs 5 -num_threads 16 -outfmt 7 
perl /work/workspace/zhouyj/script/alignKit/quick-16S2Taxonomy.pl -i /work/assembly/current/MNH/MNH051/MNH05119/genomic.fna -o MNH05119 -sq
/work/workspace/zhouyj/script/alignKit/quick-blast.pl
/work/workspace/zhouyj/script/alignKit/blast2annotation.pl

fastani: https://anaconda.org/bioconda/fastani
busco : https://anaconda.org/bioconda/busco; https://busco.ezlab.org/

fastANI -q /work/assembly/current/MNH/MNH061/MNH06145/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000007/GCA_000007525/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH06145_GCA_000007525
fastANI --ql /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/input --rl /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/input -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH06145_GCA_000007525

MNH20304
fastANI -q /work/assembly/current/MNH/MNH203/MNH20304/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000156/GCA_000156535/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH20304_intestinalis
fastANI -q /work/assembly/current/MNH/MNH203/MNH20304/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000225/GCA_000225345/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH20304_hominis
fastANI -q /work/assembly/current/MNH/MNH203/MNH20304/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_001/GCA_001405/GCA_001405615/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH20304_faecis

MNH20393
fastANI -q /work/assembly/current/MNH/MNH203/MNH20393/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000156/GCA_000156535/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH20393_intestinalis
fastANI -q /work/assembly/current/MNH/MNH203/MNH20393/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000225/GCA_000225345/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH20393_hominis
fastANI -q /work/assembly/current/MNH/MNH203/MNH20393/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_001/GCA_001405/GCA_001405615/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH20393_faecis


fastANI -q /work/assembly/current/MNH/MNH203/MNH20394/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000156/GCA_000156535/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH20394_intestinalis
fastANI -q /work/assembly/current/MNH/MNH203/MNH20394/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000225/GCA_000225345/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH20394_hominis
fastANI -q /work/assembly/current/MNH/MNH203/MNH20394/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_001/GCA_001405/GCA_001405615/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH20394_faecis


fastANI -q /work/assembly/current/MNH/MNH203/MNH20395/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000156/GCA_000156535/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH20395_intestinalis
fastANI -q /work/assembly/current/MNH/MNH203/MNH20395/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000225/GCA_000225345/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH20395_hominis
fastANI -q /work/assembly/current/MNH/MNH203/MNH20395/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_001/GCA_001405/GCA_001405615/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH20395_faecis

fastANI -q /work/assembly/current/MNH/MNH203/MNH20397/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000156/GCA_000156535/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH20397_intestinalis
fastANI -q /work/assembly/current/MNH/MNH203/MNH20397/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000225/GCA_000225345/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH20397_hominis
fastANI -q /work/assembly/current/MNH/MNH203/MNH20397/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_001/GCA_001405/GCA_001405615/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH20397_faecis

fastANI -q /work/assembly/current/MNH/MNH104/MNH10420/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000162/GCA_000162015/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH10420_prausnitzii
fastANI -q /work/assembly/current/MNH/MNH104/MNH10420/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000012/GCA_000012845/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH10420_pd


fastANI -q /work/assembly/current/MNH/MNH196/MNH19699/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000981/GCA_000981035/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH19699_Catabacter_hongkongensis

fastANI -q /work/assembly/current/MNH/MNH063/MNH06322/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000949/GCA_000949455/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH06322_Ruthenibacterium_lactatiformans

fastANI -q /work/assembly/current/MNH/MNH063/MNH06321/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000949/GCA_000949455/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH06321_Ruthenibacterium_lactatiformans


checkm ssu_finder /work/assembly/current/MNH/MNH104/MNH10420/genomic.fna /work/assembly/current/MNH/MNH104/MNH10420 /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH10420 -x fna -t 20
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH10420/ssu_filter.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH10420/MNH10420.out -max_target_seqs 5 -num_threads 16 -outfmt 7 



quast -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/quast/MNH06145 /work/assembly/current/MNH/MNH061/MNH06145/genomic.fna
checkm lineage_wf /work/assembly/current/MNH/MNH061/MNH06145 /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH06145 -x fna -t 20
python /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/run_BUSCO.py -i /work/assembly/current/MNH/MNH061/MNH06145/genomic.fna -c 16 -o MNH06145 -m geno -f -l /work/workspace/zhurj/reference/BUSCO/bacteria_odb9

/work/workspace/zhurj/reference/serverdb

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/config/config.ini
wget http://busco.ezlab.org/v2/datasets/bacteria_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/proteobacteria_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/rhizobiales_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/betaproteobacteria_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/gammaproteobacteria_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/enterobacteriales_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/deltaepsilonsub_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/actinobacteria_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/cyanobacteria_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/firmicutes_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/clostridia_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/lactobacillales_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/bacillales_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/bacteroidetes_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/spirochaetes_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/tenericutes_odb9.tar.gz

busco: only fit for several bacteria, now don't use this software for data analysis

metaxa2 -i /work/assembly/current/MNH/MNH061/MNH06145/genomic.fna -o MNH06145 --cpu 16 -t b -d /work/program/current/ncbi-blast/bin -p /work/workspace/zhurj/instal/hmmer/bin


mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH05121
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH05121/quast
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH05121/checkm
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH05121/quast /work/assembly/current/MNH/MNH051/MNH05121/genomic.fna 
checkm lineage_wf /work/assembly/current/MNH/MNH051/MNH05121 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH05121/checkm -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH051/MNH05121/genomic.fna /work/assembly/current/MNH/MNH051/MNH05121 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH05121/checkm -x fna -t 20
fastANI -q /work/assembly/current/MNH/MNH061/MNH06145/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000007/GCA_000007525/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH06145_GCA_000007525


mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MN12371
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MN12371/quast
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MN12371/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MN12371/blast
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MN12371/quast /work/assembly/current/MN/MN123/MN12371/genomic.fna 
checkm lineage_wf /work/assembly/current/MN/MN123/MN12371 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MN12371/checkm -x fna -t 20
checkm ssu_finder /work/assembly/current/MN/MN123/MN12371/genomic.fna /work/assembly/current/MN/MN123/MN12371 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MN12371/checkm -x fna -t 20
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MN12371/checkm/ssufilter.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MN12371/blast/MN12371.out -max_target_seqs 5 -num_threads 16 -outfmt 7 

mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MN111097
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MN111097/quast
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MN111097/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MN111097/blast
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MN111097/quast /work/assembly/current/MN/MN111/MN111097/genomic.fna 
checkm lineage_wf /work/assembly/current/MN/MN111/MN111097 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MN111097/checkm -x fna -t 20
checkm ssu_finder /work/assembly/current/MN/MN111/MN111097/genomic.fna /work/assembly/current/MN/MN111/MN111097 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MN111097/checkm -x fna -t 20
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MN111097/checkm/ssufilter.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MN111097/blast/MN111097.out -max_target_seqs 5 -num_threads 16 -outfmt 7 

mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH22053  Enterococcus faecium
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH22590  Enterococcus faecium
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH23386  Enterococcus faecium
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH21691  Enterococcus faecium
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH09296  Enterococcus faecium
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH23818  Enterococcus faecium
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH21870  Enterococcus faecium
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH09666  Enterococcus faecium
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH08896  Enterococcus faecium
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH23933  Enterococcus faecium
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH21491  Enterococcus faecium
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH21866  Enterococcus faecium
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH09782  Enterococcus faecium
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH10333  Enterococcus faecium
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH22871  Enterococcus faecium
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH22686  Enterococcus faecium
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH23566  Enterococcus faecium
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH26248  Enterococcus gallinarum
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH22464  Enterococcus gallinarum
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH17373  Enterococcus gallinarum
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH22599  Enterococcus faecalis
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH21146  Enterococcus hirae
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH10864  Enterococcus hirae
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH26765  Enterococcus hirae
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH10944  Enterococcus hirae
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH11355  Enterococcus faecium
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH10908  Enterococcus hirae
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH25574  Enterococcus hirae
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH25464  Enterococcus hirae
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH25945  Enterococcus hirae
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH10846  Enterococcus faecium
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH11304  Enterococcus faecium
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH25145  Enterococcus hirae
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH25848  Enterococcus hirae
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH27666  Enterococcus faecium
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH11341  Enterococcus hirae



mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH22053/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH22590/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH23386/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH21691/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH09296/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH23818/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH21870/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH09666/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH08896/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH23933/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH21491/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH21866/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH09782/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH10333/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH22871/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH22686/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH23566/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH26248/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH22464/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH17373/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH22599/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH21146/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH10864/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH26765/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH10944/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH11355/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH10908/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH25574/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH25464/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH25945/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH10846/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH11304/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH25145/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH25848/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH27666/checkm
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH11341/checkm


quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH22053/quast /work/assembly/current/MNH/MNH220/MNH22053/genomic.fna 
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH22590/quast /work/assembly/current/MNH/MNH225/MNH22590/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH23386/quast /work/assembly/current/MNH/MNH233/MNH23386/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH21691/quast /work/assembly/current/MNH/MNH216/MNH21691/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH09296/quast /work/assembly/current/MNH/MNH092/MNH09296/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH23818/quast /work/assembly/current/MNH/MNH238/MNH23818/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH21870/quast /work/assembly/current/MNH/MNH218/MNH21870/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH09666/quast /work/assembly/current/MNH/MNH096/MNH09666/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH08896/quast /work/assembly/current/MNH/MNH088/MNH08896/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH23933/quast /work/assembly/current/MNH/MNH239/MNH23933/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH21491/quast /work/assembly/current/MNH/MNH214/MNH21491/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH21866/quast /work/assembly/current/MNH/MNH218/MNH21866/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH09782/quast /work/assembly/current/MNH/MNH097/MNH09782/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH10333/quast /work/assembly/current/MNH/MNH103/MNH10333/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH22871/quast /work/assembly/current/MNH/MNH228/MNH22871/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH22686/quast /work/assembly/current/MNH/MNH226/MNH22686/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH23566/quast /work/assembly/current/MNH/MNH235/MNH23566/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH26248/quast /work/assembly/current/MNH/MNH262/MNH26248/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH22464/quast /work/assembly/current/MNH/MNH224/MNH22464/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH17373/quast /work/assembly/current/MNH/MNH173/MNH17373/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH22599/quast /work/assembly/current/MNH/MNH225/MNH22599/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH21146/quast /work/assembly/current/MNH/MNH211/MNH21146/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH10864/quast /work/assembly/current/MNH/MNH108/MNH10864/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH26765/quast /work/assembly/current/MNH/MNH267/MNH26765/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH10944/quast /work/assembly/current/MNH/MNH109/MNH10944/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH11355/quast /work/assembly/current/MNH/MNH113/MNH11355/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH10908/quast /work/assembly/current/MNH/MNH109/MNH10908/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH25574/quast /work/assembly/current/MNH/MNH255/MNH25574/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH25464/quast /work/assembly/current/MNH/MNH254/MNH25464/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH25945/quast /work/assembly/current/MNH/MNH259/MNH25945/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH10846/quast /work/assembly/current/MNH/MNH108/MNH10846/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH11304/quast /work/assembly/current/MNH/MNH113/MNH11304/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH25145/quast /work/assembly/current/MNH/MNH251/MNH25145/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH25848/quast /work/assembly/current/MNH/MNH258/MNH25848/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH27666/quast /work/assembly/current/MNH/MNH276/MNH27666/genomic.fna
quast -t 16 -o  /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH11341/quast /work/assembly/current/MNH/MNH113/MNH11341/genomic.fna

checkm ssu_finder /work/assembly/current/MNH/MNH220/MNH22053/genomic.fna /work/assembly/current/MNH/MNH220/MNH22053 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH22053/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH225/MNH22590/genomic.fna /work/assembly/current/MNH/MNH225/MNH22590 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH22590/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH233/MNH23386/genomic.fna /work/assembly/current/MNH/MNH233/MNH23386 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH23386/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH216/MNH21691/genomic.fna /work/assembly/current/MNH/MNH216/MNH21691 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH21691/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH092/MNH09296/genomic.fna /work/assembly/current/MNH/MNH092/MNH09296 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH09296/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH238/MNH23818/genomic.fna /work/assembly/current/MNH/MNH238/MNH23818 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH23818/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH218/MNH21870/genomic.fna /work/assembly/current/MNH/MNH218/MNH21870 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH21870/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH096/MNH09666/genomic.fna /work/assembly/current/MNH/MNH096/MNH09666 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH09666/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH088/MNH08896/genomic.fna /work/assembly/current/MNH/MNH088/MNH08896 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH08896/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH239/MNH23933/genomic.fna /work/assembly/current/MNH/MNH239/MNH23933 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH23933/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH214/MNH21491/genomic.fna /work/assembly/current/MNH/MNH214/MNH21491 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH21491/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH218/MNH21866/genomic.fna /work/assembly/current/MNH/MNH218/MNH21866 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH21866/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH097/MNH09782/genomic.fna /work/assembly/current/MNH/MNH097/MNH09782 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH09782/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH103/MNH10333/genomic.fna /work/assembly/current/MNH/MNH103/MNH10333 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH10333/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH228/MNH22871/genomic.fna /work/assembly/current/MNH/MNH228/MNH22871 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH22871/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH226/MNH22686/genomic.fna /work/assembly/current/MNH/MNH226/MNH22686 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH22686/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH235/MNH23566/genomic.fna /work/assembly/current/MNH/MNH235/MNH23566 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH23566/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH262/MNH26248/genomic.fna /work/assembly/current/MNH/MNH262/MNH26248 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH26248/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH224/MNH22464/genomic.fna /work/assembly/current/MNH/MNH224/MNH22464 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH22464/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH173/MNH17373/genomic.fna /work/assembly/current/MNH/MNH173/MNH17373 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH17373/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH225/MNH22599/genomic.fna /work/assembly/current/MNH/MNH225/MNH22599 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH22599/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH211/MNH21146/genomic.fna /work/assembly/current/MNH/MNH211/MNH21146 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH21146/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH108/MNH10864/genomic.fna /work/assembly/current/MNH/MNH108/MNH10864 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH10864/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH267/MNH26765/genomic.fna /work/assembly/current/MNH/MNH267/MNH26765 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH26765/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH109/MNH10944/genomic.fna /work/assembly/current/MNH/MNH109/MNH10944 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH10944/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH113/MNH11355/genomic.fna /work/assembly/current/MNH/MNH113/MNH11355 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH11355/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH109/MNH10908/genomic.fna /work/assembly/current/MNH/MNH109/MNH10908 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH10908/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH255/MNH25574/genomic.fna /work/assembly/current/MNH/MNH255/MNH25574 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH25574/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH254/MNH25464/genomic.fna /work/assembly/current/MNH/MNH254/MNH25464 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH25464/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH259/MNH25945/genomic.fna /work/assembly/current/MNH/MNH259/MNH25945 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH25945/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH108/MNH10846/genomic.fna /work/assembly/current/MNH/MNH108/MNH10846 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH10846/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH113/MNH11304/genomic.fna /work/assembly/current/MNH/MNH113/MNH11304 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH11304/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH251/MNH25145/genomic.fna /work/assembly/current/MNH/MNH251/MNH25145 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH25145/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH258/MNH25848/genomic.fna /work/assembly/current/MNH/MNH258/MNH25848 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH25848/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH276/MNH27666/genomic.fna /work/assembly/current/MNH/MNH276/MNH27666 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH27666/checkm  -x fna -t 20
checkm ssu_finder /work/assembly/current/MNH/MNH113/MNH11341/genomic.fna /work/assembly/current/MNH/MNH113/MNH11341 /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MNH11341/checkm  -x fna -t 20


/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/taxcheck/sample/MN12371/checkm/ssufilter.fna -task blastn -db /work/database/ezbio/16S/current/

perl /work/workspace/zhurj/script/git/phytreesnp/genoPhyana.pl -i /work/workspace/zhurj/project/5_tools/phytreesnp/input/group.txt -o /work/workspace/zhurj/project/5_tools/phytreesnp/output -d /work/workspace/zhurj/project/5_tools/phytreesnp/fna --ani

perl genoPhyana_nodel.pl -i /work/workspace/ruijuan/script/2_MWGS/kSNPtest/groupTest/input/group.txt -o /work/workspace/ruijuan/project/2_mwgs/run00020-25/result -d /work/workspace/ruijuan/project/2_mwgs/run00020-25/result/finalAssemble --ani

With a working anaconda installation, install the bioconda and conda-forge channels:
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forgex

my $ksnp3pro = '/work/workspace/zhurj/software/kSNP3.1_Linux_package/kSNP3/kSNP3';
my $figtreepro = 'java -Xmx10g  -jar /work/workspace/zhurj/software/FigTree_v1.4.4/lib/figtree.jar';
my $anipro = "/work/workspace/zhurj/software/pyani-master/bin/average_nucleotide_identity.py";

ln -s /work/assembly/current/MNH/MNH175/MNH17503/genomic.fna /work/workspace/zhurj/project/5_tools/phytreesnp/fna/MNH17503.fna
ln -s /work/assembly/current/MNH/MNH175/MNH17505/genomic.fna /work/workspace/zhurj/project/5_tools/phytreesnp/fna/MNH17505.fna
ln -s /work/assembly/current/MNH/MNH192/MNH19250/genomic.fna /work/workspace/zhurj/project/5_tools/phytreesnp/fna/MNH19250.fna
ln -s /work/assembly/current/MNH/MNH196/MNH19631/genomic.fna /work/workspace/zhurj/project/5_tools/phytreesnp/fna/MNH19631.fna
ln -s /work/assembly/current/MNH/MNH197/MNH19714/genomic.fna /work/workspace/zhurj/project/5_tools/phytreesnp/fna/MNH19714.fna

download human genome: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
cd /work/workspace/zhurj/reference/human
bwa index -a bwtsw hg38.fa.gz

bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH31610.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH31610.clean.2.fq.gz -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31610/MNH31610.sam
samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31610/MNH31610.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31610/MNH31610.sam
samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31610/MNH31610.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH31610.txt


mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31671
bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH31671.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH31671.clean.2.fq.gz -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31671/MNH31671.sam
samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31671/MNH31671.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31671/MNH31671.sam
samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31671/MNH31671.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH31671.txt

mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33080
bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33080.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33080.clean.2.fq.gz -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33080/MNH33080.sam
samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33080/MNH33080.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33080/MNH33080.sam
samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33080/MNH33080.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33080.txt

mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33156
bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33156.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33156.clean.2.fq.gz -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33156/MNH33156.sam
samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33156/MNH33156.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33156/MNH33156.sam
samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33156/MNH33156.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33156.txt

mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31671 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH31671.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH31671.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31671/MNH31671.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31671/MNH31671.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31671/MNH31671.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31671/MNH31671.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH31671.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31610 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH31610.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH31610.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31610/MNH31610.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31610/MNH31610.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31610/MNH31610.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31610/MNH31610.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH31610.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33080 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33080.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33080.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33080/MNH33080.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33080/MNH33080.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33080/MNH33080.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33080/MNH33080.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33080.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33156 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33156.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33156.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33156/MNH33156.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33156/MNH33156.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33156/MNH33156.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33156/MNH33156.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33156.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH16402 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH16402.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH16402.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH16402/MNH16402.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH16402/MNH16402.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH16402/MNH16402.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH16402/MNH16402.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH16402.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31703 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH31703.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH31703.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31703/MNH31703.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31703/MNH31703.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31703/MNH31703.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31703/MNH31703.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH31703.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33028 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33028.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33028.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33028/MNH33028.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33028/MNH33028.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33028/MNH33028.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33028/MNH33028.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33028.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33296 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33296.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33296.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33296/MNH33296.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33296/MNH33296.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33296/MNH33296.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33296/MNH33296.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33296.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33386 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33386.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33386.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33386/MNH33386.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33386/MNH33386.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33386/MNH33386.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33386/MNH33386.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33386.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33413 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33413.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33413.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33413/MNH33413.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33413/MNH33413.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33413/MNH33413.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33413/MNH33413.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33413.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36202 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36202.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36202.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36202/MNH36202.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36202/MNH36202.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36202/MNH36202.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36202/MNH36202.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36202.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36208 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36208.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36208.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36208/MNH36208.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36208/MNH36208.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36208/MNH36208.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36208/MNH36208.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36208.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36228 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36228.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36228.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36228/MNH36228.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36228/MNH36228.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36228/MNH36228.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36228/MNH36228.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36228.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36384 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36384.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36384.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36384/MNH36384.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36384/MNH36384.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36384/MNH36384.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36384/MNH36384.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36384.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36240 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36240.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36240.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36240/MNH36240.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36240/MNH36240.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36240/MNH36240.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36240/MNH36240.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36240.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH27992 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH27992.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH27992.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH27992/MNH27992.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH27992/MNH27992.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH27992/MNH27992.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH27992/MNH27992.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH27992.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36362 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36362.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36362.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36362/MNH36362.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36362/MNH36362.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36362/MNH36362.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36362/MNH36362.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36362.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH20648 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH20648.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH20648.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH20648/MNH20648.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH20648/MNH20648.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH20648/MNH20648.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH20648/MNH20648.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH20648.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31704 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH31704.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH31704.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31704/MNH31704.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31704/MNH31704.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31704/MNH31704.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31704/MNH31704.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH31704.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33204 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33204.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33204.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33204/MNH33204.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33204/MNH33204.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33204/MNH33204.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33204/MNH33204.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33204.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33253 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33253.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33253.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33253/MNH33253.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33253/MNH33253.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33253/MNH33253.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33253/MNH33253.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33253.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33317 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33317.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33317.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33317/MNH33317.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33317/MNH33317.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33317/MNH33317.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33317/MNH33317.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33317.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33363 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33363.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33363.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33363/MNH33363.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33363/MNH33363.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33363/MNH33363.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33363/MNH33363.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33363.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33453 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33453.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33453.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33453/MNH33453.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33453/MNH33453.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33453/MNH33453.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33453/MNH33453.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33453.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36128 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36128.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36128.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36128/MNH36128.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36128/MNH36128.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36128/MNH36128.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36128/MNH36128.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36128.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36212 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36212.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36212.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36212/MNH36212.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36212/MNH36212.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36212/MNH36212.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36212/MNH36212.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36212.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36387 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36387.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36387.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36387/MNH36387.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36387/MNH36387.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36387/MNH36387.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36387/MNH36387.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36387.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36395 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36395.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36395.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36395/MNH36395.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36395/MNH36395.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36395/MNH36395.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36395/MNH36395.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36395.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36021 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36021.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36021.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36021/MNH36021.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36021/MNH36021.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36021/MNH36021.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36021/MNH36021.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36021.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31683 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH31683.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH31683.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31683/MNH31683.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31683/MNH31683.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31683/MNH31683.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31683/MNH31683.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH31683.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33041 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33041.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33041.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33041/MNH33041.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33041/MNH33041.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33041/MNH33041.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33041/MNH33041.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33041.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36393 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36393.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36393.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36393/MNH36393.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36393/MNH36393.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36393/MNH36393.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36393/MNH36393.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36393.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33040 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33040.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33040.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33040/MNH33040.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33040/MNH33040.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33040/MNH33040.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33040/MNH33040.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33040.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33078 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33078.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33078.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33078/MNH33078.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33078/MNH33078.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33078/MNH33078.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33078/MNH33078.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33078.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33106 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33106.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33106.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33106/MNH33106.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33106/MNH33106.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33106/MNH33106.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33106/MNH33106.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33106.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33194 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33194.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33194.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33194/MNH33194.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33194/MNH33194.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33194/MNH33194.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33194/MNH33194.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33194.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33320 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33320.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33320.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33320/MNH33320.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33320/MNH33320.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33320/MNH33320.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33320/MNH33320.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33320.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36095 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36095.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36095.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36095/MNH36095.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36095/MNH36095.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36095/MNH36095.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36095/MNH36095.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36095.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36265 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36265.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36265.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36265/MNH36265.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36265/MNH36265.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36265/MNH36265.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36265/MNH36265.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36265.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33014 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33014.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33014.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33014/MNH33014.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33014/MNH33014.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33014/MNH33014.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33014/MNH33014.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33014.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33403 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33403.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33403.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33403/MNH33403.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33403/MNH33403.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33403/MNH33403.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33403/MNH33403.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33403.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36188 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36188.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36188.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36188/MNH36188.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36188/MNH36188.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36188/MNH36188.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36188/MNH36188.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36188.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36255 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36255.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36255.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36255/MNH36255.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36255/MNH36255.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36255/MNH36255.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36255/MNH36255.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36255.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH14164 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH14164.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH14164.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH14164/MNH14164.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH14164/MNH14164.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH14164/MNH14164.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH14164/MNH14164.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH14164.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH17919 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH17919.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH17919.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH17919/MNH17919.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH17919/MNH17919.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH17919/MNH17919.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH17919/MNH17919.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH17919.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33016 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33016.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33016.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33016/MNH33016.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33016/MNH33016.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33016/MNH33016.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33016/MNH33016.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33016.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33070 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33070.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33070.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33070/MNH33070.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33070/MNH33070.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33070/MNH33070.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33070/MNH33070.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33070.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33095 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33095.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33095.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33095/MNH33095.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33095/MNH33095.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33095/MNH33095.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33095/MNH33095.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33095.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33341 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33341.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33341.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33341/MNH33341.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33341/MNH33341.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33341/MNH33341.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33341/MNH33341.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33341.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH19749 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH19749.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH19749.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH19749/MNH19749.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH19749/MNH19749.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH19749/MNH19749.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH19749/MNH19749.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH19749.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33050 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33050.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33050.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33050/MNH33050.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33050/MNH33050.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33050/MNH33050.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33050/MNH33050.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33050.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33076 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33076.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33076.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33076/MNH33076.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33076/MNH33076.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33076/MNH33076.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33076/MNH33076.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33076.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33355 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33355.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33355.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33355/MNH33355.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33355/MNH33355.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33355/MNH33355.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33355/MNH33355.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33355.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31644 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH31644.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH31644.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31644/MNH31644.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31644/MNH31644.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31644/MNH31644.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31644/MNH31644.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH31644.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31681 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH31681.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH31681.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31681/MNH31681.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31681/MNH31681.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31681/MNH31681.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH31681/MNH31681.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH31681.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33149 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33149.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33149.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33149/MNH33149.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33149/MNH33149.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33149/MNH33149.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33149/MNH33149.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33149.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33229 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33229.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33229.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33229/MNH33229.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33229/MNH33229.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33229/MNH33229.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33229/MNH33229.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33229.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33293 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33293.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33293.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33293/MNH33293.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33293/MNH33293.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33293/MNH33293.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33293/MNH33293.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33293.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33313 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33313.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33313.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33313/MNH33313.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33313/MNH33313.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33313/MNH33313.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33313/MNH33313.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33313.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33404 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33404.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33404.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33404/MNH33404.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33404/MNH33404.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33404/MNH33404.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33404/MNH33404.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33404.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33454 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33454.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33454.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33454/MNH33454.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33454/MNH33454.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33454/MNH33454.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33454/MNH33454.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33454.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36002 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36002.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36002.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36002/MNH36002.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36002/MNH36002.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36002/MNH36002.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36002/MNH36002.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36002.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36025 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36025.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36025.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36025/MNH36025.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36025/MNH36025.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36025/MNH36025.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36025/MNH36025.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36025.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36070 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36070.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36070.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36070/MNH36070.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36070/MNH36070.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36070/MNH36070.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36070/MNH36070.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36070.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36211 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36211.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36211.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36211/MNH36211.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36211/MNH36211.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36211/MNH36211.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36211/MNH36211.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36211.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36311 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36311.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36311.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36311/MNH36311.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36311/MNH36311.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36311/MNH36311.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36311/MNH36311.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36311.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36323 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36323.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36323.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36323/MNH36323.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36323/MNH36323.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36323/MNH36323.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36323/MNH36323.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36323.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36345 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36345.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36345.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36345/MNH36345.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36345/MNH36345.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36345/MNH36345.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36345/MNH36345.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36345.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36401 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36401.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36401.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36401/MNH36401.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36401/MNH36401.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36401/MNH36401.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36401/MNH36401.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36401.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH08404 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH08404.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH08404.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH08404/MNH08404.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH08404/MNH08404.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH08404/MNH08404.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH08404/MNH08404.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH08404.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH09154 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH09154.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH09154.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH09154/MNH09154.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH09154/MNH09154.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH09154/MNH09154.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH09154/MNH09154.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH09154.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36148 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36148.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36148.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36148/MNH36148.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36148/MNH36148.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36148/MNH36148.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36148/MNH36148.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36148.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36440 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36440.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36440.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36440/MNH36440.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36440/MNH36440.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36440/MNH36440.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36440/MNH36440.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36440.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33015 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33015.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33015.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33015/MNH33015.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33015/MNH33015.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33015/MNH33015.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33015/MNH33015.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33015.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33162 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33162.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33162.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33162/MNH33162.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33162/MNH33162.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33162/MNH33162.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33162/MNH33162.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33162.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33199 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33199.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33199.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33199/MNH33199.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33199/MNH33199.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33199/MNH33199.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33199/MNH33199.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33199.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33257 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33257.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33257.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33257/MNH33257.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33257/MNH33257.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33257/MNH33257.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33257/MNH33257.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33257.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33382 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33382.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33382.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33382/MNH33382.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33382/MNH33382.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33382/MNH33382.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33382/MNH33382.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33382.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36181 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36181.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36181.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36181/MNH36181.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36181/MNH36181.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36181/MNH36181.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36181/MNH36181.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36181.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36205 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36205.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36205.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36205/MNH36205.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36205/MNH36205.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36205/MNH36205.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36205/MNH36205.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36205.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33371 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33371.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH33371.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33371/MNH33371.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33371/MNH33371.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33371/MNH33371.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH33371/MNH33371.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH33371.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36049 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36049.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36049.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36049/MNH36049.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36049/MNH36049.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36049/MNH36049.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36049/MNH36049.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36049.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36299 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36299.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36299.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36299/MNH36299.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36299/MNH36299.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36299/MNH36299.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36299/MNH36299.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36299.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36333 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36333.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36333.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36333/MNH36333.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36333/MNH36333.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36333/MNH36333.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36333/MNH36333.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36333.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36067 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36067.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36067.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36067/MNH36067.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36067/MNH36067.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36067/MNH36067.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36067/MNH36067.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36067.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36171 && bwa mem -t 16 /work/workspace/zhurj/reference/human/hg38.fa.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36171.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/taxcheck/output/clean/MNH36171.clean.2.fq.gz  -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36171/MNH36171.sam && samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36171/MNH36171.sorted.sam /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36171/MNH36171.sam && samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/sam/MNH36171/MNH36171.sorted.sam > /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat/MNH36171.txt


cd /work/workspace/zhurj/project/2_swgs/taxcheck/output/human/stat
if [ -f "summary.txt" ]; then
  rm summary.txt
fi
cat MNH31671.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH31671\t"$2}' >> summary.txt
cat MNH31610.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH31610\t"$2}' >> summary.txt
cat MNH33080.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33080\t"$2}' >> summary.txt
cat MNH33156.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33156\t"$2}' >> summary.txt
cat MNH16402.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH16402\t"$2}' >> summary.txt
cat MNH31703.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH31703\t"$2}' >> summary.txt
cat MNH33028.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33028\t"$2}' >> summary.txt
cat MNH33296.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33296\t"$2}' >> summary.txt
cat MNH33386.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33386\t"$2}' >> summary.txt
cat MNH33413.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33413\t"$2}' >> summary.txt
cat MNH36202.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36202\t"$2}' >> summary.txt
cat MNH36208.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36208\t"$2}' >> summary.txt
cat MNH36228.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36228\t"$2}' >> summary.txt
cat MNH36384.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36384\t"$2}' >> summary.txt
cat MNH36240.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36240\t"$2}' >> summary.txt
cat MNH27992.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH27992\t"$2}' >> summary.txt
cat MNH36362.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36362\t"$2}' >> summary.txt
cat MNH20648.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH20648\t"$2}' >> summary.txt
cat MNH31704.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH31704\t"$2}' >> summary.txt
cat MNH33204.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33204\t"$2}' >> summary.txt
cat MNH33253.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33253\t"$2}' >> summary.txt
cat MNH33317.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33317\t"$2}' >> summary.txt
cat MNH33363.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33363\t"$2}' >> summary.txt
cat MNH33453.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33453\t"$2}' >> summary.txt
cat MNH36128.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36128\t"$2}' >> summary.txt
cat MNH36212.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36212\t"$2}' >> summary.txt
cat MNH36387.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36387\t"$2}' >> summary.txt
cat MNH36395.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36395\t"$2}' >> summary.txt
cat MNH36021.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36021\t"$2}' >> summary.txt
cat MNH31683.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH31683\t"$2}' >> summary.txt
cat MNH33041.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33041\t"$2}' >> summary.txt
cat MNH36393.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36393\t"$2}' >> summary.txt
cat MNH33040.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33040\t"$2}' >> summary.txt
cat MNH33078.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33078\t"$2}' >> summary.txt
cat MNH33106.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33106\t"$2}' >> summary.txt
cat MNH33194.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33194\t"$2}' >> summary.txt
cat MNH33320.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33320\t"$2}' >> summary.txt
cat MNH36095.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36095\t"$2}' >> summary.txt
cat MNH36265.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36265\t"$2}' >> summary.txt
cat MNH33014.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33014\t"$2}' >> summary.txt
cat MNH33403.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33403\t"$2}' >> summary.txt
cat MNH36188.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36188\t"$2}' >> summary.txt
cat MNH36255.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36255\t"$2}' >> summary.txt
cat MNH14164.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH14164\t"$2}' >> summary.txt
cat MNH17919.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH17919\t"$2}' >> summary.txt
cat MNH33016.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33016\t"$2}' >> summary.txt
cat MNH33070.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33070\t"$2}' >> summary.txt
cat MNH33095.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33095\t"$2}' >> summary.txt
cat MNH33341.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33341\t"$2}' >> summary.txt
cat MNH19749.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH19749\t"$2}' >> summary.txt
cat MNH33050.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33050\t"$2}' >> summary.txt
cat MNH33076.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33076\t"$2}' >> summary.txt
cat MNH33355.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33355\t"$2}' >> summary.txt
cat MNH31644.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH31644\t"$2}' >> summary.txt
cat MNH31681.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH31681\t"$2}' >> summary.txt
cat MNH33149.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33149\t"$2}' >> summary.txt
cat MNH33229.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33229\t"$2}' >> summary.txt
cat MNH33293.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33293\t"$2}' >> summary.txt
cat MNH33313.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33313\t"$2}' >> summary.txt
cat MNH33404.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33404\t"$2}' >> summary.txt
cat MNH33454.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33454\t"$2}' >> summary.txt
cat MNH36002.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36002\t"$2}' >> summary.txt
cat MNH36025.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36025\t"$2}' >> summary.txt
cat MNH36070.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36070\t"$2}' >> summary.txt
cat MNH36211.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36211\t"$2}' >> summary.txt
cat MNH36311.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36311\t"$2}' >> summary.txt
cat MNH36323.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36323\t"$2}' >> summary.txt
cat MNH36345.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36345\t"$2}' >> summary.txt
cat MNH36401.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36401\t"$2}' >> summary.txt
cat MNH08404.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH08404\t"$2}' >> summary.txt
cat MNH09154.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH09154\t"$2}' >> summary.txt
cat MNH36148.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36148\t"$2}' >> summary.txt
cat MNH36440.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36440\t"$2}' >> summary.txt
cat MNH33015.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33015\t"$2}' >> summary.txt
cat MNH33162.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33162\t"$2}' >> summary.txt
cat MNH33199.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33199\t"$2}' >> summary.txt
cat MNH33257.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33257\t"$2}' >> summary.txt
cat MNH33382.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33382\t"$2}' >> summary.txt
cat MNH36181.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36181\t"$2}' >> summary.txt
cat MNH36205.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36205\t"$2}' >> summary.txt
cat MNH33371.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH33371\t"$2}' >> summary.txt
cat MNH36049.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36049\t"$2}' >> summary.txt
cat MNH36299.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36299\t"$2}' >> summary.txt
cat MNH36333.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36333\t"$2}' >> summary.txt
cat MNH36067.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36067\t"$2}' >> summary.txt
cat MNH36171.txt | grep "mapped (" | awk -F "[(:]" '{print "MNH36171\t"$2}' >> summary.txt


nohup perl /work/workspace/zhurj/script/git/genoassemble/genoassemble.pl  -i /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/clean -o /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/output -t 16 > /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/log/sample4_20200115.txt 2>&1 &

mkdir -p /work/workspace/zhurj/reference/NCBI/genome/GCA_001640865/bwa && ln -s /work/database/ncbi/genome/all/GCA/GCA_001/GCA_001640/GCA_001640865/genomic.fna.gz /work/workspace/zhurj/reference/NCBI/genome/GCA_001640865/bwa/GCA_001640865.fna.gz && bwa index /work/workspace/zhurj/reference/NCBI/genome/GCA_001640865/bwa/GCA_001640865.fna.gz 


mkdir -p /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38555/bwa
bwa mem -t 16 /work/workspace/zhurj/reference/NCBI/genome/GCA_001640865/bwa/GCA_001640865.fna.gz /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/clean/MNH38555.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/clean/MNH38555.clean.2.fq.gz -o /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38555/bwa/MNH38555.sam

samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38555/bwa/MNH38555.sorted.sam /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38555/bwa/MNH38555.sam
samtools flagstat --threads 16 //work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38555/bwa/MNH38555.sorted.sam > /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38555/bwa/MNH38555.stat.sam

mkdir /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38555/quast
mkdir /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38555/checkm
mkdir /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38555/fastANI
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38555/quast /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/output/MNH38555/MNH38555.fna

checkm ssu_finder /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/output/MNH38555/MNH38555.fna /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/output/MNH38555 /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38555/checkm -x fna -t 20 

checkm lineage_wf /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/output/MNH38555 /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38555/checkm -x fna -t 20 --tab_table -f /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38555/checkm/MNH38555_checkm.txt

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38555/checkm/ssu_filter.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38555/checkm/MNH38555.out -max_target_seqs 5 -num_threads 16 -outfmt 7 

MNH38555&&NODE_60_length_5248_cov_2286.874492   ABWZ01000093    100.000 1454    0       0       33      1486    1       1454    0.0     2623
MNH38555&&NODE_60_length_5248_cov_2286.874492   CP000139        97.111  1454    38      3       33      1486    1       1450    0.0     2424
MNH38555&&NODE_60_length_5248_cov_2286.874492   PAC001204       96.217  1454    53      1       33      1486    1       1452    0.0     2372
MNH38555&&NODE_60_length_5248_cov_2286.874492   HQ778885        96.217  1454    54      1       34      1486    1       1454    0.0     2371
MNH38555&&NODE_60_length_5248_cov_2286.874492   HQ767830        96.583  1434    47      2       55      1486    1       1434    0.0     2359

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "ABWZ01000093|CP000139|PAC001204|HQ778885|HQ767830"

ABWZ01000093    Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;Bacteroides dorei;DSM 17855(T)
CP000139        Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;Bacteroides vulgatus;ATCC 8482(T)
HQ767830        Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;ELU0062-T425-S-NIPCRAMgANa_000105
HQ778885        Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;ELU0086-T395-S-NI_000069
PAC001204       Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;None

fastANI -q -r /work/workspace/zhurj/reference/NCBI/genome/GCA_001640865/bwa/GCA_001640865.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38555/fastANI/MNH38555_bacreroides_dorei


mkdir -p /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38946/bwa
bwa mem -t 16 /work/workspace/zhurj/reference/NCBI/genome/GCA_001640865/bwa/GCA_001640865.fna.gz /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/clean/MNH38946.clean.1.fq.gz /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/clean/MNH38946.clean.2.fq.gz -o /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38946/bwa/MNH38946.sam

samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38946/bwa/MNH38946.sorted.sam /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38946/bwa/MNH38946.sam
samtools flagstat --threads 16 //work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38946/bwa/MNH38946.sorted.sam > /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38946/bwa/MNH38946.stat.sam

mkdir /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38946/quast
mkdir /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38946/checkm
mkdir /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38946/fastANI
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38946/quast /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/output/MNH38946/MNH38946.fna

checkm ssu_finder /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/output/MNH38946/MNH38946.fna /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/output/MNH38946 /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38946/checkm -x fna -t 20 

checkm lineage_wf /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/output/MNH38946 /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38946/checkm -x fna -t 20 --tab_table -f /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38946/checkm/MNH38946_checkm.txt

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38946/checkm/ssu.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38946/checkm/MNH38946.out -max_target_seqs 5 -num_threads 16 -outfmt 7 

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "ABWZ01000093|CP000139|PAC001204|HQ778885|HQ767830"
ABWZ01000093    Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;Bacteroides dorei;DSM 17855(T)
CP000139        Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;Bacteroides vulgatus;ATCC 8482(T)
HQ767830        Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;ELU0062-T425-S-NIPCRAMgANa_000105
HQ778885        Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;ELU0086-T395-S-NI_000069
PAC001204       Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;None

MNH38946&&NODE_40_length_5358_cov_1925.246734   ABWZ01000093    100.000 1454    0       0       34      1487    1       1454    0.0     2623
MNH38946&&NODE_40_length_5358_cov_1925.246734   CP000139        97.111  1454    38      3       34      1487    1       1450    0.0     2424
MNH38946&&NODE_40_length_5358_cov_1925.246734   PAC001204       96.217  1454    53      1       34      1487    1       1452    0.0     2372
MNH38946&&NODE_40_length_5358_cov_1925.246734   HQ778885        96.217  1454    54      1       35      1487    1       1454    0.0     2371
MNH38946&&NODE_40_length_5358_cov_1925.246734   HQ767830        96.583  1434    47      2       56      1487    1       1434    0.0     2359


fastANI -q /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/output/MNH38946/MNH38946.fna -r /work/workspace/zhurj/reference/NCBI/genome/GCA_001640865/bwa/GCA_001640865.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38946/fastANI/MNH38946_bacreroides_dorei
fastANI -q /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/output/MNH38946/MNH38946.fna -r /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/output/MNH38555/MNH38555.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38946/fastANI/MNH38555_MNH38946

diff /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38555/checkm/ssu_filter.fna  /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38946/checkm/ssu.fna

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38555/checkm/ssu_filter.fna -subject /work/workspace/zhurj/project/2_swgs/run00046/sample4_20200115/tax/MNH38946/checkm/ssu.fna

-----------------------------------------
mkdir -p /work/workspace/zhurj/reference/NCBI/genome/GCA_002874775/bwa && ln -s /work/database/ncbi/genome/all/GCA/GCA_002/GCA_002874/GCA_002874775/genomic.fna.gz /work/workspace/zhurj/reference/NCBI/genome/GCA_002874775/bwa/GCA_002874775.fna.gz && bwa index /work/workspace/zhurj/reference/NCBI/genome/GCA_002874775/bwa/GCA_002874775.fna.gz 

mkdir -p /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/bwa
bwa mem -t 16 /work/workspace/zhurj/reference/NCBI/genome/GCA_002874775/bwa/GCA_002874775.fna.gz /work/rawdata/fastq/genome/MNH/MNH193/MNH19366/E1/L1/S1/MNH19366.1.fq.gz /work/rawdata/fastq/genome/MNH/MNH193/MNH19366/E1/L1/S1/MNH19366.2.fq.gz -o /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/bwa/MNH19366.sam

samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/bwa/MNH19366.sorted.sam /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/bwa/MNH19366.sam
samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/bwa/MNH19366.sorted.sam > /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/bwa/MNH19366.stat

mkdir -p /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/quast
mkdir -p /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/checkm
mkdir -p /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/fastANI
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/quast /work/workspace/zhurj/project/2_swgs/run00020/sample/output/MNH19366/MNH19366.fna
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/GCA_002874775 /work/workspace/zhurj/reference/NCBI/genome/GCA_002874775/bwa/GCA_002874775.fna.gz

checkm ssu_finder /work/workspace/zhurj/project/2_swgs/run00020/sample/output/MNH19366/MNH19366.fna /work/workspace/zhurj/project/2_swgs/run00020/sample/output/MNH19366 /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/checkm -x fna -t 20 

checkm lineage_wf /work/workspace/zhurj/project/2_swgs/run00020/sample/output/MNH19366 /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/checkm -x fna -t 20 --tab_table -f /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/checkm/MNH19366_checkm.txt

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/checkm/ssu_filter.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/checkm/MNH19366.out -max_target_seqs 5 -num_threads 16 -outfmt 7 

MNH19366&&MNH19366_33   GQ451293        96.998  1266    31      4       55      1317    1457    196     0.0     2101
MNH19366&&MNH19366_33   PAC000672       96.840  1266    33      6       55      1317    1457    196     0.0     2083
MNH19366&&MNH19366_33   PAC003100       96.443  1265    41      4       55      1317    1458    196     0.0     2065
MNH19366&&MNH19366_33   PAC001560       96.364  1265    41      5       55      1317    1457    196     0.0     2057
MNH19366&&MNH19366_33   EF445167        95.656  1266    48      6       55      1317    1457    196     0.0     2015

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "GQ451293|PAC000672|PAC003100|PAC001560|EF445167"
EF445167        Bacteria;Firmicutes;Clostridia;Clostridiales;Oscillospiraceae;PAC000672_g;NED2E3
GQ451293        Bacteria;Firmicutes;Clostridia;Clostridiales;Oscillospiraceae;PAC000672_g;J299
PAC000672       Bacteria;Firmicutes;Clostridia;Clostridiales;Oscillospiraceae;PAC000672_g;
PAC001560       Bacteria;Firmicutes;Clostridia;Clostridiales;Oscillospiraceae;PAC000672_g;None
PAC003100       Bacteria;Firmicutes;Clostridia;Clostridiales;Oscillospiraceae;PAC000672_g;None

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/checkm/ssu_filter.fna -task blastn -db /work/database/ncbi/16S/current/16SMicrobial -out /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/checkm/ncbi/MNH19366.out -max_target_seqs 5 -num_threads 16 -outfmt 7 
MNH19366&&MNH19366_33   NR_024829.1     89.768  1251    121     6       54      1302    1447    202     0.0     1659
MNH19366&&MNH19366_33   NR_026099.2     89.026  1294    135     6       6       1298    1482    195     0.0     1675
MNH19366&&MNH19366_33   NR_074629.1     89.486  1303    130     6       4       1302    1521    222     0.0     1712
MNH19366&&MNH19366_33   NR_102987.1     89.775  1291    125     5       9       1298    1634    350     0.0     1718
MNH19366&&MNH19366_33   NR_159227.1     90.406  1282    115     6       36      1316    1493    219     0.0     1738

cat /work/database/ncbi/16S/current/16SMicrobial.tab | grep -E "NR_159227|NR_102987|NR_074629|NR_026099|NR_024829"
NR_024829.1     Hungateiclostridium straminisolvens strain CSK1 16S ribosomal RNA, partial sequence
NR_026099.2     Hungateiclostridium aldrichii strain P-1 16S ribosomal RNA, partial sequence
NR_074629.1     Hungateiclostridium thermocellum strain ATCC 27405 16S ribosomal RNA, partial sequence
NR_102987.1     Hungateiclostridium clariflavum DSM 19732 16S ribosomal RNA, partial sequence
NR_159227.1     Monoglobus pectinilyticus strain 14 16S ribosomal RNA, partial sequence


fastANI -q /work/workspace/zhurj/project/2_swgs/run00020/sample/output/MNH19366/MNH19366.fna -r /work/workspace/zhurj/reference/NCBI/genome/GCA_002874775/bwa/GCA_002874775.fna.gz -t 16 --visualize -o /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/fastANI/MNH19366_Monoglobus_pectinilyticus
无结果，相似比例太小


mkdir -p /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/GCA_002874775/genome
cp /work/workspace/zhurj/reference/NCBI/genome/GCA_002874775/bwa/GCA_002874775.fna.gz /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/GCA_002874775/genome/GCA_002874775.fna.gz
mkdir -p /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/GCA_002874775/checkm
checkm ssu_finder /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/GCA_002874775/genome/GCA_002874775.fna /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/GCA_002874775/genome /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/GCA_002874775/checkm -x fna -t 20 


blastn  -query /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/checkm/ssu_filter.fna  -subject /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/GCA_002874775/checkm/ssu_f1.fna  -outfmt 6 > /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/GCA_002874775/checkm/ssu_f1.out

blastn  -query /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/checkm/ssu_filter.fna  -subject /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/GCA_002874775/checkm/ssu_f2.fna  -outfmt 6 > /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/GCA_002874775/checkm/ssu_f2.out

--------------------------------------------------------------------------

fastANI -q /work/workspace/zhurj/project/2_swgs/run00020/sample/output/MNH19366/MNH19366.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000237/GCA_000237085/genomic.fna.gz -t 16 --visualize -o /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/fastANI/MNH19366_Hungateiclostridium_clariflavum
if [[ ! -s /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/fastANI/MNH19366_Hungateiclostridium_clariflavum ]] 
then
  rm /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/fastANI/MNH19366_Hungateiclostridium_clariflavum
fi

没有结果

mkdir -p /work/workspace/zhurj/reference/NCBI/genome/GCA_000237085/bwa && ln -s /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000237/GCA_000237085/genomic.fna.gz /work/workspace/zhurj/reference/NCBI/genome/GCA_000237085/bwa/GCA_000237085.fna.gz && bwa index /work/workspace/zhurj/reference/NCBI/genome/GCA_000237085/bwa/GCA_000237085.fna.gz 

mkdir -p /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/bwa
bwa mem -t 16 /work/workspace/zhurj/reference/NCBI/genome/GCA_000237085/bwa/GCA_000237085.fna.gz /work/rawdata/fastq/genome/MNH/MNH193/MNH19366/E1/L1/S1/MNH19366.1.fq.gz /work/rawdata/fastq/genome/MNH/MNH193/MNH19366/E1/L1/S1/MNH19366.2.fq.gz -o /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/bwa/MNH19366_GCA_000237085.sam
samtools sort -@ 16 -o /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/bwa/MNH19366_GCA_000237085.sorted.sam /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/bwa/MNH19366_GCA_000237085.sam
samtools flagstat --threads 16 /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/bwa/MNH19366_GCA_000237085.sorted.sam > /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/bwa/MNH19366_GCA_000237085.stat

quast -t 16 -o /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/GCA_000237085 /work/workspace/zhurj/reference/NCBI/genome/GCA_000237085/bwa/GCA_000237085.fna.gz
mkdir -p /work/workspace/zhurj/project/2_swgs/run00020/sample/tax/MNH19366/GCA_000237085

taxonkit list --indent "-" --ids 2 --data-dir /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/taxonkitdb -n  -r --line-buffered -o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/bacteriaTaxid20200220.txt

#----------------------------------------------------------------------------------------------------------------------------------------------------------------
# 2020.Feb24-14:49
# further sequencing strains selection
# 1. obtain MNid of SWGS sequenced strains
#    obtain MNid of Sanger sequenced strains

bash /work/workspace/zhurj/project/2_swgs/strainscreen20200217/p/createinput.sh

# 
D:\moonbio\work\12.tmp


/work/workspace/zhurj/script/1_tools/9_strain2seq/strainscreen.py \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/straindb_20200220.txt \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/output 

/work/workspace/zhurj/script/1_tools/9_strain2seq/strainscreen.py \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/test.txt \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/output 

# test -d , -d > 0, natual number
/work/workspace/zhurj/script/1_tools/9_strain2seq/strainscreen.py \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/straindb_20200220.txt \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/output \
-d -1 

# test -n, -n start with char, contain only number, char, _
/work/workspace/zhurj/script/1_tools/9_strain2seq/strainscreen.py \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/straindb_20200220.txt \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/output \
-n 123

/work/workspace/zhurj/script/1_tools/9_strain2seq/strainscreen.py \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/straindb_20200220.txt \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/output \
-n a123?  

《We Are All Fighters》 -- a Speech in English - about Chinese fight the 2019-nCoV epidemic

a beautiful girl with excellent grades =
Curve Wrecker =
walking encyclopedia =
ivy league =
the triangle =
a beautiful summa cum laude girl

/work/workspace/zhurj/script/1_tools/9_strain2seq/strainscreen.py \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/straindb_20200220.txt \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/output \
--sangered /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/sangereds_strain_20200224.txt \
-n test  

/work/workspace/zhurj/script/1_tools/9_strain2seq/strainscreen.py \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/straindb_20200220.txt \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/output \
--sangered /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/sangereds_strain_20200224.txt \
-s /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/swgsed_strain_20200224.txt \
-i /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan.txt \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/straindb_20200227.txt \
-n test  

/work/workspace/zhurj/script/1_tools/9_strain2seq/strainscreen.py \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/straindb_20200220.txt \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/output \
--sangered /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/sangereds_strain_20200224.txt \
-s /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/swgsed_strain_20200224.txt \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/straindb_20200227.txt \
-d 1 \
-n test  

/work/workspace/zhurj/script/1_tools/9_strain2seq/strainscreen.py \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/straindb_20200220.txt \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/output \
--sangered /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/sangereds_strain_20200224.txt \
-s /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/swgsed_strain_20200224.txt \
-i /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan.txt \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/straindb_20200227.txt \
-n strainselect20200227 


# 测序计划：每月成功测序200个菌株（成功提取50个菌株/周）
# 1. 由于此列表中的菌种活化成功率未知，建议3月份第一周，第二周，活化70个/周，评估菌株活化比例后，根据实际情况对以后每周处理菌株数进行调整。
# 共798个菌株，719个 已知菌种,亚种+已属+可能新属（456个菌株 有Sanger 16S测序数据，342个菌株无 16S测序数据；761个菌株taxid，37个无taxid）

#菌种鉴定
#2020.Mar02-14:54
ln -s -r /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/.taxonkit /home/zhurj
cd /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/input
cat test.txt | taxonkit name2taxid --show-rank 
cat /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/input/gramNegative_20200228.txt | taxonkit name2taxid | taxonkit lineage --taxid-field 2 > /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/output/gramNegative_info.txt
cat /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/input/gramNegative_20200228.txt | taxonkit name2taxid --show-rank > /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/output/gramNegative_rank.txt
cat /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/input/gramPositive_20200228.txt | taxonkit name2taxid | taxonkit lineage --taxid-field 2 > /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/output/gramPositive_info.txt
cat /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/input/gramPositive_20200228.txt | taxonkit name2taxid --show-rank > /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/output/gramPositive_rank.txt

cat /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/input/gramNegative_adj_20200228.txt | taxonkit name2taxid | taxonkit lineage --taxid-field 2 > /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/output/gramNegative_adj_info.txt
cat /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/input/gramNegative_adj_20200228.txt | taxonkit name2taxid --show-rank > /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/output/gramNegative_adj_rank.txt

cat /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/input/gramPositive_adj_20200228.txt | taxonkit name2taxid | taxonkit lineage --taxid-field 2 > /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/output/gramPositive_adj_info.txt
cat /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/input/gramPositive_adj_20200228.txt | taxonkit name2taxid --show-rank > /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/output/gramPositive_adj_rank.txt

cat /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/input/stainSubmit20200302.txt | taxonkit name2taxid | taxonkit lineage --taxid-field 2 > /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/output/stainSubmit20200302_info.txt

2020.Mar03-15:10
/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/test/genome -b /work/workspace/zhurj/project/2_swgs/test/output/test -t 16

The second argument to each function is the target format. Currently, the following formats are supported:
newick
nexus
nexml
phyloxml
cdao

from Bio import Phylo
tree = Phylo.read('test.nwk', 'newick')
print(tree)
Phylo.convert('test.nwk', 'newick', 'test.xml', 'phyloxml')
tree = Phylo.parse('test.xml', 'phyloxml')
tree = Phylo.read('test.xml', 'phyloxml')
tree.ladderize()   # Flip branches so deeper clades are displayed at top
Phylo.draw(tree)

from Bio import Phylo
import pylab
file = '/work/workspace/zhurj/project/2_swgs/test/output/test.nwk'
tree = Phylo.read(file, 'newick')
tree.ladderize()
Phylo.draw(tree, do_show=False)
pylab.axis("off")
pylab.savefig("tree2.svg",format='svg', bbox_inches='tight', dpi=300)

清楚比那辆
from IPython import get_ipython
get_ipython().magic('reset -sf') 

%reset -f
del x # delete single variable

from ete3 import Tree, TreeStyle
t = Tree()
t.populate(10, random_branches=True)
ts = TreeStyle()
ts.show_leaf_name = True
ts.scale =  120 # 120 pixels per branch length unit
t.render("mytree0.png", w=183, units="mm", tree_style=ts)

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/test/output/test.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True
t.render("mytree1.png", w=183, units="mm", tree_style=ts)

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/test/output/test.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.scale =  120 # 120 pixels per branch length unit
t.render("mytree2.png", w=183, units="mm", tree_style=ts)

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/test/output/test.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.branch_vertical_margin = 10 # 10 pixels between adjacent branches
t.render("mytree3.png", w=183, units="mm", tree_style=ts)

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/test/output/test.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.rotation = 90
t.render("mytree4.png", w=183, units="mm", tree_style=ts)

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/test/output/test.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.mode = "c"
ts.arc_start = -180 # 0 degrees = 3 o'clock
ts.arc_span = 180
t.render("mytree5.png", w=183, units="mm", tree_style=ts)

from ete3 import Tree, TreeStyle, TextFace
file = '/work/workspace/zhurj/project/2_swgs/test/output/test.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.title.add_face(TextFace("Hello ETE", fsize=20), column=0)
t.render("mytree6.png", w=183, units="mm", tree_style=ts)


#----------------------------------------------------------------
from ete3 import Tree, NodeStyle, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/test/output/test.nwk'
t = Tree(file)

# Basic tree style
ts = TreeStyle()
ts.show_leaf_name = True

# Draws nodes as small red spheres of diameter equal to 10 pixels
nstyle = NodeStyle()
nstyle["shape"] = "sphere"
nstyle["size"] = 10
nstyle["fgcolor"] = "darkred"

# Gray dashed branch lines
nstyle["hz_line_type"] = 1
nstyle["hz_line_color"] = "#cccccc"

# Applies the same static style to all nodes in the tree. Note that,
# if "nstyle" is modified, changes will affect to all nodes
for n in t.traverse():
   n.set_style(nstyle)

t.render("mytree7.png", w=183, units="mm", tree_style=ts)

#----------------------------------------------------------------
from ete3 import Tree, TreeStyle, TextFace
file = '/work/workspace/zhurj/project/2_swgs/test/output/test.nwk'
t = Tree(file)
# Basic tree style
ts = TreeStyle()
ts.show_leaf_name = True

# Add two text faces to different columns
t.add_face(TextFace("hola "), column=0, position = "branch-right")
t.add_face(TextFace("mundo!"), column=1, position = "branch-right")
t.render("mytree8.png", w=183, units="mm", tree_style=ts)

#----------------------------------------------------------------
from ete3 import Tree, TextFace, NodeStyle, TreeStyle
t = Tree("((a,b),c);")

right_c0_r0 = TextFace("right_col0_row0")
right_c0_r1 = TextFace("right_col0_row1")
right_c1_r0 = TextFace("right_col1_row0")
right_c1_r1 = TextFace("right_col1_row1")
right_c1_r2 = TextFace("right_col1_row2")

top_c0_r0 = TextFace("top_col0_row0")
top_c0_r1 = TextFace("top_col0_row1")

bottom_c0_r0 = TextFace("bottom_col0_row0")
bottom_c0_r1 = TextFace("bottom_col0_row1")

aligned_c0_r0 = TextFace("aligned_col0_row0")
aligned_c0_r1 = TextFace("aligned_col0_row1")

aligned_c1_r0 = TextFace("aligned_col1_row0")
aligned_c1_r1 = TextFace("aligned_col1_row1")

all_faces = [right_c0_r0, right_c0_r1, right_c1_r0, right_c1_r1, right_c1_r2, top_c0_r0, \
     top_c0_r1, bottom_c0_r0, bottom_c0_r1, aligned_c0_r0, aligned_c0_r1,\
     aligned_c1_r0, aligned_c1_r1]

# set a border in all faces
for f in all_faces:
    f.margin_left = 1
    f.margin_bottom = 5
    f.margin_top = 5
    f.margin_right = 10


t.add_face(right_c0_r0, column=0, position="branch-right")
t.add_face(right_c0_r1, column=0, position="branch-right")

t.add_face(right_c1_r0, column=1, position="branch-right")
t.add_face(right_c1_r1, column=1, position="branch-right")
t.add_face(right_c1_r2, column=1, position="branch-right")

t.add_face(top_c0_r0, column=0, position="branch-top")
t.add_face(top_c0_r1, column=0, position="branch-top")

t.add_face(bottom_c0_r0, column=0, position="branch-bottom")
t.add_face(bottom_c0_r1, column=0, position="branch-bottom")

a = t&"a"
a.set_style(NodeStyle())
a.img_style["bgcolor"] = "lightgreen"

b = t&"b"
b.set_style(NodeStyle())
b.img_style["bgcolor"] = "indianred"

c = t&"c"
c.set_style(NodeStyle())
c.img_style["bgcolor"] = "lightblue"

t.set_style(NodeStyle())
t.img_style["bgcolor"] = "lavender"
t.img_style["size"] = 12

for leaf in t.iter_leaves():
    leaf.img_style["size"] = 12
    leaf.add_face(right_c0_r0, 0, "branch-right")
    leaf.add_face(aligned_c0_r1, 0, "aligned")
    leaf.add_face(aligned_c0_r0, 0, "aligned")
    leaf.add_face(aligned_c1_r1, 1, "aligned")
    leaf.add_face(aligned_c1_r0, 1, "aligned")

ts = TreeStyle()
ts.show_scale = False
t.render("face_positions.png", w=800, tree_style=ts)

#----------------------------------------------------------------
file = '/work/workspace/zhurj/project/2_swgs/test/output/test.nwk'
t = Tree(file)

from ete3 import Tree, TreeStyle, TextFace

t = Tree( "(a,b);" )

# Basic tree style
ts = TreeStyle()
ts.show_leaf_name = True

# Creates two faces
hola = TextFace("hola")
mundo = TextFace("mundo")

# Set some attributes
hola.margin_top = 10
hola.margin_right = 10
hola.margin_left = 10
hola.margin_bottom = 10
hola.opacity = 0.5 # from 0 to 1
hola.inner_border.width = 1 # 1 pixel border
hola.inner_border.type = 1  # dashed line
hola.border.width = 1
hola.background.color = "LightGreen"

t.add_face(hola, column=0, position = "branch-top")
t.add_face(mundo, column=1, position = "branch-bottom")

t.render("mytree9.png", w=183, units="mm", tree_style=ts)

# test
ln -s /work/assembly/current/MNH/MNH285/MNH28506/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH28506_va.fna
ln -s /work/assembly/current/MNH/MNH297/MNH29773/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH29773_va.fna
ln -s /work/assembly/current/MNH/MNH326/MNH32607/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH32607_va.fna
ln -s /work/assembly/current/MNH/MNH326/MNH32642/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH32642_va.fna
ln -s /work/assembly/current/MNH/MNH332/MNH33285/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH33285_va.fna
ln -s /work/assembly/current/MNH/MNH129/MNH12950/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH12950_va.fna
ln -s /work/assembly/current/MNH/MNH276/MNH27649/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH27649_va.fna
ln -s /work/assembly/current/MNH/MNH369/MNH36991/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH36991_vde.fna
ln -s /work/assembly/current/MNH/MNH094/MNH09418/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH09418_vdi.fna
ln -s /work/assembly/current/MNH/MNH158/MNH15895/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH15895_vdi.fna
ln -s /work/assembly/current/MNH/MNH158/MNH15898/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH15898_vdi.fna
ln -s /work/assembly/current/MNH/MNH189/MNH18969/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH18969_vdi.fna
ln -s /work/assembly/current/MNH/MNH190/MNH19071/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH19071_vdi.fna
ln -s /work/assembly/current/MNH/MNH231/MNH23156/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH23156_vdi.fna
ln -s /work/assembly/current/MNH/MNH257/MNH25706/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH25706_vdi.fna
ln -s /work/assembly/current/MNH/MNH261/MNH26111/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH26111_vdi.fna
ln -s /work/assembly/current/MNH/MNH263/MNH26312/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH26312_vdi.fna
ln -s /work/assembly/current/MNH/MNH312/MNH31247/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH31247_vdi.fna
ln -s /work/assembly/current/MNH/MNH359/MNH35903/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH35903_vdi.fna
ln -s /work/assembly/current/MNH/MNH363/MNH36315/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH36315_vdi.fna
ln -s /work/assembly/current/MNH/MNH378/MNH37837/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH37837_vm.fna
ln -s /work/assembly/current/MNH/MNH326/MNH32688/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH32688_vp.fna
ln -s /work/assembly/current/MNH/MNH341/MNH34110/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH34110_vp.fna
ln -s /work/assembly/current/MNH/MNH094/MNH09468/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH09468_vp.fna
ln -s /work/assembly/current/MNH/MNH259/MNH25911/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH25911_vp.fna
ln -s /work/assembly/current/MNH/MNH349/MNH34931/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH34931_vp.fna
ln -s /work/assembly/current/MNH/MNH283/MNH28327/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH28327_vr.fna
ln -s /work/assembly/current/MNH/MNH341/MNH34134/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH34134_vr.fna
ln -s /work/assembly/current/MNH/MNH363/MNH36331/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH36331_vr.fna
ln -s /work/assembly/current/MNH/MNH400/MNH40010/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH40010_vsp.fna
ln -s /work/assembly/current/MNH/MNH252/MNH25215/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH25215_vt.fna
ln -s /work/assembly/current/MNH/MNH254/MNH25495/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH25495_vt.fna
ln -s /work/assembly/current/MNH/MNH270/MNH27064/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/input/MNH27064_vt.fna

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/phylotree20200304/input -b /work/workspace/zhurj/project/2_swgs/phylotree20200304/output/veillonella -t 16

file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/output/veillonella.nwk'
t = Tree(file)
rectangular (A), circular (B) or radial (C) modes.

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/output/veillonella.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True
t.render("mytree1.png", w=300, units="mm", tree_style=ts)

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/output/veillonella.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.mode = "c"
ts.arc_start = -180 # 0 degrees = 3 o'clock
ts.arc_span = 270
t.render("mytree2.svg", w=300, h=300, units="mm", tree_style=ts)

# different genus phylogenetic tree
ln -s /work/assembly/current/MNH/MNH358/MNH35854/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH35854_Ap.fna
ln -s /work/assembly/current/MNH/MNH397/MNH39764/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH39764_Ap.fna
ln -s /work/assembly/current/MNH/MNH359/MNH35909/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH35909_Ap.fna
ln -s /work/assembly/current/MNH/MNH396/MNH39613/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH39613_Ap.fna
ln -s /work/assembly/current/MNH/MNH329/MNH32902/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH32902_Ap.fna
ln -s /work/assembly/current/MNH/MNH349/MNH34998/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH34998_Bc.fna
ln -s /work/assembly/current/MNH/MNH385/MNH38555/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH38555_Bd.fna
ln -s /work/assembly/current/MNH/MNH379/MNH37902/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH37902_Bd.fna
ln -s /work/assembly/current/MNH/MNH389/MNH38946/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH38946_Bd.fna
ln -s /work/assembly/current/MNH/MNH328/MNH32831/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH32831_Bf.fna
ln -s /work/assembly/current/MNH/MNH381/MNH38125/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH38125_Bf.fna
ln -s /work/assembly/current/MNH/MNH378/MNH37836/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH37836_Bf.fna
ln -s /work/assembly/current/MNH/MNH383/MNH38362/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH38362_Bf.fna
ln -s /work/assembly/current/MNH/MNH382/MNH38247/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH38247_Bf.fna
ln -s /work/assembly/current/MNH/MNH384/MNH38407/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH38407_Bf.fna
ln -s /work/assembly/current/MNH/MNH345/MNH34562/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH34562_Bt.fna
ln -s /work/assembly/current/MNH/MNH346/MNH34661/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH34661_Bt.fna
ln -s /work/assembly/current/MNH/MNH333/MNH33341/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH33341_Ba.fna
ln -s /work/assembly/current/MNH/MNH332/MNH33296/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH33296_Bl.fna
ln -s /work/assembly/current/MNH/MNH333/MNH33355/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH33355_Bl.fna
ln -s /work/assembly/current/MNH/MNH333/MNH33313/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH33313_Ca.fna
ln -s /work/assembly/current/MNH/MNH333/MNH33386/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH33386_Ca.fna
ln -s /work/assembly/current/MNH/MNH334/MNH33404/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH33404_Ca.fna
ln -s /work/assembly/current/MNH/MNH334/MNH33454/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH33454_Ca.fna
ln -s /work/assembly/current/MNH/MNH360/MNH36002/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH36002_Ca.fna
ln -s /work/assembly/current/MNH/MNH360/MNH36025/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH36025_Ca.fna
ln -s /work/assembly/current/MNH/MNH360/MNH36070/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH36070_Ca.fna
ln -s /work/assembly/current/MNH/MNH362/MNH36211/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH36211_Ca.fna
ln -s /work/assembly/current/MNH/MNH363/MNH36311/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH36311_Ca.fna
ln -s /work/assembly/current/MNH/MNH363/MNH36323/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH36323_Ca.fna
ln -s /work/assembly/current/MNH/MNH363/MNH36345/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH36345_Ca.fna
ln -s /work/assembly/current/MNH/MNH364/MNH36401/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH36401_Ca.fna
ln -s /work/assembly/current/MNH/MNH362/MNH36208/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH36208_Ef.fna
ln -s /work/assembly/current/MNH/MNH364/MNH36424/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH36424_Ef.fna
ln -s /work/assembly/current/MNH/MNH361/MNH36148/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH36148_Eh.fna
ln -s /work/assembly/current/MNH/MNH364/MNH36440/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/MNH36440_Eh.fna

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_in/ -b /work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_out/difgenus -t 16

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_out/difgenus.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_out/mytree1.png", w=300, units="mm", tree_style=ts)

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_out/difgenus.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.mode = "c"
ts.arc_start = -180 # 0 degrees = 3 o'clock
ts.arc_span = 270
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/difgenus_out/mytree2.svg", w=300, h=300, units="mm", tree_style=ts)

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/phylotree20200304/cig_in -b /work/workspace/zhurj/project/2_swgs/phylotree20200304/cig_out/cig -t 16

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/cig_out/cig.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/cig_out/mytree1.png", w=300, units="mm", tree_style=ts)

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/cig_out/cig.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.mode = "c"
ts.arc_start = -180 # 0 degrees = 3 o'clock
ts.arc_span = 270
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/cig_out/mytree2.svg", w=300, h=300, units="mm", tree_style=ts)
      
family: Christensenellaceae: 990719  

1. Genus: Christensenella: 990721
species: 626937  |       Christensenella minuta 
GCA_001652705.1 2.93715 51.5 ref
GCA_003628755.1 2.96929 51.4 complete

species: 1805714 |       Christensenella massiliensis    
GCA_900155415.1 2.56019 50.4 Complete

species: 1816678 |       Christensenella timonensis      
GCA_900087015.1 2.65085 51.7 Contig
GCA_902376065.1 2.65085 51.7 Contig

2. Genus: Beduinibacterium : 1987009
species: 1917875  Beduinibacterium massiliense
无

3. no rank: environmental samples: 1229254
species: 1229255 : uncultured Christensenellaceae bacterium

4. no rank: unclassified Christensenellaceae: 1917874
species: 1930012 Christensenellaceae bacterium Phil1
GCA_001940855.1 2.84802 52.2 Scaffold
species: 1930013 Christensenellaceae bacterium Phil7
GCA_001940845.1 2.21367 46.2 Scaffold
species: 1930014 Christensenellaceae bacterium Phil10
GCA_001940895.1 1.01531 46.7 Scaffold
species: 2054177 Christensenellaceae bacterium
GCA_902386865.1    Scaffold

MNH04745_299
MNH04863_305
MNH04759_266
MNH04798_299
MNH04933_273
MNH04949_280
MNH19961_291


/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenella20200305/input -b /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenella20200305/output/Christensenellaceae -t 16

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenella20200305/output/Christensenellaceae.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenella20200305/output/mytree1.png", w=300, units="mm", tree_style=ts)

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenella20200305/output/Christensenellaceae.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.mode = "c"
ts.arc_start = -180 # 0 degrees = 3 o'clock
ts.arc_span = 270
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenella20200305/output/mytree2.svg", w=300, h=300, units="mm", tree_style=ts)

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenella20200305/MNH04863/input -b /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenella20200305/MNH04863/output/MNH04863 -t 16

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenella20200305/MNH04863/output/MNH04863.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenella20200305/MNH04863/output/mytree1.png", w=300, units="mm", tree_style=ts)

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenella20200305/MNH04863/output/MNH04863.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.mode = "c"
ts.arc_start = -180 # 0 degrees = 3 o'clock
ts.arc_span = 270
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenella20200305/MNH04863/output/mytree2.svg", w=300, h=300, units="mm", tree_style=ts)

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenella20200305/allseq/input -b /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenella20200305/allseq/output/MNH04863 -t 16

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenella20200305/allseq/output/MNH04863.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenella20200305/allseq/output/mytree1.png", w=300, units="mm", tree_style=ts)

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenella20200305/allseq/output/MNH04863.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.mode = "c"
ts.arc_start = -180 # 0 degrees = 3 o'clock
ts.arc_span = 270
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenella20200305/allseq/output/mytree2.svg", w=300, h=300, units="mm", tree_style=ts)

cat prokaryotes.txt | grep -P "Bacteroides thetaiotaomicron" | grep "REFR" && cat prokaryotes.txt | grep -P "Bacteroides thetaiotaomicron" | grep -Pi "Complete" | awk -F "\t" '{print $7}' | sort -u
cat prokaryotes.txt | grep -P "Enterococcus gallinarum" | grep -Pi "(REFR|REPR)" && cat prokaryotes.txt | grep -P "Enterococcus gallinarum" | grep -Pi "Complete" | awk -F "\t" '{print $7}' | sort -u

2020.Mar06-10:04
checkm RESULT illustration
1. SOURCE https://github.com/Ecogenomics/CheckM/issues/65
Completeness and contamination are a percentage. Contamination >100% indicates the recovered bin likely contains multiple organisms. For example, contamination of 800% indicates, that on average, each single copy marker gene was observed 8 times! The 0, 1, 2, ..., 5+ columns reported by CheckM indicate the number of times each marker gene was observed.
Strain heterogeneity is more easily viewed as an index between 0 (no strain heterogeneity) and 100 (all markers present >1 appear to be from closely related organisms).
The strain heterogeneity (SH) index indicates the proportion of the contamination that appears to be from the same or similar strains (as determined with an AAI threshold). As such, the primary concern is the amount of contamination and the SH index gives an indication of the source of the contamination (i.e., highly similar or more divergent organisms).

Q: Hello Dparks. Do we need to consider the SH index while doing the downstream analysis after Binning? For example, if we get a Bin with completeness over 90% and contamination below 10%, what should we choose if the SH over 50%?
A: The SH index is worth considering, but isn't nearly as critical as the estimated percentage of contamination. If the SH index is high (ideally 100%), it suggest the majority of contamination is from very similar species and thus any contamination is likely from the pangenome of the species being considered. Alternatively, if the SH index is very low (ideally 0%) this indicates all the contamination is likely from other species (perhaps very divergent species) and thus you may be able to identify it and remove it from the genome. If you wish to try and remove contamination you can look at my companion tool RefineM (https://github.com/dparks1134/RefineM).

Devadas07 commented on 29 Apr 2019
Q: what is the minimum contamination to consider?
A: I generally consider MAGs with contamination<10% or where completeness - 5*contamination>50. See the following:
https://www.nature.com/articles/s41564-017-0012-7
https://www.nature.com/articles/nbt.3893

https://www.nature.com/articles/s41564-017-0012-7
Paper: Recovery of nearly 8,000 metagenome-assembled genomes substantially expands the tree of life
All genomes are estimated to be ≥50% complete and nearly half are ≥90% complete with ≤5% contamination. 

https://www.nature.com/articles/nbt.3893
Paper: Minimum information about a single amplified genome (MISAG) and a metagenome-assembled genome (MIMAG) of bacteria and archaea
but not limited to, assembly quality, and estimates of genome completeness and contamination. 
High-quality draft (SAG/MAG)
Assembly qualitya Multiple fragments where gaps span repetitive regions. Presence of the 23S, 16S, and 5S rRNA genes and at least 18 tRNAs.
Completionb >90%
Contaminationc <5%
Medium-quality draft (SAG/MAG)
Assembly qualitya Many fragments with little to no review of assembly other than reporting of standard assembly statistics.
Completionb ≥50%
Contaminationc <10%
Assembly statistics include but are not limited to: N50, L50, largest contig, number of contigs, assembly size, percentage of reads that map back to the assembly, and number of predicted
genes per genome. bCompletion: ratio of observed single-copy marker genes to total single-copy marker genes in chosen marker gene set. cContamination: ratio of observed single-copy marker
genes in ≥2 copies to total single-copy marker genes in chosen marker gene set.

Reported Statistics
Donovan Parks edited this page on 2 May 2019 · 5 revisions
bin id: unique identifier of genome bin (derived from input fasta file)
marker lineage: indicates the taxonomic rank of the lineage-specific marker set used to estimated genome completeness, contamination, and strain heterogeneity. More detailed information about the placement of a genome within the reference genome tree can be obtained with the tree_qa command. The UID indicates the branch within the reference tree used to infer the marker set applied to estimate the bins quality.
# genomes: number of reference genomes used to infer the lineage-specific marker set
markers: number of marker genes within the inferred lineage-specific marker set
marker sets: number of co-located marker sets within the inferred lineage-specific marker set
0-5+: number of times each marker gene is identified
completeness: estimated completeness of genome as determined from the presence/absence of marker genes and the expected collocalization of these genes (see Methods in the PeerJ preprint for details)
contamination: estimated contamination of genome as determined by the presence of multi-copy marker genes and the expected collocalization of these genes (see Methods in the PeerJ preprint for details)
strain heterogeneity: estimated strain heterogeneity as determined from the number of multi-copy marker pairs which exceed a specified amino acid identity threshold (default = 90%). High strain heterogeneity suggests the majority of reported contamination is from one or more closely related organisms (i.e. potentially the same species), while low strain heterogeneity suggests the majority of contamination is from more phylogenetically diverse sources (see Methods in the CheckM manuscript for more details).
genome size: number of nucleotides (including unknowns specified by N's) in the genome
# ambiguous bases: number of ambiguous (N's) bases in the genome
# scaffolds: number of scaffolds within the genome
# contigs: number of contigs within the genome as determined by splitting scaffolds at any position consisting of more than 10 consecutive ambiguous bases
N50 (scaffolds): N50 statistics as calculated over all scaffolds
N50 (contigs): N50 statistics as calculated over all contigs
longest scaffold: the longest scaffold within the genome
longest contig: the longest contig within the genome
GC: number of G/C nucleotides relative to all A,C,G, and T nucleotides in the genome
coding density: the number of nucleotides within a coding sequence (CDS) relative to all nucleotides in the genome
translation table: indicates which genetic code was used to translate nucleotides into amino acids
# predicted genes: number of predicted coding sequences (CDS) within the genome as determined using Prodigal


checkm tree <bin folder> <output folder>
checkm tree /work/assembly/current/MNH/MNH327/MNH32731 /work/workspace/zhurj/project/2_swgs/run00050/test -t 16 
checkm tree_qa /work/workspace/zhurj/project/2_swgs/run00050/checkM/MNH40010 -f /work/workspace/zhurj/project/2_swgs/run00050/checkM/MNH40010/treeqa.txt --tab_table 

fastANI -q /work/assembly/current/MNH/MNH260/MNH26050/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000020/GCA_000020225/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH26050_AKK

Genus Veillonella 29465
 cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/nodes.dmp | grep -P "\t29465\t" | grep "species" | awk -F "\t" '{print $1}'
29466 Veillonella parvula
39777 Veillonella atypica
39778 Veillonella dispar 
103891 Veillonella criceti
103892 Veillonella ratti
187328 Veillonella montpellierensis 
248315 Veillonella rodentium
248316 Veillonella caviae 
419208 Veillonella denticariosi
423477 Veillonella rogosae
464322 Veillonella magna
1110546 Veillonella tobetsuensis
1502943 Veillonella seminalis
1911679 Veillonella infantium
1936062 Veillonella massiliensis
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^29466\t" | grep -P "scientific name"
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^39777\t" | grep -P "scientific name"
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^39778\t" | grep -P "scientific name"
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^103891\t" | grep -P "scientific name"
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^103892\t" | grep -P "scientific name"
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^187328\t" | grep -P "scientific name"
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^248315\t" | grep -P "scientific name"
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^248316\t" | grep -P "scientific name"
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^419208\t" | grep -P "scientific name"
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^423477\t" | grep -P "scientific name"
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^464322\t" | grep -P "scientific name"
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1110546\t" | grep -P "scientific name"
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1502943\t" | grep -P "scientific name"
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1911679\t" | grep -P "scientific name"
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1936062\t" | grep -P "scientific name"

cat prokaryotes.txt | grep "Veillonella atypica" | grep "Complete"
cat prokaryotes.txt | grep "Veillonella dispar " | grep "REFR"
cat prokaryotes.txt | grep "Veillonella criceti" | grep "Contig"
cat prokaryotes.txt | grep "Veillonella ratti" | grep "Contig"
cat prokaryotes.txt | grep "Veillonella montpellierensis " | grep "REFR"
cat prokaryotes.txt | grep "Veillonella rodentium" | grep "Complete"
cat prokaryotes.txt | grep "Veillonella caviae " | grep "Contig"
cat prokaryotes.txt | grep "Veillonella denticariosi" | grep "Contig"
cat prokaryotes.txt | grep "Veillonella rogosae" | grep "Contig"
cat prokaryotes.txt | grep "Veillonella magna" | grep "REFR"
cat prokaryotes.txt | grep "Veillonella tobetsuensis" | grep "Scaffold"
cat prokaryotes.txt | grep "Veillonella seminalis" | grep "REFR"
cat prokaryotes.txt | grep "Veillonella infantium" | grep "Contig"

cat prokaryotes.txt | grep "Veillonella massiliensis" | grep "REFR"

REFR
Complete
Scaffold
Contig

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input -b /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/output/Veillonella -t 16

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/output/Veillonella.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/output/Veillonella_rec.png", w=300, units="mm", tree_style=ts)

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/output/Veillonella.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.mode = "c"
ts.arc_start = -180 # 0 degrees = 3 o'clock
ts.arc_span = 270
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/output/Veillonella_cir.svg", w=300, h=300, units="mm", tree_style=ts)

fastANI -q /work/assembly/current/MNH/MNH260/MNH26050/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000020/GCA_000020225/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH26050_AKK

/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000315505_Vse232.fna
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000426745_Vmo187.fna
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000428745_Vma220.fna
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900637515_Vdi211.fna
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_001078375_Vto216.fna
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002082765_Vat207.fna
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959775_Vro218.fna
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959855_Vde198.fna
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_003992115_Vra223.fna
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_003992125_Vca196.fna
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900187285_Vro204.fna
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900460315_Vcr224.fna

fastANI -q /work/assembly/current/MNH/MNH400/MNH40010/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959775_Vro218.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH40010_GCA_002959775_Vro218
fastANI -q  -r  -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/
fastANI -q /work/assembly/current/MNH/MNH378/MNH37837/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000428745_Vma220.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH37837_GCA_000428745_Vma220

fastANI -q /work/assembly/current/MNH/MNH053/MNH05351/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH05351_GCA_000024945_Vpa213.txt
fastANI -q /work/assembly/current/MNH/MNH094/MNH09418/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900637515_Vdi211.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH09418_GCA_900637515_Vdi211.txt
fastANI -q /work/assembly/current/MNH/MNH094/MNH09468/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH09468_GCA_000024945_Vpa213.txt
fastANI -q /work/assembly/current/MNH/MNH129/MNH12950/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002082765_Vat207.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH12950_GCA_002082765_Vat207.txt
fastANI -q /work/assembly/current/MNH/MNH158/MNH15895/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH15895_GCA_002959895_Vin202.txt
fastANI -q /work/assembly/current/MNH/MNH158/MNH15898/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH15898_GCA_002959895_Vin202.txt
fastANI -q /work/assembly/current/MNH/MNH189/MNH18969/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900637515_Vdi211.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH18969_GCA_900637515_Vdi211.txt
fastANI -q /work/assembly/current/MNH/MNH190/MNH19071/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900637515_Vdi211.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH19071_GCA_900637515_Vdi211.txt
fastANI -q /work/assembly/current/MNH/MNH198/MNH19841/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900637515_Vdi211.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH19841_GCA_900637515_Vdi211.txt
fastANI -q /work/assembly/current/MNH/MNH200/MNH20015/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900637515_Vdi211.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH20015_GCA_900637515_Vdi211.txt
fastANI -q /work/assembly/current/MNH/MNH202/MNH20256/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH20256_GCA_002959895_Vin202.txt
fastANI -q /work/assembly/current/MNH/MNH208/MNH20893/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900637515_Vdi211.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH20893_GCA_900637515_Vdi211.txt
fastANI -q /work/assembly/current/MNH/MNH226/MNH22663/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900637515_Vdi211.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH22663_GCA_900637515_Vdi211.txt
fastANI -q /work/assembly/current/MNH/MNH227/MNH22792/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH22792_GCA_002959895_Vin202.txt
fastANI -q /work/assembly/current/MNH/MNH231/MNH23156/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900637515_Vdi211.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH23156_GCA_900637515_Vdi211.txt
fastANI -q /work/assembly/current/MNH/MNH233/MNH23366/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002082765_Vat207.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH23366_GCA_002082765_Vat207.txt
fastANI -q /work/assembly/current/MNH/MNH252/MNH25215/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_001078375_Vto216.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH25215_GCA_001078375_Vto216.txt
fastANI -q /work/assembly/current/MNH/MNH254/MNH25495/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_001078375_Vto216.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH25495_GCA_001078375_Vto216.txt
fastANI -q /work/assembly/current/MNH/MNH257/MNH25706/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH25706_GCA_002959895_Vin202.txt
fastANI -q /work/assembly/current/MNH/MNH257/MNH25776/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002082765_Vat207.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH25776_GCA_002082765_Vat207.txt
fastANI -q /work/assembly/current/MNH/MNH259/MNH25911/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900637515_Vdi211.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH25911_GCA_900637515_Vdi211.txt
fastANI -q /work/assembly/current/MNH/MNH261/MNH26111/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002082765_Vat207.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH26111_GCA_002082765_Vat207.txt
fastANI -q /work/assembly/current/MNH/MNH263/MNH26312/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900637515_Vdi211.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH26312_GCA_900637515_Vdi211.txt
fastANI -q /work/assembly/current/MNH/MNH270/MNH27064/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH27064_GCA_000024945_Vpa213.txt
fastANI -q /work/assembly/current/MNH/MNH272/MNH27205/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH27205_GCA_002959895_Vin202.txt
fastANI -q /work/assembly/current/MNH/MNH274/MNH27431/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002082765_Vat207.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH27431_GCA_002082765_Vat207.txt
fastANI -q /work/assembly/current/MNH/MNH276/MNH27649/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002082765_Vat207.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH27649_GCA_002082765_Vat207.txt
fastANI -q /work/assembly/current/MNH/MNH276/MNH27652/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959775_Vro218.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH27652_GCA_002959775_Vro218.txt
fastANI -q /work/assembly/current/MNH/MNH299/MNH29953/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH29953_GCA_002959895_Vin202.txt
fastANI -q /work/assembly/current/MNH/MNH300/MNH30038/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH30038_GCA_002959895_Vin202.txt
fastANI -q /work/assembly/current/MNH/MNH300/MNH30076/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH30076_GCA_002959895_Vin202.txt
fastANI -q /work/assembly/current/MNH/MNH308/MNH30885/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900637515_Vdi211.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH30885_GCA_900637515_Vdi211.txt
fastANI -q /work/assembly/current/MNH/MNH312/MNH31247/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900637515_Vdi211.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH31247_GCA_900637515_Vdi211.txt
fastANI -q /work/assembly/current/MNH/MNH326/MNH32610/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900637515_Vdi211.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH32610_GCA_900637515_Vdi211.txt
fastANI -q /work/assembly/current/MNH/MNH329/MNH32936/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002082765_Vat207.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH32936_GCA_002082765_Vat207.txt
fastANI -q /work/assembly/current/MNH/MNH341/MNH34124/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959775_Vro218.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH34124_GCA_002959775_Vro218.txt
fastANI -q /work/assembly/current/MNH/MNH349/MNH34931/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH34931_GCA_000024945_Vpa213.txt
fastANI -q /work/assembly/current/MNH/MNH359/MNH35903/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH35903_GCA_002959895_Vin202.txt
fastANI -q /work/assembly/current/MNH/MNH363/MNH36315/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH36315_GCA_002959895_Vin202.txt
fastANI -q /work/assembly/current/MNH/MNH369/MNH36991/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959855_Vde198.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH36991_GCA_002959855_Vde198.txt


awk -F "\t" '{print "MNH31247\tVeillonella dispar\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI/MNH31247_GCA_900637515.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI.txt

if [ -f "/work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt" ]; then
  rm /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt -f
fi

awk -F "\t" '{print "MNH05351\tGCA_000024945_Vpa213\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH05351_GCA_000024945_Vpa213.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH09418\tGCA_900637515_Vdi211\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH09418_GCA_900637515_Vdi211.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH09468\tGCA_000024945_Vpa213\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH09468_GCA_000024945_Vpa213.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH12950\tGCA_002082765_Vat207\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH12950_GCA_002082765_Vat207.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH15895\tGCA_002959895_Vin202\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH15895_GCA_002959895_Vin202.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH15898\tGCA_002959895_Vin202\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH15898_GCA_002959895_Vin202.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH18969\tGCA_900637515_Vdi211\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH18969_GCA_900637515_Vdi211.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH19071\tGCA_900637515_Vdi211\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH19071_GCA_900637515_Vdi211.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH19841\tGCA_900637515_Vdi211\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH19841_GCA_900637515_Vdi211.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH20015\tGCA_900637515_Vdi211\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH20015_GCA_900637515_Vdi211.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH20256\tGCA_002959895_Vin202\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH20256_GCA_002959895_Vin202.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH20893\tGCA_900637515_Vdi211\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH20893_GCA_900637515_Vdi211.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH22663\tGCA_900637515_Vdi211\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH22663_GCA_900637515_Vdi211.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH22792\tGCA_002959895_Vin202\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH22792_GCA_002959895_Vin202.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH23156\tGCA_900637515_Vdi211\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH23156_GCA_900637515_Vdi211.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH23366\tGCA_002082765_Vat207\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH23366_GCA_002082765_Vat207.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH25215\tGCA_001078375_Vto216\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH25215_GCA_001078375_Vto216.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH25495\tGCA_001078375_Vto216\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH25495_GCA_001078375_Vto216.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH25706\tGCA_002959895_Vin202\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH25706_GCA_002959895_Vin202.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH25776\tGCA_002082765_Vat207\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH25776_GCA_002082765_Vat207.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH25911\tGCA_900637515_Vdi211\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH25911_GCA_900637515_Vdi211.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH26111\tGCA_002082765_Vat207\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH26111_GCA_002082765_Vat207.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH26312\tGCA_900637515_Vdi211\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH26312_GCA_900637515_Vdi211.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH27064\tGCA_000024945_Vpa213\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH27064_GCA_000024945_Vpa213.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH27205\tGCA_002959895_Vin202\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH27205_GCA_002959895_Vin202.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH27431\tGCA_002082765_Vat207\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH27431_GCA_002082765_Vat207.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH27649\tGCA_002082765_Vat207\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH27649_GCA_002082765_Vat207.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH27652\tGCA_002959775_Vro218\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH27652_GCA_002959775_Vro218.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH29953\tGCA_002959895_Vin202\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH29953_GCA_002959895_Vin202.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH30038\tGCA_002959895_Vin202\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH30038_GCA_002959895_Vin202.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH30076\tGCA_002959895_Vin202\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH30076_GCA_002959895_Vin202.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH30885\tGCA_900637515_Vdi211\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH30885_GCA_900637515_Vdi211.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH31247\tGCA_900637515_Vdi211\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH31247_GCA_900637515_Vdi211.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH32610\tGCA_900637515_Vdi211\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH32610_GCA_900637515_Vdi211.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH32936\tGCA_002082765_Vat207\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH32936_GCA_002082765_Vat207.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH34124\tGCA_002959775_Vro218\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH34124_GCA_002959775_Vro218.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH34931\tGCA_000024945_Vpa213\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH34931_GCA_000024945_Vpa213.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH35903\tGCA_002959895_Vin202\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH35903_GCA_002959895_Vin202.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH36315\tGCA_002959895_Vin202\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH36315_GCA_002959895_Vin202.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH36991\tGCA_002959855_Vde198\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH36991_GCA_002959855_Vde198.txt >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH37837\tGCA_000428745_Vma220\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH37837_GCA_000428745_Vma220 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt
awk -F "\t" '{print "MNH40010\tGCA_002959775_Vro218\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_modify/MNH40010_GCA_002959775_Vro218 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_modify.txt

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input -b /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/output/Veillonella -t 16

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/output/Veillonella.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/output/Veillonella_rec.png", w=300, units="mm", tree_style=ts)

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/output/Veillonella.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.mode = "c"
ts.arc_start = -180 # 0 degrees = 3 o'clock
ts.arc_span = 270
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/output/Veillonella_cir.svg", w=300, h=300, units="mm", tree_style=ts)

fastANI -q  -r  -t 16 -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/

fastANI -t 16 -q /work/assembly/current/MNH/MNH053/MNH05351/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH05351_GCA_000024945_Vpa213
fastANI -t 16 -q /work/assembly/current/MNH/MNH094/MNH09418/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH09418_GCA_000024945_Vpa213
fastANI -t 16 -q /work/assembly/current/MNH/MNH094/MNH09468/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH09468_GCA_000024945_Vpa213
fastANI -t 16 -q /work/assembly/current/MNH/MNH129/MNH12950/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002082765_Vat207.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH12950_GCA_002082765_Vat207
fastANI -t 16 -q /work/assembly/current/MNH/MNH158/MNH15895/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH15895_GCA_002959895_Vin202
fastANI -t 16 -q /work/assembly/current/MNH/MNH158/MNH15898/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH15898_GCA_002959895_Vin202
fastANI -t 16 -q /work/assembly/current/MNH/MNH189/MNH18969/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH18969_GCA_000024945_Vpa213
fastANI -t 16 -q /work/assembly/current/MNH/MNH190/MNH19071/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH19071_GCA_000024945_Vpa213
fastANI -t 16 -q /work/assembly/current/MNH/MNH198/MNH19841/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH19841_GCA_000024945_Vpa213
fastANI -t 16 -q /work/assembly/current/MNH/MNH200/MNH20015/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH20015_GCA_000024945_Vpa213
fastANI -t 16 -q /work/assembly/current/MNH/MNH202/MNH20256/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH20256_GCA_002959895_Vin202
fastANI -t 16 -q /work/assembly/current/MNH/MNH208/MNH20893/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH20893_GCA_000024945_Vpa213
fastANI -t 16 -q /work/assembly/current/MNH/MNH226/MNH22663/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH22663_GCA_000024945_Vpa213
fastANI -t 16 -q /work/assembly/current/MNH/MNH227/MNH22792/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900637515_Vdi211.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH22792_GCA_900637515_Vdi211
fastANI -t 16 -q /work/assembly/current/MNH/MNH231/MNH23156/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH23156_GCA_000024945_Vpa213
fastANI -t 16 -q /work/assembly/current/MNH/MNH233/MNH23366/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002082765_Vat207.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH23366_GCA_002082765_Vat207
fastANI -t 16 -q /work/assembly/current/MNH/MNH252/MNH25215/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_001078375_Vto216.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH25215_GCA_001078375_Vto216
fastANI -t 16 -q /work/assembly/current/MNH/MNH254/MNH25495/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_001078375_Vto216.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH25495_GCA_001078375_Vto216
fastANI -t 16 -q /work/assembly/current/MNH/MNH257/MNH25706/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH25706_GCA_002959895_Vin202
fastANI -t 16 -q /work/assembly/current/MNH/MNH257/MNH25776/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002082765_Vat207.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH25776_GCA_002082765_Vat207
fastANI -t 16 -q /work/assembly/current/MNH/MNH259/MNH25911/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH25911_GCA_000024945_Vpa213
fastANI -t 16 -q /work/assembly/current/MNH/MNH261/MNH26111/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002082765_Vat207.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH26111_GCA_002082765_Vat207
fastANI -t 16 -q /work/assembly/current/MNH/MNH263/MNH26312/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH26312_GCA_000024945_Vpa213
fastANI -t 16 -q /work/assembly/current/MNH/MNH270/MNH27064/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH27064_GCA_000024945_Vpa213
fastANI -t 16 -q /work/assembly/current/MNH/MNH272/MNH27205/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH27205_GCA_002959895_Vin202
fastANI -t 16 -q /work/assembly/current/MNH/MNH274/MNH27431/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002082765_Vat207.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH27431_GCA_002082765_Vat207
fastANI -t 16 -q /work/assembly/current/MNH/MNH276/MNH27649/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002082765_Vat207.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH27649_GCA_002082765_Vat207
fastANI -t 16 -q /work/assembly/current/MNH/MNH276/MNH27652/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900637515_Vdi211.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH27652_GCA_900637515_Vdi211
fastANI -t 16 -q /work/assembly/current/MNH/MNH299/MNH29953/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002082765_Vat207.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH29953_GCA_002082765_Vat207
fastANI -t 16 -q /work/assembly/current/MNH/MNH300/MNH30038/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH30038_GCA_002959895_Vin202
fastANI -t 16 -q /work/assembly/current/MNH/MNH300/MNH30076/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900637515_Vdi211.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH30076_GCA_900637515_Vdi211
fastANI -t 16 -q /work/assembly/current/MNH/MNH308/MNH30885/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH30885_GCA_000024945_Vpa213
fastANI -t 16 -q /work/assembly/current/MNH/MNH312/MNH31247/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH31247_GCA_000024945_Vpa213
fastANI -t 16 -q /work/assembly/current/MNH/MNH326/MNH32610/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH32610_GCA_000024945_Vpa213
fastANI -t 16 -q /work/assembly/current/MNH/MNH329/MNH32936/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002082765_Vat207.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH32936_GCA_002082765_Vat207
fastANI -t 16 -q /work/assembly/current/MNH/MNH341/MNH34124/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959775_Vro218.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH34124_GCA_002959775_Vro218
fastANI -t 16 -q /work/assembly/current/MNH/MNH349/MNH34931/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH34931_GCA_000024945_Vpa213
fastANI -t 16 -q /work/assembly/current/MNH/MNH359/MNH35903/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH35903_GCA_002959895_Vin202
fastANI -t 16 -q /work/assembly/current/MNH/MNH363/MNH36315/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900637515_Vdi211.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH36315_GCA_900637515_Vdi211
fastANI -t 16 -q /work/assembly/current/MNH/MNH369/MNH36991/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900187285_Vro204.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH36991_GCA_900187285_Vro204
fastANI -t 16 -q /work/assembly/current/MNH/MNH378/MNH37837/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000428745_Vma220.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH37837_GCA_000428745_Vma220
fastANI -t 16 -q /work/assembly/current/MNH/MNH400/MNH40010/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH40010_GCA_002959895_Vin202


/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH05351_GCA_000024945_Vpa213
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH09418_GCA_000024945_Vpa213
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH09468_GCA_000024945_Vpa213
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002082765_Vat207.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH12950_GCA_002082765_Vat207
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH15895_GCA_002959895_Vin202
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH15898_GCA_002959895_Vin202
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH18969_GCA_000024945_Vpa213
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH19071_GCA_000024945_Vpa213
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH19841_GCA_000024945_Vpa213
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH20015_GCA_000024945_Vpa213
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH20256_GCA_002959895_Vin202
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH20893_GCA_000024945_Vpa213
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH22663_GCA_000024945_Vpa213
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900637515_Vdi211.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH22792_GCA_900637515_Vdi211
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH23156_GCA_000024945_Vpa213
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002082765_Vat207.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH23366_GCA_002082765_Vat207
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_001078375_Vto216.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH25215_GCA_001078375_Vto216
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_001078375_Vto216.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH25495_GCA_001078375_Vto216
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH25706_GCA_002959895_Vin202
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002082765_Vat207.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH25776_GCA_002082765_Vat207
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH25911_GCA_000024945_Vpa213
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002082765_Vat207.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH26111_GCA_002082765_Vat207
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH26312_GCA_000024945_Vpa213
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH27064_GCA_000024945_Vpa213
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH27205_GCA_002959895_Vin202
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002082765_Vat207.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH27431_GCA_002082765_Vat207
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002082765_Vat207.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH27649_GCA_002082765_Vat207
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900637515_Vdi211.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH27652_GCA_900637515_Vdi211
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002082765_Vat207.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH29953_GCA_002082765_Vat207
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH30038_GCA_002959895_Vin202
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900637515_Vdi211.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH30076_GCA_900637515_Vdi211
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH30885_GCA_000024945_Vpa213
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH31247_GCA_000024945_Vpa213
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH32610_GCA_000024945_Vpa213
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002082765_Vat207.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH32936_GCA_002082765_Vat207
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959775_Vro218.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH34124_GCA_002959775_Vro218
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000024945_Vpa213.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH34931_GCA_000024945_Vpa213
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH35903_GCA_002959895_Vin202
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900637515_Vdi211.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH36315_GCA_900637515_Vdi211
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_900187285_Vro204.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH36991_GCA_900187285_Vro204
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_000428745_Vma220.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH37837_GCA_000428745_Vma220
/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959895_Vin202.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH40010_GCA_002959895_Vin202

if [ -f "/work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt" ];then
rm -f
fi

awk -F "\t" '{print "MNH05351\tGCA_000024945_Vpa213\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH05351_GCA_000024945_Vpa213 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH09418\tGCA_000024945_Vpa213\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH09418_GCA_000024945_Vpa213 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH09468\tGCA_000024945_Vpa213\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH09468_GCA_000024945_Vpa213 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH12950\tGCA_002082765_Vat207\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH12950_GCA_002082765_Vat207 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH15895\tGCA_002959895_Vin202\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH15895_GCA_002959895_Vin202 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH15898\tGCA_002959895_Vin202\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH15898_GCA_002959895_Vin202 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH18969\tGCA_000024945_Vpa213\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH18969_GCA_000024945_Vpa213 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH19071\tGCA_000024945_Vpa213\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH19071_GCA_000024945_Vpa213 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH19841\tGCA_000024945_Vpa213\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH19841_GCA_000024945_Vpa213 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH20015\tGCA_000024945_Vpa213\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH20015_GCA_000024945_Vpa213 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH20256\tGCA_002959895_Vin202\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH20256_GCA_002959895_Vin202 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH20893\tGCA_000024945_Vpa213\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH20893_GCA_000024945_Vpa213 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH22663\tGCA_000024945_Vpa213\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH22663_GCA_000024945_Vpa213 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH22792\tGCA_900637515_Vdi211\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH22792_GCA_900637515_Vdi211 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH23156\tGCA_000024945_Vpa213\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH23156_GCA_000024945_Vpa213 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH23366\tGCA_002082765_Vat207\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH23366_GCA_002082765_Vat207 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH25215\tGCA_001078375_Vto216\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH25215_GCA_001078375_Vto216 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH25495\tGCA_001078375_Vto216\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH25495_GCA_001078375_Vto216 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH25706\tGCA_002959895_Vin202\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH25706_GCA_002959895_Vin202 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH25776\tGCA_002082765_Vat207\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH25776_GCA_002082765_Vat207 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH25911\tGCA_000024945_Vpa213\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH25911_GCA_000024945_Vpa213 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH26111\tGCA_002082765_Vat207\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH26111_GCA_002082765_Vat207 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH26312\tGCA_000024945_Vpa213\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH26312_GCA_000024945_Vpa213 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH27064\tGCA_000024945_Vpa213\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH27064_GCA_000024945_Vpa213 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH27205\tGCA_002959895_Vin202\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH27205_GCA_002959895_Vin202 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH27431\tGCA_002082765_Vat207\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH27431_GCA_002082765_Vat207 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH27649\tGCA_002082765_Vat207\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH27649_GCA_002082765_Vat207 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH27652\tGCA_900637515_Vdi211\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH27652_GCA_900637515_Vdi211 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH29953\tGCA_002082765_Vat207\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH29953_GCA_002082765_Vat207 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH30038\tGCA_002959895_Vin202\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH30038_GCA_002959895_Vin202 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH30076\tGCA_900637515_Vdi211\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH30076_GCA_900637515_Vdi211 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH30885\tGCA_000024945_Vpa213\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH30885_GCA_000024945_Vpa213 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH31247\tGCA_000024945_Vpa213\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH31247_GCA_000024945_Vpa213 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH32610\tGCA_000024945_Vpa213\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH32610_GCA_000024945_Vpa213 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH32936\tGCA_002082765_Vat207\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH32936_GCA_002082765_Vat207 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH34124\tGCA_002959775_Vro218\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH34124_GCA_002959775_Vro218 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH34931\tGCA_000024945_Vpa213\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH34931_GCA_000024945_Vpa213 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH35903\tGCA_002959895_Vin202\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH35903_GCA_002959895_Vin202 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH36315\tGCA_900637515_Vdi211\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH36315_GCA_900637515_Vdi211 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH36991\tGCA_900187285_Vro204\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH36991_GCA_900187285_Vro204 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH37837\tGCA_000428745_Vma220\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH37837_GCA_000428745_Vma220 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt
awk -F "\t" '{print "MNH40010\tGCA_002959895_Vin202\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH40010_GCA_002959895_Vin202 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_202003060616.txt

fastANI -t 16 -q /work/assembly/current/MNH/MNH400/MNH40010/genomic.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Veillonella20200306/input/GCA_002959775_Vro218.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastANI_202003060616/MNH40010_GCA_002959775_Vro218

cat names.dmp | grep -P "\tCollinsella\t"
Collinsella taxid 102106

cat nodes.dmp | grep -P "\t102106\t" | grep "species" | awk -F "\t" '{print $1}'
74426 Collinsella aerofaciens 
147206
147207
626935
1720204
1870987
1871016
1907654
1937461

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Collinsella20200306/input -b /work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Collinsella20200306/output/Collinsella -t 16

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Collinsella20200306/output/Collinsella.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Collinsella20200306/output/Collinsella_rec.png", w=300, units="mm", tree_style=ts)

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Collinsella20200306/output/Collinsella.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.mode = "c"
ts.arc_start = -180 # 0 degrees = 3 o'clock
ts.arc_span = 270
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/run00050Collinsella20200306/output/Collinsella_cir.svg", w=300, h=300, units="mm", tree_style=ts)

fastANI -t 16 -q /work/assembly/current/MNH/MNH326/MNH32695/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000169/GCA_000169035/genomic.fna.gz -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANICollinsella/MNH32695_GCA_000169035
fastANI -t 16 -q /work/assembly/current/MNH/MNH327/MNH32721/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000169/GCA_000169035/genomic.fna.gz -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANICollinsella/MNH32721_GCA_000169035
fastANI -t 16 -q /work/assembly/current/MNH/MNH327/MNH32724/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000169/GCA_000169035/genomic.fna.gz -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANICollinsella/MNH32724_GCA_000169035
fastANI -t 16 -q /work/assembly/current/MNH/MNH327/MNH32746/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000169/GCA_000169035/genomic.fna.gz -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANICollinsella/MNH32746_GCA_000169035
fastANI -t 16 -q /work/assembly/current/MNH/MNH327/MNH32751/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000169/GCA_000169035/genomic.fna.gz -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANICollinsella/MNH32751_GCA_000169035
fastANI -t 16 -q /work/assembly/current/MNH/MNH327/MNH32753/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000169/GCA_000169035/genomic.fna.gz -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANICollinsella/MNH32753_GCA_000169035
fastANI -t 16 -q /work/assembly/current/MNH/MNH327/MNH32778/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000169/GCA_000169035/genomic.fna.gz -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANICollinsella/MNH32778_GCA_000169035


if [ -f "/work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_Collinsella_20200307.txt" ];then
rm -f /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_Collinsella_20200307.txt
fi

awk -F "\t" '{print "MNH32695\tGCA_000169035\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANICollinsella/MNH32695_GCA_000169035 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_Collinsella_20200307.txt
awk -F "\t" '{print "MNH32721\tGCA_000169035\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANICollinsella/MNH32721_GCA_000169035 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_Collinsella_20200307.txt
awk -F "\t" '{print "MNH32724\tGCA_000169035\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANICollinsella/MNH32724_GCA_000169035 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_Collinsella_20200307.txt
awk -F "\t" '{print "MNH32746\tGCA_000169035\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANICollinsella/MNH32746_GCA_000169035 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_Collinsella_20200307.txt
awk -F "\t" '{print "MNH32751\tGCA_000169035\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANICollinsella/MNH32751_GCA_000169035 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_Collinsella_20200307.txt
awk -F "\t" '{print "MNH32753\tGCA_000169035\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANICollinsella/MNH32753_GCA_000169035 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_Collinsella_20200307.txt
awk -F "\t" '{print "MNH32778\tGCA_000169035\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANICollinsella/MNH32778_GCA_000169035 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/fastANI_Collinsella_20200307.txt

FAMILY Coriobacteriaceae 84107
33870
102106
1427376
1472762
1473205
1979843
="cat names.dmp | grep -P "&""""&"^"&A66&"\t"&""""&" | grep -P "&""""&"scientific name"&""""


ln -s /work/assembly/current/MNH/MNH326/MNH32695/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/MNH32695_245.fna
ln -s /work/assembly/current/MNH/MNH327/MNH32721/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/MNH32721_236.fna
ln -s /work/assembly/current/MNH/MNH327/MNH32724/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/MNH32724_236.fna
ln -s /work/assembly/current/MNH/MNH327/MNH32746/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/MNH32746_238.fna
ln -s /work/assembly/current/MNH/MNH327/MNH32751/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/MNH32751_236.fna
ln -s /work/assembly/current/MNH/MNH327/MNH32753/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/MNH32753_238.fna
ln -s /work/assembly/current/MNH/MNH327/MNH32778/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/MNH32778_230.fna
cp /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000195/GCA_000195315/genomic.fna.gz /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_000195315_Cgl211.fna.gz && gunzip /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_000195315_Cgl211.fna.gz
cp /work/database/ncbi/genome/all/GCA/GCA_003/GCA_003725/GCA_003725335/genomic.fna.gz /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_003725335_Pca247.fna.gz && gunzip /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_003725335_Pca247.fna.gz
cp /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000321/GCA_000321165/genomic.fna.gz /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_000321165_Eti233.fna.gz && gunzip /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_000321165_Eti233.fna.gz
cp /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000311/GCA_000311845/genomic.fna.gz /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_000311845_Ema228.fna.gz && gunzip /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_000311845_Ema228.fna.gz
cp /work/database/ncbi/genome/all/GCA/GCA_900/GCA_900186/GCA_900186505/genomic.fna.gz /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_900186505_Eph237.fna.gz && gunzip /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_900186505_Eph237.fna.gz
cp /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000236/GCA_000236865/genomic.fna.gz /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_000236865_San239.fna.gz && gunzip /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_000236865_San239.fna.gz
cp /work/database/ncbi/genome/all/GCA/GCA_004/GCA_004135/GCA_004135645/genomic.fna.gz /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_004135645_Sfa274.fna.gz && gunzip /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_004135645_Sfa274.fna.gz
cp /work/database/ncbi/genome/all/GCA/GCA_900/GCA_900078/GCA_900078545/genomic.fna.gz /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_900078545_Oma180.fna.gz && gunzip /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_900078545_Oma180.fna.gz
cp /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000169/GCA_000169035/genomic.fna.gz /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_000169035_Cae243.fna.gz && gunzip /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_000169035_Cae243.fna.gz
cp /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000156/GCA_000156215/genomic.fna.gz /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_000156215_Cst247.fna.gz && gunzip /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_000156215_Cst247.fna.gz
cp /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000156/GCA_000156175/genomic.fna.gz /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_000156175_Cin180.fna.gz && gunzip /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_000156175_Cin180.fna.gz
cp /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000225/GCA_000225705/genomic.fna.gz /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_000225705_Cta249.fna.gz && gunzip /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_000225705_Cta249.fna.gz
cp /work/database/ncbi/genome/all/GCA/GCA_900/GCA_900046/GCA_900046475/genomic.fna.gz /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_900046475_Cih283.fna.gz && gunzip /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_900046475_Cih283.fna.gz
cp /work/database/ncbi/genome/all/GCA/GCA_900/GCA_900176/GCA_900176655/genomic.fna.gz /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_900176655_Cva216.fna.gz && gunzip /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_900176655_Cva216.fna.gz
cp /work/database/ncbi/genome/all/GCA/GCA_900/GCA_900119/GCA_900119895/genomic.fna.gz /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_900119895_Cph210.fna.gz && gunzip /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_900119895_Cph210.fna.gz
cp /work/database/ncbi/genome/all/GCA/GCA_900/GCA_900155/GCA_900155365/genomic.fna.gz /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_900155365_Cbo187.fna.gz && gunzip /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_900155365_Cbo187.fna.gz
cp /work/database/ncbi/genome/all/GCA/GCA_900/GCA_900199/GCA_900199705/genomic.fna.gz /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_900199705_Cpr173.fna.gz && gunzip /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input/GCA_900199705_Cpr173.fna.gz

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/input -b /work/workspace/zhurj/project/2_swgs/phylotree20200304/Coriobacteriaceae20200307/output/f_Coriobacteriaceae -t 16

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input -b /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/output/Enterococcusgallinarum -t 16

-o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/

fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna

fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_001558875_Com382.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_001558875_Com382
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_001886155_Sca373.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_001886155_Sca373
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_001544275_Con377.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_001544275_Con377
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_902363355_Sca315.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_902363355_Sca315
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_902162015_Sca333.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_902162015_Sca333
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_902160415_Sca340.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_902160415_Sca340
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_902160665_Sca370.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_902160665_Sca370
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_902159065_Sca354.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_902159065_Sca354
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_902159025_Sca354.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_902159025_Sca354
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_902159265_Con370.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_902159265_Con370
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_902159105_Con340.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_902159105_Con340
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_003320875_Sca361.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_003320875_Sca361
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_003434285_Sca315.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_003434285_Sca315
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_003796985_Con352.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_003796985_Con352
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_003797625_Con352.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_003797625_Con352
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_003797465_Con352.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_003797465_Con352
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_003797555_Con352.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_003797555_Con352
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_003797485_Con330.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_003797485_Con330
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_002140495_Sca339.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_002140495_Sca339
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_002383155_Sca290.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_002383155_Sca290
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_002360215_Con331.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_002360215_Con331
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_002901245_Con363.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_002901245_Con363
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_002077495_Con338.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_002077495_Con338
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_008016345_Sca237.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_008016345_Sca237
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_008082695_Con326.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_008082695_Con326
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_008081895_Con316.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_008081895_Con316
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_008083145_Con322.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_008083145_Con322
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_000157255_Sca316.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_000157255_Sca316
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_000703185_Con347.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_000703185_Con347
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_000829545_Con298.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_000829545_Con298
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_901543425_Con318.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_901543425_Con318
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_900447935_Con393.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_900447935_Con393
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_900461195_Con389.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_900461195_Con389
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_004570995_Con366.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_004570995_Con366

awk -F "\t" '{print "MNH32695\tGCA_000169035\t"$3"\t"$4"\t"$5}' 

if [ -f "/work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt" ]; then
rm /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt -f
fi

awk -F "\t" '{print "MNH39848\tGCA_001558875_Com382\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_001558875_Com382 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_001886155_Sca373\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_001886155_Sca373 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_001544275_Con377\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_001544275_Con377 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_902363355_Sca315\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_902363355_Sca315 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_902162015_Sca333\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_902162015_Sca333 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_902160415_Sca340\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_902160415_Sca340 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_902160665_Sca370\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_902160665_Sca370 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_902159065_Sca354\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_902159065_Sca354 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_902159025_Sca354\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_902159025_Sca354 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_902159265_Con370\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_902159265_Con370 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_902159105_Con340\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_902159105_Con340 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_003320875_Sca361\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_003320875_Sca361 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_003434285_Sca315\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_003434285_Sca315 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_003796985_Con352\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_003796985_Con352 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_003797625_Con352\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_003797625_Con352 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_003797465_Con352\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_003797465_Con352 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_003797555_Con352\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_003797555_Con352 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_003797485_Con330\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_003797485_Con330 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_002140495_Sca339\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_002140495_Sca339 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_002383155_Sca290\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_002383155_Sca290 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_002360215_Con331\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_002360215_Con331 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_002901245_Con363\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_002901245_Con363 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_002077495_Con338\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_002077495_Con338 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_008016345_Sca237\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_008016345_Sca237 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_008082695_Con326\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_008082695_Con326 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_008081895_Con316\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_008081895_Con316 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_008083145_Con322\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_008083145_Con322 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_000157255_Sca316\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_000157255_Sca316 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_000703185_Con347\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_000703185_Con347 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_000829545_Con298\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_000829545_Con298 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_901543425_Con318\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_901543425_Con318 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_900447935_Con393\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_900447935_Con393 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_900461195_Con389\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_900461195_Con389 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt
awk -F "\t" '{print "MNH39848\tGCA_004570995_Con366\t"$3"\t"$4"\t"$5}' /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_004570995_Con366 >> /work/workspace/zhurj/project/2_swgs/run00050/merge/Enterococcusgallinarum200309.txt

cat nodes.dmp | grep -P "\t906\t" | grep species | awk -F "\t" '{print $1}'
cp: cannot stat ‘/work/database/ncbi/genome/all/GCA/GCA_009/GCA_009730/GCA_009730155/genomic.fna.gz’: No such file or directory
cp: cannot stat ‘/work/database/ncbi/genome/all/GCA/GCA_009/GCA_009830/GCA_009830395/genomic.fna.gz’: No such file or directory
cp: cannot stat ‘/work/database/ncbi/genome/all/GCA/GCA_009/GCA_009832/GCA_009832495/genomic.fna.gz’: No such file or directory
cp: cannot stat ‘/work/database/ncbi/genome/all/GCA/GCA_901/GCA_901212/GCA_901212575/genomic.fna.gz’: No such file or directory
cp: cannot stat ‘/work/database/ncbi/genome/all/GCA/GCA_009/GCA_009734/GCA_009734765/genomic.fna.gz’: No such file or directory
cp: cannot stat ‘/work/database/ncbi/genome/all/GCA/GCA_009/GCA_009734/GCA_009734745/genomic.fna.gz’: No such file or directory
cp: cannot stat ‘/work/database/ncbi/genome/all/GCA/GCA_009/GCA_009930/GCA_009930755/genomic.fna.gz’: No such file or directory

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/phylotree20200304/Megasphaera20200309/input -b /work/workspace/zhurj/project/2_swgs/phylotree20200304/Megasphaera20200309/output/Megasphaera -t 16

/work/workspace/zhurj/project/2_swgs/run00050/checkM/MNH33226

>genomic&&NODE_4_length_224896_cov_581.289895
ATATATCATGGAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGAGAAGAGATGAGAAGCTTGCTTCTTATCGATTCGAGTGGCAAACGGGTGAGTAACGCGTAAGCAACCTGCCCTCTAGATGGGGACAACAGCTGGAAACGGCTGCTAATACCGAATACGTTCTTTCTGTCGCATGGCAGAGGGAAGAAAGGGAGGCTCTTCGGAGCTTTCGCTGGAGGAGGGGCTTGCGTCTGATTAGCTAGTTGGAGGGGTAACGGCCCACCAAGGCGACGATCAGTAGCCGGTCTGAGAGGATGAACGGCCACATTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATCTTCCGCAATGGACGAAAGTCTGACGGAGCAACGCCGCGTGAACGATGACGGCCTTCGGGTTGTAAAGTTCTGTTATACGGGACGAATGGTACGACGGTCAATACCCGTCGTAAGTGACGGTACCGTAAGAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGGGCGCGCAGGCGGCGTCGTAAGTCGGTCTTAAAAGTGCGGGGCTTAACCCCGTGAGGGGACCGAAACTGCGATGCTAGAGTATCGGAGAGGAAAGCGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTGGCGAAAGCGGCTTTCTGGACGACAACTGACGCTGAGGCGCGAAAGCCAGGGGAGCAAACGGGATTAGATACCCCGGTAGTCCTGGCCGTAAACGATGGATACTAGGTGTAGGAGGTATCGACCCCTTCTGTGCCGGAGTTAACGCAATAAGTATCCCGCCTGGGGAGTACGGCCGCAAGGCTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGTATGTGGTTTAATTCGACGCAACGCGAAGAACCTTACCAAGCCTTGACAT
Description
Max Score
Total Score
Query Cover
E value
Per. Ident
Accession
Select seq NR_134080.1  Megasphaera indica strain NMBHI-10 16S ribosomal RNA, partial sequence  1832  1832  98% 0.0 99.41%  NR_134080.1
Select seq NR_102980.1  Megasphaera elsdenii DSM 20460 16S ribosomal RNA, partial sequence  1792  1792  98% 0.0 98.71%  NR_102980.1
Select seq NR_029207.1  Megasphaera elsdenii strain LC 1 16S ribosomal RNA, partial sequence  1786  1786  98% 0.0 98.62%  NR_029207.1
Select seq NR_133027.1  Megasphaera massiliensis strain NP3 16S ribosomal RNA, partial sequence 1685  1685  98% 0.0 96.83%  NR_133027.1

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/run00050/checkM/MNH33226/ssu_filter.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/run00050/checkM/MNH33226/MNH33226.out -max_target_seqs 5 -num_threads 16 -outfmt 7 
# BLASTN 2.9.0+
# Query: genomic&&NODE_4_length_224896_cov_581.289895
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
genomic&&NODE_4_length_224896_cov_581.289895    HM990965        99.395  992     3       3       32      1020    1       992     0.0     1752
genomic&&NODE_4_length_224896_cov_581.289895    HE576794        98.584  989     14      0       32      1020    1       989     0.0     1721
genomic&&NODE_4_length_224896_cov_581.289895    JX424772        96.764  989     30      2       32      1020    1       987     0.0     1633
genomic&&NODE_4_length_224896_cov_581.289895    GU366015        93.857  993     57      2       32      1020    1       993     0.0     1511
genomic&&NODE_4_length_224896_cov_581.289895    CP029462        93.737  990     61      1       32      1020    1       990     0.0     1503


cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "HM990965|HE576794"
HE576794        Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Megasphaera;Megasphaera elsdenii;DSM 20460(T)
HM990965        Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Megasphaera;Megasphaera indica;NMBHI-10(T)

>genomic&&NODE_27_length_816_cov_3376.280108
AATGTCAAGGCTTGGTAAGGTTCTTCGCGTTGCGTCGAATTAAACCACATACTCCACCGCTTGTGCGGGCCCCCGTCAATTCCTTTGAGTTTCAGCCTTGCGGCCGTACTCCCCAGGCGGGATACTTATTGCGTTAACTCCGGCACAGAAGGGGTCGATACCTCCTACACCTAGTATCCATCGTTTACGGCCAGGACTACCGGGGTATCTAATCCCGTTTGCTCCCCTGGCTTTCGCGCCTCAGCGTCAGTTGTCGTCCAGAAAGCCGCTTTCGCCACTGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCGCTTTCCTCTCCGATACTCTAGCATCGCAGTTTCGGTCCCCTCACGGGGTTAAGCCCCGCACTTTTAAGACCGACTTACGACGCCGCCTGCGCGCCCTTTACGCCCAATAATTCCGGACAACGCTTGCCACCTACGTATTACCGCGGCTGCTGGCACGTAGTTAGCCGTGGCTTTCTCTTACGGTACCGTCACTTACGACGGGTATTGACCGTCGTACCATTCGTCCCGTATAACAGAACTTTACAACCCGAAGGCCGTCATCGTTCACGCGGCGTTGCTCCGTCAGACTTTCGTCCATTGCGGAAGATTCCCCACTGCTGCCTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCATCCTCTCAGACCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCTCCAACTAGCTAATCAGACGCAAGCCCCTCCTCCAGCGAAA
Select seq NR_113306.1  Megasphaera elsdenii DSM 20460 16S ribosomal RNA, partial sequence  1445  1445  100%  0.0 99.87%  NR_113306.1
Select seq NR_113305.1  Megasphaera elsdenii DSM 20460 16S ribosomal RNA, partial sequence  1445  1445  100%  0.0 99.87%  NR_113305.1
Select seq NR_134080.1  Megasphaera indica strain NMBHI-10 16S ribosomal RNA, partial sequence  1428  1428  100%  0.0 99.49%  NR_134080.1
Select seq NR_102980.1  Megasphaera elsdenii DSM 20460 16S ribosomal RNA, partial sequence  1399  1399  100%  0.0 98.85%  NR_102980.1
Select seq NR_029207.1  Megasphaera elsdenii strain LC 1 16S ribosomal RNA, partial sequence  1387  1387  100%  0.0 98.60%  NR_029207.1


/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/run00050/checkM/MNH19072/ssu_filter1.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/run00050/checkM/MNH19072/MNH19072_1.out -max_target_seqs 5 -num_threads 16 -outfmt 7 
# BLASTN 2.9.0+
# Query: genomic&&NODE_27_length_816_cov_3376.280108
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
genomic&&NODE_27_length_816_cov_3376.280108     HM990965        99.492  787     1       3       1       784     993     207     0.0     1391
genomic&&NODE_27_length_816_cov_3376.280108     HE576794        98.724  784     10      0       1       784     990     207     0.0     1370
genomic&&NODE_27_length_816_cov_3376.280108     JX424772        97.449  784     20      0       1       784     988     205     0.0     1324

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "HM990965|HE576794|JX424772"
HE576794        Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Megasphaera;Megasphaera elsdenii;DSM 20460(T)
HM990965        Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Megasphaera;Megasphaera indica;NMBHI-10(T)
JX424772        Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Megasphaera;Megasphaera massiliensis;NP3(T)

>genomic&&NODE_24_length_956_cov_2154.617747
CTAGAAAGGAGGTGATCCAGCCGCACCTTCCGATACGGCTACCTTGTTACGACTTCACCCCAATCATCGCCCCCACCTTCGACGGCTGGCTCCTTGCGGTTACCTCACCGGCTTCGGGTGTGAATGACTTTCGTGGTGTGACGGGCGGTGTGTACAAGGCCCGGGAACGTATTCACCGCAGTATGCTGACCTGCGATTACTAGCGATTCCTGCTTCATGCAGGCGGGTTGCAGCCTGCAATCCGAACTGGGACTCTGTTTTTGGGGTTTGCTCCGGATCGCTCCTTCGCTTCCCTCTATTAAGAGCCATTGTAGTACGTGTGTAGCCCAAGCCATAAGGGGCATGATGACTTGACGTCATCCCCGCCTTCCTCCGCATTGTCTGCGGCAGTCTCTCCTGAGTCCCCGACTTAACTCGCTGGTAACAGAAGATAGGGGTTGCGCTCGTTGCGGGACTTAACCCAACATCTCACGACACGAGCTGACGACAGCCGTGCACCACCTGTTTTCTTGTCCTCCGAAGAGGAACGGGACATCTCTGTCCCTAGCAATCAATGTCAAGGCTTGGTAAGGTTCTTCGCGTTGCGTCGAATTAAACCACATACTCCACCGCTTGTGCGG
Select seq NR_134080.1  Megasphaera indica strain NMBHI-10 16S ribosomal RNA, partial sequence  1064  1064  100%  0.0 97.74%  NR_134080.1
Select seq NR_102980.1  Megasphaera elsdenii DSM 20460 16S ribosomal RNA, partial sequence  1055  1055  98% 0.0 97.87%  NR_102980.1

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/run00050/checkM/MNH19072/ssu_filter2.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/run00050/checkM/MNH19072/MNH19072_2.out -max_target_seqs 5 -num_threads 16 -outfmt 7 
# BLASTN 2.9.0+
# Query: genomic&&NODE_24_length_956_cov_2154.617747
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
genomic&&NODE_24_length_956_cov_2154.617747     HM990965        98.230  565     9       1       56      620     1489    926     0.0     971
genomic&&NODE_24_length_956_cov_2154.617747     HE576794        97.699  565     12      1       56      620     1486    923     0.0     957
genomic&&NODE_24_length_956_cov_2154.617747     AWXA01000023    97.522  565     14      0       56      620     1490    926     0.0     957
genomic&&NODE_24_length_956_cov_2154.617747     GQ480004        97.345  565     15      0       56      620     1491    927     0.0     952

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "HM990965|HE576794|AWXA01000023|GQ480004"
AWXA01000023    Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Megasphaera;AWXA_s;BV3C16-1
GQ480004        Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Megasphaera;BXHA50
HE576794        Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Megasphaera;Megasphaera elsdenii;DSM 20460(T)
HM990965        Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Megasphaera;Megasphaera indica;NMBHI-10(T)

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/Megasphaera20200309/output/Megasphaera.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/Megasphaera20200309/output/Megasphaera_rec.png", w=300, units="mm", tree_style=ts)

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/Megasphaera20200309/output/Megasphaera.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.mode = "c"
ts.arc_start = -180 # 0 degrees = 3 o'clock
ts.arc_span = 270
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/Megasphaera20200309/output/Megasphaera_cir.svg", w=300, h=300, units="mm", tree_style=ts)


2020.Mar10-11:16
checkM
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/quast/MNH27614 /work/assembly/current/MNH/MNH276/MNH27614/genomic.fna
checkm tree <bin folder> <output folder>
checkm tree /work/assembly/current/MNH/MNH327/MNH32731 /work/workspace/zhurj/project/2_swgs/run00050/test -t 16 
checkm ssu_finder /work/assembly/current/MNH/MNH276/MNH27614/genomic.fna /work/assembly/current/MNH/MNH276/MNH27614 /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH27614 -x fna -t 16 
checkm lineage_wf /work/assembly/current/MNH/MNH276/MNH27614 /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH27614 -x fna -t 16 --tab_table -f /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH27614/MNH27614_checkm.txt
checkm tree_qa /work/workspace/zhurj/project/2_swgs/run00050/checkM/MNH40010 -f /work/workspace/zhurj/project/2_swgs/run00050/checkM/MNH40010/treeqa.txt --tab_table 

checkm ssu_finder /work/assembly/current/MNH/MNH056/MNH05676/genomic.fna /work/assembly/current/MNH/MNH056/MNH05676 /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/test -x fna -t 16
checkm lineage_wf /work/assembly/current/MNH/MNH056/MNH05676 /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/test1 -x fna -t 16 --tab_table -f /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/test1/test1_checkm.txt
checkm tree_qa /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/test1 -f /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/test1/treeqa.txt --tab_table 
# /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/test can be created automatically
# /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/test1 can be created automatically

# names.dmp obtain taxid from species name
cat names.dmp | grep -P "\tVerrucomicrobiaceae\t"
Verrucomicrobiaceae 203557
#genus of Verrucomicrobiaceae
cat nodes.dmp | grep -P "\t203557\t" | grep -P "genus" | awk -F "\t" '{print $1}'
2735
48463
177411
518753
518754
518755
574899
1348508
1383070
2576070

# taxid to specific name
cat names.dmp | grep -P "^2735\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' 
cat names.dmp | grep -P "^48463\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}'
cat names.dmp | grep -P "^177411\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}'
cat names.dmp | grep -P "^518753\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}'
cat names.dmp | grep -P "^518754\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}'
cat names.dmp | grep -P "^518755\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}'
cat names.dmp | grep -P "^574899\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}'
cat names.dmp | grep -P "^1348508\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}'
cat names.dmp | grep -P "^1383070\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}'
cat names.dmp | grep -P "^2576070\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}'
2735    Verrucomicrobium
48463   Prosthecobacter
177411  Fucophilus
518753  Luteolibacter
518754  Persicirhabdus
518755  Roseibacillus
574899  Haloferula
1348508 Roseimicrobium
1383070 Brevifollis
2576070 Verrucobacter

#------------------------------------------------------------------------------------------
family: 
1647988 Akkermansiaceae
cat nodes.dmp | grep -P "\t1647988\t" | grep -P "genus" | awk -F "\t" '{print $1}'
# 
1951308 no rank
239934 Akkermansia
cat nodes.dmp | grep -P "\t239934\t" | awk -F "\t" '{print $1"\t"$5}'
239935  species
512293  no rank
1679444 species
2608915 no rank
cat names.dmp | grep -P "^239935\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}'
cat names.dmp | grep -P "^512293\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}'
cat names.dmp | grep -P "^1679444\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}'
cat names.dmp | grep -P "^2608915\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}'
239935  Akkermansia muciniphila
512293  environmental samples no genome
1679444 Akkermansia glycaniphila
2608915 unclassified Akkermansia no genome

  239935 [species] Akkermansia muciniphila
    349741 [no rank] Akkermansia muciniphila ATCC BAA-835
  512293 [no rank] environmental samples
    512294 [species] uncultured Akkermansia sp.
    1131822 [species] uncultured Akkermansia sp. SMG25
    1262691 [species] Akkermansia sp. CAG:344
    1263034 [species] Akkermansia muciniphila CAG:154
  1679444 [species] Akkermansia glycaniphila
  2608915 [no rank] unclassified Akkermansia
    1131336 [species] Akkermansia sp. KLE1605
    1574264 [species] Akkermansia sp. KLE1797
    1574265 [species] Akkermansia sp. KLE1798
    1638783 [species] Akkermansia sp. UNK.MGS-1
    1755639 [species] Akkermansia sp. MC_55
    1872421 [species] Akkermansia sp.
    1896967 [species] Akkermansia sp. 54_46
    1929996 [species] Akkermansia sp. Phil8
    1945963 [species] Akkermansia sp. UBA3271
    1945964 [species] Akkermansia sp. UBA5128
    1945965 [species] Akkermansia sp. UBA7059
    1945966 [species] Akkermansia sp. UBA7090
    2478952 [species] Akkermansia sp. aa_0143
    2584557 [species] Akkermansia sp. BIOML-A1
    2584558 [species] Akkermansia sp. BIOML-A2
    2584559 [species] Akkermansia sp. BIOML-A3
    2584560 [species] Akkermansia sp. BIOML-A4
    2584561 [species] Akkermansia sp. BIOML-A5
    2584562 [species] Akkermansia sp. BIOML-A6
    2584563 [species] Akkermansia sp. BIOML-A7
    2584564 [species] Akkermansia sp. BIOML-A8
    2584565 [species] Akkermansia sp. BIOML-A9
    2584566 [species] Akkermansia sp. BIOML-A10
    2584567 [species] Akkermansia sp. BIOML-A11
    2584568 [species] Akkermansia sp. BIOML-A12
    2584569 [species] Akkermansia sp. BIOML-A13
    2584570 [species] Akkermansia sp. BIOML-A14
    2584571 [species] Akkermansia sp. BIOML-A15
    2584572 [species] Akkermansia sp. BIOML-A16
    2584573 [species] Akkermansia sp. BIOML-A17
    2584574 [species] Akkermansia sp. BIOML-A18
    2584575 [species] Akkermansia sp. BIOML-A19
    2584576 [species] Akkermansia sp. BIOML-A20
    2584577 [species] Akkermansia sp. BIOML-A21
    2584578 [species] Akkermansia sp. BIOML-A22
    2584579 [species] Akkermansia sp. BIOML-A23
    2584580 [species] Akkermansia sp. BIOML-A24
    2584581 [species] Akkermansia sp. BIOML-A25
    2584582 [species] Akkermansia sp. BIOML-A26
    2584583 [species] Akkermansia sp. BIOML-A27
    2584584 [species] Akkermansia sp. BIOML-A28
    2584585 [species] Akkermansia sp. BIOML-A29
    2584586 [species] Akkermansia sp. BIOML-A30
    2584587 [species] Akkermansia sp. BIOML-A31
    2584588 [species] Akkermansia sp. BIOML-A32
    2584589 [species] Akkermansia sp. BIOML-A33
    2584590 [species] Akkermansia sp. BIOML-A34
    2584591 [species] Akkermansia sp. BIOML-A35
    2584592 [species] Akkermansia sp. BIOML-A36
    2584593 [species] Akkermansia sp. BIOML-A37
    2584594 [species] Akkermansia sp. BIOML-A38
    2584595 [species] Akkermansia sp. BIOML-A39
    2584596 [species] Akkermansia sp. BIOML-A40
    2584597 [species] Akkermansia sp. BIOML-A41
    2584598 [species] Akkermansia sp. BIOML-A42
    2584599 [species] Akkermansia sp. BIOML-A43
    2584600 [species] Akkermansia sp. BIOML-A44
    2584601 [species] Akkermansia sp. BIOML-A45
    2584602 [species] Akkermansia sp. BIOML-A46
    2584603 [species] Akkermansia sp. BIOML-A47
    2584604 [species] Akkermansia sp. BIOML-A48
    2584605 [species] Akkermansia sp. BIOML-A49
    2584606 [species] Akkermansia sp. BIOML-A50
    2584607 [species] Akkermansia sp. BIOML-A51
    2584608 [species] Akkermansia sp. BIOML-A52
    2584609 [species] Akkermansia sp. BIOML-A53
    2584610 [species] Akkermansia sp. BIOML-A54
    2584611 [species] Akkermansia sp. BIOML-A55
    2584612 [species] Akkermansia sp. BIOML-A56
    2584613 [species] Akkermansia sp. BIOML-A57
    2584614 [species] Akkermansia sp. BIOML-A58
    2584615 [species] Akkermansia sp. BIOML-A59
    2584616 [species] Akkermansia sp. BIOML-A60
    2584617 [species] Akkermansia sp. BIOML-A61
    2584618 [species] Akkermansia sp. BIOML-A62
    2584619 [species] Akkermansia sp. BIOML-A63
    2584620 [species] Akkermansia sp. BIOML-A64
    2584621 [species] Akkermansia sp. BIOML-A65
    2584622 [species] Akkermansia sp. BIOML-A66
    2584623 [species] Akkermansia sp. BIOML-A67
    2608916 [species] Akkermansia sp. H38


1951308 no rank
cat nodes.dmp | grep -P "\t1951308\t" | awk -F "\t" '{print $1"\t"$5}'

1951347 species 
1951348 species 
1951349 species 
1951350 species 
1951351 species 
1951352 species 
1951353 species 
1951354 species 
1951355 species 
1951356 species 
1951357 species 
1951358 species 
1951359 species 
1951360 species 
1951361 species 
1951362 species 
1951363 species 
1951364 species 
1951365 species 
1951366 species 
1951367 species 
1951368 species 
1951369 species 
1951370 species 
1951371 species 
1951372 species 
1951373 species 
1951374 species 
1951375 species 
1951376 species 
2562705 species 

if [ -f "/work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt" ]; then
rm /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt -f
fi
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951347\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951348\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951349\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951350\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951351\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951352\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951353\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951354\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951355\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951356\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951357\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951358\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951359\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951360\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951361\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951362\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951363\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951364\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951365\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951366\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951367\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951368\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951369\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951370\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951371\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951372\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951373\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951374\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951375\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^1951376\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt
cat /work/workspace/zhurj/reference/NCBI/taxonomy/20200219/names.dmp | grep -P "^2562705\t" | grep -P "scientific name" | awk -F "\t" '{print $1"\t"$3}' >> /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/genusSpecies/output/species_norank_1951308.txt

1951347 Akkermansiaceae bacterium UBA1009
1951348 Akkermansiaceae bacterium UBA1011
1951349 Akkermansiaceae bacterium UBA1090
1951350 Akkermansiaceae bacterium UBA1315
1951351 Akkermansiaceae bacterium UBA1506
1951352 Akkermansiaceae bacterium UBA1507
1951353 Akkermansiaceae bacterium UBA1518
1951354 Akkermansiaceae bacterium UBA1527
1951355 Akkermansiaceae bacterium UBA1955
1951356 Akkermansiaceae bacterium UBA1982
1951357 Akkermansiaceae bacterium UBA2367
1951358 Akkermansiaceae bacterium UBA2381
1951359 Akkermansiaceae bacterium UBA3157
1951360 Akkermansiaceae bacterium UBA3385
1951361 Akkermansiaceae bacterium UBA4145
1951362 Akkermansiaceae bacterium UBA4581
1951363 Akkermansiaceae bacterium UBA4589
1951364 Akkermansiaceae bacterium UBA5020
1951365 Akkermansiaceae bacterium UBA5688
1951366 Akkermansiaceae bacterium UBA5689
1951367 Akkermansiaceae bacterium UBA6138
1951368 Akkermansiaceae bacterium UBA6541
1951369 Akkermansiaceae bacterium UBA6946
1951370 Akkermansiaceae bacterium UBA7329
1951371 Akkermansiaceae bacterium UBA7331
1951372 Akkermansiaceae bacterium UBA7421
1951373 Akkermansiaceae bacterium UBA7811
1951374 Akkermansiaceae bacterium UBA7890
1951375 Akkermansiaceae bacterium UBA956
1951376 Akkermansiaceae bacterium UBA985
2562705 Akkermansiaceae bacterium

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/phylotree20200304/Akk20200310/input -b /work/workspace/zhurj/project/2_swgs/phylotree20200304/Akk20200310/output/Akkermansia -t 16

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/Akk20200310/output/Akkermansia.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/Akk20200310/output/Akkermansia_rec.png", w=500, units="mm", tree_style=ts)

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/Akk20200310/output/Akkermansia.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.mode = "c"
ts.arc_start = -180 # 0 degrees = 3 o'clock
ts.arc_span = 270
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/Akk20200310/output/Akkermansia_cir.svg", w=500, h=500, units="mm", tree_style=ts)


/work/database/ncbi/genome/all/GCA/GCA_000/GCA_000020/GCA_000020225/genomic.fna.gz
/work/database/ncbi/genome/all/GCA/GCA_900/GCA_900184/GCA_900184975/genomic.fna.gz

/work/assembly/current/MNH/MNH192/MNH19250/genomic.fna
/work/assembly/current/MNH/MNH206/MNH20651/genomic.fna

fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/MNH39848_317.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Enterococcusgallinarum200309/input/GCA_000157255_Sca316.fna -o /work/workspace/zhurj/project/2_swgs/run00050/fastqANIEnterococcusgallinarum200309/MNH39848_GCA_000157255_Sca316
fastANI -t 16 -q /work/assembly/current/MNH/MNH192/MNH19250/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000020/GCA_000020225/genomic.fna.gz -o /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/fastANI/MNH19250_282_GCA_000020225_Amu266

fastANI -t 16 -q /work/assembly/current/MNH/MNH206/MNH20651/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000020/GCA_000020225/genomic.fna.gz -o /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/fastANI/MNH20651_313_GCA_000020225_Amu266
2020.Mar11-12:00
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/checkM/MNH20651/ssu.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/checkM/MNH20651/MNH20651.out -max_target_seqs 5 -num_threads 16 -outfmt 7 

# BLASTN 2.9.0+
# BLASTN 2.9.0+
# Query: genomic&&NODE_18_length_5761_cov_1362.738037
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
genomic&&NODE_18_length_5761_cov_1362.738037    PJKF01000002    100.000 1434    0       0       37      1470    1       1434    0.0     2587
genomic&&NODE_18_length_5761_cov_1362.738037    PJKB01000002    99.582  1434    6       0       37      1470    1       1434    0.0     2560
genomic&&NODE_18_length_5761_cov_1362.738037    CP001071        99.024  1434    12      2       37      1470    1       1432    0.0     2516
genomic&&NODE_18_length_5761_cov_1362.738037    KQ968618        98.047  1434    28      0       37      1470    1       1434    0.0     2461
genomic&&NODE_18_length_5761_cov_1362.738037    KT254068        94.024  1439    78      6       37      1470    1       1436    0.0     2188
# BLAST processed 1 queries


cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "KQ968618|CP001071|PJKB01000002|PJKF01000002"
CP001071        Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;Akkermansia muciniphila;ATCC BAA-835(T)
KQ968618        Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;KLE1797
PJKB01000002    Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;PJKB_s;GP22
PJKF01000002    Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;PJKF_s;GP15

echo "# GCA2genome" >> README.md
git init
git add README.md
git commit -m "first commit"
git remote add origin git@github.com:server2-2019/GCA2genome.git
git push -u origin master

2020.Mar12-13:35
cat names.txt | taxonkit name2taxid | csvtk pretty -t
cat names.txt | taxonkit name2taxid --show-rank | csvtk pretty -t
cat names.txt | taxonkit name2taxid | taxonkit lineage --taxid-field 2
taxonkit list --show-rank --show-name --ids 239935

echo "Christensenella minuta" | taxonkit name2taxid
echo "Christensenella minuta" | taxonkit name2taxid | taxonkit lineage --taxid-field 2
Christensenella minuta  626937  cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;Christensenella;Christensenella minuta
Clostridia
Christensenellaceae, Clostridiaceae, Dehalobacteriaceae
Christensenellaceae     990719  cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Clostridia;Clostridiales;Christensenellaceae
Clostridiaceae  31979   cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Clostridia;Clostridiales;Clostridiaceae

echo "Oscillospira" | taxonkit name2taxid | taxonkit lineage --taxid-field 2
Oscillospira    119852  cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Oscillospira


echo "Clostridiales" | taxonkit name2taxid | taxonkit lineage --taxid-field 2
Clostridiales   186802  cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Clostridia;Clostridiales
Clostridia      186801  cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Clostridia
taxonkit list --show-rank --show-name --ids 186802 | grep -P "\[family\]" 
  31979 [family] Clostridiaceae
  31984 [family] Heliobacteriaceae
  68298 [family] Syntrophomonadaceae
  186803 [family] Lachnospiraceae ************************
  186804 [family] Peptostreptococcaceae
  186806 [family] Eubacteriaceae
  186807 [family] Peptococcaceae
  216572 [family] Oscillospiraceae
  424536 [family] Catabacteriaceae
    539000 [family] Clostridiales Family XVII. Incertae Sedis  ************************
    543313 [family] Clostridiales Family XII. Incertae Sedis ************************
    543314 [family] Clostridiales Family XIII. Incertae Sedis ************************
    543317 [family] Clostridiales Family XIV. Incertae Sedis ************************
    543347 [family] Clostridiales Family XVI. Incertae Sedis ************************
    543350 [family] Clostridiales Family XIX. Incertae Sedis ************************
    1689146 [family] Clostridiales Family III. Incertae Sedis ************************
    1689151 [family] Clostridiales Family IV. Incertae Sedis ************************
  541000 [family] Ruminococcaceae ************************
  541019 [family] Gracilibacteraceae
  543349 [family] Symbiobacteriaceae
  715221 [family] Caldicoprobacteraceae
  990719 [family] Christensenellaceae
  1185407 [family] Defluviitaleaceae
  1491775 [family] Proteinivoraceae
  2304686 [family] Hungateiclostridiaceae
  2603322 [family] Vallitaleaceae

echo "Enterobacteriaceae" | taxonkit name2taxid | taxonkit lineage --taxid-field 2
Methanobrevibacter      2172    cellular organisms;Archaea;Euryarchaeota;Methanomada group;Methanobacteria;Methanobacteriales;Methanobacteriaceae;Methanobrevibacter
echo "Methanobrevibacter" | taxonkit name2taxid | taxonkit lineage --taxid-field 2
family
Enterobacteriaceae      543     cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae



域(Domain)、界（Kingdom）、门（Phylum）、纲（Class）、目（Order）、科（Family）、属（Genus）、种（Species）
taxonkit list --show-rank --show-name --ids 186801
taxonkit list --show-rank --show-name --ids 186801 | grep -P "order" 
echo "Christensenellaceae" | taxonkit name2taxid | taxonkit lineage --taxid-field 2
echo "Clostridiaceae" | taxonkit name2taxid | taxonkit lineage --taxid-field 2
echo "Dehalobacteriaceae" | taxonkit name2taxid | taxonkit lineage --taxid-field 2

Bacteroidetes, Firmicutes, Staphylococcus spp., Akkermansia muciniphila and methanogens were significantly higher in lean volunteers
Ruminococcus spp., Christensenella minuta, γ-Proteobacteria and A. muciniphila were more abundant in overweight as opposed to obese individuals.

echo "Clostridiales" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t
echo "Christensenellaceae" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t
echo "Christensenella" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t
taxonkit list --show-rank --show-name --ids 990719 | grep -P "species" 
小蔡:
特别是与BMI呈负相关的，Firmicutes, clostridia,clostridiales
clostridia这个是界门纲目clostridiales
626937 [species] Christensenella minuta
  1545742 [species] uncultured Christensenella sp.
1805714 [species] Christensenella massiliensis
1816678 [species] Christensenella timonensis
  1851429 [species] Christensenella sp. AF73-05CM02
  1935934 [species] Christensenella sp.
  2086585 [species] Christensenella sp. Marseille-P3954
1229255 [species] uncultured Christensenellaceae bacterium
1930012 [species] Christensenellaceae bacterium Phil1
1930013 [species] Christensenellaceae bacterium Phil7
1930014 [species] Christensenellaceae bacterium Phil10
2054177 [species] Christensenellaceae bacterium
1917875 [species] Beduinibacterium massiliense

if [ -f "/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/mid/gca20200316.txt" ];then
rm /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/mid/gca20200316.txt -f
fi

cat /work/workspace/zhurj/reference/NCBI/prokaryotes/20200305/prokaryotes.txt | grep -P "\t626937\t" | awk -F "\t" '{print $1"\t"$2"\t"$19"\t"$7"\t"$8"\t"$16"\t"$20}'  >> /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/mid/gca20200316.txt
cat /work/workspace/zhurj/reference/NCBI/prokaryotes/20200305/prokaryotes.txt | grep -P "\t1545742\t" | awk -F "\t" '{print $1"\t"$2"\t"$19"\t"$7"\t"$8"\t"$16"\t"$20}' >> /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/mid/gca20200316.txt
cat /work/workspace/zhurj/reference/NCBI/prokaryotes/20200305/prokaryotes.txt | grep -P "\t1805714\t" | awk -F "\t" '{print $1"\t"$2"\t"$19"\t"$7"\t"$8"\t"$16"\t"$20}' >> /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/mid/gca20200316.txt
cat /work/workspace/zhurj/reference/NCBI/prokaryotes/20200305/prokaryotes.txt | grep -P "\t1816678\t" | awk -F "\t" '{print $1"\t"$2"\t"$19"\t"$7"\t"$8"\t"$16"\t"$20}' >> /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/mid/gca20200316.txt
cat /work/workspace/zhurj/reference/NCBI/prokaryotes/20200305/prokaryotes.txt | grep -P "\t1851429\t" | awk -F "\t" '{print $1"\t"$2"\t"$19"\t"$7"\t"$8"\t"$16"\t"$20}' >> /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/mid/gca20200316.txt
cat /work/workspace/zhurj/reference/NCBI/prokaryotes/20200305/prokaryotes.txt | grep -P "\t1935934\t" | awk -F "\t" '{print $1"\t"$2"\t"$19"\t"$7"\t"$8"\t"$16"\t"$20}' >> /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/mid/gca20200316.txt
cat /work/workspace/zhurj/reference/NCBI/prokaryotes/20200305/prokaryotes.txt | grep -P "\t2086585\t" | awk -F "\t" '{print $1"\t"$2"\t"$19"\t"$7"\t"$8"\t"$16"\t"$20}' >> /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/mid/gca20200316.txt
cat /work/workspace/zhurj/reference/NCBI/prokaryotes/20200305/prokaryotes.txt | grep -P "\t1229255\t" | awk -F "\t" '{print $1"\t"$2"\t"$19"\t"$7"\t"$8"\t"$16"\t"$20}' >> /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/mid/gca20200316.txt
cat /work/workspace/zhurj/reference/NCBI/prokaryotes/20200305/prokaryotes.txt | grep -P "\t1930012\t" | awk -F "\t" '{print $1"\t"$2"\t"$19"\t"$7"\t"$8"\t"$16"\t"$20}' >> /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/mid/gca20200316.txt
cat /work/workspace/zhurj/reference/NCBI/prokaryotes/20200305/prokaryotes.txt | grep -P "\t1930013\t" | awk -F "\t" '{print $1"\t"$2"\t"$19"\t"$7"\t"$8"\t"$16"\t"$20}' >> /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/mid/gca20200316.txt
cat /work/workspace/zhurj/reference/NCBI/prokaryotes/20200305/prokaryotes.txt | grep -P "\t1930014\t" | awk -F "\t" '{print $1"\t"$2"\t"$19"\t"$7"\t"$8"\t"$16"\t"$20}' >> /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/mid/gca20200316.txt
cat /work/workspace/zhurj/reference/NCBI/prokaryotes/20200305/prokaryotes.txt | grep -P "\t2054177\t" | awk -F "\t" '{print $1"\t"$2"\t"$19"\t"$7"\t"$8"\t"$16"\t"$20}' >> /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/mid/gca20200316.txt
cat /work/workspace/zhurj/reference/NCBI/prokaryotes/20200305/prokaryotes.txt | grep -P "\t1917875\t" | awk -F "\t" '{print $1"\t"$2"\t"$19"\t"$7"\t"$8"\t"$16"\t"$20}' >> /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/mid/gca20200316.txt


/work/workspace/zhurj/script/git/GCA2genome/gca2genome.py -i /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/mid/gcainput.txt -o /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/input -f /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/mid/gcawithgenome.txt


/work/workspace/zhurj/script/git/mnid2genome/mnid2genome.py -i /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/mid/cmsample.txt -o /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/input -f /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/mid/mnidwithgenome.txt


/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/input -b /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/output/Christensenellaceae -t 16

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/output/Christensenellaceae.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/output/Christensenellaceae_rec.png", w=500, units="mm", tree_style=ts)

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/output/Christensenellaceae.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.mode = "c"
ts.arc_start = -180 # 0 degrees = 3 o'clock
ts.arc_span = 270
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/output/Christensenellaceae_cir.svg", w=500, h=500, units="mm", tree_style=ts)


from ete3 import Tree, TextFace, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/output/Christensenellaceae.nwk'
t = Tree(file)
D_leaf_color = {"MNH04863":"red","GCA_001652705":"Lime","GCA_003628755":"Lime","GCA_001571425":"Lime","GCA_001678855":"Lime","GCA_902388075":"Lime","GCA_900155415":"Indigo","GCA_900087015":"Cyan","GCA_902376065":"Cyan","GCA_001678845":"DarkMagenta","GCA_900604345":"DarkMagenta","GCA_001940855":"SteelBlue","GCA_001940845":"SteelBlue","GCA_001940895":"SteelBlue","GCA_902386865":"SteelBlue"}
for node in t.traverse():
    # Hide node circles
    node.img_style['size'] = 0
    if node.is_leaf():
        color = D_leaf_color.get(node.name, None)
        if color:
            name_face = TextFace(node.name, fgcolor=color, fsize=8)
            node.add_face(name_face, column=0, position='branch-right')           
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/output/Christensenellaceae_rec.png", w=500, units="mm", tree_style=ts)

from ete3 import Tree, TextFace, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/output/Christensenellaceae.nwk'
t = Tree(file)
D_leaf_color = {"MNH04863":"red","GCA_001652705":"Lime","GCA_003628755":"Lime","GCA_001571425":"Lime","GCA_001678855":"Lime","GCA_902388075":"Lime","GCA_900155415":"Indigo","GCA_900087015":"Cyan","GCA_902376065":"Cyan","GCA_001678845":"DarkMagenta","GCA_900604345":"DarkMagenta","GCA_001940855":"SteelBlue","GCA_001940845":"SteelBlue","GCA_001940895":"SteelBlue","GCA_902386865":"SteelBlue"}
for node in t.traverse():
    # Hide node circles
    node.img_style['size'] = 0
    if node.is_leaf():
        color = D_leaf_color.get(node.name, None)
        if color:
            name_face = TextFace(node.name, fgcolor=color, fsize=8)
            node.add_face(name_face, column=0, position='branch-right')       
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True
ts.mode = "c"
ts.arc_start = -180 # 0 degrees = 3 o'clock
ts.arc_span = 270
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/output/Christensenellaceae_cir.svg", w=500, h=500, units="mm", tree_style=ts)


fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/input/MNH04863.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/input/GCA_001678845.fna -o /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/fastANI/MNH04863_GCA_001678845
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/input/MNH04863.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/input/GCA_902386865.fna -o /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/fastANI/MNH04863_GCA_902386865

fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/input/GCA_001678845.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/input/GCA_902386865.fna -o /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/fastANI/GCA_001678845_GCA_902386865

checkm ssu_finder /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/input/MNH04863.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/input /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/checkM/MNH04863 -x fna -t 20
checkm ssu_finder /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/input/GCA_001678845.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/input /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/checkM/GCA_001678845 -x fna -t 20
checkm ssu_finder /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/input/GCA_902386865.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/input /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/checkM/GCA_902386865 -x fna -t 20

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/sanger/MNH04863_all.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/blastn/MNH04863/MNH04863.out -max_target_seqs 5 -num_threads 16 -outfmt 7 
# BLASTN 2.9.0+
# Query: MNH04863
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
MNH04863        MAIQ01000023    100.000 1370    0       0       1       1370    35      1404    0.0     2471
MNH04863        KQ965520        98.686  1370    18      0       1       1370    35      1404    0.0     2390
MNH04863        LAYJ01000068    97.299  1370    37      0       1       1370    35      1404    0.0     2305
MNH04863        FLKP01000002    96.934  1370    42      0       1       1370    35      1404    0.0     2282
MNH04863        LT700187        96.861  1370    43      0       1       1370    35      1404    0.0     2278
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "MAIQ01000023|KQ965520"
KQ965520        Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;Christensenella;Christensenella minuta;DSM 22607(T)
MAIQ01000023    Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;Christensenella;MAIQ_s;AF73-05CM02


/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/checkM/MNH04863/ssu.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/checkM/MNH04863/MNH04863.out -max_target_seqs 5 -num_threads 16 -outfmt 7 
# BLASTN 2.9.0+
# Query: MNH04863&&MNH04863_15
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
MNH04863&&MNH04863_15   MAIQ01000023    100.000 1458    0       0       32      1489    1       1458    0.0     2630
MNH04863&&MNH04863_15   KQ965520        98.765  1458    18      0       32      1489    1       1458    0.0     2549
MNH04863&&MNH04863_15   LAYJ01000068    97.462  1458    37      0       32      1489    1       1458    0.0     2463
MNH04863&&MNH04863_15   FLKP01000002    97.119  1458    42      0       32      1489    1       1458    0.0     2441
MNH04863&&MNH04863_15   LT700187        96.982  1458    44      0       32      1489    1       1458    0.0     2432

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "MAIQ01000023|KQ965520"
KQ965520        Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;Christensenella;Christensenella minuta;DSM 22607(T)
MAIQ01000023    Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;Christensenella;MAIQ_s;AF73-05CM02

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/checkM/GCA_001678845/ssu.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/checkM/GCA_001678845/GCA_001678845.out -max_target_seqs 5 -num_threads 16 -outfmt 7 
# BLASTN 2.9.0+
# Query: GCA_001678845&&MAIQ01000023.1
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
GCA_001678845&&MAIQ01000023.1   MAIQ01000023    100.000 1458    0       0       32      1489    1       1458    0.0     2630
GCA_001678845&&MAIQ01000023.1   KQ965520        98.765  1458    18      0       32      1489    1       1458    0.0     2549
GCA_001678845&&MAIQ01000023.1   LAYJ01000068    97.462  1458    37      0       32      1489    1       1458    0.0     2463
GCA_001678845&&MAIQ01000023.1   FLKP01000002    97.119  1458    42      0       32      1489    1       1458    0.0     2441
GCA_001678845&&MAIQ01000023.1   LT700187        96.982  1458    44      0       32      1489    1       1458    0.0     2432

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "MAIQ01000023|KQ965520"
KQ965520        Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;Christensenella;Christensenella minuta;DSM 22607(T)
MAIQ01000023    Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;Christensenella;MAIQ_s;AF73-05CM02

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/checkM/GCA_902386865/ssu.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200316/checkM/GCA_902386865/GCA_902386865.out -max_target_seqs 5 -num_threads 16 -outfmt 7
# BLASTN 2.9.0+
# Query: GCA_902386865&&CABMKF010000023.1
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
GCA_902386865&&CABMKF010000023.1        MAIQ01000023    100.000 1458    0       0       32      1489    1       1458    0.0     2630
GCA_902386865&&CABMKF010000023.1        KQ965520        98.765  1458    18      0       32      1489    1       1458    0.0     2549
GCA_902386865&&CABMKF010000023.1        LAYJ01000068    97.462  1458    37      0       32      1489    1       1458    0.0     2463
GCA_902386865&&CABMKF010000023.1        FLKP01000002    97.119  1458    42      0       32      1489    1       1458    0.0     2441
GCA_902386865&&CABMKF010000023.1        LT700187        96.982  1458    44      0       32      1489    1       1458    0.0     2432
# BLAST processed 1 queries

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "MAIQ01000023|KQ965520"
KQ965520        Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;Christensenella;Christensenella minuta;DSM 22607(T)
MAIQ01000023    Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;Christensenella;MAIQ_s;AF73-05CM02



echo "Actinobacteria" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t
echo "Arcanobacterium" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t
taxonkit list --show-rank --show-name --ids 201174 | grep -P "species" 
taxonkit list --show-rank --show-name --ids 28263 | grep -P "species" 
taxonkit list --ids 28263 | grep -P "species"
echo "Alicyclobacillus ferrooxydans" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t 
echo "Bacillus" | taxonkit name2taxid | taxonkit lineage --taxid-field 2

taxonkit list --show-rank --show-name --indent "" --ids 214904 | grep -P "species"

/work/workspace/zhurj/script/git/bacteriaGramdb/bacteriaGramdb.py 

import os
import re

var=os.popen('taxonkit list --show-rank --show-name --indent "" --ids 214904 | grep -P "species"').read()
list_ = []
list_ = var.split("\n")
_SPECIES = re.compile(r"\s\[species\]\s")

for line in list_:
  if _SPECIES.search(line):
    print(line + "***")
    [taxid,name] = line.strip().split(r" [species] ")
    print('{}'.format(taxid))

import os
path1=os.path.abspath('.')   #表示当前所处的文件夹的绝对路径
path2=os.path.abspath('..')  #表示当前所处的文件夹上一级文件夹的绝对路径
os.path.abspath(os.path.dirname(__file__))

nohup /work/workspace/zhurj/script/git/bacteriaGramdb/bacteriaGramdb.py -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaGramdb -n bacteriaGramResult > /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaGramdb/bacteriaGram.log &


echo "Clostridium amygdalinum" | taxonkit name2taxid

/work/workspace/zhurj/script/git/mnid2genome/mnid2genome.py -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/tmp -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/output/mnid2genome46.txt -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/input/mnid.txt


checkm ssu_finder /work/assembly/current/MNH/MNH056/MNH05676/genomic.fna /work/assembly/current/MNH/MNH056/MNH05676 /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/test -x fna -t 16
checkm lineage_wf /work/assembly/current/MNH/MNH056/MNH05676 /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/test1 -x fna -t 16 --tab_table -f /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/test1/test1_checkm.txt
checkm tree_qa /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/test1 -f /work/workspace/zhurj/project/2_swgs/Akkermansia20200310/test1/treeqa.txt --tab_table



echo "Clostridiales" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 | cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t

echo "Clostridiales" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 | cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t
cat /work/workspace/zhurj/project/4_tmp/taxonkit/name.txt | taxonkit name2taxid | taxonkit lineage --taxid-field 2 | cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t

(python3.6) [zhurj@mnhead taxonkit]$ cat /work/workspace/zhurj/project/4_tmp/taxonkit/name.txt | taxonkit name2taxid | taxonkit lineage --taxid-field 2 | cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t
taxid    kindom     phylum           class                 order             family             genus           species
1263     Bacteria   Firmicutes       Clostridia            Clostridiales     Ruminococcaceae    Ruminococcus    
33042    Bacteria   Firmicutes       Clostridia            Clostridiales     Lachnospiraceae    Coprococcus     
1730     Bacteria   Firmicutes       Clostridia            Clostridiales     Eubacteriaceae     Eubacterium     
189330   Bacteria   Firmicutes       Clostridia            Clostridiales     Lachnospiraceae    Dorea           
976      Bacteria   Bacteroidetes                                                                               
200643   Bacteria   Bacteroidetes    Bacteroidia                                                                
816      Bacteria   Bacteroidetes    Bacteroidia           Bacteroidales     Bacteroidaceae     Bacteroides     
1654     Bacteria   Actinobacteria   Actinobacteria        Actinomycetales   Actinomycetaceae   Actinomyces     
729      Bacteria   Proteobacteria   Gammaproteobacteria   Pasteurellales    Pasteurellaceae    Haemophilus     Haemophilus parainfluenzae
39491    Bacteria   Firmicutes       Clostridia            Clostridiales     Lachnospiraceae                    [Eubacterium] rectale
1301     Bacteria   Firmicutes       Bacilli               Lactobacillales   Streptococcaceae   Streptococcus   

f__Coriobacteriaceae
echo "Coriobacteriaceae" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t 
taxonkit list --show-rank --show-name --indent "" --ids 84107 | grep -P "species"

/work/workspace/zhurj/script/git/mnid2genome/mnid2genome.py -l MNH36311 -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/f__Coriobacteriaceae_MNH36311.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f__Coriobacteriaceae/input

/work/workspace/zhurj/script/git/GCA2genome/gca2genome.py -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/input/f__Coriobacteriaceae_150GCA.txt -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/f__Coriobacteriaceae_GCAgenome.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f__Coriobacteriaceae/input

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f__Coriobacteriaceae/input -b /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f__Coriobacteriaceae/output/Christensenellaceae -t 16

fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f__Coriobacteriaceae/input/MNH36311.fna -r /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f__Coriobacteriaceae/input/GCA_002320125.fna -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/MNH36311/MNH36311_GCA_002320125_201
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f__Coriobacteriaceae/input/MNH36311.fna -r /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f__Coriobacteriaceae/input/GCA_002391315.fna -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/MNH36311/MNH36311_GCA_002391315_209
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f__Coriobacteriaceae/input/MNH36311.fna -r /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f__Coriobacteriaceae/input/GCA_002396135.fna -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/MNH36311/MNH36311_GCA_002396135_174
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f__Coriobacteriaceae/input/MNH36311.fna -r /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f__Coriobacteriaceae/input/GCA_003452675.fna -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/MNH36311/MNH36311_GCA_003452675_160
if [ -f "/work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/MNH36311/merge.txt" ];then
rm -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/MNH36311/merge.txt
fi
cat /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/MNH36311/MNH36311_GCA_002320125_201 | awk -F "\t" '{print "MNH36311\tGCA_002320125_201\t"$(NF-2)"\t"$(NF-1)"\t"$NF}' >> /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/MNH36311/merge.txt
cat /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/MNH36311/MNH36311_GCA_002391315_209 | awk -F "\t" '{print "MNH36311\tGCA_002391315_209\t"$(NF-2)"\t"$(NF-1)"\t"$NF}' >> /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/MNH36311/merge.txt
cat /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/MNH36311/MNH36311_GCA_002396135_174 | awk -F "\t" '{print "MNH36311\tGCA_002396135_174\t"$(NF-2)"\t"$(NF-1)"\t"$NF}' >> /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/MNH36311/merge.txt
cat /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/MNH36311/MNH36311_GCA_003452675_160 | awk -F "\t" '{print "MNH36311\tGCA_003452675_160\t"$(NF-2)"\t"$(NF-1)"\t"$NF}' >> /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/MNH36311/merge.txt
check GCA_002320125
mkdir /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_002320125
checkm ssu_finder /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f__Coriobacteriaceae/input/GCA_002320125.fna /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f__Coriobacteriaceae/input /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_002320125 -x fna -t 16
checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/all/GCA/GCA_002/GCA_002320/GCA_002320125 /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_002320125 -x fna -t 16 --tab_table -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_002320125/GCA_002320125_checkm.txt

checkm tree_qa /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_002320125 -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_002320125/treeqa.txt --tab_table




g__Butyricimonas
echo "Butyricimonas" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t 

taxid    kindom     phylum          class         order           family             genus           species
574697   Bacteria   Bacteroidetes   Bacteroidia   Bacteroidales   Odoribacteraceae   Butyricimonas

taxonkit list --show-rank --show-name --indent "" --ids 574697 | grep -P "species" 

/work/workspace/zhurj/script/git/mnid2genome/mnid2genome.py -l MNH19394 -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/g__Butyricimonas_MNH19394.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Butyricimonas/input

/work/workspace/zhurj/script/git/GCA2genome/gca2genome.py -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/input/g__Butyricimonas_23GCA.txt -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/g__Butyricimonas_GCAgenome.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Butyricimonas/input

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Butyricimonas/input -b /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Butyricimonas/output/Butyricimonas -t 16

mkdir -p /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/MNH19394
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Butyricimonas/input/MNH19394.fna -r /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Butyricimonas/input/GCA_900547335.fna -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/MNH19394/MNH19394_GCA_900547335_465

mkdir -p /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_900547335
checkm ssu_finder /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Butyricimonas/input/GCA_900547335.fna /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Butyricimonas/input /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_900547335 -x fna -t 16

checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/all/GCA/GCA_900/GCA_900547/GCA_900547335 /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_900547335 -x fna -t 16 --tab_table -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_900547335/GCA_900547335_checkm.txt


g__Parabacteroides
echo "Parabacteroides" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t 

taxid    kindom     phylum          class         order           family           genus             species
375288   Bacteria   Bacteroidetes   Bacteroidia   Bacteroidales   Tannerellaceae   Parabacteroides

taxonkit list --show-rank --show-name --indent "" --ids 375288 | grep -P "species"
cat /work/workspace/zhurj/reference/NCBI/prokaryotes/20200305/prokaryotes.txt | grep -P "\t2211362\t" | awk -F "\t" '{print $1"\t"$2"\t"$16"\t"$19"\t"$20"\t"$7"\t"$8}' >> /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/taxonkit/g__Butyricimonas.txt

/work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Parabacteroides

/work/workspace/zhurj/script/git/mnid2genome/mnid2genome.py -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/input/g__Parabacteroides_4MNH.txt -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/g__Parabacteroides_4MNH.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Parabacteroides/input

/work/workspace/zhurj/script/git/GCA2genome/gca2genome.py -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/input/g__Parabacteroides_50GCA.txt -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/g__Parabacteroides_GCAgenome.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Parabacteroides/input

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Parabacteroides/input -b /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Parabacteroides/output/Parabacteroides -t 16

mkdir /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/g__Parabacteroides
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Parabacteroides/input/MNH09207.fna -r /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Parabacteroides/input/GCA_000969845.fna -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/g__Parabacteroides/MNH09207_GCA_000969845
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Parabacteroides/input/MNH08509.fna -r /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Parabacteroides/input/GCA_003439895.fna -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/g__Parabacteroides/MNH08509_GCA_003439895
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Parabacteroides/input/MNH09691.fna -r /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Parabacteroides/input/GCA_003439895.fna -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/g__Parabacteroides/MNH09691_GCA_003439895

cat /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/g__Parabacteroides/MNH09207_GCA_000969845 | awk '{print "MNH09207\tGCA_000969845\t"$(NF-2)"\t"$(NF-1)"\t"$NF}'
cat /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/g__Parabacteroides/MNH08509_GCA_003439895 | awk '{print "MNH08509\tGCA_003439895\t"$(NF-2)"\t"$(NF-1)"\t"$NF}'
cat /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/g__Parabacteroides/MNH09691_GCA_003439895 | awk '{print "MNH09691\tGCA_003439895\t"$(NF-2)"\t"$(NF-1)"\t"$NF}'

checkm ssu_finder /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f__Coriobacteriaceae/input/GCA_002320125.fna /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f__Coriobacteriaceae/input /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_002320125 -x fna -t 16
mkdir /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_000969845
checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/all/GCA/GCA_000/GCA_000969/GCA_000969845 /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_000969845 -x fna -t 16 --tab_table -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_000969845/GCA_000969845_checkm.txt
mkdir /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_003439895
checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/all/GCA/GCA_003/GCA_003439/GCA_003439895 /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_003439895 -x fna -t 16 --tab_table -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_003439895/GCA_003439895_checkm.txt

g__Alistipes
echo "Alistipes" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t 

taxid    kindom     phylum          class         order           family          genus       species
239759   Bacteria   Bacteroidetes   Bacteroidia   Bacteroidales   Rikenellaceae   Alistipes

taxonkit list --show-rank --show-name --indent "" --ids 239759 | grep -P "species"
cat /work/workspace/zhurj/reference/NCBI/prokaryotes/20200305/prokaryotes.txt | grep -P "\t2211362\t" | awk -F "\t" '{print $1"\t"$2"\t"$16"\t"$19"\t"$20"\t"$7"\t"$8}' >> /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/taxonkit/g__Alistipes.txt


/work/workspace/zhurj/script/git/mnid2genome/mnid2genome.py -l MNH05020 -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/g__Alistipes_MNH.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Alistipes/input

/work/workspace/zhurj/script/git/GCA2genome/gca2genome.py -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/input/g__Alistipes_65GCA.txt -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/g__Alistipes_GCAgenome.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Alistipes/input

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Alistipes/input -b /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Alistipes/output/Alistipes -t 16

mkdir /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/g__Alistipes
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Alistipes/input/MNH05020.fna -r /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g__Alistipes/input/GCA_003477565.fna -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/g__Alistipes/MNH05020_GCA_003477565

mkdir /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_003477565
checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/all/GCA/GCA_003/GCA_003477/GCA_003477565 /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_003477565 -x fna -t 16 --tab_table -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_003477565/GCA_003477565_checkm.txt


s__Enterococcus_faecium
echo "Enterococcus faecium" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t 

taxid   kindom     phylum       class     order             family            genus          species
1352    Bacteria   Firmicutes   Bacilli   Lactobacillales   Enterococcaceae   Enterococcus   Enterococcus faecium

taxonkit list --show-rank --show-name --indent "" --ids 1352 | grep -P "species"
cat /work/workspace/zhurj/reference/NCBI/prokaryotes/20200305/prokaryotes.txt | grep -P "\t1352\t" | awk -F "\t" '{print $1"\t"$2"\t"$16"\t"$19"\t"$20"\t"$7"\t"$8}' >> /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/taxonkit/s__Enterococcus_faecium.txt


/work/workspace/zhurj/script/git/mnid2genome/mnid2genome.py -l MNH31116 -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/s__EnterococcusFaecium_MNH.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/s__EnterococcusFaecium/input

/work/workspace/zhurj/script/git/GCA2genome/gca2genome.py -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/input/s__EnterococcusFaecium_38GCA.txt -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/s__EnterococcusFaecium_GCAgenome.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/s__EnterococcusFaecium/input

GCA 36
/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/s__EnterococcusFaecium/input -b /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/s__EnterococcusFaecium/output/EnterococcusFaecium -t 16

mkdir /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/s__EnterococcusFaecium
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/s__EnterococcusFaecium/input/MNH31116.fna -r /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/s__EnterococcusFaecium/input/GCA_002140455.fna -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/s__EnterococcusFaecium/MNH31116_GCA_002140455

mkdir /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_002140455
checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/all/GCA/GCA_002/GCA_002140/GCA_002140455 /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_002140455 -x fna -t 16 --tab_table -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_002140455/GCA_002140455_checkm.txt

o__Clostridiales
echo "Clostridiales" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t 

taxid    kindom     phylum       class        order           family   genus   species
186802   Bacteria   Firmicutes   Clostridia   Clostridiales


taxonkit list --show-rank --show-name --indent "" --ids 186802 | grep -P "species"  > /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/taxonkit/species/o__Clostridiales_species.txt
cat /work/workspace/zhurj/reference/NCBI/prokaryotes/20200305/prokaryotes.txt | grep -P "\t1352\t" | awk -F "\t" '{print $1"\t"$2"\t"$16"\t"$19"\t"$20"\t"$7"\t"$8}' >> /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/taxonkit/s__Enterococcus_faecium.txt


/work/workspace/zhurj/script/git/mnid2genome/mnid2genome.py -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/input/o__Clostridiales_b1_4MNH.txt -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/o__Clostridiales_MNH.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/o__Clostridiales/sample4GC54len31/input

/work/workspace/zhurj/script/git/GCA2genome/gca2genome.py -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/input/o__Clostridiales_b1_48GCA.txt -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/o__Clostridiales_GCAgenome.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/o__Clostridiales/sample4GC54len31/input

GCA30


MNH04831 已确定
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH04831/ssu.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH04831/MNH04831.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&MNH04831_17    MAIQ01000023    100.000 1458    0       0       56      1513    1458    1       0.0     2630
genomic&&MNH04831_17    KQ965520        98.765  1458    18      0       56      1513    1458    1       0.0     2549
genomic&&MNH04831_17    LAYJ01000068    97.462  1458    37      0       56      1513    1458    1       0.0     2463
genomic&&MNH04831_17    FLKP01000002    97.119  1458    42      0       56      1513    1458    1       0.0     2441
genomic&&MNH04831_17    LT700187        96.982  1458    44      0       56      1513    1458    1       0.0     2432
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "KQ965520|MAIQ01000023"
KQ965520        Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;Christensenella;Christensenella minuta;DSM 22607(T)
MAIQ01000023    Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;Christensenella;MAIQ_s;AF73-05CM02

mkdir /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_001678845
checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/all/GCA/GCA_001/GCA_001678/GCA_001678845 /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_001678845 -x fna -t 16 --tab_table -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_001678845/GCA_001678845_checkm.txt

MNH06365 - g_Christensenella
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH06365/ssu.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH06365/MNH06365.out -max_target_seqs 5 -num_threads 16 -outfmt 7
# 5 hits found
genomic&&MNH06365_19    KQ965520        97.668  1458    34      0       56      1513    1458    1       0.0     2477
genomic&&MNH06365_19    LT700187        97.531  1458    36      0       56      1513    1458    1       0.0     2468
genomic&&MNH06365_19    MAIQ01000023    97.462  1458    37      0       56      1513    1458    1       0.0     2463
genomic&&MNH06365_19    FLKP01000002    96.708  1458    48      0       56      1513    1458    1       0.0     2414
genomic&&MNH06365_19    LAYJ01000068    95.610  1458    64      0       56      1513    1458    1       0.0     2342
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "KQ965520|LT700187"
KQ965520        Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;Christensenella;Christensenella minuta;DSM 22607(T)
LT700187        Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;Christensenella;Christensenella massiliensis;Marseille-P2438(T)

MNH19729 - Christensenella
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH19729/ssu_filter.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH19729/MNH19729.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&NODE_1_length_765916_cov_501.972117    FLKP01000002    97.668  1458    32      2       34      1489    1       1458    0.0     2470
genomic&&NODE_1_length_765916_cov_501.972117    KQ965520        96.982  1458    42      2       34      1489    1       1458    0.0     2425
genomic&&NODE_1_length_765916_cov_501.972117    MAIQ01000023    96.710  1459    44      4       34      1489    1       1458    0.0     2401
genomic&&NODE_1_length_765916_cov_501.972117    LT700187        95.819  1459    57      4       34      1489    1       1458    0.0     2342
genomic&&NODE_1_length_765916_cov_501.972117    LAYJ01000068    95.000  1460    67      4       34      1489    1       1458    0.0     2292
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "FLKP01000002"
FLKP01000002    Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;Christensenella;Christensenella timonensis;Marseille-P2437(T)

MNH05089 - Christensenella
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH05089/ssu.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH05089/MNH05089.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&MNH05089_9     FLKP01000002    98.217  1458    26      0       56      1513    1458    1       0.0     2513
genomic&&MNH05089_9     KQ965520        97.942  1458    30      0       56      1513    1458    1       0.0     2495
genomic&&MNH05089_9     MAIQ01000023    97.805  1458    32      0       56      1513    1458    1       0.0     2486
genomic&&MNH05089_9     LT700187        96.776  1458    47      0       56      1513    1458    1       0.0     2418
genomic&&MNH05089_9     LAYJ01000068    96.776  1458    47      0       56      1513    1458    1       0.0     2418
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "FLKP01000002|KQ965520|MAIQ01000023"
FLKP01000002    Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;Christensenella;Christensenella timonensis;Marseille-P2437(T)
KQ965520        Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;Christensenella;Christensenella minuta;DSM 22607(T)
MAIQ01000023    Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;Christensenella;MAIQ_s;AF73-05CM02

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/o__Clostridiales/sample4GC54len31/input -b /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/o__Clostridiales/sample4GC54len31/output/Clostridiales -t 16

g_Christensenella
echo "Christensenella" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t 

taxid    kindom     phylum       class        order           family                genus             species
990721   Bacteria   Firmicutes   Clostridia   Clostridiales   Christensenellaceae   Christensenella

taxonkit list --show-rank --show-name --indent "" --ids 990721 | grep -P "species"

/work/workspace/zhurj/script/git/mnid2genome/mnid2genome.py -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/input/g_Christensenella_8MNH.txt -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/g_christensenella_MNH.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Christensenella/input

/work/workspace/zhurj/script/git/GCA2genome/gca2genome.py -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/input/g_Christensenella_GCA.txt -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/g_Christensenella_GCAgenome.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Christensenella/input

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Christensenella/input -b /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Christensenella/output/Christensenella -t 16


f_Christensenellaceae
echo "Christensenellaceae" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t 

taxid    kindom     phylum       class        order           family                genus             species
990719   Bacteria   Firmicutes   Clostridia   Clostridiales   Christensenellaceae   Christensenella

taxonkit list --show-rank --show-name --indent "" --ids 990719 | grep -P "species"

/work/workspace/zhurj/script/git/mnid2genome/mnid2genome.py -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/input/f_Christensenellaceae_9MNH.txt -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/f_Christensenellaceae_MNH.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f_Christensenellaceae/input

/work/workspace/zhurj/script/git/GCA2genome/gca2genome.py -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/input/f_Christensenellaceae_GCA.txt -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/f_Christensenellaceae_GCAgenome.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f_Christensenellaceae/input

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f_Christensenellaceae/input -b /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f_Christensenellaceae/output/Christensenellaceae -t 16



s__Sutterella_wadsworthensis
echo "Sutterella" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t 

taxid   kindom     phylum           class                order             family           genus        species
40544   Bacteria   Proteobacteria   Betaproteobacteria   Burkholderiales   Sutterellaceae   Sutterella

taxonkit list --show-rank --show-name --indent "" --ids 40544 | grep -P "species"
cat /work/workspace/zhurj/reference/NCBI/prokaryotes/20200305/prokaryotes.txt | grep -P "\t40545\t" | awk -F "\t" '{print $1"\t"$2"\t"$16"\t"$19"\t"$20"\t"$7"\t"$8}' >> /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/taxonkit/g__Parabacteroides.txt

/work/workspace/zhurj/script/git/mnid2genome/mnid2genome.py -l MNH36384 -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/s__Sutterella_MNH.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/s__Sutterella/input

/work/workspace/zhurj/script/git/GCA2genome/gca2genome.py -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/input/s__Sutterella_9GCA.txt -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/s__Sutterella_GCAgenome.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/s__Sutterella/input

GCA30
/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/s__Sutterella/input -b /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/s__Sutterella/output/Sutterella -t 16

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH36384/ssu_filter1.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH36384/MNH36384.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&NODE_77_length_1099_cov_3253.789628    AP018786        100.000 640     0       0       1       640     818     1457    0.0     1155
genomic&&NODE_77_length_1099_cov_3253.789628    LT223579        99.219  640     5       0       1       640     818     1457    0.0     1132
genomic&&NODE_77_length_1099_cov_3253.789628    AJ566849        97.816  641     13      1       1       640     808     1448    0.0     1090
genomic&&NODE_77_length_1099_cov_3253.789628    LT623892        97.812  640     13      1       2       640     819     1458    0.0     1088
genomic&&NODE_77_length_1099_cov_3253.789628    AB300989        97.500  640     16      0       1       640     820     1459    0.0     1083


cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "AP018786|LT223579|AJ566849|LT623892|AB300989"
AB300989        Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Sutterellaceae;Sutterella;Sutterella parvirubra;YIT 11816(T)
AJ566849        Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Sutterellaceae;Sutterella;Sutterella stercoricanis;CCUG 47620(T)
AP018786        Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Sutterellaceae;Sutterella;Sutterella megalosphaeroides;6FBBBH3(T)
LT223579        Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Sutterellaceae;Sutterella;Sutterella massiliensis;Marseille-P2435(T)
LT623892        Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Sutterellaceae;Sutterella;Sutterella timonensis;Marseille-P3282(T)

fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/s__Sutterella/input/MNH36384.fna -r /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/s__Sutterella/input/GCA_003609995.fna -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/s__Sutterella/MNH36384_GCA_003609995
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/assemble/MNH36384/MNH36384.fna -r /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/s__Sutterella/input/GCA_003609995.fna -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/s__Sutterella/MNH36384_GCA_003609995_all
mkdir /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_003609995
checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/all/GCA/GCA_003/GCA_003609/GCA_003609995 /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_003609995 -x fna -t 16 --tab_table -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_003609995/GCA_003609995_checkm.txt
MNH36384 基因组过滤后所有序列完整性，和污染评估结果一致

/work/assembly/current/MNH/MNH363/MNH36384/genomic.fna

g__Megasphaera
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH05117/ssu_filter.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH05117/MNH05117.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&MNH05117_80    KX021300        97.580  1281    30      1       55      1334    1490    210     0.0     2168
genomic&&MNH05117_80    GU366015        97.268  1281    33      2       55      1334    1491    212     0.0     2146
genomic&&MNH05117_80    HE576794        96.013  1279    50      1       55      1333    1486    209     0.0     2074
genomic&&MNH05117_80    HM990965        95.947  1283    46      6       55      1333    1489    209     0.0     2058
genomic&&MNH05117_80    GQ480004        95.405  1284    51      2       55      1334    1491    212     0.0     2049


cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "KX021300|GU366015"
GU366015        Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Megasphaera;AB032
KX021300        Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Megasphaera;Megasphaera hexanoica;MH(T)

echo "Megasphaera" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t 

taxid   kindom     phylum       class           order            family            genus         species
906     Bacteria   Firmicutes   Negativicutes   Veillonellales   Veillonellaceae   Megasphaera

taxonkit list --show-rank --show-name --indent "" --ids 906 | grep -P "species"

/work/workspace/zhurj/script/git/mnid2genome/mnid2genome.py -l MNH05117 -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/g_Megasphaera_MNH.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Megasphaera/input

/work/workspace/zhurj/script/git/GCA2genome/gca2genome.py -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/input/g_Megasphaera_29GCA.txt -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/g_Megasphaera_GCAgenome.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Megasphaera/input

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Megasphaera/input -b /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Megasphaera/output/Megasphaera -t 16

fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Megasphaera/input/MNH05117.fna -r /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Megasphaera/input/GCA_900538925.fna -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/s__Sutterella/MNH05117_GCA_900538925
mkdir /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_900538925
checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/all/GCA/GCA_900/GCA_900538/GCA_900538925 /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_900538925 -x fna -t 16 --tab_table -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_900538925/GCA_900538925_checkm.txt

g__Blautia
echo "Blautia" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t 

taxid    kindom     phylum       class        order           family            genus     species
572511   Bacteria   Firmicutes   Clostridia   Clostridiales   Lachnospiraceae   Blautia
taxonkit list --show-rank --show-name --indent "" --ids 572511 | grep -P "species"

/work/workspace/zhurj/script/git/mnid2genome/mnid2genome.py -l MNH19801 MNH19407 -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/g_Blautia_MNH.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Blautia/input

/work/workspace/zhurj/script/git/GCA2genome/gca2genome.py -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/input/g_Blautia_200GCA.txt -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/g_Blautia_GCAgenome.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Blautia/input

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Blautia/input -b /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Blautia/output/Blautia -t 16

fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Blautia/input/MNH19407.fna -r /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Blautia/input/GCA_900066145.fna -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/s__Sutterella/MNH19407_GCA_900066145
mkdir /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_900066145
checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/all/GCA/GCA_900/GCA_900066/GCA_900066145 /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_900066145 -x fna -t 16 --tab_table -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_900066145/GCA_900066145_checkm.txt

MNH19801
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH19801/ssu_filter1.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH19801/MNH19801.out -max_target_seqs 5 -num_threads 16 -outfmt 7

genomic&&MNH19801_45    PAC002897       96.166  913     35      0       52      964     1448    536     0.0     1489
genomic&&MNH19801_45    QUCK01000052    95.728  913     39      0       52      964     1456    544     0.0     1471
genomic&&MNH19801_45    LT860099        94.137  904     47      6       52      953     1446    547     0.0     1370
genomic&&MNH19801_45    EF529620        93.319  913     59      2       52      964     1456    546     0.0     1370
genomic&&MNH19801_45    AY442823        93.428  913     57      3       52      964     1453    544     0.0     1366


cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "PAC002897"
PAC002897       Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;QUCK_g;None


g__Veillonella
MNH34134
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH34134/ssu_filter.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH34134/MNH34134.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&NODE_10_length_1343_cov_2970.832543    AENU01000007    99.686  1273    4       0       1       1273    214     1486    0.0     2278
genomic&&NODE_10_length_1343_cov_2970.832543    CP001820        99.450  1273    6       1       1       1273    214     1485    0.0     2261
genomic&&NODE_10_length_1343_cov_2970.832543    AB679109        99.372  1273    7       1       1       1273    213     1484    0.0     2257
genomic&&NODE_10_length_1343_cov_2970.832543    EF185167        98.743  1273    15      1       1       1273    207     1478    0.0     2223
genomic&&NODE_10_length_1343_cov_2970.832543    KQ960728        98.743  1273    15      1       1       1273    214     1485    0.0     2221
# BLAST processed 1 queries

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "AENU01000007|CP001820|AB679109|EF185167|KQ960728"
AB679109        Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Veillonella;Veillonella tobetsuensis;B16(T)
AENU01000007    Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Veillonella;Veillonella rogosae;F0412
CP001820        Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Veillonella;Veillonella parvula;DSM 2008(T)
EF185167        Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Veillonella;Veillonella denticariosi;RBV106(T)
KQ960728        Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Veillonella;DNF00926

MNH36331
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH36331/ssu_filter1.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH36331/MNH36331.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&NODE_28_length_965_cov_1870.565315     AENU01000007    100.000 896     0       0       55      950     1486    591     0.0     1617
genomic&&NODE_28_length_965_cov_1870.565315     AB679109        99.554  896     3       1       55      950     1484    590     0.0     1595
genomic&&NODE_28_length_965_cov_1870.565315     CP001820        99.330  896     5       1       55      950     1485    591     0.0     1586
genomic&&NODE_28_length_965_cov_1870.565315     RQUY01000036    98.996  896     8       1       55      950     1485    591     0.0     1572
genomic&&NODE_28_length_965_cov_1870.565315     AMEX01000022    98.996  896     8       1       55      950     1485    591     0.0     1572
# BLAST processed 1 queries
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "AENU01000007|CP001820|AB679109|RQUY01000036|AMEX01000022"
AB679109        Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Veillonella;Veillonella tobetsuensis;B16(T)
AENU01000007    Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Veillonella;Veillonella rogosae;F0412
AMEX01000022    Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Veillonella;Veillonella atypica;KON(T)
CP001820        Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Veillonella;Veillonella parvula;DSM 2008(T)
RQUY01000036    Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Veillonella;Veillonella caviae;DSM 20738(T)

echo "Veillonella" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t 

taxid   kindom     phylum       class           order            family            genus         species
29465   Bacteria   Firmicutes   Negativicutes   Veillonellales   Veillonellaceae   Veillonella

taxonkit list --show-rank --show-name --indent "" --ids 29465 | grep -P "species"

/work/workspace/zhurj/script/git/mnid2genome/mnid2genome.py -l MNH36331 MNH34134 -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/g_Veillonella_MNH.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Veillonella/input

/work/workspace/zhurj/script/git/GCA2genome/gca2genome.py -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/input/g_Veillonella_GCA.txt -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/g_Veillonella_GCAgenome.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Veillonella/input

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Veillonella/input -b /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Veillonella/output/Veillonella -t 16

fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Veillonella/input/MNH34134.fna -r /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Veillonella/input/GCA_900552445.fna -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/g_Veillonella/MNH34134_GCA_900552445
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Veillonella/input/MNH36331.fna -r /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Veillonella/input/GCA_900552445.fna -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/g_Veillonella/MNH36331_GCA_900552445

mkdir /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_900552445
checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/all/GCA/GCA_900/GCA_900552/GCA_900552445 /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_900552445 -x fna -t 16 --tab_table -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_900552445/GCA_900552445_checkm.txt

fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Veillonella/input/MNH34134.fna -r /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Veillonella/input/GCA_002959835.fna -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/g_Veillonella/MNH34134_GCA_002959835
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Veillonella/input/MNH36331.fna -r /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Veillonella/input/GCA_002959835.fna -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/g_Veillonella/MNH36331_GCA_002959835
mkdir /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_002959835
checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/all/GCA/GCA_002/GCA_002959/GCA_002959835 /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_002959835 -x fna -t 16 --tab_table -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_002959835/GCA_002959835_checkm.txt

f__Verrucomicrobiaceae
MNH20648
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH20648/ssu.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH20648/MNH20648.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&NODE_18_length_5761_cov_1617.580577    PJKF01000002    100.000 1434    0       0       37      1470    1       1434    0.0     2587
genomic&&NODE_18_length_5761_cov_1617.580577    PJKB01000002    99.582  1434    6       0       37      1470    1       1434    0.0     2560
genomic&&NODE_18_length_5761_cov_1617.580577    CP001071        99.024  1434    12      2       37      1470    1       1432    0.0     2516
genomic&&NODE_18_length_5761_cov_1617.580577    KQ968618        98.047  1434    28      0       37      1470    1       1434    0.0     2461
genomic&&NODE_18_length_5761_cov_1617.580577    KT254068        94.024  1439    78      6       37      1470    1       1436    0.0     2188
# BLAST processed 1 queries

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "PJKF01000002|PJKB01000002|CP001071|KQ968618"
CP001071        Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;Akkermansia muciniphila;ATCC BAA-835(T)
KQ968618        Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;KLE1797
PJKB01000002    Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;PJKB_s;GP22
PJKF01000002    Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;PJKF_s;GP15

echo "Akkermansia" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t 

taxid    kindom     phylum            class              order                family            genus         species
239934   Bacteria   Verrucomicrobia   Verrucomicrobiae   Verrucomicrobiales   Akkermansiaceae   Akkermansia

taxonkit list --show-rank --show-name --indent "" --ids 239934 | grep -P "species"

/work/workspace/zhurj/script/git/mnid2genome/mnid2genome.py -l MNH20648 MNH20651 -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/g_Akkermansia_MNH.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Akkermansia/input

/work/workspace/zhurj/script/git/GCA2genome/gca2genome.py -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/input/g_Akkermansia_GCA.txt -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/g_Akkermansia_GCAgenome.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Akkermansia/input

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Akkermansia/input -b /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Akkermansia/output/Akkermansia -t 16

fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Akkermansia/input/MNH20648.fna -r /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Akkermansia/input/GCA_002885535.fna -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/g_Akkermansia/MNH20648_GCA_002885535
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Akkermansia/input/MNH20651.fna -r /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Akkermansia/input/GCA_002885535.fna -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/g_Akkermansia/MNH20651_GCA_002885535
mkdir /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_002885535
checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/all/GCA/GCA_002/GCA_002885/GCA_002885535 /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_002885535 -x fna -t 16 --tab_table -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_002885535/GCA_002885535_checkm.txt

f__Lachnospiraceae
echo "Lachnospiraceae" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t 

taxid    kindom     phylum       class        order           family            genus   species
186803   Bacteria   Firmicutes   Clostridia   Clostridiales   Lachnospiraceae

taxonkit list --show-rank --show-name --indent "" --ids 186803 | grep -P "species"

MNH20130
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH20130/ssu_filter1.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH20130/MNH20130.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&MNH20130_80    JRFU01000184    98.524  1016    15      0       1       1016    1098    83      0.0     1765
genomic&&MNH20130_80    L34623  97.933  1016    21      0       1       1016    1098    83      0.0     1757
genomic&&MNH20130_80    AJ312385        95.074  1015    49      1       3       1016    1095    81      0.0     1602
genomic&&MNH20130_80    CP003040        94.877  1015    51      1       3       1016    1095    81      0.0     1593
genomic&&MNH20130_80    JNIN01000001    94.778  1015    53      0       2       1016    1096    82      0.0     1592

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "JRFU01000184|L34623"
JRFU01000184    Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Eubacterium_g21;JRFU_s;21
L34623  Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Eubacterium_g21;Eubacterium ramulus;ATCC 29099(T)

MNH20128
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH20128/ssu_filter.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH20128/MNH20128.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&MNH20128_96    JRFU01000184    99.117  453     4       0       55      507     1455    1003    0.0     800
genomic&&MNH20128_96    L34623  97.092  447     13      0       62      508     1448    1002    0.0     759
genomic&&MNH20128_96    PAC001269       96.256  454     16      1       55      508     1457    1005    0.0     739
genomic&&MNH20128_96    EU124830        95.796  452     19      0       55      506     1457    1006    0.0     732
genomic&&MNH20128_96    AF371632        95.344  451     21      0       57      507     1452    1002    0.0     719

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "JRFU01000184|L34623"
JRFU01000184    Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Eubacterium_g21;JRFU_s;21
L34623  Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Eubacterium_g21;Eubacterium ramulus;ATCC 29099(T)

MNH20160
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH20160/ssu_filter.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH20160/MNH20160.out -max_target_seqs 5 -num_threads 16 -outfmt 7
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
genomic&&MNH20160_58    JRFU01000184    98.485  858     13      0       1       858     1098    241     0.0     1489
genomic&&MNH20160_58    L34623  97.786  858     19      0       1       858     1098    241     0.0     1480
genomic&&MNH20160_58    AJ312385        95.683  857     36      1       3       858     1095    239     0.0     1376
genomic&&MNH20160_58    PAC002348       95.216  857     41      0       2       858     1097    241     0.0     1361
genomic&&MNH20160_58    CP003040        95.099  857     41      1       3       858     1095    239     0.0     1353

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "JRFU01000184|L34623"
JRFU01000184    Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Eubacterium_g21;JRFU_s;21
L34623  Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Eubacterium_g21;Eubacterium ramulus;ATCC 29099(T)

g__Eubacterium
echo "Eubacterium" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t 

taxid   kindom     phylum       class        order           family           genus         species
1730    Bacteria   Firmicutes   Clostridia   Clostridiales   Eubacteriaceae   Eubacterium

taxonkit list --show-rank --show-name --indent "" --ids 1730 | grep -P "species"

/work/workspace/zhurj/script/git/mnid2genome/mnid2genome.py -l MNH20130 MNH20128 MNH20160 -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/g_Eubacterium_MNH.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Eubacterium/input

/work/workspace/zhurj/script/git/GCA2genome/gca2genome.py -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/input/gg_Eubacterium_GCA.txt -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/g_Eubacterium_GCAgenome.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Eubacterium/input

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Eubacterium/input -b /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Eubacterium/output/Eubacterium -t 16

MNH20130
MNH20128
MNH20160

fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Eubacterium/input/MNH20130.fna -r /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Eubacterium/input/GCA_004167445.fna -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/g_Eubacterium/MNH20130_GCA_004167445
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Eubacterium/input/MNH20128.fna -r /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Eubacterium/input/GCA_004167445.fna -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/g_Eubacterium/MNH20128_GCA_004167445
fastANI -t 16 -q /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Eubacterium/input/MNH20160.fna -r /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/g_Eubacterium/input/GCA_004167445.fna -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/fastANI/g_Eubacterium/MNH20160_GCA_004167445
mkdir /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_004167445
checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/all/GCA/GCA_004/GCA_004167/GCA_004167445 /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_004167445 -x fna -t 16 --tab_table -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_004167445/GCA_004167445_checkm.txt
=IF(or(and(F7533>=2.08,F7533<=3.68,G7533>=37.3,G7533<=38.9),and(F7533>=1.82,F7533<=3.49,G7533>=47.6,G7533<=49.8),and(F7533>=2.34,F7533<=3.94,G7533>=48.3,G7533<=49.9),and(F7533>=1.46,F7533<=3.06,G7533>=49.3,G7533<=50.9),and(F7533>=2.22,F7533<=3.86,G7533>=51.27,G7533<=52.9),and(F7533>=2.35,F7533<=4,G7533>=53.2,G7533<=55)),1,"")

=IF(AND(F7533>=2.66,F7533<=3.7,G7533>=53.5,G7533<54.6),1,"")
MNH06154
MNH06163
MNH06165
MNH05119
/work/workspace/zhurj/script/git/mnid2genome/mnid2genome.py -l MNH06154 MNH06163 MNH06165 MNH05119 -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/o__Clostridiales_g1_4MNH.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/o__Clostridiales_g1/input

/work/workspace/zhurj/script/git/GCA2genome/gca2genome.py -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/input/o__Clostridiales_g1_GCA.txt -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/o__Clostridiales_g1_GCAgenome.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/o__Clostridiales_g1/input

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/o__Clostridiales_g1/input -b /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/o__Clostridiales_g1/output/Eubacterium -t 16

MNH06154
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH06154/ssu.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH06154/MNH06154.out -max_target_seqs 5 -num_threads 16 -outfmt 7
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 10 hits found
genomic&&MNH06154_32    KQ965520        98.767  1379    17      0       56      1434    1458    80      0.0     2411     80.6
genomic&&MNH06154_32    MAIQ01000023    98.615  1372    19      0       56      1427    1458    87      0.0     2389
genomic&&MNH06154_32    LT700187        97.741  1372    31      0       56      1427    1458    87      0.0     2335
genomic&&MNH06154_32    FLKP01000002    97.317  1379    37      0       56      1434    1458    80      0.0     2321
genomic&&MNH06154_32    LAYJ01000068    96.793  1372    44      0       56      1427    1458    87      0.0     2277

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "KQ965520|MAIQ01000023"
KQ965520        Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;Christensenella;Christensenella minuta;DSM 22607(T)
MAIQ01000023    Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;Christensenella;MAIQ_s;AF73-05CM02

MNH04997
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH04997/ssu.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH04997/MNH04997.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&MNH04997_35    MAIQ01000023    98.151  1460    25      2       33      1492    1       1458    0.0     2505
genomic&&MNH04997_35    KQ965520        97.603  1460    33      2       33      1492    1       1458    0.0     2469
genomic&&MNH04997_35    FLKP01000002    97.192  1460    39      2       33      1492    1       1458    0.0     2442
genomic&&MNH04997_35    LAYJ01000068    96.986  1460    42      2       33      1492    1       1458    0.0     2428
genomic&&MNH04997_35    LT700187        96.164  1460    54      1       33      1492    1       1458    0.0     2379

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "KQ965520|MAIQ01000023"
KQ965520        Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;Christensenella;Christensenella minuta;DSM 22607(T)
MAIQ01000023    Bacteria;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;Christensenella;MAIQ_s;AF73-05CM02

MNH19366 g_Oscillospiraceae
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH19366/ssu_filter.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/MNH19366/MNH19366.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&MNH19366_33    GQ451293        96.998  1266    31      4       55      1317    1457    196     0.0     2101
genomic&&MNH19366_33    PAC000672       96.840  1266    33      6       55      1317    1457    196     0.0     2083
genomic&&MNH19366_33    PAC003100       96.443  1265    41      4       55      1317    1458    196     0.0     2065
genomic&&MNH19366_33    PAC001560       96.364  1265    41      5       55      1317    1457    196     0.0     2057
genomic&&MNH19366_33    EF445167        95.656  1266    48      6       55      1317    1457    196     0.0     2015
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "GQ451293|PAC000672"
GQ451293        Bacteria;Firmicutes;Clostridia;Clostridiales;Oscillospiraceae;PAC000672_g;J299
PAC000672       Bacteria;Firmicutes;Clostridia;Clostridiales;Oscillospiraceae;PAC000672_g;

echo "Oscillospiraceae" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t 

taxid    kindom     phylum       class        order           family             genus   species
216572   Bacteria   Firmicutes   Clostridia   Clostridiales   Oscillospiraceae

taxonkit list --show-rank --show-name --indent "" --ids 216572 | grep -P "species"

/work/workspace/zhurj/script/git/mnid2genome/mnid2genome.py -l MNH19366  -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/f_Oscillospiraceae_MNH.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f_Oscillospiraceae/input

/work/workspace/zhurj/script/git/GCA2genome/gca2genome.py -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/input/f_Oscillospiraceae_GCA.txt -f /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/ofile/f_Oscillospiraceae_GCAgenome.txt -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f_Oscillospiraceae/input

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f_Oscillospiraceae/input -b /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f_Oscillospiraceae/output/Oscillospiraceae -t 16

echo "Phascolarctobacterium faecium" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t

错误
Carnobacteriaceae Phascolarctobacterium faecium 1 1429  AACCATGCAGTCGAACGGAGAATTTTATTTCGGTAGAATTCTTAGTGGCGAACGGGTGAGTAACGCGTAGGCAACCTGCCCTTTAGACGGGGACAACATTCCGAAAGGAGTGCTAATACCGGATGTGATCATCGTGCCGCATGGCAGGATGAAGAAAGATGGCCTCTACAAGTAAGCTATCGCTAAAGGATGGGCCTGCGTCTGATTAGCTAGTTGGTAGTGTAACGGACTACCAAGGCGATGATCAGTAGCCGGTCTGAGAGGATGAACGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTGGGGAATCTTCCGCAATGGACGAAAGTCTGACGGAGCAACGCCGCGTGAGTGATGAAGGATTTCGGTCTGTAAAGCTCTGTTGTTTATGACGAACGTGCAGTGTGTGAACAATGCATTGCAATGACGGTAGTAAACGAGGAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCGAGCGTTGTCCGGAATTATTGGGCGTAAAGAGCATGTAGGCGGCTTAATAAGTCGAGCGTGAAAATGCGGGGCTCAACCCCGTATGGCGCTGGAAACTGTTAGGCTTGAGTGCAGGAGAGGAAAGGGGAATTCCCAGTGTAGCGGTGAAATGCGTAGATATTGGGAGGAACACCAGTGGCGAAGGCGCCTTTCTGGACTGTGTCTGACGCTGAGATGCGAAAGCCAGGGTAGCGAACGGGATTAGATACCCCGGTAGTCCTGGCCGTAAACGATGGGTACTAGGTGTAGGAGGTATCGACCCCTTCTGTGCCGGAGTTAACGCAATAAGTACCCCGCCTGGGGAGTACGGCCGCAAGGTTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGTATGTGGTTTAATTCGACGCAACGCGAAGAACCTTACCAAGGCTTGACATTGATTGAACGCTCTAGAGATAGAGATTTCCCTTCGGGGACAAGAAAACAGGTGGTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCCTATCCTATGTTACCAGCAAGTAAAGTTGGGGACTCATGGGAGACTGCCAGGGACAACCTGGAGGAAGGCGGGGATGACGTCAAGTCATCATGCCCCTTATGTCTTGGGCTACACACGTACTACAATGGTCGGAAACAGAGGGAAGCGAAGCCGCGAGGCAGAGCAAACCCCAGAAACCCGATCTCAGTTCGGATCGCAGGCTGCAACCCGCCTGCGTGAAGTCGGAATCGCTAGTAATCGCAGGTCAGCATACTGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGGTAACCTATTAGGAGCCAGCC

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/3_16S/phylogenetic16S20200330/input/MNH04317.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/3_16S/phylogenetic16S20200330/input/MNH04317.out -max_target_seqs 5 -num_threads 16 -outfmt 7
# BLASTN 2.9.0+
# Query: MNH04317
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
MNH04317        X72865  99.720  1427    3       1       4       1429    27      1453    0.0     2554
MNH04317        PAC002124       94.398  1428    77      3       4       1429    27      1453    0.0     2205
MNH04317        FJ825568        93.886  1423    79      4       9       1429    5       1421    0.0     2164
MNH04317        HQ775656        93.627  1428    81      6       4       1429    27      1446    0.0     2148
MNH04317        EU009791        93.702  1429    83      7       4       1429    27      1451    0.0     2147
# BLAST processed 1 queries


cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "X72865"
X72865  Bacteria;Firmicutes;Negativicutes;Acidaminococcales;Acidaminococcaceae;Phascolarctobacterium;Phascolarctobacterium faecium;ACM 3679(T)

taxid   kindom     phylum       class           order               family               genus                   species
33025   Bacteria   Firmicutes   Negativicutes   Acidaminococcales   Acidaminococcaceae   Phascolarctobacterium   Phascolarctobacterium faecium

checkm ssu_finder /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f__Coriobacteriaceae/input/GCA_002320125.fna /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/ete3/f__Coriobacteriaceae/input /work/workspace/zhurj/project/2_swgs/taxonomy_20200110/bacteriaReclass20200323/checkm/GCA_002320125 -x fna -t 16

echo "Barnesiella intestinihominis" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH02743/ssu1038.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH02743/ssu_1038.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&MNH02743_48    ABYH01000014    97.769  986     22      0       53      1038    1453    468     0.0     1680
genomic&&MNH02743_48    LT700189        97.262  986     27      0       53      1038    1453    468     0.0     1657
genomic&&MNH02743_48    QRMP01000005    97.059  986     29      0       53      1038    1451    466     0.0     1648
genomic&&MNH02743_48    AAXE02000112    97.059  986     29      0       53      1038    1453    468     0.0     1648
genomic&&MNH02743_48    JN680579        95.951  988     36      4       54      1038    1448    462     0.0     1590

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "ABYH01000014|LT700189"
ABYH01000014    Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Parabacteroides;Parabacteroides johnsonii;DSM 18315(T)
LT700189        Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Parabacteroides;SN4

echo "Parabacteroides johnsonii" | taxonkit name2taxid
Parabacteroides johnsonii       387661

fastANI -t 16 -q /work/assembly/current/MNH/MNH027/MNH02743/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_900/GCA_900624/GCA_900624885/genomic.fna.gz -o /work/workspace/zhurj/project/2_swgs/checkmvsGCdepth20200401/fastANI/GCA_900624885

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH26050/ssu_1568.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH26050/ssu_1568.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&NODE_147_length_4998_cov_32.961796     BCVU01000117    99.933  1485    1       0       31      1515    1       1485    0.0     2674
genomic&&NODE_147_length_4998_cov_32.961796     AYYP01000002    94.418  1487    79      4       31      1515    1       1485    0.0     2294
genomic&&NODE_147_length_4998_cov_32.961796     AEOF01000010    94.217  1487    83      3       31      1515    1       1486    0.0     2284
genomic&&NODE_147_length_4998_cov_32.961796     BCVJ01000104    94.149  1487    84      3       31      1515    1       1486    0.0     2279
genomic&&NODE_147_length_4998_cov_32.961796     AB812750        93.880  1487    89      2       31      1515    1       1487    0.0     2267
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "BCVU01000117"
BCVU01000117    Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;Lactobacillus ruminis;NBRC 102161(T)
echo "Lactobacillus ruminis" | taxonkit name2taxid
Lactobacillus ruminis   1623


/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH26050/ssu_1523.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH26050/ssu_1523.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&NODE_123_length_5782_cov_2306.254338   CP001071        99.861  1434    0       2       56      1489    1432    1       0.0     2571
genomic&&NODE_123_length_5782_cov_2306.254338   PJKF01000002    99.163  1434    12      0       56      1489    1434    1       0.0     2533
genomic&&NODE_123_length_5782_cov_2306.254338   PJKB01000002    98.954  1434    15      0       56      1489    1434    1       0.0     2519
genomic&&NODE_123_length_5782_cov_2306.254338   KQ968618        98.047  1434    28      0       56      1489    1434    1       0.0     2461
genomic&&NODE_123_length_5782_cov_2306.254338   KT254068        94.024  1439    78      6       56      1489    1436    1       0.0     2188

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "CP001071|PJKF01000002|PJKB01000002"
CP001071        Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;Akkermansia muciniphila;ATCC BAA-835(T)
PJKB01000002    Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;PJKB_s;GP22
PJKF01000002    Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;PJKF_s;GP15
echo "Akkermansia muciniphila" | taxonkit name2taxid
Akkermansia muciniphila 239935

/work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH02793/ssu_832.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH02793/ssu_832.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&MNH02793_256   CP000140        97.772  808     16      2       1       806     808     1       0.0     1370
genomic&&MNH02793_256   HQ807067        97.436  780     20      0       1       780     780     1       0.0     1317
genomic&&MNH02793_256   DQ824372        95.916  808     30      3       1       806     807     1       0.0     1298
genomic&&MNH02793_256   HM124306        95.533  806     36      0       1       806     806     1       0.0     1292
genomic&&MNH02793_256   LARM01000059    98.246  741     13      0       1       741     741     1       0.0     1278

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH02793/ssu_632.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH02793/ssu_632.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&MNH02793_300   CP000140        99.205  629     5       0       54      682     1449    821     0.0     1113
genomic&&MNH02793_300   LARM01000059    97.440  625     16      0       54      678     1382    758     0.0     1056
genomic&&MNH02793_300   HQ760014        97.134  628     18      0       55      682     1421    794     0.0     1052
genomic&&MNH02793_300   HQ786179        96.820  629     20      0       54      682     1425    797     0.0     1045
genomic&&MNH02793_300   AMCI01005881    96.343  629     23      0       54      682     1453    825     0.0     1031
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "CP000140|LARM01000059"

CP000140        Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Parabacteroides;Parabacteroides distasonis;ATCC 8503(T)
LARM01000059    Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Parabacteroides;LARM_s;N54.MGS-20

GCA_000012845.1
GCA_900683725.1

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH12192/ssu_649.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH12192/ssu_649.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&MNH12192_52    CP000140        99.846  649     1       0       1       649     450     1098    0.0     1167
genomic&&MNH12192_52    LARM01000059    98.459  649     10      0       1       649     383     1031    0.0     1126
genomic&&MNH12192_52    HQ807067        97.385  650     16      1       1       649     422     1071    0.0     1093
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "CP000140|LARM01000059"

/work/database/ncbi/genome/all/GCA/GCA_000/GCA_000012/GCA_000012845/genomic.fna.gz
/work/database/ncbi/genome/all/GCA/GCA_900/GCA_900683/GCA_900683725/genomic.fna.gz
fastANI -t 16 -q /work/database/ncbi/genome/all/GCA/GCA_900/GCA_900683/GCA_900683725/genomic.fna.gz -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000012/GCA_000012845/genomic.fna.gz -o /work/workspace/zhurj/project/2_swgs/checkmvsGCdepth20200401/fastANI/GCA_900683725_GCA_000012845

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH04850/ssu.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH04850/ssu.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&NODE_28_length_5071_cov_1890.105727    KE350284        99.933  1482    1       0       55      1536    1482    1       0.0     2669
genomic&&NODE_28_length_5071_cov_1890.105727    KE350284        78.261  46      10      0       1453    1498    39      84      0.23    39.2
genomic&&NODE_28_length_5071_cov_1890.105727    BCQE01000074    99.933  1482    1       0       55      1536    1482    1       0.0     2669
genomic&&NODE_28_length_5071_cov_1890.105727    BCQE01000074    78.261  46      10      0       1453    1498    39      84      0.23    39.2
genomic&&NODE_28_length_5071_cov_1890.105727    LHOX01000013    99.865  1482    2       0       55      1536    1482    1       0.0     2664
genomic&&NODE_28_length_5071_cov_1890.105727    LHOX01000013    78.261  46      10      0       1453    1498    39      84      0.23    39.2
genomic&&NODE_28_length_5071_cov_1890.105727    AF039903        99.798  1482    3       0       55      1536    1482    1       0.0     2660
genomic&&NODE_28_length_5071_cov_1890.105727    AF039903        78.261  46      10      0       1453    1498    39      84      0.23    39.2
genomic&&NODE_28_length_5071_cov_1890.105727    NGKU01000001    99.730  1482    4       0       55      1536    1482    1       0.0     2655
genomic&&NODE_28_length_5071_cov_1890.105727    NGKU01000001    76.087  46      11      0       1453    1498    39      84      9.7     34.6

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "KE350284|BCQE01000074|LHOX01000013|AF039903|NGKU01000001"
AF039903        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus casseliflavus;MUTK 20(T)
BCQE01000074    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus gallinarum;NBRC 100675(T)
KE350284        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;13.SD.W.09
LHOX01000013    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;LHOX_s;RIT-PI-f
NGKU01000001    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;NGKU_s;8G7_MSG3316
GCA_000157355.2
GCA_001544275.1
GCA_000414965.1
GCA_001297065.1
GCA_002140915.1

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH33931/ssu_1325.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH33931/ssu_1325.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&NODE_52_length_4676_cov_935.318113     AAYH02000049    100.000 1263    0       0       1       1263    187     1449    0.0     2278
genomic&&NODE_52_length_4676_cov_935.318113     BAKJ01000105    98.654  1263    16      1       1       1263    187     1448    0.0     2198
genomic&&NODE_52_length_4676_cov_935.318113     OEST01000016    97.432  1246    32      0       18      1263    204     1449    0.0     2104
genomic&&NODE_52_length_4676_cov_935.318113     QROH01000022    96.675  1263    42      0       1       1263    183     1445    0.0     2089
genomic&&NODE_52_length_4676_cov_935.318113     PAC002443       96.437  1263    45      0       1       1263    187     1449    0.0     2076
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "AAYH02000049|BAKJ01000105"
AAYH02000049    Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;Bacteroides uniformis;ATCC 8492(T)
BAKJ01000105    Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;Bacteroides rodentium;JCM 16496(T)
GCA_000154205.1
GCA_000614125.1

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH15751/ssu_832.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH15751/ssu_832.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&MNH15751_44    CP000140        97.772  808     16      2       1       806     808     1       0.0     1370
genomic&&MNH15751_44    HQ807067        97.436  780     20      0       1       780     780     1       0.0     1317
genomic&&MNH15751_44    DQ824372        95.916  808     30      3       1       806     807     1       0.0     1298
genomic&&MNH15751_44    HM124306        95.533  806     36      0       1       806     806     1       0.0     1292
genomic&&MNH15751_44    LARM01000059    98.246  741     13      0       1       741     741     1       0.0     1278

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "CP000140"

GCA_000012845.1

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH20105/ssu_630.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH20105/ssu_630.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&MNH20105_147   LARM01000059    99.048  630     6       0       1       630     741     112     0.0     1110
genomic&&MNH20105_147   CP000140        99.048  630     6       0       1       630     808     179     0.0     1110
genomic&&MNH20105_147   DQ805826        98.254  630     11      0       1       630     814     185     0.0     1087
genomic&&MNH20105_147   HQ807067        98.095  630     12      0       1       630     780     151     0.0     1083
genomic&&MNH20105_147   HQ760262        96.667  630     21      0       1       630     783     154     0.0     1042
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "LARM01000059|CP000140|DQ805826|HQ807067"
LARM01000059    Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Parabacteroides;LARM_s;N54.MGS-20
CP000140        Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Parabacteroides;Parabacteroides distasonis;ATCC 8503(T)
DQ805826        Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Parabacteroides;RL306aal91e10
HQ807067        Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Parabacteroides;ELU0158-T300-S-NIPCRAMgANa_000283

GCA_000980475.1 LARM01000059
GCA_000012845.1 ATCC 8503(T)


/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH20069/ssu_897.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH20069/ssu_897.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&MNH20069_114   CP000140        98.765  891     11      0       1       891     963     73      0.0     1558
genomic&&MNH20069_114   LARM01000059    97.980  891     18      0       1       891     896     6       0.0     1526
genomic&&MNH20069_114   HQ807067        97.868  891     19      0       1       891     935     45      0.0     1522
genomic&&MNH20069_114   HQ760262        96.409  891     32      0       1       891     938     48      0.0     1463
genomic&&MNH20069_114   DQ824372        96.629  890     27      3       4       891     959     71      0.0     1460


/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH20069/ssu_612.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH20069/ssu_612.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&MNH20069_116   CP000140        99.104  558     5       0       1       558     892     1449    0.0     985
genomic&&MNH20069_116   LARM01000059    97.312  558     15      0       1       558     825     1382    0.0     939
genomic&&MNH20069_116   HQ760014        96.948  557     17      0       1       557     865     1421    0.0     929
genomic&&MNH20069_116   HQ786179        96.595  558     19      0       1       558     868     1425    0.0     921
genomic&&MNH20069_116   HQ760262        96.237  558     21      0       1       558     867     1424    0.0     912

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "CP000140|LARM01000059"
CP000140        Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Parabacteroides;Parabacteroides distasonis;ATCC 8503(T)
LARM01000059    Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Parabacteroides;LARM_s;N54.MGS-20

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH12776/ssu.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH12776/ssu.out -max_target_seqs 5 -num_threads 16 -outfmt 7
MNH12776
genomic&&NODE_21_length_5076_cov_3050.583717    AP010888        99.793  1448    2       1       55      1502    1447    1       0.0     2595
genomic&&NODE_21_length_5076_cov_3050.583717    JGZA01000002    99.448  1448    7       1       55      1502    1447    1       0.0     2572
genomic&&NODE_21_length_5076_cov_3050.583717    CP001095        99.309  1448    10      0       55      1502    1448    1       0.0     2567
genomic&&NODE_21_length_5076_cov_3050.583717    AB924514        99.309  1448    10      0       55      1502    1448    1       0.0     2567
genomic&&NODE_21_length_5076_cov_3050.583717    ACCG01000002    96.903  1453    40      2       55      1502    1453    1       0.0     2414
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "AP010888|JGZA01000002|CP001095|AB924514"
AB924514        Bacteria;Actinobacteria;Actinobacteria_c;Bifidobacteriales;Bifidobacteriaceae;Bifidobacterium;Bifidobacterium longum;Bifidobacterium longum subsp. suillum;Su 851(T)
AP010888        Bacteria;Actinobacteria;Actinobacteria_c;Bifidobacteriales;Bifidobacteriaceae;Bifidobacterium;Bifidobacterium longum;Bifidobacterium longum subsp. longum;JCM 1217(T)
CP001095        Bacteria;Actinobacteria;Actinobacteria_c;Bifidobacteriales;Bifidobacteriaceae;Bifidobacterium;Bifidobacterium longum;Bifidobacterium longum subsp. infantis;ATCC 15697(T)
JGZA01000002    Bacteria;Actinobacteria;Actinobacteria_c;Bifidobacteriales;Bifidobacteriaceae;Bifidobacterium;Bifidobacterium longum;Bifidobacterium longum subsp. suis;LMG 21814(T)

cat prokaryotes.txt | grep -P "Bifidobacterium longum" | grep REFR
GCA_000007525.1

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH31221/ssu_1527.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH31221/ssu_1527.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&NODE_20_length_5361_cov_2489.836109    AP010888        99.793  1448    2       1       55      1502    1447    1       0.0     2595
genomic&&NODE_20_length_5361_cov_2489.836109    JGZA01000002    99.448  1448    7       1       55      1502    1447    1       0.0     2572
genomic&&NODE_20_length_5361_cov_2489.836109    CP001095        99.309  1448    10      0       55      1502    1448    1       0.0     2567
genomic&&NODE_20_length_5361_cov_2489.836109    AB924514        99.309  1448    10      0       55      1502    1448    1       0.0     2567
genomic&&NODE_20_length_5361_cov_2489.836109    ACCG01000002    96.903  1453    40      2       55      1502    1453    1       0.0     2414

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "AP010888|JGZA01000002|CP001095|AB924514"
cat prokaryotes.txt | grep -P "Bifidobacterium longum" | grep REFR
GCA_000007525.1

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH32716/ssu.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH32716/ssu.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&NODE_18_length_5644_cov_2174.778157    AP010888        99.793  1448    2       1       55      1502    1447    1       0.0     2595
genomic&&NODE_18_length_5644_cov_2174.778157    JGZA01000002    99.448  1448    7       1       55      1502    1447    1       0.0     2572
genomic&&NODE_18_length_5644_cov_2174.778157    CP001095        99.309  1448    10      0       55      1502    1448    1       0.0     2567
genomic&&NODE_18_length_5644_cov_2174.778157    AB924514        99.309  1448    10      0       55      1502    1448    1       0.0     2567
genomic&&NODE_18_length_5644_cov_2174.778157    ACCG01000002    96.903  1453    40      2       55      1502    1453    1       0.0     2414

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH36171/ssu.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH36171/ssu.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&NODE_46_length_5523_cov_1087.156996    JNJP01000178    99.796  1473    3       0       55      1527    1473    1       0.0     2644
genomic&&NODE_46_length_5523_cov_1087.156996    HM123944        94.297  1473    83      1       55      1527    1472    1       0.0     2275
genomic&&NODE_46_length_5523_cov_1087.156996    PAC002477       93.890  1473    89      1       55      1527    1472    1       0.0     2248
genomic&&NODE_46_length_5523_cov_1087.156996    AGFG01000025    91.859  1474    115     5       55      1527    1470    1       0.0     2100
genomic&&NODE_46_length_5523_cov_1087.156996    KI632512        91.588  1474    119     5       55      1527    1470    1       0.0     2082

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "JNJP01000178"
JNJP01000178    Bacteria;Proteobacteria;Deltaproteobacteria;Desulfovibrionales;Desulfovibrionaceae;Bilophila;Bilophila wadsworthia;ATCC 49260(T)
cat prokaryotes.txt | grep -P "ATCC 49260"
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH31374/ssu.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH31374/ssu.out -max_target_seqs 5 -num_threads 16 -outfmt 7

genomic&&NODE_37_length_5234_cov_1425.001745    AAXE02000112    99.036  1453    10      2       27      1475    1       1453    0.0     2553
genomic&&NODE_37_length_5234_cov_1425.001745    LT700189        97.660  1453    30      2       27      1475    1       1453    0.0     2462
genomic&&NODE_37_length_5234_cov_1425.001745    ABYH01000014    97.591  1453    31      2       27      1475    1       1453    0.0     2458
genomic&&NODE_37_length_5234_cov_1425.001745    JN680579        96.416  1451    46      6       27      1474    1       1448    0.0     2363
genomic&&NODE_37_length_5234_cov_1425.001745    QRMP01000005    95.865  1451    58      1       27      1475    1       1451    0.0     2344
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "AAXE02000112"
AAXE02000112    Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Parabacteroides;Parabacteroides merdae;ATCC 43184(T)

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH11355/ssu_1276.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH11355/ssu_1276.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&NODE_7_length_92179_cov_549.867321     CP003504        99.679  1246    4       0       1       1246    1246    1       0.0     2230
genomic&&NODE_7_length_92179_cov_549.867321     BCQB01000108    99.599  1246    5       0       1       1246    1246    1       0.0     2225
genomic&&NODE_7_length_92179_cov_549.867321     JXLE01000039    99.117  1246    11      0       1       1246    1246    1       0.0     2198
genomic&&NODE_7_length_92179_cov_549.867321     AJ301830        99.359  1248    3       5       1       1246    1245    1       0.0     2197
genomic&&NODE_7_length_92179_cov_549.867321     AY321376        99.037  1246    12      0       1       1246    1246    1       0.0     2195
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "CP003504|BCQB01000108|JXLE01000039|AJ301830|AY321376"
CP003504        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus hirae;ATCC 9790(T)
BCQB01000108    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus durans;NBRC 100479(T)
JXLE01000039    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus thailandicus;DSM 21767(T)
AJ301830        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus faecium;LMG 11423(T)
AY321376        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus sanguinicola;SS-1729(T)
GCA_000271405.2
GCA_001544215.1
GCA_001886265.1
GCA_000174395.2

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH22080/ssu_1565.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH22080/ssu_1565.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&NODE_103_length_1699_cov_1259.857583   AF039903        99.798  1482    3       0       55      1536    1482    1       0.0     2660
genomic&&NODE_103_length_1699_cov_1259.857583   BCQE01000074    99.933  1482    1       0       55      1536    1482    1       0.0     2669
genomic&&NODE_103_length_1699_cov_1259.857583   KE350284        99.933  1482    1       0       55      1536    1482    1       0.0     2669
genomic&&NODE_103_length_1699_cov_1259.857583   LHOX01000013    99.865  1482    2       0       55      1536    1482    1       0.0     2664
genomic&&NODE_103_length_1699_cov_1259.857583   NGKU01000001    99.730  1482    4       0       55      1536    1482    1       0.0     2655
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "KE350284|BCQE01000074|LHOX01000013|AF039903|NGKU01000001"
AF039903        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus casseliflavus;MUTK 20(T)
BCQE01000074    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus gallinarum;NBRC 100675(T)
KE350284        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;13.SD.W.09
LHOX01000013    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;LHOX_s;RIT-PI-f
NGKU01000001    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;NGKU_s;8G7_MSG3316
GCA_000157355.2
GCA_001544275.1
GCA_000414965.1
GCA_001297065.1
GCA_002140915.1

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH22080/ssu_1457.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH22080/ssu_1457.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&NODE_104_length_1546_cov_2059.185160   AP007281        99.715  1403    4       0       55      1457    1498    96      0.0     2513
genomic&&NODE_104_length_1546_cov_2059.185160   KT343143        98.574  1403    19      1       55      1457    1469    68      0.0     2437
genomic&&NODE_104_length_1546_cov_2059.185160   FTOY01000011    98.147  1403    24      2       55      1457    1496    96      0.0     2406
genomic&&NODE_104_length_1546_cov_2059.185160   PNGN01000029    97.648  1403    32      1       55      1457    1497    96      0.0     2379
genomic&&NODE_104_length_1546_cov_2059.185160   PNFV01000020    97.577  1403    34      0       55      1457    1498    96      0.0     2378
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "AP007281|KT343143|FTOY01000011"
AP007281        Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;Lactobacillus reuteri;JCM 1112(T)
KT343143        Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;Lactobacillus caviae;MOZM2(T) 无基因组
FTOY01000011    Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;FTOY_s;Marseille-P3519
GCA_000010005.1
GCA_900156885.1

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH02743/ssu_1038.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH02743/ssu_1038.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&MNH02743_48    ABYH01000014    97.769  986     22      0       53      1038    1453    468     0.0     1680
genomic&&MNH02743_48    QRMP01000005    97.059  986     29      0       53      1038    1451    466     0.0     1648
genomic&&MNH02743_48    AAXE02000112    97.059  986     29      0       53      1038    1453    468     0.0     1648
genomic&&MNH02743_48    JN680579        95.951  988     36      4       54      1038    1448    462     0.0     1590
genomic&&MNH02743_48    LT700189        97.262  986     27      0       53      1038    1453    468     0.0     1657
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "ABYH01000014|LT700189|QRMP01000005|AAXE02000112"
ABYH01000014    Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Parabacteroides;Parabacteroides johnsonii;DSM 18315(T)
LT700189        Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Parabacteroides;SN4
QRMP01000005    Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Parabacteroides;QRMP_s;AM08-6
AAXE02000112    Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Parabacteroides;Parabacteroides merdae;ATCC 43184(T)
GCA_000156495.1
GCA_003473295.1
GCA_000154105.1


/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH02743/ssu_564.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH02743/ssu_564.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&MNH02743_5     QRMP01000005    96.840  538     17      0       27      564     1       538     0.0     894
genomic&&MNH02743_5     MH697664        95.508  512     22      1       53      564     1       511     0.0     817
genomic&&MNH02743_5     ABYH01000014    93.727  542     28      3       27      564     1       540     0.0     817
genomic&&MNH02743_5     JN680579        93.866  538     29      3       27      564     1       534     0.0     812
genomic&&MNH02743_5     LT700189        92.804  542     33      3       27      564     1       540     0.0     794

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH06247/ssu.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH06247/ssu.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&MNH06247_14    AB081585        98.619  1376    16      2       103     1478    1373    1       0.0     2394
genomic&&MNH06247_14    PAC002236       95.439  1425    64      1       55      1478    1425    1       0.0     2274
genomic&&MNH06247_14    KE993469        94.881  1426    69      3       55      1478    1424    1       0.0     2233
genomic&&MNH06247_14    CP040924        94.807  1425    71      3       55      1478    1423    1       0.0     2226
genomic&&MNH06247_14    LN879450        94.530  1426    75      3       55      1478    1425    1       0.0     2210

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "AB081585"
AB081585        Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae;Clostridium;A4d
GCA_000013285.1
GCA_000016505.1
GCA_000017045.1
GCA_000154805.1
GCA_000154345.1
GCA_000154505.1
GCA_000189595.1
GCA_000371705.1
GCA_000008765.1
GCA_000007625.1
GCA_000833105.2
GCA_000014125.1
GCA_001456065.2
GCA_001038625.1
GCA_000145275.1
GCA_000158075.1
GCA_000156515.1
GCA_000156055.1
GCA_000158655.1
GCA_000371425.1
GCA_000144625.1
GCA_000233455.1
GCA_000807255.1
GCA_000320405.1
GCA_000511955.1
GCA_000340885.1
GCA_000424205.1
GCA_001593985.1
GCA_000246895.2
GCA_000401215.1
GCA_001642655.1
GCA_000647895.1
GCA_000469625.2
GCA_000383295.1
GCA_000424025.1
GCA_000484505.1
GCA_000473995.1
GCA_900168365.1
GCA_000953215.1
GCA_000285575.1
GCA_000612845.1
GCA_000619945.1
GCA_000620945.1
GCA_000686665.1
GCA_000686705.1
GCA_000686725.1
GCA_000687555.1
GCA_000703125.1
GCA_000711825.1
GCA_000732635.1
GCA_000789395.1
GCA_002074155.1
GCA_001042715.1
GCA_000577895.1
GCA_001047375.1
GCA_001263795.1
GCA_000499525.1
GCA_001405015.1
GCA_001404895.1
GCA_000820705.1
GCA_001458595.1
GCA_001584565.1
GCA_001594005.1
GCA_001623875.1
GCA_001735765.2
GCA_001758365.1
GCA_001854185.1
GCA_900103025.1
GCA_900104115.1
GCA_900102365.1
GCA_900111235.1
GCA_900111595.1
GCA_900112485.1
GCA_900111985.1
GCA_900112775.1
GCA_001877035.1
GCA_900142075.1
GCA_900130005.1
GCA_900129365.1
GCA_900141845.1
GCA_900129965.1
GCA_001244495.1
GCA_002006345.1
GCA_002006355.1
GCA_002029295.1
GCA_002029235.1
GCA_002029255.1
GCA_002050515.1
GCA_900176305.1
GCA_900176635.1
GCA_002008345.1
GCA_000389635.1
GCA_000020165.1
GCA_000063585.1
GCA_001705235.1
GCA_000204565.1

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH26240/ssu_892.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH26240/ssu_892.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&NODE_148_length_5539_cov_1070.614061   CP001071        99.883  858     0       1       1       858     857     1       0.0     1540
genomic&&NODE_148_length_5539_cov_1070.614061   PJKF01000002    98.718  858     11      0       1       858     858     1       0.0     1498
genomic&&NODE_148_length_5539_cov_1070.614061   PJKB01000002    98.601  858     12      0       1       858     858     1       0.0     1494
genomic&&NODE_148_length_5539_cov_1070.614061   KQ968618        97.552  858     21      0       1       858     858     1       0.0     1453
genomic&&NODE_148_length_5539_cov_1070.614061   KT254068        92.343  862     59      5       1       858     859     1       0.0     1242
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "CP001071|PJKF01000002|PJKB01000002|KQ968618"

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/result/checkm/MNH/MNH26240/ssu_698.fna -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/result/checkm/MNH/MNH26240/ssu_698.out -max_target_seqs 5 -num_threads 16 -outfmt 7
genomic&&NODE_461_length_1058_cov_3378.401631   PJKF01000002    99.844  643     1       0       56      698     1434    792     0.0     1156
genomic&&NODE_461_length_1058_cov_3378.401631   CP001071        99.844  643     0       1       56      698     1432    791     0.0     1152
genomic&&NODE_461_length_1058_cov_3378.401631   PJKB01000002    99.533  643     3       0       56      698     1434    792     0.0     1147
genomic&&NODE_461_length_1058_cov_3378.401631   KQ968618        98.911  643     7       0       56      698     1434    792     0.0     1129
genomic&&NODE_461_length_1058_cov_3378.401631   KT254068        96.729  642     20      1       56      696     1436    795     0.0     1060

CP001071        Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;Akkermansia muciniphila;ATCC BAA-835(T)
PJKF01000002    Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;PJKF_s;GP15
PJKB01000002    Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;PJKB_s;GP22
KQ968618        Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;KLE1797

GCA_000020225.1
GCA_002885025.1
GCA_002884975.1
GCA_001578645.1

perl /work/workspace/zhouyj/script/alignKit/quick-bowtie2.pl -i /work/rawdata/fastq/genome/MNH/MNH050/MNH05020/E1/L1/S1/MNH05020.1.fq.gz,/work/rawdata/fastq/genome/MNH/MNH050/MNH05020/E1/L1/S1/MNH05020.2.fq.gz -d MNH05020.fas -o MNH05020.Self -n 100

https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc

perl /work/workspace/zhouyj/script/alignKit/quick-bowtie2.pl -i /work/workspace/zhurj/download/SRA/SRR3496270.1.fastq -d /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000701/GCA_000701705/genomic.fna.gz -o /work/workspace/zhurj/project/2_swgs/checkmvsGCdepth20200401/GCdepth/GCA_000701705.Self -n 100

/work/assembly/current/MNH/MNH262/MNH26240/genomic.fna

fastANI -t 16 -q /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000020/GCA_000020225/genomic.fna.gz -r /work/assembly/fasta/genome/MNH/MNH262/MNH26240/E1/L1/S1/MNH26240.fa.gz -o /work/workspace/zhurj/project/2_swgs/checkmvsGCdepth20200401/fastANI/GCA_000020225_MNH26240all

/work/workspace/zhurj/script/1_tools/9_strain2seq/strainscreen_lstr.py \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/straindb.txt \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/output \
-s /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/swgseq_strain.txt \
--ongoing /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/ongoing_strain.txt \
--question /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/questionstrain.txt \
-i /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan.txt \
-n strainselect20200407

/work/workspace/zhurj/script/1_tools/9_strain2seq/strainscreen_lstr.py \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/straindb.txt \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/output \
-s /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/swgseq_strain.txt \
--ongoing /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/ongoing_strain.txt \
--question /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/questionstrain.txt \
-d 1 \
-n strainselect20200407

/work/workspace/zhurj/script/1_tools/9_strain2seq/renewPlanstrain.py \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/output \
-i /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/processing_20200408.txt \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan_20200220.txt

/work/workspace/zhurj/script/1_tools/9_strain2seq/renewPlanstrain.py \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input \
-i /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/processing_20200408.txt \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan_20200220.txt

/work/workspace/zhurj/script/1_tools/9_strain2seq/strainscreen_lstr.py \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/straindb.txt \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/output \
-s /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/swgseq_strain.txt \
--ongoing /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/ongoing_strain.txt \
--question /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/questionstrain.txt \
-i /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan.txt \
-n strainselect20200414

/work/workspace/zhurj/script/1_tools/9_strain2seq/strainscreen_lstr.py \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/straindb.txt \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/output \
-s /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/swgseq_strain.txt \
-d 1 \
-n strainMaxMNID20200414

# 1. seqed: straindb['taxid'][taxid]['seqed'][donor] = set()
  # 2. ongoing: straindb['taxid'][taxid]['ongoing'][donor] = set()
  # 3. not seqed: straindb['taxid'][taxid]['unseq'][mnid] = donor
  # 4. seqed: straindb['name'][name]['seqed'][donor] = set()
  # 5. ongoing: straindb['taxid'][taxid]['ongoing'][donor] = set()
  # 6. not seqed: straindb['name'][name]['unseq'][mnid] = donor
  # if specific_name, straindb['id2n'][taxid] = specific_name
  # else straindb['id2n'][taxid] = input_name 

cd  /work/workspace/zhurj/project/2_swgs/sanger16S/test
blastdbcmd -db /work/workspace/zhurj/reference/NCBI/16S/2020.Feb28/16S_ribosomal_RNA -entry 'NR_026331.1' > query.fa
blastdbcmd -db /work/workspace/zhurj/reference/NCBI/16S/2020.Feb28/16S_ribosomal_RNA -entry_batch list.txt > target.fa
content of list.txt
NR_024570.1
NR_112558.1
NR_114042.1
cat 16S_ribosomal_RNA.fas | grep "Streptococcus lutetiensis" | awk -F "[> ]" '{print $2}' > /work/workspace/zhurj/project/2_swgs/sanger16S/test/Streptococcus_lutetiensis.txt
blastdbcmd -db /work/workspace/zhurj/reference/NCBI/16S/2020.Feb28/16S_ribosomal_RNA -entry_batch /work/workspace/zhurj/project/2_swgs/sanger16S/test/Streptococcus_lutetiensis.txt > /work/workspace/zhurj/project/2_swgs/sanger16S/test/Streptococcus_lutetiensis.fa
cat /work/workspace/zhurj/reference/NCBI/16S/2020.Feb28/16S_ribosomal_RNA.fas | grep "Streptococcus equinus" | awk -F "[> ]" '{print $2}' > /work/workspace/zhurj/project/2_swgs/sanger16S/test/Streptococcus_equinus.txt
blastdbcmd -db /work/workspace/zhurj/reference/NCBI/16S/2020.Feb28/16S_ribosomal_RNA -entry_batch /work/workspace/zhurj/project/2_swgs/sanger16S/test/Streptococcus_equinus.txt > /work/workspace/zhurj/project/2_swgs/sanger16S/test/Streptococcus_equinus.fa
blastn -query /work/workspace/zhurj/project/2_swgs/sanger16S/test/Streptococcus_equinus.fa -subject /work/workspace/zhurj/project/2_swgs/sanger16S/test/Streptococcus_lutetiensis.fa -outfmt 6 -max_target_seqs 1

/work/workspace/moonis/bin/genoassemble/soft/SPAdes-3.13.1-Linux/bin/spades.py --plasmid \
-1 /work/rawdata/fastq/genome/MNH/MNH192/MNH19250/E1/L1/S1/MNH19250.1.fq.gz \
-2 /work/rawdata/fastq/genome/MNH/MNH192/MNH19250/E1/L1/S1/MNH19250.2.fq.gz \
-t 20 \
--careful \
--cov-cutoff 5 \
-o /work/workspace/zhurj/project/2_swgs/plasmidSPAdes/MNH19250

/work/rawdata/fastq/genome/MNH/MNH192/MNH19250/E1/L1/S1/MNH19250.2.fq.gz \
plasmidfinder.py \
-i /work/workspace/zhurj/project/2_swgs/plasmidfinder/MNH19250/input/MNH19250.fasta \
-mp /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/blastn \
-t 0.6 \
-p /work/workspace/zhurj/bin/miniconda3/envs/python3.6/share/plasmidfinder-2.1-0/database \
-o /work/workspace/zhurj/project/2_swgs/plasmidfinder/MNH19250/fasta

plasmidfinder.py \
-i /work/workspace/zhurj/project/2_swgs/plasmidfinder/MNH19250/input/MNH19250.1.fastq /work/workspace/zhurj/project/2_swgs/plasmidfinder/MNH19250/input/MNH19250.2.fastq \
-mp /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/kma \
-t 0.9 \
-p /work/workspace/zhurj/bin/miniconda3/envs/python3.6/share/plasmidfinder-2.1-0/database \
-o /work/workspace/zhurj/project/2_swgs/plasmidfinder/MNH19250/fastq

plasmidfinder.py \
-i /work/rawdata/fastq/genome/MNH/MNH192/MNH19250/E1/L1/S1/MNH19250.1.fq.gz /work/rawdata/fastq/genome/MNH/MNH192/MNH19250/E1/L1/S1/MNH19250.2.fq.gz \
-mp /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/kma \
-t 0.9 \
-p /work/workspace/zhurj/bin/miniconda3/envs/python3.6/share/plasmidfinder-2.1-0/database \
-o /work/workspace/zhurj/project/2_swgs/plasmidfinder/MNH19250/fastq

plasmidfinder.py \
-i /work/rawdata/fastq/genome/MNH/MNH192/MNH19250/E1/L1/S1/MNH19250.1.fq.gz /work/rawdata/fastq/genome/MNH/MNH192/MNH19250/E1/L1/S1/MNH19250.2.fq.gz \
-mp /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/kma \
-t 0.9 \
-p /work/workspace/zhurj/bin/miniconda3/envs/python3.6/share/plasmidfinder-2.1-0/database/acfd0096c01a \
-o /work/workspace/zhurj/project/2_swgs/plasmidfinder/MNH19250/fastq

echo "Prevotella copri" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t
cat nodes.dmp | grep -P "\t838\t" | awk -F "\t" '{print $1"\t"$5}'
taxonkit list --show-rank --show-name --indent "" --ids 838 | grep -P "species"

plasmidfinder.py \
-i /work/assembly/fasta/genome/MNH/MNH192/MNH19250/E1/L1/S1/MNH19250.fa.gz \
-mp /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/blastn \
-t 0.8 \
-p /work/workspace/zhurj/bin/miniconda3/envs/python3.6/share/plasmidfinder-2.1-0/database \
-o /work/workspace/zhurj/project/2_swgs/plasmidfinder/MNH19250/fasta

/work/workspace/zhurj/script/1_tools/9_strain2seq/renewPlanstrain.py \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input \
-i /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/processing_20200509.txt \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan_20200220.txt

/work/workspace/zhurj/script/1_tools/9_strain2seq/strainscreen_lstr.py \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/straindb.txt \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/output \
-s /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/swgseq_strain.txt \
--ongoing /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/ongoing_strain.txt \
--question /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/questionstrain.txt \
-i /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan.txt \
-n strainselect20200509

a. copy sheet“需测序菌株B-I列”
#MNID SpeciesName Taxid 提交时间  Submit_batch  菌株出库  菌株入库  currentStatus
b. 序列 212-315 行，删除 281-315 行， 将 212-280 行中标黄色的行，在currentStatus 列标注为 “yes”，备注正在活化菌株 currentStatus列标注为 “uncertain”

echo "Megasphaera massiliensis" | taxonkit name2taxid | taxonkit lineage --taxid-field 2

GGACTAC[ATC]][GAC]GGGT[AT]TCTAAT
ATTAGA[AT]ACCC[GTC][GAT]GTAGTCC

zcat MNH27992.fa.gz | grep "GGACTAC[ATC][GAC]GGGT[AT]TCTAAT"

/work/assembly/fasta/genome/MNH/MNH279/MNH27992/E1/L1/S1

MNH27359

/work/assembly/fasta/genome/MNH/MNH279/MNH27992/E1/L1/S1

checkm ssu_finder /work/assembly/current/MNH/MNH047/MNH04744/genomic.fna /work/assembly/current/MNH/MNH047/MNH04744 ./extract16S -x fna -t 20
/work/assembly/fasta/genome/MNH/MNH279/MNH27992/E1/L1/S1/MNH27992.fa.gz

checkm ssu_finder /work/workspace/zhurj/project/2_swgs/extract16S/checkM/MNH27992/input/MNH27992.fna /work/workspace/zhurj/project/2_swgs/extract16S/checkM/MNH27992/input /work/workspace/zhurj/project/2_swgs/extract16S/checkM/MNH27992/16S -x fna -t 20

/work/workspace/zhurj/script/1_tools/9_strain2seq/renewPlanstrain.py \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input \
-i /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/processing_20200526.txt \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan_20200220.txt

rm /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/questionstrain.txt -f
rm /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/ongoing_strain.txt -f
rm /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan.txt -f
ln -s /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/questionstrain_20200526.txt /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/questionstrain.txt
ln -s /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/ongoing_20200526.txt /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/ongoing_strain.txt
ln -s /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan_20200526.txt /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan.txt

/work/workspace/zhurj/script/1_tools/9_strain2seq/strainscreen_lstr.py \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/straindb.txt \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/output \
-s /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/swgseq_strain.txt \
--ongoing /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/ongoing_strain.txt \
--question /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/questionstrain.txt \
-i /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan.txt \
-n strainselect20200526

/work/workspace/zhurj/script/1_tools/9_strain2seq/renewPlanstrain.py \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input \
-i /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/processing_20200605.txt \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan_20200220.txt


rm  /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/ongoing_strain.txt -f
ln -s /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/ongoing_20200605.txt  /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/ongoing_strain.txt
rm  /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/questionstrain.txt -f
ln -s /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/questionstrain_20200605.txt  /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/questionstrain.txt
rm  /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan.txt -f
ln -s /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan_20200605.txt  /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan.txt


/work/workspace/zhurj/script/1_tools/9_strain2seq/strainscreen_lstr.py \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/straindb.txt \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/output \
-s /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/swgseq_strain.txt \
--ongoing /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/ongoing_strain.txt \
--question /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/questionstrain.txt \
-i /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan.txt \
-n strainselect20200605



/work/workspace/zhurj/script/1_tools/9_strain2seq/renewPlanstrain.py \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input \
-i /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/processing_20200610.txt \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan_20200220.txt


rm  /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/ongoing_strain.txt -f
ln -s /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/ongoing_20200610.txt  /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/ongoing_strain.txt
rm  /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/questionstrain.txt -f
ln -s /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/questionstrain_20200610.txt  /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/questionstrain.txt
rm  /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan.txt -f
ln -s /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan_20200610.txt  /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan.txt


/work/workspace/zhurj/script/1_tools/9_strain2seq/strainscreen20200610.py \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/straindb.txt \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/output \
-s /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/swgseq_strain.txt \
--ongoing /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/ongoing_strain.txt \
--question /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/questionstrain.txt \
-i /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan.txt \
-n strainselect20200610


/work/workspace/zhurj/script/1_tools/15_selectStrain2store/input
肥胖所需数据： /work/workspace/liangzj/project/metagenomics/Obesity/PRJEB12123/data
样本信息：/work/workspace/liangzj/project/metagenomics/Obesity/PRJEB12123/map/final_mapping
/work/classify/species/classify/


/work/workspace/zhurj/script/1_tools/16_sumtaxo/sumtaxo.py \
-r /work/classify/species/classify \
-o /work/workspace/zhurj/script/1_tools/16_sumtaxo/example/output \
-i /work/workspace/zhurj/script/1_tools/16_sumtaxo/example/input/mnlist \
-n run00053_54_taxo 

/work/workspace/zhurj/script/1_tools/16_sumtaxo/sumcheckm.py \
-r  /work/cleandata/fastq/genome/evalue \
-o /work/workspace/zhurj/script/1_tools/16_sumtaxo/example/output \
-i /work/workspace/zhurj/script/1_tools/16_sumtaxo/example/input/mnlist \
-n run00053_54_checm


fastANI -q /work/assembly/current/MNH/MNH005/MNH00519/genomic.fna -r /work/workspace/zhurj/reference/NCBI/genome/GCA_000152245/GCF_000152245.2_ASM15224v2_genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxonomy_20200619/fastANI/MNH00519/GCF_000152245
/work/workspace/zhurj/script/1_tools/17_aniMNspename

C:\Users\MoonBiotech\AppData\Local\Programs\Python\Python36\Scritp\pyinstaller -F strain2store2_v2.py

=IF(AND(NOT(EXACT(C4,D4)),NOT(ISBLANK(D4)),IFERROR(NOT(SEARCH("no ",C2,1)=1),"TRUE"),"TRUE")),"1","")

添加pyinstaller到系统环境变量
2020.06.30 添加环境变量
windows桌面搜索框中输入：“系统” => 点击“系统”（蓝色小电脑，非系统信息）=> 高级系统设置 => 高级 => 环境变量 => 系统环境变量 => Path => 新建  输入 需要添加的环境变量
C:\Users\MoonBiotech\AppData\Local\Programs\Python\Python36\Scripts
D:
cd D:\moonbio\work\项目\2.软件开发\14.tools\保藏菌株筛选\strain2store_v2\

/work/workspace/zhurj/script/1_tools/17_aniMNspename/aniMNspename.py 


fastANI -q /work/assembly/current/MNH/MNH260/MNH26050/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000020/GCA_000020225/genomic.fna.gz -t 16 -o /work/workspace/zhurj/script/1_tools/17_aniMNspename/output/test1
fastANI -q /work/assembly/current/MNH/MNH260/MNH26050/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_003/GCA_003096/GCA_003096055/genomic.fna.gz  -t 16 -o /work/workspace/zhurj/script/1_tools/17_aniMNspename/output/test2 --minFrag 0
fastANI -q /work/assembly/current/MNH/MNH260/MNH26050/genomic.fna --rl /work/workspace/zhurj/script/1_tools/17_aniMNspename/output/igeno.txt  -t 16 -o /work/workspace/zhurj/script/1_tools/17_aniMNspename/output/test3 --minFrag 0


/work/workspace/zhurj/script/1_tools/17_aniMNspename/aniMNspename.py \
-i /work/workspace/zhurj/script/1_tools/17_aniMNspename/input/input.txt \
-o /work/workspace/zhurj/script/1_tools/17_aniMNspename/output \
-n output

conda install -c bioconda prodigal
conda install -c bioconda hmmer

gtdbtk classify_wf --genome_dir /work/workspace/zhurj/project/2_swgs/run00053-54/GTDB20200703/ingenodir/ --out_dir /work/workspace/zhurj/project/2_swgs/run00053-54/GTDB20200703/output/test --cpus 2


nohup gtdbtk classify_wf --genome_dir /work/workspace/zhurj/project/2_swgs/run00053-54/GTDB20200703/sample72/inputdir/ --out_dir /work/workspace/zhurj/project/2_swgs/run00053-54/GTDB20200703/sample72/output --cpus 2 > log.txt &
cd /work/workspace/zhurj/project/2_swgs/run00056/taxonomy/tmp
nohup gtdbtk classify_wf --genome_dir /work/workspace/zhurj/project/2_swgs/run00056/taxonomy/inputdir/ --out_dir /work/workspace/zhurj/project/2_swgs/run00056/taxonomy/output --cpus 2 > log.txt &

/work/workspace/zhurj/script/1_tools/16_sumtaxo/sumtaxo.py \
-r /work/classify/species/classify \
-o /work/workspace/zhurj/project/2_swgs/run00056/taxoserver/output \
-i /work/workspace/zhurj/project/2_swgs/run00056/taxoserver/input/mnlist \
-n run00056_taxo

/work/workspace/zhurj/script/1_tools/16_sumtaxo/sumcheckm.py \
-r /work/cleandata/fastq/genome/evalue \
-o /work/workspace/zhurj/project/2_swgs/run00056/taxoserver/output \
-i /work/workspace/zhurj/project/2_swgs/run00056/taxoserver/input/mnlist \
-n run00056_checkm


/work/assembly/current/MNH/MNH260/MNH26050/genomic.fna MNH26050
gtdbtk classify_wf --batchfile /work/workspace/zhurj/project/2_swgs/GTDB/run00053-56/akk20200707/input/input.txt --out_dir /work/workspace/zhurj/project/2_swgs/GTDB/run00053-56/akk20200707/output --cpus 20
export GTDBTK_DATA_PATH=/work/workspace/zhurj/reference/GTDB/release89
gtdbtk classify_wf --batchfile /work/workspace/zhurj/project/2_swgs/GTDB/run00053-56/run00053_54_sam171/input/input.txt --out_dir /work/workspace/zhurj/project/2_swgs/GTDB/run00053-56/run00053_54_sam171/output --cpus 20

java -Xmx10g  -jar /work/workspace/zhurj/software/FigTree_v1.4.4/lib/figtree.jar -graphic PDF /work/workspace/zhurj/project/2_swgs/GTDB/run00053-56/akk20200707/output/classify/gtdbtk.bac120.classify.tree /work/workspace/zhurj/project/2_swgs/GTDB/run00053-56/akk20200707/output/classify/gtdbtk.bac120.classify.pdf

echo "[Clostridium] lavalense" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t


/work/workspace/zhurj/script/1_tools/10_mnid2assgenome/mnid2assgenome.py  \
-o /work/workspace/zhurj/script/1_tools/10_mnid2assgenome/output \
--batchfile /work/workspace/zhurj/script/1_tools/10_mnid2assgenome/input/mnlist \
--type 2

/work/workspace/zhurj/script/1_tools/10_mnid2assgenome/mnid2assgenome.py  \
--type 2 \
--print \
-i MNH41088 MNH41089


WT_List = [ var for var in ["KRAS","NRAS","BRAF"] if any(re.search(var,x) for x in hash_clinical.keys())]
/work/workspace/zhurj/script/git/seqedCheck
perl /work/workspace/zhurj/script/git/seqedCheck/seqedCheck.pl -i /work/workspace/zhurj/script/git/seqedCheck/checkPro/seq20200709/input.txt -o /work/workspace/zhurj/script/git/seqedCheck/checkPro/seq20200709/

python /work/workspace/zhurj/script/1_tools/19_gca2taxidsciname/gca2taxidsciname.py \
--batchfile /work/workspace/zhurj/script/1_tools/19_gca2taxidsciname/input/input \
-o /work/workspace/zhurj/script/1_tools/19_gca2taxidsciname/output \
--type 2


python /work/workspace/zhurj/script/1_tools/19_gca2taxidsciname/gca2taxidsciname.py \
-i  GCF_000154565.1 GCW_000406945.1 \
-o /work/workspace/zhurj/script/1_tools/19_gca2taxidsciname/output

python /work/workspace/zhurj/script/1_tools/19_gca2taxidsciname/gca2taxidsciname.py \
-i  GCF_000154565.1 GCA_000406945.1 \
--print --type 1

python /work/workspace/zhurj/script/1_tools/19_gca2taxidsciname/gca2taxidsciname.py \
-i  GCF_000154565.1 GCA_000406945.1 \
--print --type 2

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/input.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
--test


gtdbtk classify_wf --batchfile /work/workspace/zhurj/project/2_swgs/GTDB/SGBzj20200714/input/input.txt --out_dir /work/workspace/zhurj/project/2_swgs/GTDB/SGBzj20200714/output --cpus 20

echo "1697784"| taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/input.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
--test

fastANI -q /work/assembly/current/MNH/MNH383/MNH38360/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000239/GCA_000239295/genomic.fna.gz -t 16 -o /work/workspace/zhurj/script/1_tools/17_aniMNspename/output/test1 -k 21
/work/workspace/zhouyj/script/bin/fastANI -q /work/assembly/current/MNH/MNH383/MNH38360/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000239/GCA_000239295/genomic.fna.gz -t 16 -o /work/workspace/zhurj/script/1_tools/17_aniMNspename/output/test1 

echo "Tilletia controversa" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t
Tilletia controversa

fastANI -q /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_003/GCA_003269/GCA_003269465/genomic.fna.gz -t 16 -o /work/workspace/zhurj/script/1_tools/17_aniMNspename/output/test1
fastANI -q /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_006/GCA_006337/GCA_006337145/genomic.fna.gz -t 16 -o /work/workspace/zhurj/script/1_tools/17_aniMNspename/output/test1
fastANI -q /work/assembly/current/MNH/MNH228/MNH22871/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_006/GCA_006337/GCA_006337145/genomic.fna.gz -t 16 -o /work/workspace/zhurj/script/1_tools/17_aniMNspename/output/MNH22871_GCA_006337145
fastANI -q /work/assembly/current/MNH/MNH228/MNH22871/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000299/GCA_000299455/genomic.fna.gz -t 16 -o /work/workspace/zhurj/script/1_tools/17_aniMNspename/output/MNH22871_GCA_000299455
fastANI -q /work/assembly/current/MNH/MNH228/MNH22871/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_902/GCA_902165/GCA_902165245/genomic.fna.gz -t 16 -o /work/workspace/zhurj/script/1_tools/17_aniMNspename/output/MNH22871_GCA_902165245
fastANI -q /work/database/ncbi/genome/all/GCA/GCA_902/GCA_902165/GCA_902165245/genomic.fna.gz -r  /work/assembly/current/MNH/MNH228/MNH22871/genomic.fna -t 16 -o /work/workspace/zhurj/script/1_tools/17_aniMNspename/output/GCA_902165245_MNH22871

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run00056_s18.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run00056_s18

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniGCAmn.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run00056_s18.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run00056_s18_rev

fastANI -q /work/assembly/current/MNH/MNH050/MNH05050/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000179/GCA_000179335/genomic.fna.gz -t 16 -o /work/workspace/zhurj/script/1_tools/17_aniMNspename/output/test1
fastANI -q /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000179/GCA_000179335/genomic.fna.gz -r /work/assembly/current/MNH/MNH050/MNH05050/genomic.fna -t 16 -o /work/workspace/zhurj/script/1_tools/17_aniMNspename/output/test2


python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run00056_s33.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run00056_s33

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniGCAmn.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run00056_s33.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run00056_s33_rev

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run00056_64s_self.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run00056_64s_self

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniGCAmn.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run00056_64s_self.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run00056_64s_self_rev


python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run00056_64s_GDTB.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run00056_64s_GDTB


python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniGCAmn.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run00056_64s_GDTB.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run00056_64s_GDTB_rev


python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run00056_s18_self.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run00056_s18_self

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniGCAmn.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run00056_s18_self.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run00056_s18_self_rev


python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run00056_s18_GDTB.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run00056_s18_GDTB


python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniGCAmn.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run00056_s18_GDTB.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run00056_s18_GDTB_rev

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run00056_33s_final.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run00056_33s_final

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run00056_17s_final.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run00056_17s_final

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run00056_64s_final.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run00056_64s_final

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run00056_17s_new_final.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run00056_17s_new_final

/work/workspace/zhurj/script/1_tools/16_sumtaxo/sumtaxo.py \
-r /work/classify/species/classify \
-o /work/workspace/zhurj/script/1_tools/16_sumtaxo/example/output \
-i /work/workspace/zhurj/script/1_tools/16_sumtaxo/example/input/mnlist248 \
-n run00053_54_taxo 


python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run0005354_45s_self.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run0005354_45s_self

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniGCAmn.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run0005354_45s_self.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run0005354_45s_self_rev


python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run0005354_45s_GTDB.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run0005354_45s_GTDB


python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniGCAmn.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run0005354_45s_GTDB.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run0005354_45s_GTDB_rev

python /work/workspace/zhurj/script/1_tools/19_gca2taxidsciname/gca2taxidsciname.py \
--batchfile /work/workspace/zhurj/script/1_tools/19_gca2taxidsciname/input/run0005354_GTDB.txt \
-o /work/workspace/zhurj/script/1_tools/19_gca2taxidsciname/output \
--type 2

python /work/workspace/zhurj/script/1_tools/19_gca2taxidsciname/gca2taxidsciname.py \
--batchfile /work/workspace/zhurj/script/1_tools/19_gca2taxidsciname/input/run0005354_GTDB.txt \
--print \
--type 2

python /work/workspace/zhurj/script/1_tools/19_gca2taxidsciname/gca2taxidsciname.py \
-i GCA_000210015.1 \
--print \
--type 2
GCA_000210015.1


python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run0005354_22s16s_self.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run0005354_22s16s_self

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniGCAmn.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run0005354_22s16s_self.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run0005354_22s16s_self_rev


python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run0005354_22s16s_GTDB.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run0005354_22s16s_GTDB

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniGCAmn.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run0005354_22s16s_GTDB.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run0005354_22s16s_GTDB_rev

GCA_003340345.1
python /work/workspace/zhurj/script/1_tools/19_gca2taxidsciname/gca2taxidsciname.py \
-i GCA_003340345.1 \
--print \
--type 2

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run0005354_8s_GTDB.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run0005354_8s_GTDB

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniGCAmn.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run0005354_8s_GTDB.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run0005354_8s_GTDB_rev

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run0005354_38s_GTDB.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run0005354_38s_GTDB

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniGCAmn.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run0005354_38s_GTDB.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run0005354_38s_GTDB_rev


python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run0005354_125s_self.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run0005354_125s_self

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniGCAmn.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run0005354_125s_self.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run0005354_125s_self_rev


python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run0005354_125s_GTDB.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run0005354_125s_GTDB

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniGCAmn.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run0005354_125s_GTDB.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run0005354_125s_GTDB_rev

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run0005354_238s_final.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run0005354_238s_final

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniGCAmn.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run0005354_238s_final.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run0005354_238s_final_rev

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run00056_131s_final.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run00056_131s_final

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniGCAmn.py \
-i /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/input/run00056_131s_final.txt \
-o /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output \
-n run00056_131s_final_rev


echo "Clostridiales" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t

D:\moonbio\work\项目\3.研发项目\8.菌种鉴定\2. 数据分析\run00053-56\SelfAna\
/work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/output

fastANI -q /work/assembly/current/MNH/MNH258/MNH25867/genomic.fna -r /work/workspace/zhurj/reference/GTDB/release89/fastani/database/GCA_900315625.1_genomic.fna.gz -t 16 -o /work/workspace/zhurj/script/1_tools/17_aniMNspename/output/test1
fastANI -q /work/workspace/zhurj/reference/GTDB/release89/fastani/database/GCA_900315625.1_genomic.fna.gz -r /work/assembly/current/MNH/MNH258/MNH25867/genomic.fna -t 16 -o /work/workspace/zhurj/script/1_tools/17_aniMNspename/output/test1

fastANI -q /work/assembly/current/MNH/MNH258/MNH25893/genomic.fna -r /work/workspace/zhurj/reference/GTDB/release89/fastani/database/GCA_900315625.1_genomic.fna.gz -t 16 -o /work/workspace/zhurj/script/1_tools/17_aniMNspename/output/test1
fastANI -q /work/workspace/zhurj/reference/GTDB/release89/fastani/database/GCA_900315625.1_genomic.fna.gz -r /work/assembly/current/MNH/MNH258/MNH25893/genomic.fna -t 16 -o /work/workspace/zhurj/script/1_tools/17_aniMNspename/output/test1

/work/workspace/zhurj/reference/NCBI/prokaryotes/current

python /work/workspace/zhurj/script/1_tools/20_fastANIgca2gca/anigca2gca.py \
-i /work/workspace/zhurj/script/1_tools/20_fastANIgca2gca/input/run53_56_29s.txt \
-o /work/workspace/zhurj/script/1_tools/20_fastANIgca2gca/output \
-n run53_56_29s

python /work/workspace/zhurj/script/1_tools/20_fastANIgca2gca/anigca2gcarev.py \
-i /work/workspace/zhurj/script/1_tools/20_fastANIgca2gca/input/run53_56_29s.txt \
-o /work/workspace/zhurj/script/1_tools/20_fastANIgca2gca/output \
-n run53_56_29s_rev

python /work/workspace/zhurj/script/1_tools/19_gca2taxidsciname/gca2taxidsciname.py \
--batchfile /work/workspace/zhurj/script/1_tools/19_gca2taxidsciname/input/run0005356_gcx.txt \
-o  /work/workspace/zhurj/script/1_tools/19_gca2taxidsciname/output \
--type 2


fastANI -q /work/workspace/zhurj/reference/GTDB/release89/fastani/database/GCF_000311925.1_genomic.fna.gz -r /work/database/ncbi/genome/all/GCA/GCA_006/GCA_006542/GCA_006542665/genomic.fna.gz -t 16 -o /work/workspace/zhurj/script/1_tools/17_aniMNspename/output/test1
fastANI -q /work/database/ncbi/genome/all/GCA/GCA_006/GCA_006542/GCA_006542665/genomic.fna.gz -r /work/workspace/zhurj/reference/GTDB/release89/fastani/database/GCF_000311925.1_genomic.fna.gz -t 16 -o /work/workspace/zhurj/script/1_tools/17_aniMNspename/output/test1


gtdbtk classify_wf --batchfile /work/workspace/zhurj/project/2_swgs/GTDB/SGBzj20200714/input/input.txt --out_dir /work/workspace/zhurj/project/2_swgs/GTDB/SGBzj20200714/output --cpus 20

export GTDBTK_DATA_PATH=/work/workspace/zhurj/reference/GTDB/release89
nohup gtdbtk classify_wf --batchfile /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/input/input732s --out_dir /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/output --cpus 20 > /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/log.txt &



ln -s /work/database/ncbi/genome/all/GCA/GCA_001/GCA_001678/GCA_001678845/genomic.fna.gz /work/workspace/zhurj/project/2_swgs/phylotree20200304/CMbaojia20200728/input/AF73_05CM02.fna
ln -s /work/database/ncbi/genome/all/GCA/GCA_900/GCA_900155/GCA_900155415/genomic.fna.gz /work/workspace/zhurj/project/2_swgs/phylotree20200304/CMbaojia20200728/input/DSM_102344.fna
ln -s /work/database/ncbi/genome/all/GCA/GCA_900/GCA_900087/GCA_900087015/genomic.fna.gz /work/workspace/zhurj/project/2_swgs/phylotree20200304/CMbaojia20200728/input/DSM_102800.fna
ln -s /work/workspace/zhurj/project/2_swgs/phylotree20200304/CMbaojia20200728/NZ_CP029256.1.fasta /work/workspace/zhurj/project/2_swgs/phylotree20200304/CMbaojia20200728/input/DSM_22607.fna
ln -s /work/assembly/current/MNH/MNH048/MNH04863/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/CMbaojia20200728/input/MNH04863.fna
ln -s /work/assembly/current/MNH/MNH061/MNH06163/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/CMbaojia20200728/input/MNH06163.fna

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/phylotree20200304/CMbaojia20200728/input -b /work/workspace/zhurj/project/2_swgs/phylotree20200304/CMbaojia20200728/output/cm -t 16

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/CMbaojia20200728/output/cm.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/CMbaojia20200728/output/cm_rec.png", w=500, units="mm", tree_style=ts)

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/CMbaojia20200728/output/cm.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.mode = "c"
ts.arc_start = -180 # 0 degrees = 3 o'clock
ts.arc_span = 270
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/CMbaojia20200728/output/cm_cir.svg", w=500, h=500, units="mm", tree_style=ts)

/work/workspace/zhurj/script/1_tools/9_strain2seq/renewPlanstrain.py \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input \
-i /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/processing_20200803.txt \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan_20200220.txt


rm  /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/ongoing_strain.txt -f
ln -s /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/ongoing_20200803.txt  /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/ongoing_strain.txt
rm  /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/questionstrain.txt -f
ln -s /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/questionstrain_20200803.txt  /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/questionstrain.txt
rm  /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan.txt -f
ln -s /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan_20200803.txt  /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan.txt


/work/workspace/zhurj/script/1_tools/9_strain2seq/strainscreen20200610.py \
-r /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/straindb.txt \
-o /work/workspace/zhurj/project/2_swgs/strainscreen20200217/output \
-s /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/swgseq_strain.txt \
--ongoing /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/ongoing_strain.txt \
--question /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/questionstrain.txt \
-i /work/workspace/zhurj/project/2_swgs/strainscreen20200217/input/seqplan.txt \
-n strainselect20200803


/work/workspace/zhurj/script/1_tools/16_sumtaxo/sumcheckm.py \
-r /work/cleandata/fastq/genome/evalue \
-o /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/checkm \
-i /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/input/mnlist732 \
-n run00022_52_checm

/work/workspace/zhurj/script/1_tools/16_sumtaxo/sumtaxo.py \
-r /work/classify/species/classify \
-o /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/taxonomy \
-i /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/input/mnid36 \
-n run00022_52_taxo 

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/input/mngcalist693 \
-o /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/mngca \
-n run0002252_693s_GDTB


python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniGCAmn.py \
-i /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/input/mngcalist693 \
-o /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/mngca \
-n run0002252_693s_GDTB_rev

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/input/mngcalist12 \
-o /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/mngca \
-n run0002252_12s_GDTB

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniGCAmn.py \
-i /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/input/mngcalist12 \
-o /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/mngca \
-n run0002252_12s_GDTB_rev

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/input/mngcalist8 \
-o /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/mngca \
-n run0002252_8s_self

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniGCAmn.py \
-i /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/input/mngcalist8 \
-o /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/mngca \
-n run0002252_8s_self_rev


/work/classify/current/MNH/MNH316/MNH31644/MNH31644.all.report

/work/workspace/zhurj/script/1_tools/16_sumtaxo/sumtaxo.py \
-r /work/classify/species/classify \
-o /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/taxonomy \
-i /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/input/baoxia4 \
-n baoxia_4s_taxo 

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniGCAmn.py \
-i /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/input/mngcabaoxi4 \
-o /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/mngca \
-n baoxia_4s_self_rev

cdsnum=$(cat /work/predict/current/MNH/MNH149/MNH14923/MNH14923.tsv | grep -P "\tCDS" |wc -l ) && avlen=$(cat /work/predict/current/MNH/MNH149/MNH14923/MNH14923.tsv | grep -P "\tCDS" | awk -F "\t" '{print $3}' | awk '{sum+=$1} END {print sum/NR}') && echo -e "MNH14923\t$cdsnum\t$avlen"

/work/workspace/zhurj/script/1_tools/16_sumtaxo/sumtaxo.py \
-r /work/classify/species/classify \
-o /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/taxonomy \
-i /work/workspace/zhurj/project/2_swgs/GTDB/run00022-50/input/mnlist732 \
-n run0002252_732s_taxo 


cat /work/workspace/zhurj/project/2_swgs/phylotree20200304/CM20200805/inputf/cmlist77 | awk '{print "ln -s /work/assembly/current/MNH/"substr($1,1,6)"/"$1"/genomic.fna /work/workspace/zhurj/project/2_swgs/phylotree20200304/CM20200805/input/"$1".fna"}' > /work/workspace/zhurj/project/2_swgs/phylotree20200304/CM20200805/P/softlink.sh

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/phylotree20200304/CM20200805/input -b /work/workspace/zhurj/project/2_swgs/phylotree20200304/CM20200805/output/cm -t 16

1. GTDB 
export GTDBTK_DATA_PATH=/work/workspace/zhurj/reference/GTDB/release89
nohup gtdbtk classify_wf --batchfile /work/workspace/zhurj/project/2_swgs/GTDB/run00059-60/input/gtdbin139 --out_dir /work/workspace/zhurj/project/2_swgs/GTDB/run00059-60/output --cpus 20 > /work/workspace/zhurj/project/2_swgs/GTDB/run00059-60/log.txt &

2. taxo & checkm
/work/workspace/zhurj/script/1_tools/16_sumtaxo/sumcheckm.py \
-r /work/cleandata/fastq/genome/evalue \
-o /work/workspace/zhurj/project/2_swgs/GTDB/run00059-60/checkm \
-i /work/workspace/zhurj/project/2_swgs/GTDB/run00059-60/input/mnlist139 \
'-n run00059_60_checm

/work/workspace/zhurj/script/1_tools/16_sumtaxo/sumtaxo.py \
-r /work/classify/species/classify \
-o /work/workspace/zhurj/project/2_swgs/GTDB/run00059-60/taxonomy \
-i /work/workspace/zhurj/project/2_swgs/GTDB/run00059-60/input/mnlist139 \
'-n run00059_60_taxo

3. 核查GDTB&self分析结果
python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniMNgca.py \
-i /work/workspace/zhurj/project/2_swgs/GTDB/run00059-60/input/mngca134 \
-o /work/workspace/zhurj/project/2_swgs/GTDB/run00059-60/mngca \
'-n run0005960_134s

python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniGCAmn.py \
-i /work/workspace/zhurj/project/2_swgs/GTDB/run00059-60/input/mngca134 \
-o /work/workspace/zhurj/project/2_swgs/GTDB/run00059-60/mngca \
'-n run0005960_134s_rev

/work/workspace/liangzj/download/database/UHGG/uhgg_catalogue/1_100/MGYG-HGUT-00039.fna

fastANI -q /work/workspace/zhurj/reference/GTDB/release89/fastani/database/GCF_000311925.1_genomic.fna.gz -r /work/database/ncbi/genome/all/GCA/GCA_006/GCA_006542/GCA_006542665/genomic.fna.gz -t 16 -o /work/workspace/zhurj/script/1_tools/17_aniMNspename/output/test1
/work/assembly/current/MNH/MNH048/MNH04863/genomic.fna

/work/workspace/liangzj/download/database/UHGG/uhgg_catalogue/1_100/MGYG-HGUT-00031.fna

cat /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/input/mgygmnlist97 | awk -F "\t" '{print "ln -s /work/assembly/current/MNH/"substr($1,1,6)"/"$1"/genomic.fna /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/softlink/"$1".fna"}' > /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/P/softlink.sh
bash /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/P/softlink.sh
cat /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/input/mgygmnlist97 | awk -F "\t" '{print "/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/fastANI -q "$3" -r /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/softlink/"$1".fna -t 16 -o /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/output/"$1"_"$2}' > /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/P/fastAni.sh
bash /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/P/fastAni.sh

if [ -f /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/merge/animerge.txt ];
then
  rm /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/merge/animerge.txt -f
fi
find /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/output -type f  | awk '{print "awk 'NR==1' "$0" >> /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/merge/animerge.txt"}' > /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/P/mergeANI.sh
bash /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/P/mergeANI.sh


cat /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/input/mgygmnlist1105 | awk -F "\t" '{print "ln -s /work/assembly/current/MNH/"substr($1,1,6)"/"$1"/genomic.fna /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/softlink/"$1".fna"}' > /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/P/softlink.sh
bash /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/P/softlink.sh
cat /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/input/mgygmnlist1105 | awk -F "\t" '{print "/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/fastANI -q "$3" -r /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/softlink/"$1".fna -t 16 -o /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/output/"$1"_"$2}' > /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/P/fastAni.sh
bash /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/P/fastAni.sh

if [ -f /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/merge/animerge1105.txt ];
then
  rm /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/merge/animerge1105.txt -f
fi
find /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/output -type f  | awk '{print "awk 'NR==1' "$0" >> /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/merge/animerge1105.txt"}' > /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/P/mergeANI.sh
bash /work/workspace/zhurj/project/2_swgs/MGYG-HGUT/run0002260_20200806/P/mergeANI.sh


# https://github.com/tseemann/prokka#installation
# prokka
# under no environment
#[zhurj@mnhead sh]$
# /work/workspace/zhurj/software/prokka-master/bin/prokka

/work/workspace/zhurj/lib/sh
bash mnsoftlink MNH04863 /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/fna
bash mnsoftlink MNH06163 /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/fna

export PATH=$PATH:`pwd`/bin
prokka --cpus 20  --outdir --prefix 
prodigal -i PROKKA_08072020\/PROKKA_08072020\.fna -c -m -g 11 -p single -f sco -q
# options ‘-c’ (predict proteins with closed ends only), 
#‘-m’ (prevent genes from being built across stretches of sequence marked as Ns) and 
# ‘-p single’ (single mode for genome assemblies containing a single species)

cat PROKKA_08072020\/PROKKA_08072020\.IS\.tmp\.135961\.faa | parallel --gnu --plain -j 8 --block 58189 --recstart '>' --pipe blastp -query - -db /work/workspace/zhurj/software/prokka-master/db/kingdom/Bacteria/IS -evalue 1e-30 -qcov_hsp_perc 90 -num_threads 1 -num_descriptions 1 -num_alignments 1 -seg no > PROKKA_08072020\/PROKKA_08072020\.IS\.tmp\.135961\.blast 2> /dev/null
#Generating Genbank and Sequin files
tbl2asn -V b -a r10k -l paired-ends -M n -N 1 -y 'Annotated using prokka 1.14.6 from https://github.com/tseemann/prokka' -Z PROKKA_08072020\/PROKKA_08072020\.err -i PROKKA_08072020\/PROKKA_08072020\.fsa 2> /dev/null

表1. Prokka 结果说明

Extension Description
.gff  基因注释文件，包括gff和序列，可用igv直接查看
.gbk  Genebank格式，来自gff
.fna  输入contig核酸文件
.faa  翻译CDS的AA序列
.ffn  所有转录本核酸序列
.sqn  用于提交的序列
.fsa  输入序列，但有sqn的描述，用于tbl2asn生成sqn文件
.tbl  特征表，用于tbl2asn生成sqn文件
.err  错误报告
.log  日志
.txt  统计结果
.tsv  所有注释基因特征表格


export GTDBTK_DATA_PATH=/work/workspace/zhurj/reference/GTDB/release89
nohup gtdbtk classify_wf --batchfile /work/workspace/zhurj/project/2_swgs/GTDB/MGYG43-20200807/input/input43 --out_dir /work/workspace/zhurj/project/2_swgs/GTDB/MGYG43-20200807/output --cpus 20 > /work/workspace/zhurj/project/2_swgs/GTDB/MGYG43-20200807/log.txt &

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/GTDB/MGYG43-20200807/softlink -b /work/workspace/zhurj/project/2_swgs/GTDB/MGYG43-20200807/phylogeneticTree/Eubacterium -t 16
/work/workspace/zhurj/project/2_swgs/GTDB/MGYG43-20200807/phylogeneticTree/Eubacterium.nwk

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/GTDB/MGYG43-20200807/phylogeneticTree/Eubacterium.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True
t.render("/work/workspace/zhurj/project/2_swgs/GTDB/MGYG43-20200807/phylogeneticTree/Eubacterium_rec.png", w=500, units="mm", tree_style=ts)

from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/GTDB/MGYG43-20200807/phylogeneticTree/Eubacterium.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.mode = "c"
ts.arc_start = -180 # 0 degrees = 3 o'clock
ts.arc_span = 270
t.render("/work/workspace/zhurj/project/2_swgs/GTDB/MGYG43-20200807/phylogeneticTree/Eubacterium_cir.svg", w=500, h=500, units="mm", tree_style=ts)

##mmseqs protein cluster
cd /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/fna/PROKKA_08072020
mmseqs createdb /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/fna/PROKKA_08072020/PROKKA_08072020.faa /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/fna/PROKKA_08072020/inDB/test_proteins
mmseqs cluster --cov-mode 0 -c 0.8 --min-seq-id 0.9 --threads 20 /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/fna/PROKKA_08072020/PROKKA_08072020.faa /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/fna/PROKKA_08072020/oDB/test /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/fna/PROKKA_08072020/tmpdir
mmseqs easy-cluster /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/fna/PROKKA_08072020/PROKKA_08072020.faa  /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/fna/PROKKA_08072020/oDB/test /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/fna/PROKKA_08072020/tmpdir


source activate py2
python emapper.py -m diamond -i /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/fna/PROKKA_08072020/oDB_rep_seq.fasta --output /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/fna/PROKKA_08072020/emapper -d bact --cpu 16


/work/workspace/zhurj/software/prokka-master/bin/prokka --outdir /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/prokka/MNH04863 --cpus 36 --prefix MNH04863 --kingdom Bacteria /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/fna/MNH04863.fna 
Enterococcus faecium



/work/workspace/zhurj/project/2_swgs/prokka/cm20200806/signalp/MNH05168/MNH05168_short_summary.signalp5

cat /work/workspace/zhurj/members/LZJ/downmetastat20200814/check_ByRun_sralist | awk -F "\t" '{print "du -sh "$2" >> /work/workspace/zhurj/members/LZJ/downmetastat20200814/result_check_filesize.txt "}' > /work/workspace/zhurj/members/LZJ/downmetastat20200814/p/sizestat.sh



export PATH="/work/program/instal/miniconda/bin:$PATH";
source activate HUMAnN
humann2 --input ERR1190619_1.fq -o 2020.May18.HUMAnN/ --threads 12 --metaphlan-options "--mpa_pkl /work/workspace/zhouyj/script/MGWAS/metaphlan2/metaphlan_databases/mpa_v20_m200.pkl --bowtie2db /work/workspace/zhouyj/script/MGWAS/metaphlan2/metaphlan_databases/"


/work/program/bin/megahit
python /work/program/instal/metaphlan2/metaphlan2.py -h

mnid='MNH04863' && cat /work/classify/current/MNH/${mnid:0:6}/${mnid}/${mnid}.all.report | awk -F "\t" 'NR==2{print "taxid:"$1"\nSpecies:"$2"\nRef: "$3"\nANI: "$4"\nCoverage: "$5"\n"}'

$fastpro -i read1 -o $gr_oread1 -I read2 -O $gr_oread2 --poly_g_min_len 10 --poly_x_min_len 10 -q 15 -u 40 -n 5 -l 50 -w $thread 


echo "Christensenella minuta" | taxonkit name2taxid | taxonkit lineage --taxid-field 2
Christensenella minuta  626937  cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;Christensenella;Christensenella minuta
Christensenellaceae     990719  cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Clostridia;Clostridiales;Christensenellaceae
Clostridiaceae  31979   cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Clostridia;Clostridiales;Clostridiaceae

echo "Christensenella" | taxonkit name2taxid | taxonkit lineage --taxid-field 2
Christensenella 990721  cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Clostridia;Clostridiales;Christensenellaceae;Christensenella
taxonkit list --show-rank --show-name --ids 990721 | grep -P "\[species\]" 
626937 [species] Christensenella minuta
    1545742 [species] uncultured Christensenella sp.
  1805714 [species] Christensenella massiliensis
  1816678 [species] Christensenella timonensis
    1851429 [species] Christensenella sp. AF73-05CM02
    1935934 [species] Christensenella sp.
    2086585 [species] Christensenella sp. Marseille-P3954


echo "Christensenellaceae" | taxonkit name2taxid | taxonkit lineage --taxid-field 2
Christensenellaceae     990719  cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Clostridia;Clostridiales;Christensenellaceae
taxonkit list --show-rank --show-name --ids 990719 | grep -P "\[species\]" 
5   626937 [species] Christensenella minuta 
      1545742 [species] uncultured Christensenella sp.
1   1805714 [species] Christensenella massiliensis
2    1816678 [species] Christensenella timonensis
1      1851429 [species] Christensenella sp. AF73-05CM02
      1935934 [species] Christensenella sp.
1      2086585 [species] Christensenella sp. Marseille-P3954
    1229255 [species] uncultured Christensenellaceae bacterium
1    1930012 [species] Christensenellaceae bacterium Phil1
1    1930013 [species] Christensenellaceae bacterium Phil7
1    1930014 [species] Christensenellaceae bacterium Phil10
1    2054177 [species] Christensenellaceae bacterium
    1917875 [species] Beduinibacterium massiliense


/work/workspace/zhurj/lib/sh/gcasoftgz GCA_900155415.1 /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/input
/work/workspace/zhurj/lib/sh/gcasoftgz GCA_900087015.1 /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/input
/work/workspace/zhurj/lib/sh/gcasoftgz GCA_001678845.1 /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/input
/work/workspace/zhurj/lib/sh/gcasoftgz GCA_900604345.1 /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/input
/work/workspace/zhurj/lib/sh/gcasoftgz GCA_001940855.1 /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/input
/work/workspace/zhurj/lib/sh/gcasoftgz GCA_001940845.1 /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/input
/work/workspace/zhurj/lib/sh/gcasoftgz GCA_001940895.1 /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/input
/work/workspace/zhurj/lib/sh/gcasoftgz GCA_902386865.1 /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/input

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/input -b /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/output/Christensenellaceae -t 16


from ete3 import Tree, TreeStyle
file = '/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/output/Christensenellaceae.nwk'
t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/output/Christensenellaceae_rec.png", w=300, units="mm", tree_style=ts)

t = Tree(file)
ts = TreeStyle()
ts.show_leaf_name = True
ts.mode = "c"
ts.arc_start = -180 # 0 degrees = 3 o'clock
ts.arc_span = 270
t.render("/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/output/Christensenellaceae_cir.svg", w=300, h=300, units="mm", tree_style=ts)

fastANI -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/input/GCA_902386865.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/input/GCA_001678845.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/fastANI/GCA_902386865_GCA_001678845
fastANI -q /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/input/GCA_001678845.fna -r /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/input/GCA_902386865.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/fastANI/GCA_001678845_GCA_902386865

export GTDBTK_DATA_PATH=/work/workspace/zhurj/reference/GTDB/release89
gtdbtk classify_wf --batchfile /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/inputf/input10 --out_dir /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/gtdb --cpus 20

iqtree -s /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/gtdb/gtdbtk.bac120.user_msa.fasta -m MFP -b 1000  -T 20  -cmax 15 --prefix /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/iqtree/iqtreecm  --redo --bnni 

raxmlHPC-PTHREADS-SSE3 -s /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/raxml/input/markeralign.phylip -n tree -m PROTGAMMAJTT -T 30 -N 1000 -p 20170808 -f a -x 20170808

from Bio import SeqIO
records = SeqIO.parse("/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/raxml/input/gtdbtk.bac120.user_msa.fasta", "fasta")
count = SeqIO.write(records, "/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/raxml/input/markeralign.phylip", "phylip")
print("Converted %i records" % count)

GTDB-Tk v1.2.0 (database release89)
RAxML version 8.2.12
MUSCLE v3.8.1551

muscle -in /work/workspace/zhurj/project/3_16S/MNH05168_20200820/phylo/input/Efaecium.fa -out /work/workspace/zhurj/project/3_16S/MNH05168_20200820/phylo/input/Efaecium_align.fa -maxiters 16 

muscle -in /work/workspace/zhurj/project/3_16S/MNH05168_20200820/phylo/input/Egallinarum.fa -out /work/workspace/zhurj/project/3_16S/MNH05168_20200820/phylo/input/Egallinarum_align.fa -maxiters 16   


from Bio import SeqIO
records = SeqIO.parse("/work/workspace/zhurj/project/3_16S/MNH05168_20200820/phylo/input/Efaecium_align.fa", "fasta")
count = SeqIO.write(records, "/work/workspace/zhurj/project/3_16S/MNH05168_20200820/phylo/input/Efaecium.phylip", "phylip")
print("Converted %i records" % count)

records = SeqIO.parse("/work/workspace/zhurj/project/3_16S/MNH05168_20200820/phylo/input/Egallinarum_align.fa", "fasta")
count = SeqIO.write(records, "/work/workspace/zhurj/project/3_16S/MNH05168_20200820/phylo/input/Egallinarum.phylip", "phylip")
print("Converted %i records" % count)


raxmlHPC-PTHREADS-SSE3 -s /work/workspace/zhurj/project/3_16S/MNH05168_20200820/phylo/input/Efaecium.phylip -n tree -m GTRGAMMAI -T 30 -N 1000 -p 20170808 -f a -x 20170808

raxmlHPC-PTHREADS-SSE3 -s /work/workspace/zhurj/project/3_16S/MNH05168_20200820/phylo/input/Egallinarum.phylip -n tree -m GTRGAMMAI -T 30 -N 1000 -p 20170808 -f a -x 20170808

muscle -in /work/workspace/zhurj/project/3_16S/MNH05168_20200820/phylo/input/Efaecium_noMRX010.fa -out /work/workspace/zhurj/project/3_16S/MNH05168_20200820/phylo/input/Efaecium_noMRX010_align.fa -maxiters 16  

from Bio import SeqIO
records = SeqIO.parse("/work/workspace/zhurj/project/3_16S/MNH05168_20200820/phylo/input/Efaecium_noMRX010_align.fa", "fasta")
count = SeqIO.write(records, "/work/workspace/zhurj/project/3_16S/MNH05168_20200820/phylo/input/Efaecium_noMRX010.phylip", "phylip")
print("Converted %i records" % count)

raxmlHPC-PTHREADS-SSE3 -s /work/workspace/zhurj/project/3_16S/MNH05168_20200820/phylo/input/Efaecium_noMRX010.phylip -n tree -m GTRGAMMAI -T 30 -N 1000 -p 20170808 -f a -x 20170808

muscle -in /work/workspace/zhurj/project/3_16S/MNH05168_20200820/phylo/input/Efaecium_MRX010rc.fa -out /work/workspace/zhurj/project/3_16S/MNH05168_20200820/phylo/input/Efaecium_MRX010rc_align.fa -maxiters 16  

from Bio import SeqIO
records = SeqIO.parse("/work/workspace/zhurj/project/3_16S/MNH05168_20200820/phylo/input/Efaecium_MRX010rc_align.fa", "fasta")
count = SeqIO.write(records, "/work/workspace/zhurj/project/3_16S/MNH05168_20200820/phylo/input/Efaecium_MRX010rc.phylip", "phylip")
print("Converted %i records" % count)

nohup raxmlHPC-PTHREADS-SSE3 -s /work/workspace/zhurj/project/3_16S/MNH05168_20200820/phylo/input/Efaecium_MRX010rc.phylip -n tree -m GTRGAMMAI -T 30 -N 1000 -p 20170808 -f a -x 20170808 > Efaecium_MRX010rc_log.txt &


/work/workspace/zhurj/software/prokka-master/bin/prokka --outdir /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/prokka/MNH05168 --cpus 36 --prefix MNH05168 --kingdom Bacteria /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/fna/MNH05168.fna 

经典分泌蛋白预测
综合SignalP和TMHMM的分析结果（存在信号肽但无跨膜结构域）
signalp -org gram+ -fasta /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/prokka/MNH05168/MNH05168.faa -format short -prefix /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/signalp/MNH05168/MNH05168_short
/work/workspace/zhurj/software/tmhmm-2.0c/bin/tmhmm /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/prokka/MNH05168/MNH05168.faa > /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/tmhmm/MNH05168_tmhmm.txt 

islandpath /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/prokka/MNH05168/MNH05168.gbk /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/islandpath/gi_test.txt

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/3_16S/MNH05168_20200820/phylo/input/MRXO10.fa -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/3_16S/MNH05168_20200820/phylo/input/MRXO10.out -max_target_seqs 5 -num_threads 16 -outfmt 7 


# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 6 hits found
MRX010  GU983697        99.651  1431    5       0       1       1431    1437    7       0.0     2559
MRX010  CP003504        99.510  1428    7       0       1       1428    1460    33      0.0     2544
MRX010  BCQB01000108    99.510  1428    7       0       1       1428    1460    33      0.0     2544
MRX010  BCQB01000108    78.049  41      9       0       1382    1422    44      84      8.8     34.6
MRX010  AJ301830        99.510  1430    2       5       1       1428    1458    32      0.0     2530
MRX010  AJAN01000023    99.020  1428    14      0       1       1428    1460    33      0.0     2513

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "GU983697|CP003504|BCQB01000108|AJ301830|AJAN01000023"
AJ301830        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus faecium;LMG 11423(T)
AJAN01000023    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus villorum;ATCC 700913(T)
BCQB01000108    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus durans;NBRC 100479(T)
CP003504        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus hirae;ATCC 9790(T)
GU983697        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus lactis;BT159(T)

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/3_16S/MNH05168_20200820/phylo/input/MRX518.fa -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/3_16S/MNH05168_20200820/phylo/input/MRX518.out -max_target_seqs 5 -num_threads 16 -outfmt 7
# Query: AF039900.1 Enterococcus gallinarum 16S ribosomal RNA gene, partial sequence
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 10 hits found
AF039900.1      KE350284        99.932  1461    1       0       1       1461    22      1482    0.0     2631
AF039900.1      BCQE01000074    99.932  1461    1       0       1       1461    22      1482    0.0     2631
AF039900.1      LHOX01000013    99.863  1461    2       0       1       1461    22      1482    0.0     2626
AF039900.1      AF039903        99.795  1461    3       0       1       1461    22      1482    0.0     2622
AF039900.1      NGKU01000001    99.726  1461    4       0       1       1461    22      1482    0.0     2617
# BLASTN 2.9.0+
# Query: MRX518
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 10 hits found
MRX518  KE350284        99.861  1441    1       1       5       1444    24      1464    0.0     2587
MRX518  BCQE01000074    99.861  1441    1       1       5       1444    24      1464    0.0     2587
MRX518  LHOX01000013    99.792  1441    2       1       5       1444    24      1464    0.0     2582
MRX518  AF039903        99.722  1441    3       1       5       1444    24      1464    0.0     2578
MRX518  NGKU01000001    99.653  1441    4       1       5       1444    24      1464    0.0     2573
# BLAST processed 2 queries
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "KE350284|BCQE01000074|LHOX01000013|AF039903|NGKU01000001"
AF039903        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus casseliflavus;MUTK 20(T)
BCQE01000074    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus gallinarum;NBRC 100675(T)
KE350284        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;13.SD.W.09
LHOX01000013    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;LHOX_s;RIT-PI-f
NGKU01000001    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;NGKU_s;8G7_MSG3316

python /work/workspace/zhurj/script/1_tools/10_mnid2assgenome/mnid2assgenome.py --batchfile /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/genome/mnid -o /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/genome/dir -t 2

export GTDBTK_DATA_PATH=/work/workspace/zhurj/reference/GTDB/release89
nohup gtdbtk classify_wf --batchfile /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/genome/dir/mn2assemble.genome.txt --out_dir /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/gtdbtk --cpus 20 > /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/gtdbtk_96_log.txt &

muscle -in /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/phygenus/input/gtdbtk.bac120.user_msa.fasta -out /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/phygenus/input/Efaecium96_align.fa -maxiters 16  

from Bio import SeqIO
records = SeqIO.parse("/work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/phygenus/input/Efaecium96_align.fa", "fasta")
count = SeqIO.write(records, "/work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/phygenus/input/Efaecium96.phylip", "phylip")
print("Converted %i records" % count)

nohup raxmlHPC-PTHREADS-SSE3 -s /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/phygenus/input/Efaecium96.phylip -n tree -m PROTGAMMAJTT -T 30 -N 1000 -p 20170808 -f a -x 20170808 >  Efaecium96.log &


seqkit grep -f /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/secretomes/inids -o /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/secretomes/MNH05168_secretome.faa /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/prokka/MNH05168_IDadj.faa 

from Bio import SeqIO
records = SeqIO.parse("/work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/secretomes/MNH05168_secretome.faa", "fasta")
count = SeqIO.write(records, "/work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/secretomes/MNH05168_secretome.tab", "tab")

source activate py2
python /work/workspace/zhurj/software/eggnog-mapper-master/emapper.py -m diamond -d bact -i /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/secretomes/MNH05168_secretome.faa --output /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/eggNOG/MNH05168_secretomes --cpu 20

export PATH=/work/workspace/zhurj/software/eggnog-mapper-master/bin/:$PATH
python /work/workspace/zhurj/software/eggnog-mapper-master/emapper.py -i /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/secretomes/MNH05168_secretome.faa --output /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/eggNOG/MNH05168_secretomes -d bact --usemem --cpu 20 --override

source activate py2
interproscan
/work/workspace/zhurj/software/my_interproscan/interproscan-5.36-75.0/interproscan.sh -i /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/secretomes/MNH05168_secretome.faa -f tsv -dp -b /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/interpro/MNH05168_interpro



from Bio import SeqIO
records = SeqIO.parse("/work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/prokka/MNH05168_IDadj.faa", "fasta")
count = SeqIO.write(records, "/work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/prokka/MNH05168_IDadj.tab", "tab")

export GTDBTK_DATA_PATH=/work/workspace/zhurj/reference/GTDB/release89
gtdbtk classify_wf --batchfile /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/genome/dir/Efaecium29 --out_dir /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/GTDB_Efaecium --cpus 24

muscle -in /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/GTDB_Efaecium/gtdbtk.bac120.user_msa.fasta -out /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/phylo_Efaecium/input/Efaecium29_align.fa -maxiters 16  

from Bio import SeqIO
records = SeqIO.parse("/work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/phylo_Efaecium/input/Efaecium29_align.fa", "fasta")
count = SeqIO.write(records, "/work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/phylo_Efaecium/input/Efaecium29.phylip", "phylip")
print("Converted %i records" % count)

nohup raxmlHPC-PTHREADS-SSE3 -s /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/phylo_Efaecium/input/Efaecium29.phylip -n tree -m PROTGAMMAJTT -T 30 -N 1000 -p 20170808 -f a -x 20170808 >  Efaecium29.log &

source activate /work/program/instal/miniconda/envs/Prokka && prokka -v
/work/program/instal/miniconda/envs/Prokka/bin/prokka --outdir /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/prokka/NCTC7171 --cpus 36 --prefix NCTC7171 --kingdom Bacteria /work/workspace/zhurj/reference/NCBI/genome/GCA_900447735/genomic.fna 
/work/program/instal/miniconda/envs/Prokka/bin/prokka --outdir /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/prokka20200823 --cpus 36 --prefix MNH05168 --kingdom Bacteria /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna --force

ref: http://blog.sciencenet.cn/blog-2970729-1174911.html
mummer /work/database/ncbi/genome/all/GCA/GCA_900/GCA_900447/GCA_900447735/genomic.fna.gz /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna
# step1: run nucmer for alignment
nucmer --prefix=ref_qry /work/workspace/zhurj/reference/NCBI/genome/GCA_900447735/genomic.fna /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna
nucmer --mum -p /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_NCTC7171/nucmer/MNH05168_NCTC7171 /work/workspace/zhurj/reference/NCBI/genome/GCA_900447735/genomic.fna /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna
delta-filter -i 85 -l 8000 -o 85 /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_NCTC7171/nucmer/MNH05168_NCTC7171.delta > /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_NCTC7171/nucmer/MNH05168_NCTC7171.best_delta
mummerplot -p /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_NCTC7171/nucmer/MNH05168_NCTC7171 /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_NCTC7171/nucmer/MNH05168_NCTC7171.best_delta -t postscript
ps2pdf /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_NCTC7171/nucmer/MNH05168_NCTC7171.ps /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_NCTC7171/nucmer/MNH05168_NCTC7171.pdf
convert -density 300 /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_NCTC7171/nucmer/MNH05168_NCTC7171.pdf /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_NCTC7171/nucmer/MNH05168_NCTC7171.png


PhiSpy.py -o /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/phispy /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/prokka/NCTC7171/NCTC7171.gbk --threads 30 
PhiSpy.py -o /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/phispy /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/prokka20200823/MNH05168_adj.gbk --threads 30
PhiSpy.py -o /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/phispy /work/predict/current/MNH/MNH051/MNH05168/MNH05168.gbf --threads 30

islandpath /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/prokka/NCTC7171/NCTC7171.gbk /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/islandpath/gi_NCTC7171.txt


source activate py2
python /work/workspace/zhurj/software/eggnog-mapper-master/emapper.py -m diamond -d bact -i /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/prokka/MNH05168.faa --output /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/eggNOG_all/MNH05168 --cpu 20

python /work/workspace/zhurj/software/eggnog-mapper-master/emapper.py -m diamond -d bact -i /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/prokka/NCTC7171/NCTC7171.faa --output /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/eggNOG_all/NCTC7171 --cpu 20

signalp -org gram+ -fasta /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/prokka/NCTC7171/NCTC7171.faa -format short -prefix /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/signalp/NCTC7171_sec
/work/workspace/zhurj/software/tmhmm-2.0c/bin/tmhmm /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/prokka/NCTC7171/NCTC7171.faa > /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/tmhmm/NCTC7171_tmhmm.txt 

cdsnum=$(cat /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/prokka/NCTC7171/NCTC7171.tsv | grep -P "\tCDS" |wc -l ) && avlen=$(cat /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/prokka/NCTC7171/NCTC7171.tsv | grep -P "\tCDS" | awk -F "\t" '{print $3}' | awk '{sum+=$1} END {print sum/NR}') && echo -e "NCTC7171\t$cdsnum\t$avlen"

cdsnum=$(cat /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/prokka/MNH05168.tsv | grep -P "\tCDS" |wc -l ) && avlen=$(cat /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/prokka/MNH05168.tsv | grep -P "\tCDS" | awk -F "\t" '{print $3}' | awk '{sum+=$1} END {print sum/NR}') && echo -e "MNH05168\t$cdsnum\t$avlen"

cdsnum=$(cat /work/predict/current/MNH/MNH051/MNH05168/MNH05168.tsv | grep -P "\tCDS" |wc -l ) && avlen=$(cat /work/predict/current/MNH/MNH051/MNH05168/MNH05168.tsv | grep -P "\tCDS" | awk -F "\t" '{print $3}' | awk '{sum+=$1} END {print sum/NR}') && echo -e "MNH05168\t$cdsnum\t$avlen"


quast -t 16 -o /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/quast/NCTC7171 /work/workspace/zhurj/reference/NCBI/genome/GCA_900447735/genomic.fna
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/quast/MNH05168 /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/quast/MNH09296 /work/assembly/current/MNH/MNH092/MNH09296/genomic.fna
cdsnum=$(cat /work/predict/current/MNH/MNH092/MNH09296/MNH09296.tsv | grep -P "\tCDS" |wc -l ) && avlen=$(cat /work/predict/current/MNH/MNH092/MNH09296/MNH09296.tsv | grep -P "\tCDS" | awk -F "\t" '{print $3}' | awk '{sum+=$1} END {print sum/NR}') && echo -e "MNH09296\t$cdsnum\t$avlen"


plasmid
plasmidfinder.py \
-i /work/workspace/zhurj/reference/NCBI/genome/GCA_900447735/genomic.fna \
-mp /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/blastn \
-t 0.6 \
-p /work/workspace/zhurj/bin/miniconda3/envs/python3.6/share/plasmidfinder-2.1-0/database \
-o /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/plasmidfinder


plasmidfinder.py \
-i /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna \
-mp /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/blastn \
-t 0.6 \
-p /work/workspace/zhurj/bin/miniconda3/envs/python3.6/share/plasmidfinder-2.1-0/database \
-o /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/plasmidfinder

vfdb
perl /work/workspace/zhouyj/script/alignKit/quick-blast.pl -q /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/prokka/NCTC7171/NCTC7171.faa -s /work/database/vfdb/current/VFDB_setA_pro.fas -o /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/VFDB/NCTC7171 -t /work/database/vfdb/current/VFDB_setA_pro.fas.tab -m blastp -d 60 
perl /work/workspace/zhouyj/script/alignKit/quick-blast.pl -q /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/prokka/MNH05168.faa -s /work/database/vfdb/current/VFDB_setA_pro.fas -o /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/VFDB/MNH05168 -t /work/database/vfdb/current/VFDB_setA_pro.fas.tab -m blastp -d 60 

perl /work/workspace/zhouyj/script/isolate-strain2analysis/AMR-annotation.pl -i /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/prokka/NCTC7171/NCTC7171.faa -p /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/ANTI/NCTC7171 -t protein
bash /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/ANTI/NCTC7171.sh

perl /work/workspace/zhouyj/script/isolate-strain2analysis/AMR-annotation.pl -i /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/prokka/MNH05168.faa -p /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/ANTI/MNH05168 -t protein
bash /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/ANTI/MNH05168.sh

# seqkit grep -f /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/secretomes/inids -o /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/secretomes/MNH05168_secretome.faa /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/prokka/MNH05168_IDadj.faa 

from Bio import SeqIO
records = SeqIO.parse("/work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/secretomes/MNH05168_secretome.faa", "fasta")
count = SeqIO.write(records, "/work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/secretomes/MNH05168_secretome.tab", "tab")

s

from Bio import SeqIO
records = SeqIO.parse("/work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/prokka/NCTC7171/NCTC7171.faa", "fasta")
count = SeqIO.write(records, "/work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/prokka/NCTC7171/NCTC7171.tab", "tab")
from Bio import SeqIO
records = SeqIO.parse("/work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/prokka/NCTC7171/NCTC7171.tab", "tab")
count = SeqIO.write(records, "/work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/prokka/NCTC7171/NCTC7171_adj.faa", "fasta")

seqkit grep -f /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/secretome/inids -o /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/secretome/NCTC7171_secretome.faa /work/workspace/zhurj/project/2_swgs/prokka/NCTC7171/prokka/NCTC7171/NCTC7171_adj.faa

Kegg 基因标注颜色
https://www.genome.jp/kegg-bin/show_pathway?ko01503+K12973

call snps
ref: http://mummer.sourceforge.net/manual/
#nucmer --prefix=ref_qry /work/workspace/zhurj/reference/NCBI/genome/GCA_900447735/genomic.fna /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna
delta-filter -r -q  /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_NCTC7171/nucmer/MNH05168_NCTC7171.delta > /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_NCTC7171/nucmer/MNH05168_NCTC7171.filter
show-snps -Clr /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_NCTC7171/nucmer/MNH05168_NCTC7171.filter > /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_NCTC7171/nucmer/MNH05168_NCTC7171.snps

/work/workspace/zhurj/software/my_interproscan/interproscan-5.36-75.0/interproscan.sh -i /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/prokka/MNH05168.faa -goterms -pa -f tsv -b /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/interpro/MNH05168_go



/work/assembly/current/MNH/MNH092/MNH09296/genomic.fna
source activate /work/program/instal/miniconda/envs/Prokka 
/work/program/instal/miniconda/envs/Prokka/bin/prokka --outdir /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/prokka/MNH09296 --cpus 36 --prefix MN09296 --kingdom Bacteria /work/assembly/current/MNH/MNH092/MNH09296/genomic.fna

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
plasmidfinder.py \
-i /work/assembly/current/MNH/MNH092/MNH09296/genomic.fna \
-mp /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/blastn \
-t 0.6 \
-p /work/workspace/zhurj/bin/miniconda3/envs/python3.6/share/plasmidfinder-2.1-0/database \
-o /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/plasmidfinder


perl /work/workspace/zhouyj/script/alignKit/quick-blast.pl -q /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/prokka/MNH09296/MNH09296.faa -s /work/database/vfdb/current/VFDB_setA_pro.fas -o /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/VFDB/MNH09296 -t /work/database/vfdb/current/VFDB_setA_pro.fas.tab -m blastp -d 60 

perl /work/workspace/zhouyj/script/isolate-strain2analysis/AMR-annotation.pl -i /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/prokka/MNH09296/MNH09296.faa -p /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/ANTI/MNH09296 -t protein
bash /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/ANTI/MNH09296.sh

signalp -org gram+ -fasta /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/prokka/MNH09296/MNH09296.faa -format short -prefix /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/signalp/MNH09296_signalp
/work/workspace/zhurj/software/tmhmm-2.0c/bin/tmhmm /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/prokka/MNH09296/MNH09296.faa > /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/tmhmm/MNH09296_tmhmm.txt 

source activate py2
python /work/workspace/zhurj/software/eggnog-mapper-master/emapper.py -m diamond -d bact -i /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/prokka/MNH09296/MNH09296.faa --output /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/eggNOG/MNH09296 --cpu 20

source activate py2
/work/workspace/zhurj/software/my_interproscan/interproscan-5.36-75.0/interproscan.sh -i /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/prokka/MNH09296/MNH09296.faa -goterms -pa -f tsv -b /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/interproscan/MNH09296

cat /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/prokka/MNH09296/MNH09296.gbk | awk -F "_" '{if(/LOCUS/) print "LOCUS       MNH09296_"$2"   "$4" bp   DNA    linear"; else print $0} ' > /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/prokka/MNH09296/MNH09296.gbf

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
PhiSpy.py -o /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/phispy /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/prokka/MNH09296/MNH09296.gbf --threads 30
PhiSpy.py -o /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/phispy /work/predict/current/MNH/MNH092/MNH09296/MNH09296.gbf --threads 30

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
islandpath /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/prokka/MNH09296/MNH09296.gbf /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/islandpath/gi.txt

source activate roary
cd /work/workspace/zhurj/project/2_swgs/prokka/corePangene
roary -e --mafft  -p 30 /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/prokka/MNH09296/MNH09296.gff  /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/prokka/MNH05168.gff
roary -e -i 50 --mafft  -p 20 /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/prokka/MNH09296/MNH09296.gff  /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/prokka/MNH05168.gff -f /work/workspace/zhurj/project/2_swgs/prokka/corePangene/MNH05168_9296
roary -e -i 90 --mafft  -p 20 /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/prokka/MNH09296/MNH09296.gff  /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/prokka/MNH05168.gff -f /work/workspace/zhurj/project/2_swgs/prokka/corePangene/MNH05168_9296_90
roary -e -i 50 --mafft  -p 20 /work/workspace/zhurj/project/2_swgs/prokka/MNH22871/MNH22871.gff  /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_20200820/prokka/MNH05168.gff -f /work/workspace/zhurj/project/2_swgs/prokka/corePangene/MNH05168_22871_50

http://json2table.com/
plasmidfinder result to table

perl /work/workspace/moonis/bin/Prokka/prokka.pl -i "/work/workspace/liangyj/project/2020.July3.strain_select/download/OVXM01.1.fsa_nt" -o /work/workspace/liangyj/project/2020.July3.strain_select/download/prokka/ -p OVXM01 -k Bacteria -e 1e-09 -c 8 -f

source activate /work/program/instal/miniconda/envs/Prokka 
/work/program/instal/miniconda/envs/Prokka/bin/prokka --addgenes --increment 1 --gffver 2 --compliant --locustag MNH09296  --centre MoonBioTech --outdir /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/prokka/MNH09296_1 --cpus 36 --prefix MN09296 --kingdom Bacteria --genus Enterococcus --species "Enterococcus faecium" --evalue 1e-09 /work/assembly/current/MNH/MNH092/MNH09296/genomic.fna 

source activate python3.6
PhiSpy.py -o /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/phispy /work/workspace/zhurj/project/2_swgs/prokka/MNH09296/prokka/MNH09296_1/MN09296.gbk --threads 30

source activate /work/program/instal/miniconda/envs/Prokka 
/work/program/instal/miniconda/envs/Prokka/bin/prokka --addgenes --increment 1 --gffver 2 --compliant --locustag MNH08404  --centre MoonBioTech --outdir /work/workspace/zhurj/project/2_swgs/prokka/MNH08404/prokka --cpus 36 --prefix MNH08404 --kingdom Bacteria --genus Enterococcus --species faecium --evalue 1e-09 --force /work/assembly/current/MNH/MNH084/MNH08404/genomic.fna

source activate python3.6
PhiSpy.py -o /work/workspace/zhurj/project/2_swgs/prokka/MNH08404/phispy /work/workspace/zhurj/project/2_swgs/prokka/MNH08404/prokka/MNH08404.gbk --threads 30

plasmidfinder.py \
-i /work/assembly/current/MNH/MNH084/MNH08404/genomic.fna \
-mp /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/blastn \
-t 0.6 \
-p /work/workspace/zhurj/bin/miniconda3/envs/python3.6/share/plasmidfinder-2.1-0/database \
-o /work/workspace/zhurj/project/2_swgs/prokka/MNH08404/plasmidfinder

perl /work/workspace/zhouyj/script/alignKit/quick-blast.pl -q /work/workspace/zhurj/project/2_swgs/prokka/MNH08404/prokka/MNH08404.faa -s /work/database/vfdb/current/VFDB_setA_pro.fas -o /work/workspace/zhurj/project/2_swgs/prokka/MNH08404/VFDB/MNH08404 -t /work/database/vfdb/current/VFDB_setA_pro.fas.tab -m blastp -d 60 

perl /work/workspace/zhouyj/script/isolate-strain2analysis/AMR-annotation.pl -i /work/workspace/zhurj/project/2_swgs/prokka/MNH08404/prokka/MNH08404.faa -p /work/workspace/zhurj/project/2_swgs/prokka/MNH08404/ANTI/MNH08404 -t protein
bash /work/workspace/zhurj/project/2_swgs/prokka/MNH08404/ANTI/MNH08404.sh


/work/assembly/current/MNH/MNH048/MNH04863/genomic.fna
/work/assembly/current/MNH/MNH061/MNH06163/genomic.fna
/work/workspace/zhurj/reference/NCBI/genome/GCA_001678845/genomic.fna
source activate /work/program/instal/miniconda/envs/Prokka 
/work/program/instal/miniconda/envs/Prokka/bin/prokka --addgenes --increment 1 --gffver 2 --compliant --locustag MNH04863  --centre MoonBioTech --outdir /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/prokka --cpus 36 --prefix MNH04863 --kingdom Bacteria --genus Christensenella --species AF73-05CM02 --evalue 1e-09 --force /work/assembly/current/MNH/MNH048/MNH04863/genomic.fna

/work/program/instal/miniconda/envs/Prokka/bin/prokka --addgenes --increment 1 --gffver 2 --compliant --locustag MNH06163  --centre MoonBioTech --outdir /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/prokka --cpus 36 --prefix MNH06163 --kingdom Bacteria --genus Christensenella --species new_sp --evalue 1e-09 --force /work/assembly/current/MNH/MNH061/MNH06163/genomic.fna

/work/program/instal/miniconda/envs/Prokka/bin/prokka --addgenes --increment 1 --gffver 2 --compliant --locustag AF7305CM02  --centre MoonBioTech --outdir /work/workspace/zhurj/project/2_swgs/prokka/AF7305CM02/prokka --cpus 36 --prefix AF7305CM02 --kingdom Bacteria --genus Christensenella --species new_sp --evalue 1e-09 --force /work/workspace/zhurj/reference/NCBI/genome/GCA_001678845/genomic.fna

source activate python3.6
mkdir /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/plasmidfinder
mkdir /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/VFDB
mkdir /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/ANTI

plasmidfinder.py \
-i /work/assembly/current/MNH/MNH048/MNH04863/genomic.fna \
-mp /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/blastn \
-t 0.6 \
-p /work/workspace/zhurj/bin/miniconda3/envs/python3.6/share/plasmidfinder-2.1-0/database \
-o /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/plasmidfinder

perl /work/workspace/zhouyj/script/alignKit/quick-blast.pl -q /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/prokka/MNH04863.faa -s /work/database/vfdb/current/VFDB_setA_pro.fas -o /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/VFDB/MNH04863 -t /work/database/vfdb/current/VFDB_setA_pro.fas.tab -m blastp -d 60 

perl /work/workspace/zhouyj/script/isolate-strain2analysis/AMR-annotation.pl -i /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/prokka/MNH04863.faa -p /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/ANTI/MNH04863 -t protein
bash /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/ANTI/MNH04863.sh

mkdir /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/plasmidfinder
mkdir /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/VFDB
mkdir /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/ANTI

plasmidfinder.py \
-i /work/assembly/current/MNH/MNH061/MNH06163/genomic.fna \
-mp /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/blastn \
-t 0.6 \
-p /work/workspace/zhurj/bin/miniconda3/envs/python3.6/share/plasmidfinder-2.1-0/database \
-o /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/plasmidfinder

perl /work/workspace/zhouyj/script/alignKit/quick-blast.pl -q /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/prokka/MNH06163.faa -s /work/database/vfdb/current/VFDB_setA_pro.fas -o /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/VFDB/MNH06163 -t /work/database/vfdb/current/VFDB_setA_pro.fas.tab -m blastp -d 60 

perl /work/workspace/zhouyj/script/isolate-strain2analysis/AMR-annotation.pl -i /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/prokka/MNH06163.faa -p /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/ANTI/MNH06163 -t protein
bash /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/ANTI/MNH06163.sh


mkdir /work/workspace/zhurj/project/2_swgs/prokka/AF7305CM02/plasmidfinder
mkdir /work/workspace/zhurj/project/2_swgs/prokka/AF7305CM02/VFDB
mkdir /work/workspace/zhurj/project/2_swgs/prokka/AF7305CM02/ANTI

plasmidfinder.py \
-i /work/workspace/zhurj/reference/NCBI/genome/GCA_001678845/genomic.fna \
-mp /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/blastn \
-t 0.6 \
-p /work/workspace/zhurj/bin/miniconda3/envs/python3.6/share/plasmidfinder-2.1-0/database \
-o /work/workspace/zhurj/project/2_swgs/prokka/AF7305CM02/plasmidfinder

perl /work/workspace/zhouyj/script/alignKit/quick-blast.pl -q /work/workspace/zhurj/project/2_swgs/prokka/AF7305CM02/prokka/AF7305CM02.faa -s /work/database/vfdb/current/VFDB_setA_pro.fas -o /work/workspace/zhurj/project/2_swgs/prokka/AF7305CM02/VFDB/AF7305CM02 -t /work/database/vfdb/current/VFDB_setA_pro.fas.tab -m blastp -d 60 

perl /work/workspace/zhouyj/script/isolate-strain2analysis/AMR-annotation.pl -i /work/workspace/zhurj/project/2_swgs/prokka/AF7305CM02/prokka/AF7305CM02.faa -p /work/workspace/zhurj/project/2_swgs/prokka/AF7305CM02/ANTI/AF7305CM02 -t protein
bash /work/workspace/zhurj/project/2_swgs/prokka/AF7305CM02/ANTI/AF7305CM02.sh

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
mkdir /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/islandpath
mkdir /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/islandpath
mkdir /work/workspace/zhurj/project/2_swgs/prokka/AF7305CM02/islandpath
mkdir /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/signalp
mkdir /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/signalp
mkdir /work/workspace/zhurj/project/2_swgs/prokka/AF7305CM02/signalp
mkdir /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/tmhmm
mkdir /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/tmhmm
mkdir /work/workspace/zhurj/project/2_swgs/prokka/AF7305CM02/tmhmm

islandpath /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/prokka/MNH04863.gbk /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/islandpath/gi
signalp -org gram+ -fasta /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/prokka/MNH04863.faa -format short -prefix /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/signalp/signalp
/work/workspace/zhurj/software/tmhmm-2.0c/bin/tmhmm /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/prokka/MNH04863.faa > /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/tmhmm/tmhmm 

islandpath /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/prokka/MNH06163.gbk /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/islandpath/gi
signalp -org gram+ -fasta /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/prokka/MNH06163.faa -format short -prefix /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/signalp/signalp
/work/workspace/zhurj/software/tmhmm-2.0c/bin/tmhmm /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/prokka/MNH06163.faa > /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/tmhmm/tmhmm 

islandpath /work/workspace/zhurj/project/2_swgs/prokka/AF7305CM02/prokka/AF7305CM02.faa /work/workspace/zhurj/project/2_swgs/prokka/AF7305CM02/islandpath/gi
signalp -org gram+ -fasta /work/workspace/zhurj/project/2_swgs/prokka/AF7305CM02/prokka/AF7305CM02.faa -format short -prefix /work/workspace/zhurj/project/2_swgs/prokka/AF7305CM02/signalp/signalp
/work/workspace/zhurj/software/tmhmm-2.0c/bin/tmhmm /work/workspace/zhurj/project/2_swgs/prokka/AF7305CM02/prokka/AF7305CM02.faa > /work/workspace/zhurj/project/2_swgs/prokka/AF7305CM02/tmhmm/tmhmm 


source activate py2
mkdir /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/eggNOG
mkdir /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/eggNOG
mkdir /work/workspace/zhurj/project/2_swgs/prokka/AF7305CM02/eggNOG
mkdir /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/interproscan
mkdir /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/interproscan
mkdir /work/workspace/zhurj/project/2_swgs/prokka/AF7305CM02/interproscan

python /work/workspace/zhurj/software/eggnog-mapper-master/emapper.py -m diamond -d bact -i /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/prokka/MNH04863.faa --output /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/eggNOG/MNH04863 --cpu 20
/work/workspace/zhurj/software/my_interproscan/interproscan-5.36-75.0/interproscan.sh -i /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/prokka/MNH04863.faa -goterms -pa -f tsv -b /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/interproscan/MNH04863
python /work/workspace/zhurj/software/eggnog-mapper-master/emapper.py -m diamond -d bact -i /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/prokka/MNH06163.faa --output /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/eggNOG/MNH06163 --cpu 20
/work/workspace/zhurj/software/my_interproscan/interproscan-5.36-75.0/interproscan.sh -i /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/prokka/MNH06163.faa -goterms -pa -f tsv -b /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/interproscan/MNH06163


source activate /work/program/instal/miniconda/envs/Prokka 
/work/program/instal/miniconda/envs/Prokka/bin/prokka -increment 1 --locustag MNH04863  --centre MoonBioTech --outdir /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/prokka_1 --cpus 36 --prefix MNH04863 --kingdom Bacteria --genus Christensenella --species AF73-05CM02 --evalue 1e-09 --force /work/assembly/current/MNH/MNH048/MNH04863/genomic.fna

/work/program/instal/miniconda/envs/Prokka/bin/prokka --increment 1 --locustag MNH06163  --centre MoonBioTech --outdir /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/prokka_1 --cpus 36 --prefix MNH06163 --kingdom Bacteria --genus Christensenella --species new_sp --evalue 1e-09 --force /work/assembly/current/MNH/MNH061/MNH06163/genomic.fna


source activate roary
cd /work/workspace/zhurj/project/2_swgs/prokka/MNH04863_6163
roary -e --mafft -p 30 /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/prokka_1/MNH04863.gff /work/workspace/zhurj/project/2_swgs/prokka/MNH06163/prokka_1/MNH06163.gff



plasmidfinder.py -i /work/assembly/current/MNH/MNH040/MNH04077/genomic.fna -mp /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/blastn -t 0.6 -p /work/workspace/zhurj/bin/miniconda3/envs/python3.6/share/plasmidfinder-2.1-0/database -o /work/workspace/zhurj/project/2_swgs/prokka/test/plasmidfinder





/work/workspace/zhurj/lib/sh/plasmidfinder /work/workspace/zhurj/project/2_swgs/plasmidfinder/EnterococcusChristensenella_20200826/mnid173 /work/workspace/zhurj/predict/plasmid

perl /work/workspace/zhouyj/script/isolate-strain2analysis/AMR-annotation.pl -i /work/workspace/zhurj/project/2_swgs/prokka/AF7305CM02/prokka/AF7305CM02.faa -p /work/workspace/zhurj/project/2_swgs/prokka/AF7305CM02/ANTI/AF7305CM02 -t protein
bash /work/workspace/zhurj/project/2_swgs/prokka/AF7305CM02/ANTI/AF7305CM02.sh
/work/predict/current/MNH/MNH048/MNH04863/MNH04863.faa 
/work/workspace/zhurj/lib/sh/antipre /work/workspace/zhurj/project/2_swgs/plasmidfinder/EnterococcusChristensenella_20200826/mnid173  /work/workspace/zhurj/predict/ANTI

find /work/workspace/zhurj/predict/ANTI | grep tab | xargs -i cp {} /work/workspace/zhurj/project/2_swgs/ANTI/EnterococcusChristensenella173/in173
cat /work/workspace/zhurj/project/2_swgs/ANTI/EnterococcusChristensenella173/in173/*.tab > /work/workspace/zhurj/project/2_swgs/ANTI/EnterococcusChristensenella173/merge/sam173

find /work/workspace/zhurj/predict/plasmid | grep data | awk -F "/" '{print $0"\t"$(NF-1)}' > /work/workspace/zhurj/project/2_swgs/plasmidfinder/EnterococcusChristensenella_20200826/in173/in173file
python /work/workspace/zhurj/lib/python/script/mergejson.py /work/workspace/zhurj/script/1_tools/test/flist

python /work/workspace/zhurj/lib/python/script/mergejson.py /work/workspace/zhurj/project/2_swgs/plasmidfinder/EnterococcusChristensenella_20200826/in173/in173file > /work/workspace/zhurj/project/2_swgs/plasmidfinder/EnterococcusChristensenella_20200826/merge/sam173

ref: https://docs.antismash.secondarymetabolites.org/command_line/
ref: https://docs.antismash.secondarymetabolites.org/PDFmanual/antiSMASH5manual.pdf
conda activate /work/program/instal/Miniconda3-4.7.12/envs/antismash-5
antismash -c 30 --taxon bacteria --cb-general --asf --pfam2go --genefinding-tool prodigal --cb-knownclusters --cb-subclusters --smcog-trees  --output-dir  /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/antismash /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/prokka/MNH04863.fna

conda activate /work/workspace/liangyj/bin/conda_env/RNASeq
R
library(clusterProfiler)
conda activate /work/program/instal/Miniconda3-4.7.12/envs/antismash-5
antismash -c 30 --taxon bacteria --cb-general  --asf --pfam2go --genefinding-tool prodigal --cf-borders-only --clusterhmmer  --output-dir  MNH04863 /work/assembly/current/MNH/MNH048/MNH04863/genomic.fna
antismash -c 30 --taxon bacteria --cb-general  --asf --pfam2go --genefinding-tool prodigal --cf-borders-only --clusterhmmer  --output-dir  MNH04863 /work/assembly/current/MNH/MNH048/MNH04863/genomic.fna


quast -t 16 -o /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/quast/ /work/workspace/zhurj/project/2_swgs/prokka/MNH04863/prokka/MNH04863_circos.fna

source activate circos
../bin/circos -conf etc/circos.conf -debug_group summary,timer > run.out
source activate python3.6
1. input create
perl /work/workspace/zhurj/script/5_wgspro/1_filter/mnid2rawdata.pl -i /work/workspace/zhurj/project/1_metadata/metapipe/input/ifile -t 0 -r metagenome -o /work/workspace/zhurj/project/1_metadata/metapipe/input -n test

/work/rawdata/run/guangzhou/novogene/2020/07/20200714/run00057/rawdata

2. data filter
perl /work/workspace/zhurj/script/5_wgspro/1_filter/datafilter.pl -i /work/workspace/zhurj/project/1_metadata/metapipe/input/test_r1r2_2c.txt -o /work/workspace/zhurj/project/1_metadata/metapipe

srun -o mnc00230.out -e mnc00230.err -N 1 -c 20 -p slurm256 bash mnc00230.sh &
srun -o x02.out -e x02.err -N 1 -c 20 -p slurm128 bash x02_20200110 &

3. contig assemble
cd /work/workspace/zhurj/project/1_metadata/metapipe/p
srun -o mnc00230.out -e mnc00230.err -N 1 -c 20 -p slurm256 bash mnc00230.sh
content in mnc00230.sh
```
python /work/workspace/moonis/bin/genoassemble/soft/SPAdes-3.13.1-Linux/bin/spades.py --meta -1 /work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00230.clean.1.fq.gz -2 /work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00230.clean.2.fq.gz -t 36 -o /work/workspace/zhurj/project/1_metadata/metapipe/spades/MNC00230
```

选择使用bwa
4. bowtie map
ln -s /work/workspace/zhurj/project/1_metadata/metapipe/spades/MNC00230/contigs.fasta /work/workspace/zhurj/project/1_metadata/metapipe/bowtie/db/bt2.fa
bowtie2-build --threads 40 --quiet /work/workspace/zhurj/project/1_metadata/metapipe/spades/MNC00230/contigs.fasta /work/workspace/zhurj/project/1_metadata/metapipe/bowtie/db/bt2 
bowtie2 -p 40 --fast -x /work/workspace/zhurj/project/1_metadata/metapipe/bowtie/db/bt2 -1 /work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00230.clean.1.fq.gz -2 /work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00230.clean.2.fq.gz -S /work/workspace/zhurj/project/1_metadata/metapipe/bowtie/MNC00230.sam

4. bwa map
ln -s /work/workspace/zhurj/project/1_metadata/metapipe/spades/MNC00230/contigs.fasta /work/workspace/zhurj/project/1_metadata/metapipe/bwa/db/ref.fa
bwa index /work/workspace/zhurj/project/1_metadata/metapipe/bwa/db/ref.fa
bwa mem -t 40 /work/workspace/zhurj/project/1_metadata/metapipe/bwa/db/ref.fa /work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00230.clean.1.fq.gz /work/workspace/zhurj/project/1_metadata/metapipe/clean/MNC00230.clean.2.fq.gz > /work/workspace/zhurj/project/1_metadata/metapipe/bwa/MNC00230.sam
samtools view -@ 40 -S -b /work/workspace/zhurj/project/1_metadata/metapipe/bwa/MNC00230.sam -o /work/workspace/zhurj/project/1_metadata/metapipe/bwa/MNC00230.bam
samtools sort -@ 40 /work/workspace/zhurj/project/1_metadata/metapipe/bwa/MNC00230.bam -o /work/workspace/zhurj/project/1_metadata/metapipe/bwa/MNC00230.sorted.bam

5. MetaBAT2 binning
ref: https://bitbucket.org/berkeleylab/metabat/src/master/
runMetaBat.sh -o /work/workspace/zhurj/project/1_metadata/metapipe/metabat /work/workspace/zhurj/project/1_metadata/metapipe/bwa/db/ref.fa /work/workspace/zhurj/project/1_metadata/metapipe/bwa/MNC00230.sorted.bam # does not work

jgi_summarize_bam_contig_depths --outputDepth /work/workspace/zhurj/project/1_metadata/metapipe/metabat/depth.txt  /work/workspace/zhurj/project/1_metadata/metapipe/bwa/MNC00230.sorted.bam
--minContigLength   arg  The mimimum length of contig to include for mapping and shredding
--minContigDepth    arg  The minimum depth along contig at which to break the contig
--outputGC          arg  The file to print the gc coverage histogram
--gcWindow          arg  The sliding window size for GC calculations

metabat2 -i /work/workspace/zhurj/project/1_metadata/metapipe/bwa/db/ref.fa -a /work/workspace/zhurj/project/1_metadata/metapipe/metabat/depth.txt -o /work/workspace/zhurj/project/1_metadata/metapipe/metabat/bin/MNHC00230 -m 2000 
6. checkm check MAG
CheckM v1.1.1
checkm lineage_wf /work/workspace/zhurj/project/1_metadata/metapipe/metabat/bin /work/workspace/zhurj/project/1_metadata/metapipe/checkm -x fna -t 36 --pplacer_threads 


代谢的数据库
宏基因组数据分析参考文章：
https://www.jianshu.com/p/52e9b63740d9

metaphlan 参考
https://blog.csdn.net/woodcorpse/article/details/80031620
ref: chrome-extension://dagcmkpagjlhakfdhnbomgmjdpkdklff/enhanced-reader.html?openApp&pdf=https%3A%2F%2Fwww.nature.com%2Farticles%2Fs41598-019-50299-6.pdf
A multilocus sequence typing scheme of Pseudomonas putida for clinical and environmental isolates

echo "Peudomonas chlororaphis subsp. Piscium" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t

quast -t 16 -o /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/PC60/quast/ /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/rawdata/PC60.fna
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/DSM21509/quast /work/database/ncbi/genome/all/GCA/GCA_003/GCA_003850/GCA_003850345/genomic.fna.gz
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/quast/DSM19603 /work/database/ncbi/genome/all/GCA/GCA_003/GCA_003851/GCA_003851835/genomic.fna.gz
fastANI -q /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/rawdata/PC60.fna -r /work/database/ncbi/genome/all/GCA/GCA_003/GCA_003850/GCA_003850345/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/fastANI/PC60_DSM21509
fastANI -q /work/database/ncbi/genome/all/GCA/GCA_003/GCA_003850/GCA_003850345/genomic.fna.gz -r /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/rawdata/PC60.fna  -t 16 -o /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/fastANI/DSM21509_PC60
fastANI -q /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/rawdata/PC60.fna -r /work/database/ncbi/genome/all/GCA/GCA_003/GCA_003851/GCA_003851835/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/fastANI/PC60_DSM19603
fastANI -q /work/database/ncbi/genome/all/GCA/GCA_003/GCA_003851/GCA_003851835/genomic.fna.gz -r /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/rawdata/PC60.fna  -t 16 -o /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/fastANI/DSM19603_PC60

ln -s /work/workspace/zhurj/reference/NCBI/genome/GCA_003850345/genomic.fna /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/fna/DSM21509.fna
ln -s /work/workspace/zhurj/reference/NCBI/genome/GCA_003851835/genomic.fna /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/fna/DSM19630.fna

gtdbtk classify_wf --batchfile /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/input/gdtbin --out_dir /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/gtdb --cpus 36


cd /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/assemble
seqtk subseq all.fna  /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/input/gene/28847 > /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/assemble/gene/s28847.fna 
seqtk subseq all.fna  /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/input/gene/211250 > /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/assemble/gene/s211250.fna 
seqtk subseq all.fna  /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/input/gene/217420 > /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/assemble/gene/s217420.fna 
seqtk subseq all.fna  /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/input/gene/218886 > /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/assemble/gene/s218886.fna 
seqtk subseq all.fna  /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/input/gene/218932 > /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/assemble/gene/s218932.fna 
seqtk subseq all.fna  /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/input/gene/219556 > /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/assemble/gene/s219556.fna 
seqtk subseq all.fna  /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/input/gene/221574 > /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/assemble/gene/s221574.fna 
seqtk subseq all.fna  /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/input/gene/PC60 > /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/assemble/gene/sPC60.fna 


cd /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/assemble/merge
union -filter /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/assemble/gene/s28847.fna >  s28847.fna 
union -filter /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/assemble/gene/s211250.fna > s211250.fna
union -filter /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/assemble/gene/s217420.fna > s217420.fna
union -filter /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/assemble/gene/s218886.fna > s218886.fna
union -filter /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/assemble/gene/s218932.fna > s218932.fna
union -filter /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/assemble/gene/s219556.fna > s219556.fna
union -filter /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/assemble/gene/s221574.fna > s221574.fna
union -filter /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/assemble/gene/sPC60.fna >   sPC60.fna

cd /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/assemble
rm merge8sam4gene.fna -f
awk '{if( $0 ~ /^>/) print ">s28847";else print $0}' ./merge/s28847.fna  >> merge8sam4gene.fna
awk '{if( $0 ~ /^>/) print ">s211250";else print $0}' ./merge/s211250.fna >> merge8sam4gene.fna
awk '{if( $0 ~ /^>/) print ">s217420";else print $0}' ./merge/s217420.fna >> merge8sam4gene.fna
awk '{if( $0 ~ /^>/) print ">s218886";else print $0}' ./merge/s218886.fna >> merge8sam4gene.fna
awk '{if( $0 ~ /^>/) print ">s218932";else print $0}' ./merge/s218932.fna >> merge8sam4gene.fna
awk '{if( $0 ~ /^>/) print ">s219556";else print $0}' ./merge/s219556.fna >> merge8sam4gene.fna
awk '{if( $0 ~ /^>/) print ">s221574";else print $0}' ./merge/s221574.fna >> merge8sam4gene.fna
awk '{if( $0 ~ /^>/) print ">sPC60";else print $0}' ./merge/sPC60.fna >> merge8sam4gene.fna

muscle -in /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/assemble/merge8sam4gene.fna  -out /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/muscle/sam8gene4.fna -maxiters 16 

iqtree -s /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/muscle/sam8gene4.fna -m MFP -b 1000  -T 20  -cmax 15 --prefix /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/iqtree/sam8gene4  --redo --bnni 



from Bio import SeqIO
records = SeqIO.parse("/work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/muscle/sam8gene4.fna", "fasta")
count = SeqIO.write(records, "/work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/muscle/sam8gene4.phylip", "phylip")

# raxmlHPC-PTHREADS-SSE3 -s /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/muscle/sam8gene4.phylip -n tree -m GTRGAMMAI -T 30 -N 1000 -p 20170808 -f a -x 20170808 # anino acid
raxmlHPC-PTHREADS-SSE3 -s /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/muscle/sam8gene4.phylip -n tree -m PROTGAMMAJTT -T 30 -N 1000 -p 20170808 -f a -x 20170808 # nucleotide




GTDB-Tk v1.2.0 (database release89)
RAxML version 8.2.12
MUSCLE v3.8.1551  


from Bio import SeqIO
records = SeqIO.parse("/work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/gtdb/gtdbtk.bac120.user_msa.fasta", "fasta")
count = SeqIO.write(records, "/work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/gtdb/three.phylip", "phylip")

cd /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/raxml
raxmlHPC-PTHREADS-SSE3 -s /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/gtdb/three.phylip -n tree -m PROTGAMMAJTT -T 30 -N 1000 -p 20170808 -f a -x 20170808

iqtree -s /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/gtdb/gtdbtk.bac120.user_msa.fasta -m MFP -b 1000  -T 20  -cmax 15 --prefix /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/iqtree/threesam  --redo --bnni 

source activate python3.6
union /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/assemble/gene


Genomic DNA was extracted using a commercially avaliable kit (TIANamp Bacteria DNA Kit) following the manufacturer’s instructions. One microgram of genomic DNA was used for fragmentation and 350-bp library preparaion using a NEBNext®Ultra™ DNA Library Prep Kit (Illumina, NEB, USA) was sequenced on Illumina NovaSeq platform (2 x 150-bp paired-end[PE] reads). For each strain, ~1 Gbp of sequence data was required. Quality Control and data filter were made with FastQC and fastp, respectively. Clean reads from each strain were assembled with SPAdes v3.13.1. Thereafter, assembled genome was estimated with CheckM v1.1.1 using the 'lineage_wf' workflow (qualified genome with contamination < 5% and completeness > 90%). For qualified genomes, taxonomic annotation were performed with GTDB-Tk version 1.2.0. using the 'classify_wf' function and default parameters (ANI ≥ 95% and AF ≥ 60%, ANI: average nucleotide identity, AF: alignment fragment).
ref: https://mra.asm.org/content/9/18/e00334-20
ref: https://mra.asm.org/content/8/42/e00978-19
ref: https://plantmethods.biomedcentral.com/articles/10.1186/s13007-016-0152-4


GCF_001941925

source activate python3.6
mlst --label Anthrax  /work/database/ncbi/genome/all/GCA/GCA_001/GCA_001941/GCA_001941925/genomic.fna.gz
mlst 



blastn  -query /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/assemble/merge/sPC60.fna  -subject /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/rawdata/PC60.fna  -outfmt 6 

blastn  -query /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/assemble/all.fna -subject  /work/workspace/zhurj/software/mlst-master/db/blast/mlst.fa   -outfmt 6 -max_target_seqs 3

from Bio import SeqIO
records = SeqIO.parse("/work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/muscle/sam8gene4.fna", "fasta")
count = SeqIO.write(records, "//work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/muscle/sam8gene4.tab", "tab")

三代测序数据下载

/work/predict/orfs/antismash/MNH/MNH049/MNH04997/E1/L1/S1
MNH21576

/work/workspace/zhurj/lib/sh/mergeantismash /work/workspace/zhurj/temp/mergeantismash/infile /work/workspace/zhurj/temp/mergeantismash/ofile

source activate python3.6 
ln -s /work/assembly/current/MN/MN231/MN231152/genomic.fna /work/workspace/zhurj/project/2_swgs/agriculture/kejingMNH231152/input/MN231152.fna
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/agriculture/kejingMNH231152/quast /work/assembly/current/MN/MN231/MN231152/genomic.fna
checkm lineage_wf /work/assembly/current/MN/MN231/MN231152 /work/workspace/zhurj/project/2_swgs/agriculture/kejingMNH231152/checkm -x fna -t 36 --tab_table -f /work/workspace/zhurj/project/2_swgs/agriculture/kejingMNH231152/checkm/MNH231152_checkm.txt
gtdbtk classify_wf --genome_dir /work/workspace/zhurj/project/2_swgs/agriculture/kejingMNH231152/input/ --out_dir /work/workspace/zhurj/project/2_swgs/agriculture/kejingMNH231152/gtdb --cpus 36
fastANI -q /work/workspace/zhurj/project/2_swgs/agriculture/kejingMNH231152/input/MN231152.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000333/GCA_000333655/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/agriculture/kejingMNH231152/fastANI/MN231152_GCA_000333655
fastANI -q /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000333/GCA_000333655/genomic.fna.gz -r /work/workspace/zhurj/project/2_swgs/agriculture/kejingMNH231152/input/MN231152.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/agriculture/kejingMNH231152/fastANI/GCA_000333655_MN231152
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/199/985/GCA_014199985.1_ASM1419998v1/
/work/workspace/zhurj/reference/NCBI/genome/GCA_014199985

fastANI -q /work/workspace/zhurj/project/2_swgs/agriculture/kejingMNH231152/input/MN231152.fna -r /work/workspace/zhurj/reference/NCBI/genome/GCA_014199985/GCA_014199985.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/agriculture/kejingMNH231152/fastANI/MN231152_GCA_014199985
fastANI -q /work/workspace/zhurj/reference/NCBI/genome/GCA_014199985/GCA_014199985.fna.gz -r /work/workspace/zhurj/project/2_swgs/agriculture/kejingMNH231152/input/MN231152.fna -t 16 -o /work/workspace/zhurj/project/2_swgs/agriculture/kejingMNH231152/fastANI/GCA_014199985_MN231152
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/agriculture/kejingMNH231152/GCA_014199985/quast /work/workspace/zhurj/reference/NCBI/genome/GCA_014199985/GCA_014199985.fna.gz
checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/GCA_014199985 /work/workspace/zhurj/project/2_swgs/agriculture/kejingMNH231152/GCA_014199985/checkm -x fna.gz -t 36 --tab_table -f /work/workspace/zhurj/project/2_swgs/agriculture/kejingMNH231152/GCA_014199985/checkm/GCA_014199985_checkm.txt


checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/GCA_001678845 /work/workspace/zhurj/project/12_checkm/GCA/001/678/845 -x fna -t 36 --tab_table -f /work/workspace/zhurj/project/12_checkm/GCA/001/678/845/checkm.txt
echo "Christensenella" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t

fastANI -q /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/rawdata/PC60.fna -r /work/database/ncbi/genome/all/GCA/GCA_003/GCA_003850/GCA_003850345/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/zhaobo_MLST20200902/fastANI/PC60_DSM21509

GCA_009735405 /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735405/genomic.fna.gz
GCA_009735495 /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735495/genomic.fna.gz
GCA_009735475 /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735475/genomic.fna.gz
GCA_013867815 /work/database/ncbi/genome/all/GCA/GCA_013/GCA_013867/GCA_013867815/genomic.fna.gz
GCA_009735445 /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735445/genomic.fna.gz
GCF_009735435 /work/database/ncbi/genome/all/GCA/GCF_009/GCA_009735/GCA_009735435/genomic.fna.gz

fastANI -q /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna 
fastANI -o /work/workspace/zhurj/project/2_swgs/MNH05168/fastANI/GCA_009735405 -q /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735405/genomic.fna.gz -t 16 
fastANI -o /work/workspace/zhurj/project/2_swgs/MNH05168/fastANI/GCA_009735495 -q /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735495/genomic.fna.gz -t 16 
fastANI -o /work/workspace/zhurj/project/2_swgs/MNH05168/fastANI/GCA_009735475 -q /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735475/genomic.fna.gz -t 16 
fastANI -o /work/workspace/zhurj/project/2_swgs/MNH05168/fastANI/GCA_013867815 -q /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_013/GCA_013867/GCA_013867815/genomic.fna.gz -t 16 
fastANI -o /work/workspace/zhurj/project/2_swgs/MNH05168/fastANI/GCA_009735445 -q /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735445/genomic.fna.gz -t 16 
fastANI -o /work/workspace/zhurj/project/2_swgs/MNH05168/fastANI/GCA_009735435 -q /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735435/genomic.fna.gz -t 16 

fastANI -o /work/workspace/zhurj/project/2_swgs/MNH05168/fastANI/GCA_009735405_rev -r /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna -q /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735405/genomic.fna.gz -t 16 
fastANI -o /work/workspace/zhurj/project/2_swgs/MNH05168/fastANI/GCA_009735495_rev -r /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna -q /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735495/genomic.fna.gz -t 16 
fastANI -o /work/workspace/zhurj/project/2_swgs/MNH05168/fastANI/GCA_009735475_rev -r /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna -q /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735475/genomic.fna.gz -t 16 
fastANI -o /work/workspace/zhurj/project/2_swgs/MNH05168/fastANI/GCA_013867815_rev -r /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna -q /work/database/ncbi/genome/all/GCA/GCA_013/GCA_013867/GCA_013867815/genomic.fna.gz -t 16 
fastANI -o /work/workspace/zhurj/project/2_swgs/MNH05168/fastANI/GCA_009735445_rev -r /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna -q /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735445/genomic.fna.gz -t 16 
fastANI -o /work/workspace/zhurj/project/2_swgs/MNH05168/fastANI/GCA_009735435_rev -r /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna -q /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735435/genomic.fna.gz -t 16

fastANI -o /work/workspace/zhurj/project/2_swgs/MNH05168/fastANI/GCA_900447735 -q /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_900/GCA_900447/GCA_900447735/genomic.fna.gz -t 16 
fastANI -o /work/workspace/zhurj/project/2_swgs/MNH05168/fastANI/GCA_001544255 -q /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_001/GCA_001544/GCA_001544255/genomic.fna.gz -t 16 

fastANI -o /work/workspace/zhurj/project/2_swgs/MNH05168/fastANI/GCA_900447735_rev -r /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna -q /work/database/ncbi/genome/all/GCA/GCA_900/GCA_900447/GCA_900447735/genomic.fna.gz -t 16 
fastANI -o /work/workspace/zhurj/project/2_swgs/MNH05168/fastANI/GCA_001544255_rev -r /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna -q /work/database/ncbi/genome/all/GCA/GCA_001/GCA_001544/GCA_001544255/genomic.fna.gz -t 16 

fastANI -o /work/workspace/zhurj/project/2_swgs/MNH05168/fastANI/MNH09296_MNH05168 -r /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna -q /work/assembly/fasta/genome/MNH/MNH092/MNH09296/E1/L1/S1/MNH09296.fas -t 16 
fastANI -o /work/workspace/zhurj/project/2_swgs/MNH05168/fastANI/MNH05168_MNH09296 -r /work/assembly/fasta/genome/MNH/MNH092/MNH09296/E1/L1/S1/MNH09296.fas -q /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna -t 16 

fastANI -o /work/workspace/zhurj/project/2_swgs/MNH05168/fastANI/MNH22871_MNH05168 -r /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna -q /work/assembly/fasta/genome/MNH/MNH228/MNH22871/E2/L1/S1/MNH22871.fas -t 16 
fastANI -o /work/workspace/zhurj/project/2_swgs/MNH05168/fastANI/MNH05168_MNH22871 -r /work/assembly/fasta/genome/MNH/MNH228/MNH22871/E2/L1/S1/MNH22871.fas -q /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna -t 16 

quast -t 16 -o /work/workspace/zhurj/project/2_swgs/agriculture/kejingMNH231152/GCA_014199985/quast /work/workspace/zhurj/reference/NCBI/genome/GCA_014199985/GCA_014199985.fna.gz
checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/GCA_014199985 /work/workspace/zhurj/project/2_swgs/agriculture/kejingMNH231152/GCA_014199985/checkm -x fna.gz -t 36 --tab_table -f /work/workspace/zhurj/project/2_swgs/agriculture/kejingMNH231152/GCA_014199985/checkm/GCA_014199985_checkm.txt

mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/quast/GCA_009735405
mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/quast/GCA_009735495
mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/quast/GCA_009735475
mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/quast/GCA_013867815
mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/quast/GCA_009735445
mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/quast/GCA_009735435
mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/quast/GCA_006337145
mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/quast/GCA_900447735
mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/quast/GCA_001544255

mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_009735405
mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_009735495
mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_009735475
mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_013867815
mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_009735445
mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_009735435
mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_006337145
mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_900447735
mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_001544255

quast -t 16 -o /work/workspace/zhurj/project/2_swgs/MNH05168/quast/GCA_009735405  /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735405/genomic.fna.gz
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/MNH05168/quast/GCA_009735495  /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735495/genomic.fna.gz
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/MNH05168/quast/GCA_009735475  /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735475/genomic.fna.gz
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/MNH05168/quast/GCA_013867815  /work/database/ncbi/genome/all/GCA/GCA_013/GCA_013867/GCA_013867815/genomic.fna.gz
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/MNH05168/quast/GCA_009735445  /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735445/genomic.fna.gz
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/MNH05168/quast/GCA_009735435  /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735435/genomic.fna.gz
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/MNH05168/quast/GCA_006337145  /work/database/ncbi/genome/all/GCA/GCA_006/GCA_006337/GCA_006337145/genomic.fna.gz
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/MNH05168/quast/GCA_900447735  /work/database/ncbi/genome/all/GCA/GCA_900/GCA_900447/GCA_900447735/genomic.fna.gz
quast -t 16 -o /work/workspace/zhurj/project/2_swgs/MNH05168/quast/GCA_001544255  /work/database/ncbi/genome/all/GCA/GCA_001/GCA_001544/GCA_001544255/genomic.fna.gz

mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/fna/GCA_009735405
mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/fna/GCA_009735495
mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/fna/GCA_009735475
mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/fna/GCA_013867815
mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/fna/GCA_009735445
mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/fna/GCA_009735435
mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/fna/GCA_006337145
mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/fna/GCA_900447735
mkdir /work/workspace/zhurj/project/2_swgs/MNH05168/fna/GCA_001544255

rm /work/workspace/zhurj/project/2_swgs/MNH05168/fna/GCA_009735405/* -f
rm /work/workspace/zhurj/project/2_swgs/MNH05168/fna/GCA_009735495/* -f
rm /work/workspace/zhurj/project/2_swgs/MNH05168/fna/GCA_009735475/* -f
rm /work/workspace/zhurj/project/2_swgs/MNH05168/fna/GCA_013867815/* -f
rm /work/workspace/zhurj/project/2_swgs/MNH05168/fna/GCA_009735445/* -f
rm /work/workspace/zhurj/project/2_swgs/MNH05168/fna/GCA_009735435/* -f
rm /work/workspace/zhurj/project/2_swgs/MNH05168/fna/GCA_006337145/* -f
rm /work/workspace/zhurj/project/2_swgs/MNH05168/fna/GCA_900447735/* -f
rm /work/workspace/zhurj/project/2_swgs/MNH05168/fna/GCA_001544255/* -f

cp /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735405/genomic.fna.gz /work/workspace/zhurj/reference/NCBI/genome/GCA_009735405 | gunzip /work/workspace/zhurj/reference/NCBI/genome/GCA_009735405/genomic.fna.gz
cp /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735495/genomic.fna.gz /work/workspace/zhurj/reference/NCBI/genome/GCA_009735495 | gunzip /work/workspace/zhurj/reference/NCBI/genome/GCA_009735495/genomic.fna.gz
cp /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735475/genomic.fna.gz /work/workspace/zhurj/reference/NCBI/genome/GCA_009735475 | gunzip /work/workspace/zhurj/reference/NCBI/genome/GCA_009735475/genomic.fna.gz
cp /work/database/ncbi/genome/all/GCA/GCA_013/GCA_013867/GCA_013867815/genomic.fna.gz /work/workspace/zhurj/reference/NCBI/genome/GCA_013867815 | gunzip /work/workspace/zhurj/reference/NCBI/genome/GCA_013867815/genomic.fna.gz
cp /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735445/genomic.fna.gz /work/workspace/zhurj/reference/NCBI/genome/GCA_009735445 | gunzip /work/workspace/zhurj/reference/NCBI/genome/GCA_009735445/genomic.fna.gz
cp /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735435/genomic.fna.gz /work/workspace/zhurj/reference/NCBI/genome/GCA_009735435 | gunzip /work/workspace/zhurj/reference/NCBI/genome/GCA_009735435/genomic.fna.gz
cp /work/database/ncbi/genome/all/GCA/GCA_006/GCA_006337/GCA_006337145/genomic.fna.gz /work/workspace/zhurj/reference/NCBI/genome/GCA_006337145 | gunzip /work/workspace/zhurj/reference/NCBI/genome/GCA_006337145/genomic.fna.gz
cp /work/database/ncbi/genome/all/GCA/GCA_001/GCA_001544/GCA_001544255/genomic.fna.gz /work/workspace/zhurj/reference/NCBI/genome/GCA_001544255 | gunzip /work/workspace/zhurj/reference/NCBI/genome/GCA_001544255/genomic.fna.gz


checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/GCA_009735405 /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_009735405 -x fna -t 36 --tab_table -f /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_009735405/checkm.txt
checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/GCA_009735495 /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_009735495 -x fna -t 36 --tab_table -f /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_009735495/checkm.txt
checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/GCA_009735475 /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_009735475 -x fna -t 36 --tab_table -f /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_009735475/checkm.txt
checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/GCA_013867815 /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_013867815 -x fna -t 36 --tab_table -f /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_013867815/checkm.txt
checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/GCA_009735445 /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_009735445 -x fna -t 36 --tab_table -f /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_009735445/checkm.txt
checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/GCA_009735435 /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_009735435 -x fna -t 36 --tab_table -f /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_009735435/checkm.txt
checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/GCA_006337145 /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_006337145 -x fna -t 36 --tab_table -f /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_006337145/checkm.txt
checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/GCA_900447735 /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_900447735 -x fna -t 36 --tab_table -f /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_900447735/checkm.txt
checkm lineage_wf /work/workspace/zhurj/reference/NCBI/genome/GCA_001544255 /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_001544255 -x fna -t 36 --tab_table -f /work/workspace/zhurj/project/2_swgs/MNH05168/checkm/GCA_001544255/checkm.txt

mkdir /work/workspace/zhurj/reference/NCBI/genome/GCA_009735405
mkdir /work/workspace/zhurj/reference/NCBI/genome/GCA_009735495
mkdir /work/workspace/zhurj/reference/NCBI/genome/GCA_009735475
mkdir /work/workspace/zhurj/reference/NCBI/genome/GCA_013867815
mkdir /work/workspace/zhurj/reference/NCBI/genome/GCA_009735445
mkdir /work/workspace/zhurj/reference/NCBI/genome/GCA_009735435
mkdir /work/workspace/zhurj/reference/NCBI/genome/GCA_006337145
mkdir /work/workspace/zhurj/reference/NCBI/genome/GCA_001544255

RGI 
# Resistance Gene Identifier - 5.1.1
source activate rgi
rgi main --input_sequence /work/workspace/zhurj/project/2_swgs/MNH04863/assemble/MNH-4863.fna --output_file /work/workspace/zhurj/project/2_swgs/MNH04863/RGI/gri_result --input_type contig --clean


MNO-863
MNO-163
export GTDBTK_DATA_PATH=/work/workspace/zhurj/reference/GTDB/release89
gtdbtk classify_wf --batchfile /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/inputf/input10 --out_dir /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/gtdb --cpus 20

iqtree -s /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/gtdb/gtdbtk.bac120.user_msa.fasta -m MFP -b 1000  -T 20  -cmax 15 --prefix /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/iqtree/iqtreecm  --redo --bnni 

```
parallel 'mkdir -p {1}/{2}/{3}' ::: /work/workspace/zhurj/project/2_swgs/phylotree20200304/ChristensenellaceaeZP20200925 ::: MNH04863 MNH06163 ::: gtdb iqtree
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
export GTDBTK_DATA_PATH=/work/workspace/zhurj/reference/GTDB/release89
parallel --xapply -j 1 'gtdbtk classify_wf --batchfile {1}/inputf/{2}  --out_dir {1}/{2}/gtdb --cpus 20' \
::: /work/workspace/zhurj/project/2_swgs/phylotree20200304/ChristensenellaceaeZP20200925 ::: MNH04863 MNH06163

parallel --xapply -j 1 'iqtree -s {1}/{2}/gtdb/gtdbtk.bac120.user_msa.fasta -m MFP -b 1000  -T 20  -cmax 15 --prefix {1}/{2}/iqtree/iqtreecm --redo --bnni' \
::: /work/workspace/zhurj/project/2_swgs/phylotree20200304/ChristensenellaceaeZP20200925 ::: MNH04863 MNH06163
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

```
cd /work/workspace/zhurj/project/2_swgs/phylotree20200304/ChristensenellaceaeZP20200925/p
time srun -o piqtree.out -e piqtree.err -N 1 -c 20 -p slurm256 bash piqtree.sh &

'''
raxmlHPC-PTHREADS-SSE3 -s /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/raxml/input/markeralign.phylip -n tree -m PROTGAMMAJTT -T 30 -N 1000 -p 20170808 -f a -x 20170808

from Bio import SeqIO
records = SeqIO.parse("/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/raxml/input/gtdbtk.bac120.user_msa.fasta", "fasta")
count = SeqIO.write(records, "/work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenellaceae20200819/raxml/input/markeralign.phylip", "phylip")
print("Converted %i records" % count)
'''


from Bio import SeqIO
records = SeqIO.parse("uniprot-proteome_UP000030151.fasta", "fasta")
count = SeqIO.write(records, "uniprot-proteome_UP000030151.tab", "tab")

from Bio import SeqIO
records = SeqIO.parse("genelist.fa.txt", "fasta")
count = SeqIO.write(records, "genelist.tab", "tab")

from Bio import SeqIO
records = SeqIO.parse("sequence.fasta", "fasta")
count = SeqIO.write(records, "sequence.tab", "tab")

cat uniprot-proteome_UP000030151.tab | grep -P "A0A0A1UXM7|A0A014NF81|A0A0A1URT6|A0A0A1V5H9|A0A0A1UZC9|A0A0A1UVR5|A0A014P734|A0A0A1USV3"
>tr|A0A0A1UXM7|A0A0A1UXM7_9HYPO DNA-directed RNA polymerase subunit OS=Metarhizium robertsii OX=568076 GN=X797_004950 PE=3 SV=1
>tr|A0A014NF81|A0A014NF81_9HYPO DNA-directed RNA polymerase subunit beta OS=Metarhizium robertsii OX=568076 GN=X797_005628 PE=3 SV=1
>tr|A0A0A1URT6|A0A0A1URT6_9HYPO DNA-directed RNA polymerase subunit OS=Metarhizium robertsii OX=568076 GN=X797_007939 PE=3 SV=1
>tr|A0A0A1V5H9|A0A0A1V5H9_9HYPO DNA-directed RNA polymerase subunit beta OS=Metarhizium robertsii OX=568076 GN=X797_000192 PE=3 SV=1
>tr|A0A0A1UZC9|A0A0A1UZC9_9HYPO DNA-directed RNA polymerase subunit beta OS=Metarhizium robertsii OX=568076 GN=X797_003079 PE=3 SV=1
>tr|A0A0A1UVR5|A0A0A1UVR5_9HYPO DNA-directed RNA polymerase subunit OS=Metarhizium robertsii OX=568076 GN=X797_004772 PE=3 SV=1
>tr|A0A014P734|A0A014P734_9HYPO DNA-directed RNA polymerase subunit OS=Metarhizium robertsii OX=568076 GN=X797_008303 PE=3 SV=1
>tr|A0A0A1USV3|A0A0A1USV3_9HYPO DNA-directed RNA polymerase subunit OS=Metarhizium robertsii OX=568076 GN=X797_007188 PE=3 SV=1
A0A0A1UXM7
A0A014NF81
A0A0A1URT6
A0A0A1V5H9
A0A0A1UZC9
A0A0A1UVR5
A0A014P734
A0A0A1USV3
source activate python3.6
cutadapt --discard-untrimmed -g $FORWARD $INPUT 2> /dev/null | cutadapt --discard-untrimmed -a $REVERSE - 2> /dev/null > $OUTPUT
cutadapt --discard-untrimmed -g GGTCCCTTCGGTCAGCTCTTCC genelist.fa  2> /dev/null | cutadapt --discard-untrimmed -a GACCCTAAGAACATGATGGCTG - 2> /dev/null
cutadapt --discard-untrimmed -g GC[CT]CC[CT]GG[ATC]CA[CT]CGTGA[CT]TT[CT]AT genelist.fa  2> /dev/null | cutadapt --discard-untrimmed -a GACCCTAAGAACATGATGGCTG - 2> /dev/null
cutadapt --discard-untrimmed -g GGTCCCTTCGGTCAGCTCTTCC genelist.fa  2> /dev/null | cutadapt --discard-untrimmed -a GGAGATGATGCA[AG]CTTCCGCTCAA - 2> /dev/null
cutadapt --discard-untrimmed -g GGTCCCTTCGGTCAGCTCTTCC genelist.fa  2> /dev/null | cutadapt --discard-untrimmed -a GACCCTAAGAACATGATGGCTG - 2> /dev/null
cutadapt --discard-untrimmed 
cutadapt --discard-untrimmed -a CAGCCATCATGTTCTTAGGGTC genelist.fa 
cd /work/workspace/zhurj/project/14_coworker/chenjuan/Metarhizium_robertsii20200926
perl /work/workspace/zhurj/software/in_silico_pcr-master/in_silico_PCR.pl -s sequence.fasta -a GGTCCCTTCGGTCAGCTCTTCC -b GACCCTAAGAACATGATGGCTG
cutadapt --discard-untrimmed -b GGTCCCTTCGGTCAGCTCTTCC genelist.fa

makeblastdb -in sequence.fasta -dbtype nucl -input_type fasta -parse_seqids -out sequence
makeblastdb -in genelist.fa -dbtype nucl -input_type fasta -parse_seqids -out genelist
blastn -query ./input/in.fa -db sequence -out test_results
blastn -query ./input/in.fa -db genelist -out test_results

>BenA-F
GGTCCCTTCGGTCAGCTCTTCC
>BenA-R
GACCCTAAGAACATGATGGCTG
>EF-F
GC[CT]CC[CT]GG[ATC]CA[CT]CGTGA[CT]TT[CT]AT
>EF-R
CA[GA]AC[CT]GT[CT]GC[CT]GT[CT]GGTGTCAT
>RPB1-F
CG[AG]AC[AC][CT]T[AG]CC[CT]CATTTCACAA
>RPB1-R
GGAGATGATGCA[AG]CTTCCGCTCAA
>RPB2-F
GA[CT]GA[CT][AC]G[AT]GATCA[CT]TT[CT]GG
>RPB2-R
ATGGG[CT]AA[AG]CAAGC[TC]ATGGG

blastn -query ./input/assemble_gene.fa -db sequence -out output/results_4genes
blastn -query ./input/assemble_gene.fa -db sequence -out output/results_4genes

parallel -j 1 --link 'blastn -query {1}/input/gene4_all.fa -db {1}/sequence -out {1}/output/results_4genes_all' ::: /work/workspace/zhurj/project/14_coworker/chenjuan/Metarhizium_robertsii20200926

samtools faidx sequence.fasta JELW01000001.1:1360109-1360709
>JELW01000001.1:1360109-1360709
CGGAAAGAGTGAGCACCCCGGCTGGTCAGGGGGGCGAAGCCGACCATGAAGAAGTGCAGA
CGAGGGAAGGGGACCATGTTGACAGCCAGCTTACGCAGATCAGAGTTCAACTGACCGGGG
AAACGCAAGCATGTGGTGACGCCAGACATGACGGCAGAGACAAGATAGTTCAGGTCACCG
TACGAAGGGTTAGACAGCTTGAGAGTGCGCATGCAGATGTCGTACAGAGCCTCATTGTCG
ATGCAGAAAGTCTCGTCAGAGTTCTCAACGAGCTGATGGACGGAGAGGGTTGCGTTGTAG
GGCTCGACAACGGTGTCGGAAACCTTGGGAGAGGGAACGACGGAGAAAGTGGCCATCATT
CGGTCGGGAAACTCTTCACGGATCTTGGAGATCAACAGAGTACCCATACCAGCACCGGTA
CCACCACCGAGAGAGTGGGTGATCTGGAAGCCCTGGAGGCAGTCACAACCTTCCGCCTCG
CGACGGACAACATCAAGGACATTGTCGACAAGCTCAGCACCTTCAGTGTAGTGACCCTTG
GCCCAGTTGTTGCCAGCACCAGACTGACCAAAGACGAAGTTGTCGGGACGGAAAAGCTGA
C


source activate /work/workspace/zhurj/bin/miniconda3/envs/kraken2
/work/workspace/zhurj/software/linuxnd/linuxnd login -u X101SC20042206-Z01-J011 -p 2ysb6w5y
/work/workspace/zhurj/software/linuxnd/linuxnd list
/work/workspace/zhurj/software/linuxnd/linuxnd list oss://CP2019011700024 
/work/workspace/zhurj/software/linuxnd/linuxnd cp -d oss://CP2019011700024/H101SC20042206/RSCS7100/X101SC20042206-Z01/X101SC20042206-Z01-J011 /work/workspace/zhurj/rawdata/rawdata/run/guangzhou/novogene/2020/9/20200928


ln -s  /work/assembly/current/MNH/MNH048/MNH04863/genomic.fna /work/workspace/zhurj/project/2_swgs/MNH04863/Christensenella/gna/MNH04863.fna

cd /work/workspace/zhurj/project/2_swgs/MNH04863/Christensenella/input
cat mnid77 | awk '{print substr($0,1,6)}' > mnid77_6
cd /work/workspace/zhurj/project/2_swgs/MNH04863/Christensenella/input
parallel --link 'ln -s {1}/{3}/{2}/genomic.fna {4}/{2}.fna ' ::: /work/assembly/current/MNH :::: mnid77 :::: mnid77_6 ::: /work/workspace/zhurj/project/2_swgs/MNH04863/Christensenella/gna

config: /work/workspace/zhurj/bin/miniconda3/envs/phylophlan/lib/python3.8/site-packages/phylophlan/phylophlan_configs
# -d /work/workspace/zhurj/bin/miniconda3/envs/phylophlan/lib/python3.8/site-packages/phylophlan/phylophlan_databases/phylophlan
#phylophlan -i gna -o result -d /work/workspace/zhurj/bin/miniconda3/envs/phylophlan/lib/python3.8/site-packages/phylophlan/phylophlan_databases/phylophlan --diversity low --nproc 20 \
# -f supermatrix_aa.cfg \
# --genome_extension fna \
# -t a

christensenella77.sh
'''
source activate /work/workspace/zhurj/bin/miniconda3/envs/phylophlan
parallel --link 'phylophlan -i {1}/gna -o {1}/result -d phylophlan --diversity low --nproc 20 \
-f supermatrix_aa.cfg \
--genome_extension fna \
-t a' ::: /work/workspace/zhurj/project/2_swgs/MNH04863/Christensenella
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/phylophlan

'''
cd /work/workspace/zhurj/project/2_swgs/MNH04863/Christensenella/p
time srun -o christensenella77.out -e christensenella77.err -N 1 -c 20 -p slurm256 -w mnclient02 bash christensenella77.sh &

echo "Oscillospira" | taxonkit name2taxid | taxonkit lineage --taxid-field 2
echo "Christensenella" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t
Acetatifactor
Acutalibacter
Anaerotruncus
Clostridium
Dorea
Erysipelatoclostridium
Escherichia
Lachnospiraceae
Lactobacillus
Lactococcus
Mucispirillum
Oscillibacter

echo "Acetatifactor" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t
echo "Acutalibacter" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t
echo "Anaerotruncus" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t
echo "Clostridium" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t
echo "Dorea" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t
echo "Erysipelatoclostridium" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t
echo "Escherichia" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t
echo "Lachnospiraceae" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t
echo "Lactobacillus" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t
echo "Lactococcus" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t
echo "Mucispirillum" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t
echo "Oscillibacter" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t
echo "Escherichia" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t


fastANI -q /work/assembly/current/MNH/MNH048/MNH04863/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_001/GCA_001652/GCA_001652705/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/MNH04863/Christensenella_minuta/fastANI_CM_MNH04863
fastANI -r /work/assembly/current/MNH/MNH048/MNH04863/genomic.fna -q /work/database/ncbi/genome/all/GCA/GCA_001/GCA_001652/GCA_001652705/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/MNH04863/Christensenella_minuta/fastANI_MNH04863_CM
  Christensenella minuta  Similarity  Coverage
MNH04863  GCA_001652705.1 83.4453 69.82%

2020.11.16 单菌全基因基因组数据download
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
perl /work/workspace/zhurj/script/git/downloadNovogene/ddnovogene.pl -i oss://CP2019011700024/H101SC20042206/RSCS7100/X101SC20042206-Z01/X101SC20042206-Z01-J012/ -u X101SC20042206-Z01-J012 -w xagtn30c -o /work/workspace/zhurj/rawdata/rawdata/run/guangzhou/novogene/2020/11/20201116
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

cd /work/workspace/zhurj/download/rawdata
time srun -o novo_swgs_20201116.out -e novo_swgs_20201116.err -N 1 -c 20 -p slurm128 -w mnclient03 bash novo_swgs_20201116.sh &
slurm download 失败

MNH05168.type.report
MNH05168.all.report


import pandas as pd
import os
import sys
import subprocess
import linecache
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter

"""

"""

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
group = parser.add_mutually_exclusive_group()
group.add_argument("--dir",help="Input the directory",type=str)
group.add_argument("--batchfile",help="Input file",type=str)
parser.add_argument("-o",help="Output directory",required=True,type=str)
parser.add_argument("--line",help="Line selected",default=2,type=int)
parser.add_argument("--file_suffix",help="File suffix when --dir",default='type.report')
args = parser.parse_args()

"""
indir = '/work/classify/species/classify/MNH'
odir = '/work/workspace/zhurj/project/2_swgs/1_strain_classify_type_sum'
line_number = 2
"""

typereportf = con = tmpcon = ''
odir = args.o
line_number = args.line
ofile = os.path.join(odir,'typereport_summary')
search_suffix = args.file_suffix

if args.dir:
  indir = args.dir
  typereportf = os.path.join(odir,'typereport_list')
  con = "find {} | grep {} > {} ".format(indir,search_suffix,typereportf)
  subprocess.run(con,shell=True)

if args.batchfile:
  typereportf = args.batchfile


mnid = eid = lid = sid = ''
files = []
with open(ofile,'w',encoding='utf-8') as ofp, open(typereportf,'r',encoding='utf-8') as fp:
  con = 'MNID\t\tEid\tLid\tSid\tTaxID\tScientificName\tRefGCA\tANI\tCoverage\n'
  ofp.write(con)
  for line in fp.readlines():
    infile = line.strip()
    if not infile:
      continue
    [*tmplist,mnid,eid,lid,sid,_] = infile.strip().split('/')
    tmpcon = linecache.getline(infile, line_number).strip()
    con = "{}\t{}\t{}\t{}\t{}\n".format(mnid,eid,lid,sid,tmpcon)
    ofp.write(con)

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/21_multifile_defineLine_merge/merge_multifile_defineline.py --dir /work/classify/species/classify/MNH -o /work/workspace/zhurj/project/2_swgs/1_strain_classify_type_sum --line 2 --file_suffix type.report
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/21_multifile_defineLine_merge/merge_multifile_defineline.py --batchfile /work/workspace/zhurj/project/2_swgs/1_strain_classify_type_sum/typereport_list -o /work/workspace/zhurj/project/2_swgs/1_strain_classify_type_sum --line 2
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/21_multifile_defineLine_merge/merge_multifile_defineline.py --dir /work/classify/species/classify/MNH -o /work/workspace/zhurj/project/2_swgs/1_strain_classify_info --line 2 --file_suffix all.report --name allreport
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/21_multifile_defineLine_merge/merge_multifile_defineline.py --dir /work/cleandata/fastq/genome/evalue/MNH -o /work/workspace/zhurj/project/2_swgs/1_strain_classify_info --line 2 --file_suffix checkM.report --name checkm


20210112
# file:
cat /work/workspace/zhurj/project/2_swgs/summary2020/genus_list |taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t > /work/workspace/zhurj/project/2_swgs/summary2020/strain_info_from_genera_20210112.txt
cat /work/workspace/zhurj/project/2_swgs/summary2020/input/knownspecies_20210115 |taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species -o /work/workspace/zhurj/project/2_swgs/summary2020/knownspecies_levels_20210115.txt
cat /work/workspace/zhurj/project/2_swgs/summary2020/input/new_genus_20210115 |taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species -o /work/workspace/zhurj/project/2_swgs/summary2020/new_genus_levels_20210115.txt
cat /work/workspace/zhurj/project/2_swgs/summary2020/input/new_species_20210115 |taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species -o /work/workspace/zhurj/project/2_swgs/summary2020/new_species_levels_20210115.txt

# one genus
echo "[Clostridium]" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/10_mnid2assgenome/mnid2assgenome.py -o /work/workspace/zhurj/project/2_swgs/GTDB/run00068/ingenome --batchfile /work/workspace/zhurj/project/2_swgs/GTDB/run00068/input/mnlist_run00068 --type 2

ln -s /work/assembly/fasta/genome/MNH/MNH051/MNH05168/E2/L1/S1/MNH05168.fas /work/workspace/zhurj/project/2_swgs/GTDB/run00068/genome/MNH05168_30.fna
ln -s /work/assembly/fasta/genome/MNH/MNH051/MNH05168/E3/L1/S1/MNH05168.fas /work/workspace/zhurj/project/2_swgs/GTDB/run00068/genome/MNH05168_54.fna

export GTDBTK_DATA_PATH=/work/workspace/zhurj/reference/GTDB/release89
gtdbtk classify_wf --batchfile /work/workspace/zhurj/project/2_swgs/GTDB/run00068/ingenome/run00068_genome.txt --out_dir /work/workspace/zhurj/project/2_swgs/GTDB/run00068/output --cpus 20

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniELSmnGC.py -i /work/workspace/zhurj/project/2_swgs/fastANI/input/input -o /work/workspace/zhurj/project/2_swgs/fastANI --name MNH05168e2e3
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/11_fastANIfromMNgca/aniELSmnGC.py -i /work/workspace/zhurj/project/2_swgs/fastANI/input/run00068 -o /work/workspace/zhurj/project/2_swgs/fastANI --name run00068_MNH22871

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/fastANI -q /work/assembly/current/MNH/MNH260/MNH26050/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000020/GCA_000020225/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH26050_AKK

import pandas as pd
import os
import sys

def column_count(df,col_index,col_key,new_col_name):
  select_cols = [col_index,col_key]
  new_df = df[select_cols].drop_duplicates()
  select_df = new_df[col_index].value_counts().to_frame()
  select_df.columns = [new_col_name]
  return select_df

def column_sum_count(df,col_index,col_key):
  select_cols = [col_index,col_key]
  select_df = df[select_cols].groupby([col_index]).sum()
  return select_df


infile = '/work/workspace/zhurj/project/2_swgs/sum_2020/strain_classify_20210118'
odir = ''
define_col = 'phylum'
merge_df = pd.DataFrame()
df = pd.read_csv(infile,sep='\t',header=0,index_col=None)

if define_col == 'phylum':
  new_df = column_count(df,define_col,'class','class_num')
  merge_df = new_df
  new_df = column_count(df,define_col,'order','order_num')
  merge_df = merge_df.join(new_df)
  new_df = column_count(df,define_col,'family','family_num')
  merge_df = merge_df.join(new_df)
  new_df = column_count(df,define_col,'genus','genus_num')
  merge_df = merge_df.join(new_df)
  new_df = column_count(df,define_col,'species_name','species_num')
  merge_df = merge_df.join(new_df)
  tmpcount_df = column_sum_count(df,define_col,'strain_num')
  merge_df = merge_df.join(tmpcount_df)

if define_col == 'class':
  new_df = column_count(df,define_col,'order','order_num')
  merge_df = new_df
  new_df = column_count(df,define_col,'family','family_num')
  merge_df = merge_df.join(new_df)
  new_df = column_count(df,define_col,'genus','genus_num')
  merge_df = merge_df.join(new_df)
  new_df = column_count(df,define_col,'species_name','species_num')
  merge_df = merge_df.join(new_df)
  tmpcount_df = column_sum_count(df,define_col,'strain_num')
  merge_df = merge_df.join(tmpcount_df)

if define_col == 'order':
  new_df = column_count(df,define_col,'family','family_num')
  merge_df = new_df
  new_df = column_count(df,define_col,'genus','genus_num')
  merge_df = merge_df.join(new_df)
  new_df = column_count(df,define_col,'species_name','species_num')
  merge_df = merge_df.join(new_df)
  tmpcount_df = column_sum_count(df,define_col,'strain_num')
  merge_df = merge_df.join(tmpcount_df)

if define_col == 'family':
  new_df = column_count(df,define_col,'genus','genus_num')
  merge_df = new_df
  new_df = column_count(df,define_col,'species_name','species_num')
  merge_df = merge_df.join(new_df)
  tmpcount_df = column_sum_count(df,define_col,'strain_num')
  merge_df = merge_df.join(tmpcount_df)

if define_col == 'genus':
  new_df = column_count(df,define_col,'species_name','species_num')
  merge_df = new_df
  tmpcount_df = column_sum_count(df,define_col,'strain_num')
  merge_df = merge_df.join(tmpcount_df)

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/3_swgspro/1_classify_sum/classify_levels_sum.py -i /work/workspace/zhurj/project/2_swgs/sum_2020/input/strain_classify_20210118 -o /work/workspace/zhurj/project/2_swgs/sum_2020/output --level phylum
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/3_swgspro/1_classify_sum/classify_levels_sum.py -i /work/workspace/zhurj/project/2_swgs/sum_2020/input/strain_classify_20210118 -o /work/workspace/zhurj/project/2_swgs/sum_2020/output --level class
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/3_swgspro/1_classify_sum/classify_levels_sum.py -i /work/workspace/zhurj/project/2_swgs/sum_2020/input/strain_classify_20210118 -o /work/workspace/zhurj/project/2_swgs/sum_2020/output --level order
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/3_swgspro/1_classify_sum/classify_levels_sum.py -i /work/workspace/zhurj/project/2_swgs/sum_2020/input/strain_classify_20210118 -o /work/workspace/zhurj/project/2_swgs/sum_2020/output --level family
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/3_swgspro/1_classify_sum/classify_levels_sum.py -i /work/workspace/zhurj/project/2_swgs/sum_2020/input/strain_classify_20210118 -o /work/workspace/zhurj/project/2_swgs/sum_2020/output --level genus




/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/10_mnid2assgenome/mnid2assgenome.py --batchfile /work/workspace/zhurj/project/2_swgs/sum_2020/input/new_genus_species_mnid  -t 1 -o /work/workspace/zhurj/project/2_swgs/sum_2020/input 

/work/workspace/zhurj/bin/miniconda3/bin/mash sketch -o /work/workspace/zhurj/project/2_swgs/sum_2020/newgenus_species_mash/new_genus_species118 -l /work/workspace/zhurj/project/2_swgs/sum_2020/input/mn2assemble.genome.txt
/work/workspace/zhurj/bin/miniconda3/bin/mash dist -t /work/workspace/zhurj/project/2_swgs/sum_2020/newgenus_species_mash/new_genus_species118.msh /work/workspace/zhurj/project/2_swgs/sum_2020/newgenus_species_mash/new_genus_species118.msh > /work/workspace/zhurj/project/2_swgs/sum_2020/newgenus_species_mash/new_genus_species118.dist

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenella20200305/MNH04863/input -b /work/workspace/zhurj/project/2_swgs/phylotree20200304/Christensenella20200305/MNH04863/output/MNH04863 -t 16

new genus species mash dist to cluster 118 strains

library(dplyr)
library(Matrix)

infile = '/work/workspace/zhurj/project/2_swgs/sum_2020/newgenus_species_mash/new_genus_species_mnid.dist'
ofile = '/work/workspace/zhurj/project/2_swgs/sum_2020/newgenus_species_mash/new_genus_species_mnid_3cols'
df = read.table(infile,header=TRUE,sep='\t',row.names=1)
row_names = rownames(df)
mm = data.matrix(df)
mm[upper.tri(mm)] = 0
smm = Matrix(mm)
re = summary(smm)
sm = as.data.frame(re)
arrange(sm,i,j)
new_df <- mutate(sm,row=row_names[i],col=row_names[j])
select_cols <- c('row','col','x')
final_df <- new_df[,select_cols]
write.table(final_df,file=ofile,sep='\t',col.names=FALSE,row.names=FALSE)


# 毒力基因组统计
find /work/workspace/liangyj/project/2021.Jan13.VFDB/MNH | grep filter.info | xargs du -ch # 输出文件大小
find /work/workspace/liangyj/project/2021.Jan13.VFDB/MNH | grep filter.info | xargs du -ch >  /work/workspace/zhurj/project/2_swgs/sum_2020/input/vfdb_files

import pandas as pd
import os
import sys

infile = '/work/workspace/liangyj/project/2021.Jan13.VFDB/MNH/MNH049/MNH04936/E1/L1/S1/MNH04936.VFDB.blastp.filter.info'
infile = '/work/workspace/liangyj/project/2021.Jan13.VFDB/MNH/MNH323/MNH32303/E1/L1/S1/MNH32303.VFDB.blastp.filter.infoles'
df = pd.read_csv(infile,sep='\t',header=None,index_col=None)
df.columns = ['gene','vfgid','identity_info','gene_desc']
tmp_df = df['identity_info'].str.split(':',expand=True)
df['identity'] = tmp_df[3]

infile = '/work/workspace/zhurj/project/2_swgs/sum_2020/input/test'
ofile = '/work/workspace/zhurj/project/2_swgs/sum_2020/vfdb_sum/vfdb_merge'
con = ''
mnid = eid = lid = sid = ''
geneid = vfgid = tmpstr = geneinfo = identity = ''

with open(infile,'r',encoding='utf-8') as fp, open(ofile,'w',encoding='utf-8') as ofp:
  con = 'MNID\tEid\tLid\tSid\tGene\tVFGid\tIdentity\tVFG_description\n'
  ofp.write(con)
  
  for line in fp.readlines():
    [size,filepath] = line.strip().split('\t')
    [*tmp,mnid,eid,lid,sid,_] = filepath.split('/')
    if size == '0':
      con = '{}\t{}\t{}\t{}\t-\t-\t-\t-\n'.format(mnid,eid,lid,sid)
      ofp.write(con)
      continue
    with open(filepath,'r',encoding='utf-8') as fp_2:
      for line_2 in fp_2.readlines():
        [geneid,vfgid,tmpstr,geneinfo] = line_2.strip().split('\t')
        [*tmp2,identity] = tmpstr.split(':')
        con = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(mnid,eid,lid,sid,geneid,vfgid,identity,geneinfo)
        ofp.write(con)

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/23_merge_VFDB_info/merge_vfdb_info.py -i /work/workspace/zhurj/project/2_swgs/sum_2020/input//vfdb_files -o /work/workspace/zhurj/project/2_swgs/sum_2020/vfdb_sum/vfdb_merge

统计VFDB
import pandas as pd
import os

infile = '/work/workspace/zhurj/project/2_swgs/sum_2020/vfdb_sum/vfdb_merge'
odir = '/work/workspace/zhurj/project/2_swgs/sum_2020/vfdb_sum/'
sample_antigene_countf = os.path.join(odir,'sample_antigene_count_sum')
no_antigenef = os.path.join(odir,'sample_no_antigene_sum')
gene_sample_sumf = os.path.join(odir,'gene_sample_sum')

df = pd.read_csv(infile,sep='\t',header=0,index_col=None)

no_antigene_df = df[df['VFGid'] == '-']
no_antigene_df.to_csv(no_antigenef,sep='\t',index=False,header=True)
with_antigene_df = df[df['VFGid'] != '-']
select_cols = ['MNID','Eid','Lid','Sid','VFGid','VFG_description']
with_antigene_part_df = with_antigene_df[select_cols]
with_antigene_filter = with_antigene_part_df.drop_duplicates()
with_antigene_filter['count'] = 1
sample_count_cols = ['MNID','Eid','Lid','Sid','count']
gene_df = with_antigene_filter[sample_count_cols]
group_cols = ['MNID','Eid','Lid','Sid']
gene_sample_count = gene_df.groupby(group_cols).sum()
sample_gene_sum_df = gene_sample_count.sort_values(by='count',axis=0,ascending=False)
sample_gene_sum_df.to_csv(sample_antigene_countf,sep='\t',index=True,header=True)

gene_count_cols = ['VFG_description','count']
gene_df = with_antigene_filter[gene_count_cols]
gene_sample_df = gene_df.groupby('VFG_description').sum()
gene_sample_sum_df = gene_sample_df.sort_values(by='count',axis=0,ascending=False)
gene_sample_sum_df.to_csv(gene_sample_sumf,sep='\t',index=True,header=True)


/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/23_merge_VFDB_info/summry_vfdb_info.py -i /work/workspace/zhurj/project/2_swgs/sum_2020/vfdb_sum/vfdb_merge -o /work/workspace/zhurj/project/2_swgs/sum_2020/vfdb_sum


# 抗性基因组统计
find /work/workspace/liangyj/project/2021.Jan13.CARD/MNH | grep CARD.tab | xargs du -ch | grep -v total > /work/workspace/zhurj/project/2_swgs/sum_2020/input/card_files

import pandas as pd
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
import os
import sys

"""
infile = '/work/workspace/zhurj/project/2_swgs/sum_2020/input/card_files'
ofile = '/work/workspace/zhurj/project/2_swgs/sum_2020/card_sum/card_merge'
"""
parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-i',help='Input file',required=True)
parser.add_argument('-o',help='Output file',required=True)
args = parser.parse_args()

infile = args.i
ofile = args.o

con = ''
mnid = eid = lid = sid = ''
geneid = ago = identity = besthit = drug = ''

with open(infile,'r',encoding='utf-8') as fp, open(ofile,'w',encoding='utf-8') as ofp:
  con = 'MNID\tEid\tLid\tSid\tGene\tAGO\tIdentity\tBesthit\tDrug\n'
  ofp.write(con)
  
  for line in fp.readlines():
    [size,filepath] = line.strip().split('\t')
    [*tmp,mnid,eid,lid,sid,_] = filepath.split('/')
    if size == '0':
      con = '{}\t{}\t{}\t{}\t-\t-\t-\t-\t-\n'.format(mnid,eid,lid,sid)
      ofp.write(con)
      continue
    with open(filepath,'r',encoding='utf-8') as fp_2:
      lines = fp_2.readlines()
      if len(lines) == 1:
        con = '{}\t{}\t{}\t{}\t-\t-\t-\t-\t-\n'.format(mnid,eid,lid,sid)
        ofp.write(con)
        continue
      for line_2 in lines[1:]:
        [geneid,ago,identity,besthit,drug] = line_2.strip().split('\t')
        con = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(mnid,eid,lid,sid,geneid,ago,identity,besthit,drug)
        ofp.write(con)

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/24_merge_card_info/merge_card_info.py -i /work/workspace/zhurj/project/2_swgs/sum_2020/input/card_files -o /work/workspace/zhurj/project/2_swgs/sum_2020/card_sum/card_merge

统计Card
import pandas as pd
import os
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-i',help='Input file',required=True)
parser.add_argument('-o',help='Directory of output',required=True)
args = parser.parse_args()


#infile = '/work/workspace/zhurj/project/2_swgs/sum_2020/card_sum/card_merge'
#odir = '/work/workspace/zhurj/project/2_swgs/sum_2020/card_sum'
infile = args.i
odir = args.o

sample_antigene_countf = os.path.join(odir,'sample_antigene_count_sum')
no_antigenef = os.path.join(odir,'sample_no_antigene_sum')
gene_sample_sumf = os.path.join(odir,'gene_sample_sum')

df = pd.read_csv(infile,sep='\t',header=0,index_col=None)

no_antigene_df = df[df['Besthit'] == '-']
no_antigene_df.to_csv(no_antigenef,sep='\t',index=False,header=True)

with_antigene_df = df[df['Besthit'] != '-']
select_cols = ['MNID','Eid','Lid','Sid','Besthit']
with_antigene_part_df = with_antigene_df[select_cols]
with_antigene_filter = with_antigene_part_df.drop_duplicates()
with_antigene_filter['count'] = 1
sample_count_cols = ['MNID','Eid','Lid','Sid','count']
gene_df = with_antigene_filter[sample_count_cols]
group_cols = ['MNID','Eid','Lid','Sid']
gene_sample_count = gene_df.groupby(group_cols).sum()
sample_gene_sum_df = gene_sample_count.sort_values(by='count',axis=0,ascending=False)
sample_gene_sum_df.to_csv(sample_antigene_countf,sep='\t',index=True,header=True)

gene_count_cols = ['Besthit','count']
gene_df = with_antigene_filter[gene_count_cols]
gene_sample_df = gene_df.groupby('Besthit').sum()
gene_sample_sum_df = gene_sample_df.sort_values(by='count',axis=0,ascending=False)
gene_sample_sum_df.to_csv(gene_sample_sumf,sep='\t',index=True,header=True)

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/24_merge_card_info/summary_card_info.py -i /work/workspace/zhurj/project/2_swgs/sum_2020/card_sum/card_merge -o /work/workspace/zhurj/project/2_swgs/sum_2020/card_sum


MOONIS summary
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/3_swgspro/1_classify_sum/classify_levels_sum.py -i /work/workspace/zhurj/project/2_swgs/sum_2020/input/moonis_mnh -o /work/workspace/zhurj/project/q2_swgs/sum_2020/moonis/mnh --level phylum
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/3_swgspro/1_classify_sum/classify_levels_sum.py -i /work/workspace/zhurj/project/2_swgs/sum_2020/input/moonis_mnh -o /work/workspace/zhurj/project/2_swgs/sum_2020/moonis/mnh --level class
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/3_swgspro/1_classify_sum/classify_levels_sum.py -i /work/workspace/zhurj/project/2_swgs/sum_2020/input/moonis_mnh -o /work/workspace/zhurj/project/2_swgs/sum_2020/moonis/mnh --level order
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/3_swgspro/1_classify_sum/classify_levels_sum.py -i /work/workspace/zhurj/project/2_swgs/sum_2020/input/moonis_mnh -o /work/workspace/zhurj/project/2_swgs/sum_2020/moonis/mnh --level family
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/3_swgspro/1_classify_sum/classify_levels_sum.py -i /work/workspace/zhurj/project/2_swgs/sum_2020/input/moonis_mnh -o /work/workspace/zhurj/project/2_swgs/sum_2020/moonis/mnh --level genus

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/3_swgspro/1_classify_sum/classify_levels_sum.py -i /work/workspace/zhurj/project/2_swgs/sum_2020/input/mn_bacteria -o /work/workspace/zhurj/project/2_swgs/sum_2020/moonis/mn_bacteria --level phylum
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/3_swgspro/1_classify_sum/classify_levels_sum.py -i /work/workspace/zhurj/project/2_swgs/sum_2020/input/mn_bacteria -o /work/workspace/zhurj/project/2_swgs/sum_2020/moonis/mn_bacteria --level class
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/3_swgspro/1_classify_sum/classify_levels_sum.py -i /work/workspace/zhurj/project/2_swgs/sum_2020/input/mn_bacteria -o /work/workspace/zhurj/project/2_swgs/sum_2020/moonis/mn_bacteria --level order
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/3_swgspro/1_classify_sum/classify_levels_sum.py -i /work/workspace/zhurj/project/2_swgs/sum_2020/input/mn_bacteria -o /work/workspace/zhurj/project/2_swgs/sum_2020/moonis/mn_bacteria --level family
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/3_swgspro/1_classify_sum/classify_levels_sum.py -i /work/workspace/zhurj/project/2_swgs/sum_2020/input/mn_bacteria -o /work/workspace/zhurj/project/2_swgs/sum_2020/moonis/mn_bacteria --level genus

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/3_swgspro/1_classify_sum/classify_levels_sum.py -i /work/workspace/zhurj/project/2_swgs/sum_2020/input/mn_fungi -o /work/workspace/zhurj/project/2_swgs/sum_2020/moonis/mn_fungi --level phylum
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/3_swgspro/1_classify_sum/classify_levels_sum.py -i /work/workspace/zhurj/project/2_swgs/sum_2020/input/mn_fungi -o /work/workspace/zhurj/project/2_swgs/sum_2020/moonis/mn_fungi --level class
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/3_swgspro/1_classify_sum/classify_levels_sum.py -i /work/workspace/zhurj/project/2_swgs/sum_2020/input/mn_fungi -o /work/workspace/zhurj/project/2_swgs/sum_2020/moonis/mn_fungi --level order
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/3_swgspro/1_classify_sum/classify_levels_sum.py -i /work/workspace/zhurj/project/2_swgs/sum_2020/input/mn_fungi -o /work/workspace/zhurj/project/2_swgs/sum_2020/moonis/mn_fungi --level family
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/3_swgspro/1_classify_sum/classify_levels_sum.py -i /work/workspace/zhurj/project/2_swgs/sum_2020/input/mn_fungi -o /work/workspace/zhurj/project/2_swgs/sum_2020/moonis/mn_fungi --level genus

赵博新种进化树图
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/test/species_out.fasta -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/test/species_out.out -max_target_seqs 1 -num_threads 16 -outfmt 7 

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/test/strain2.fas -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/test/strain2.out -max_target_seqs 5 -num_threads 16 -outfmt 7

extract ref sequence 
seqtk subseq /work/database/ezbio/16S/current/Ezbio_16S_seqs.fa /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/ref.list > /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/ref.fasta

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/22_SeqIO_format_convert/format_converter.py -i /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/need_reverse_complement.tab --in_format tab --out_format fasta -o /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/need_reverse_complement.fasta

seqtk seq -r /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/need_reverse_complement.fasta > /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/after_reverse_complete.fas

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/22_SeqIO_format_convert/format_converter.py -i /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/after_reverse_complete.fas --in_format fasta --out_format tab -o /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/after_reverse_complete.tab
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/22_SeqIO_format_convert/format_converter.py -i /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/ref.fasta --in_format fasta --out_format tab -o /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/ref.tab
# merge ref.tab, after_reverse_complete.tab, 本来是正向序列的菌株信息（181） to merge_in.tab
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/22_SeqIO_format_convert/format_converter.py -i /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/merge_in.tab --in_format tab --out_format fasta -o /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/merge_in.fasta

muscle -in /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/merge_in.fasta -out /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/merge_muscle.fas
fasttree -nt -gtr  /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/merge_muscle.fas   >  /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/merge_fasttree.nwk

library(ggsci)
library(vegan)
mycolors <- unique(c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                      pal_simpsons("springfield")(16),pal_gsea("default")(12))) + pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
+                       pal_locuszoom("default")(7),pal_igv("default")(51),
+                       pal_uchicago("default")(9),pal_startrek("uniform")(7),
+                       pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
+                       pal_simpsons("springfield")(16),pal_gsea("default")(12)))

mycolors[11:26]

mycolors[30:64]

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/22_SeqIO_format_convert/format_converter.py -i /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/species_in.tab --in_format tab --out_format fasta -o /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/species_out.fasta
cat /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/unige_genus | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species -T -o /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/new_genus_levels

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -P "AAQK01001555|AB038361|AB081585|AB264064|AB266102|AB298774|ADDX01000083|ADLB01000035|AF338413|AM176535|AM236149|AM278590|AM404803|AM406047|AP013105|AP018533|AP018536|AQHX01000018|AY167945|CP002045|CP003040|CP014229|CP029462|CP039381|DQ777915|EU472250|EU779054|FCEY01000004|FJ363576|FJ880395|GQ451293|HF545616|HF558387|JH601096|JN559675|JN559721|JN680597|JNJP01000178|JTDA01000039|JXXK01000069|KF705042|KI440788|KR185326|KR364771|KY595971|KY992931|L16477|L34621|LC423447|LFQU01000071|LN869527|LN879468|LN908974|LN908986|LT223569|LT558846|LT623892|LT629804|LT635494|LT821227|LT896588|LT907848|LT962669|LT985808|MF322780|NFIP01000003|NKXR01000005|OKQO01000019|PAC000186|PAC000191|PAC000196|PAC000672|PAC001037|PAC001038|PAC001041|PAC001043|PAC001054|PAC001157|PAC001206|PAC001213|PAC001277|PAC001282|PAC001286|PAC001316|PAC001335|PAC001406|PAC001450|PAC001460|PAC001580|PAC001591|PAC001646|PAC001648|PAC001762|PAC001775|PAC001808|PAC002042|PAC002300|PAC002361|PAC002379|PAC002392|PAC002509|PAC002528|PAC002725|PAC002777|PAC002802|PAC002810|PAC002811|PAC002833|PAC002892|PAC002900|PAC003086|PAC003112|QRQB01000006|QTXI01000062|QUCK01000052|QUDC01000017|QUHA01000002|QUHC01000016|QUHQ01000001|QUJU01000002|QUKS01000002|RIBI01000013|GQ491924|KX025133" > /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/mid/ref_levels


/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/test/species_out.fasta -task blastn -db /work/database/ncbi/16S/2020.Sep8/16S_ribosomal_RNA -out /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/nr/species_nr_out.out -max_target_seqs 1 -num_threads 16 -outfmt 7
cat /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/nr/species_nr_out.out | grep "^MNH" > /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/nr/mnh_nr.out
bash /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/p/nr_info.sh
seqtk subseq /work/database/ncbi/16S/2020.Sep8/16S_ribosomal_RNA.fas /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/nr/ref_nr_id > /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/nr//ref_nr.fasta
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/22_SeqIO_format_convert/format_converter.py -i /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/nr/ref_nr.fasta --in_format fasta --out_format tab -o /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/nr/ref_nr.tab
cat /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/nr/ref_name | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species -T -o /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/nr/name_levels
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/22_SeqIO_format_convert/format_converter.py -i /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/nr/merge_forward_seq.tab --in_format tab --out_format fasta -o /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/nr/merge_forward_seq.fasta
muscle -in /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/nr/merge_forward_seq.fasta -out /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/nr/merge_nr_muscle.fas
fasttree -nt -gtr  /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/merge_nr_muscle.fas   >  /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/nr/merge_nr_fasttree.nwk

修改8个菌株的物种鉴定结果
MNH44546  Bacteroidetes Bacteroidales Odoribacteraceae  Butyricimonas Butyricimonas synergistica  94.5%
NR_028812.1 Anaeroglobus  geminatus Veillonellaceae Veillonellales  Firmicutes  #2ca02c #6F99ADFF #CC0000FF

MNH39949  MNH39949  Enterococcaceae Lactobacillales Firmicutes  #2ca02c #1B1919FF #79AF97FF NR_148574.1
MNH39529  MNH39529  Enterococcaceae Lactobacillales Firmicutes  #2ca02c #1B1919FF #79AF97FF NR_148574.1
MNH39884  MNH39884  Enterococcaceae Lactobacillales Firmicutes  #2ca02c #1B1919FF #79AF97FF NR_148574.1
MNH39946  MNH39946  Enterococcaceae Lactobacillales Firmicutes  #2ca02c #1B1919FF #79AF97FF NR_148574.1
MNH39931  MNH39931  Enterococcaceae Lactobacillales Firmicutes  #2ca02c #1B1919FF #79AF97FF NR_148574.1
NR_148574.1 Raoultibacter timonensis  Eggerthellaceae Eggerthellales  Actinobacteria  #1f77b4 #5F559BFF #B24745FF

MNH34135  MNH34135  Clostridiaceae  Clostridiales Firmicutes  #2ca02c #631879FF #ADB6B6FF NR_027565.1
NR_027565.1 Victivallis vadensis  Victivallaceae  Victivallales Lentisphaerae #e377c2 #FFDC91FF #FFCCCCFF

MNH41457  MNH41457  Lachnospiraceae Clostridiales Firmicutes  #2ca02c #631879FF #EFC000FF NR_147402.1
NR_147402.1 Sutturella  timonensis  Sutterellaceae  Burkholderiales Proteobacteria  #8c564b #008B45FF #CC33FFFF

1336个测序菌株进化树

find /work/assembly/fasta/genome/MNH | grep -P "fas$" > /work/workspace/zhurj/project/7_summary/2020/phylo/seq_strain/seq_fna.list
# 截止到2021.02.01 已组装的基因组，无重复，去除了有污染的基因组

/work/workspace/zhurj/software/JolyTree-master/JolyTree.sh -i /work/workspace/zhurj/project/7_summary/2020/phylo/seq_strain/softlink -b /work/workspace/zhurj/project/7_summary/2020/phylo/seq_strain/tree/seqstrains_20210201 -t 16
cat /work/workspace/zhurj/project/7_summary/2020/phylo/seq_strain/unique_strain_genus_name | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species -T -o /work/workspace/zhurj/project/7_summary/2020/phylo/seq_strain/genus_to_levels

2021.02.10
所有于已完成单菌全基因组测序菌株进化树制作
1. 提取所有已测序样品：
find /work/assembly/fasta/genome/MNH | grep fas$ | awk -F "[./]" '{print $(NF-1)"\t"$0}' > /work/workspace/zhurj/project/2_swgs/MNH05168/EnterococcusPhylo_20210210/input/seqStrains20210210
2. 挑选已测序且属于Enterococcus的菌株，去除测序有污染的菌株
cat /work/workspace/zhurj/project/2_swgs/MNH05168/EnterococcusPhylo_20210210/input/enterococcus_strain_20210210 | awk -F "\t" '{print "ln -s "$2" /work/workspace/zhurj/project/2_swgs/MNH05168/EnterococcusPhylo_20210210/softlink/"$1".fna"}' > /work/workspace/zhurj/project/2_swgs/MNH05168/EnterococcusPhylo_20210210/p/phylophlan_in.sh
bash /work/workspace/zhurj/project/2_swgs/MNH05168/EnterococcusPhylo_20210210/p/phylophlan_in.sh

#phylophlan_write_config_file -d a -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/phylophlan/config.cfg --db_aa diamond --map_dna diamond --map_aa diamond --msa mafft --trim trimal --tree1 raxml --verbose 
phylophlan -i /work/workspace/zhurj/project/2_swgs/MNH05168/EnterococcusPhylo_20210210/softlink -d phylophlan --diversity medium --accurate -f /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/phylophlan/config.cfg -o /work/workspace/zhurj/project/2_swgs/MNH05168/EnterococcusPhylo_20210210/tree --nproc 20 --verbose 

library(ggsci)
library(vegan)
#mycolors <- unique(c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                      pal_simpsons("springfield")(16),pal_gsea("default")(12)))+                       pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
+                       pal_locuszoom("default")(7),pal_igv("default")(51),
+                       pal_uchicago("default")(9),pal_startrek("uniform")(7),
+                       pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
+                       pal_simpsons("springfield")(16),pal_gsea("default")(12)))

mycolors <- pal_simpsons("springfield")(16)

################################################################################################################################################################################################################
# MERCK 项目
2021.02.18 物种鉴定
cd /work/workspace/zhurj/project/15_collaborator/merck_20210125/typestrain
cat typereport_list | grep -P "MNH34106|MNH36789|MNH32170|MNH41089|MNH40091|MNH35904|MNH01870|MNH07909|MNH31725|MNH38402|MNH17714|MNH03137|MNH19716|MNH08698|MNH36333|MNH36299|MNH36049|MNH19893|MNH20005|MNH20865|MNH34918" > merge_mnlist
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/21_multifile_defineLine_merge/merge_multifile_defineline.py --batchfile /work/workspace/zhurj/project/15_collaborator/merck_20210125/typestrain/merck_mnlist -o /work/workspace/zhurj/project/15_collaborator/merck_20210125/typestrain --line 2 --file_suffix type.report

fastp: --poly_g_min_len 10 --poly_x_min_len 10 -q 15 -u 40 -n 5 -l 50
# QC 
# pipeline

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

# /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/fastqc -o /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/rawdata/fastqc -t 10 --quiet /work/workspace/zhurj/project/1_metadata/mouse22/clean/MNA00233.1.fq.gz /work/workspace/zhurj/project/1_metadata/mouse22/clean/MNA00233.2.fq.gz
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/fastqc -o /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/rawdata/fastqc -t 10 --quiet /work/workspace/zhurj/project/1_metadata/mouse22/clean/MNA00233.2.fq.gz

/work/workspace/liangyj/bin/conda_env/RNASeq/bin/multiqc -o /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/rawdata/multiqc -n test -f /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/rawdata/fastqc

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/fastq_stat_gc_length.py -o /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/rawdata/gcsum -i /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/input/input -n testGC

perl /work/workspace/zhurj/script/1_tools/3_qcstat/qcstat.pl -i /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/input/input -o /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/rawdata/gcsum --basecutoff 0.8 --cutoffq20 90 --cutoffq30 90


perl /work/workspace/zhurj/script/5_wgspro/1_filter/mnid2rawdata.pl -i /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/input/mnid7 -t 0 -r genome -o /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/input -n rawdata7

bash /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/p/cprawfastq.sh

find /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/fastq | grep gz | sort > /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/input/raw_fastq

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/25_merge_qc/merge_qc.py -i /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/input/raw_fastq -o /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/rawdata -t 16

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/perl /work/workspace/zhurj/script/git/datafilter/datafilter.pl -i /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/input/rawdata7_r1r2_2c.txt -o /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC
cd /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/clean
rename .clean. . *.gz

find /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/clean | grep gz | sort > /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/input/clean_fastq

# clean data qc
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/25_merge_qc/merge_qc_clean.py -i /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/input/clean_fastq -o /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/cleandata -t 16 --rawqc_result /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/rawdata/qcstat/qcmerge.txt

# assembled genome
/work/workspace/zhurj/script/1_tools/10_mnid2assgenome/mnid2assgenome.py  \
-o /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/input \
--batchfile /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/input/mnid7 \
--type 1

cat /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/input/mn2assemble.genome.txt| awk -F "/" '{print "cp "$0" /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/assembly/"$(NF-1)".fna"}' > /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/p/cpassembly.sh
bash /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/p/cpassembly.sh

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
# parallel --link -j 1 'mkdir -p {1}/quast/{2}/{3} ' \
::: /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC :::: /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/input/mnid7

parallel --link -j 1 '/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/quast -l {2} -t 36 -o {1}/quast/{2} {1}/assembly/{2}.fna --silent' \
::: /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC :::: /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/input/mnid7
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

cd /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/p
time srun -o quast_7.out -e quast_7.err -N 1 -c 20 -p slurm256 -w mnclient01 bash quast_7.sh &

show_time = time.strftime("%Y-%m-%d %H:%M:%S",time.localtime())
find /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/quast | grep transposed_report.tsv | awk -F "/" '{print $(NF-1)"\t"$0}' > /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/input/quast_file

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/26_merge_quast/merge_quast.py -i /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/input/quast_file -o /work/workspace/zhurj/project/15_collaborator/merck_20210125/QC/merge/merge_quast

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/21_multifile_defineLine_merge/merge_multifile_defineline.py --dir /work/classify/species/classify/MNH -o /work/workspace/zhurj/project/2_swgs/sum_2020/type_strain_sum --line 2 --file_suffix type.report

fastANI -q /work/assembly/fasta/genome/MNH/MNH051/MNH05167/E1/L1/S1/MNH05167.fas -r /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735495/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/sum_2020/type_strain_sum/Elactis_MNH05167
fastANI -q /work/assembly/fasta/genome/MNH/MNH051/MNH05168/E6/L1/S1/MNH05168.fas -r /work/database/ncbi/genome/all/GCA/GCA_009/GCA_009735/GCA_009735495/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/sum_2020/type_strain_sum/Elactis_MNH05168E6

菌种筛选
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/21_multifile_defineLine_merge/merge_multifile_defineline.py --dir /work/classify/species/classify/MNH -o /work/workspace/zhurj/project/2_swgs/classify_type/20210421 --line 2 --file_suffix type.report --name typereport
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/21_multifile_defineLine_merge/merge_multifile_defineline.py --dir /work/classify/species/classify/MNH -o /work/workspace/zhurj/project/2_swgs/classify_type/20210421 --line 2 --file_suffix all.report --name allreport


df = pd.read_table(infile,sep='\t',index_col=None)
pre_gene_df = df.dropna(subset=['gene'])
gene_df = pre_gene_df['gene']
genes_df = gene_df.str.split('_',expand=True)
unique_gene_df = genes_df[0].to_frame().drop_duplicates()
unique_gene_df.columns = ['gene']

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/37_SCFA_prediction_from_genome/scfa_prediction_from_genome.py -i /work/workspace/zhurj/project/2_swgs/SCFA/SCFA_prediction_20210424/strain_genome_file -o /work/workspace/zhurj/project/2_swgs/SCFA/SCFA_prediction_20210424/strain_SCFA


26个Enterococcus 菌株进化树
cat /work/workspace/zhurj/project/2_swgs/MNH05168/EnterococcusPhylo_strain26_20210429/input/Enterococcus_strain26 | awk -F "\t" '{print "ln -s "$2" /work/workspace/zhurj/project/2_swgs/MNH05168/EnterococcusPhylo_strain26_20210429/softlink/"$1".fna"}' > /work/workspace/zhurj/project/2_swgs/MNH05168/EnterococcusPhylo_strain26_20210429/p/phylophlan_in.sh
bash /work/workspace/zhurj/project/2_swgs/MNH05168/EnterococcusPhylo_strain26_20210429/p/phylophlan_in.sh

#phylophlan_write_config_file -d a -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/phylophlan/config.cfg --db_aa diamond --map_dna diamond --map_aa diamond --msa mafft --trim trimal --tree1 raxml --verbose 
/work/workspace/zhurj/bin/miniconda3/envs/phylophlan/bin/phylophlan -i /work/workspace/zhurj/project/2_swgs/MNH05168/EnterococcusPhylo_strain26_20210429/softlink -d phylophlan --diversity medium --accurate -f /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/phylophlan/config.cfg -o /work/workspace/zhurj/project/2_swgs/MNH05168/EnterococcusPhylo_strain26_20210429/tree --nproc 20 --verbose 

# TEST


/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/40_orthovenn_cluster_extraction/cluster_extraction.py -i /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/test -g /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/group_info -c /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/compare_info -o /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/result.txt
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/40_orthovenn_cluster_extraction/cluster_extraction.py -i /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/orthoveen_cluster -g /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/group_info -c /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/compare_info -o /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/result.txt

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/40_orthovenn_cluster_extraction/cluster_extraction.py -i /work/workspace/zhurj/project/2_swgs/MNH05168/EnterococcusLactis5_20210514/orthovenn2_cluster_Elactis5 -g /work/workspace/zhurj/project/2_swgs/MNH05168/EnterococcusLactis5_20210514/group_info -c /work/workspace/zhurj/project/2_swgs/MNH05168/EnterococcusLactis5_20210514/compare_info -o /work/workspace/zhurj/project/2_swgs/MNH05168/EnterococcusLactis5_20210514/result.txt

find /work/predict/orfs/Prokka/MNH | grep faa > /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/all_20210511_faa

cp /work/predict/orfs/Prokka/MNH/MNH051/MNH05168/E2/L1/S1/MNH05168.ffn /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/ffn
cp /work/predict/orfs/Prokka/MNH/MNH087/MNH08779/E1/L1/S1/MNH08779.ffn /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/ffn
cp /work/predict/orfs/Prokka/MNH/MNH218/MNH21870/E2/L1/S1/MNH21870.ffn /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/ffn
cp /work/predict/orfs/Prokka/MNH/MNH214/MNH21491/E1/L1/S1/MNH21491.ffn /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/ffn

cp /work/predict/orfs/Prokka/MNH/MNH062/MNH06231/E1/L1/S1/MNH06231.ffn /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/ffn
cp /work/predict/orfs/Prokka/MNH/MNH120/MNH12018/E1/L1/S1/MNH12018.ffn /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/ffn
cp /work/predict/orfs/Prokka/MNH/MNH164/MNH16402/E1/L1/S1/MNH16402.ffn /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/ffn
cp /work/predict/orfs/Prokka/MNH/MNH258/MNH25848/E1/L1/S1/MNH25848.ffn /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/ffn
cp /work/predict/orfs/Prokka/MNH/MNH364/MNH36424/E1/L1/S1/MNH36424.ffn /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/ffn
cp /work/predict/orfs/Prokka/MNH/MNH108/MNH10846/E1/L1/S1/MNH10846.ffn /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/ffn

cat /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/ffn/* > /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/merge_ffn.ffn
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/seqtk subseq /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/merge_ffn.ffn /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/cluster2268_sortase_id > /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/cluster2268_sortase.ffn
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/seqtk subseq /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/merge_ffn.ffn /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/cluster2297_Thermophilic_serine_proteinase_id > /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/cluster2297_Thermophilic_serine_proteinase.ffn
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/seqtk subseq /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/merge_ffn.ffn /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/cluster2228_deacetylase_id > /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/cluster2228_deacetylase.ffn
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/seqtk subseq /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/merge_ffn.ffn /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/cluster1823_NIpcP60_id > /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/cluster1823_NIpcP60.ffn
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/seqtk subseq /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/merge_ffn.ffn /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/mnh05168_4candidate_id > /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/mnh05168_4candidate.ffn
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/seqtk subseq /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/merge_faa.faa /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/mnh05168_4candidate_id > /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/mnh05168_4candidate.faa

cp /work/assembly/fasta/genome/MNH/MNH051/MNH05168/E2/L1/S1/MNH05168.fas /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/blast_ref/MNH05168/MNH05168.fasta
/work/program/current/ncbi-blast/bin/makeblastdb -in /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/blast_ref/MNH05168/MNH05168.fasta -input_type fasta -dbtype 'nucl' -parse_seqids -out MNH05168 -hash_index
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/mnh05168_4candidate.ffn -task blastn -db /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/blast_ref/MNH05168/MNH05168 -out /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/blastn/mnh05168_4candidate.out -max_target_seqs 2 -num_threads 16
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/mnh05168_4candidate.ffn -task blastn -db /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/blast_ref/MNH05168/MNH05168 -out /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/blastn/mnh05168_4candidate_format7.out -max_target_seqs 1 -num_threads 16 -outfmt 7

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/samtools faidx /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/blast_ref/MNH05168/MNH05168.fasta 


/usr/bin/bash /work/workspace/zhurj/script/1_tools/41_linux_tools/extraction_ffile2_accordingfile1.sh /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/mnid /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/all_20210511_faa /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/strain8_faa
cat /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/strain8_faa | grep "E1/L1/S1" > /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/select_strain8_faa



# /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/seqtk subseq /work/database/ncbi/16S/current/../2021.Apr28/16S_ribosomal_RNA.fas /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/new/ref_dedup_id
cat /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/select_strain8_faa | awk '{print "cp "$0" /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/faa"}' > /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/p/relocate_faa.sh 
bash /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/p/relocate_faa.sh
cat /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/faa/*.faa > /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/strain8_merge.faa
cat /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/strain8_merge.faa | grep "P60 family lipoprotein" | awk -F '[> ]' '{print $2}' > /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/NlpC_P60_id
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/seqtk subseq /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/strain8_merge.faa /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/NlpC_P60_id > /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/NlpC_P60_strain8.faa
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/seqtk subseq /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/strain8_merge.faa /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/NlpC_P60_id_renew > /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/NlpC_P60_strain8gene18.faa
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/39_fasta_id_rename/id_rename_fasta.py -i /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/NlpC_P60_strain8gene18.faa -o /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/NlpC_P60_strain8gene18_rename_id.fa --delimiter space

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/22_SeqIO_format_convert/format_converter.py -i /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/NlpC_P60_strain8gene18_rename_id.fa --in_format fasta -o /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/NlpC_P60_strain8gene18_rename_id.tab --out_format tab

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/22_SeqIO_format_convert/format_converter.py -i /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/NlpC_P60_phylo_order_id.tab --in_format tab -o /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/NlpC_P60_phylo_order_id.fa --out_format fasta

NlpC/P60是否为分泌蛋白 经典分泌蛋白预测
ref: https://zhuanlan.zhihu.com/p/52744204
# 信号肽判断：1.  CS probability >=0.5 && max(SP(Sec/SPI), TAT(Tat/SPI), LIPO(Sec/SPII)) >=0.5 || 2. max(SP(Sec/SPI), TAT(Tat/SPI), LIPO(Sec/SPII)) >=0.6
综合SignalP和TMHMM的分析结果（存在信号肽但无跨膜结构域）
MNH10846_02418
MNH12018_02252
MNH16402_00081
MNH25848_02487
MNH36424_01138
MNH05168_00452
MNH06231_02904

signalp -org gram+ -fasta /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/NlpC_P60_phylo_order_id.fa -format short -prefix /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/secretion_prediction/signalp/NlpC_P60
/work/workspace/zhurj/software/tmhmm-2.0c/bin/tmhmm /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/NlpC_P60_phylo_order_id.fa > /work/workspace/zhurj/project/2_swgs/MNH05168/SagA_strain8_20210511/secretion_prediction/tmhmm/NlpC_P60_tmhmm.txt 

checkm lineage_wf /work/assembly/current/MNH/MNH276/MNH27614 /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH27614 -x fna -t 20 --tab_table -f /work/workspace/zhurj/project/2_swgs/taxcheck/checkM/MNH27614/MNH27614_checkm.txt
checkm lineage_wf /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/fna /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/checkm -x fa -t 20 --tab_table -f /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/checkm/checm_result.txt

cat /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/input/input | awk -F "[./]" '{print $(NF-1)}' > /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/input/sample
source activate python3.6
parallel -j 1 'mkdir -p /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/quast/{1} ' :::: /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/input/sample
parallel -j 1 'quast -t 16 -o {1}/quast/{2} {1}/fna/{2}.fa' ::: /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/ :::: /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/input/sample
source deactivate python3.6
cd /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/p
time srun -o quast.out -e quast.err -N 1 -c 20 -p slurm256 -w mnclient02 bash quast.sh &
find /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/quast | grep transposed_report.tsv | awk -F "[./]" '{print $(NF-2)"\t"$0}' > /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/input/quast_file
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/26_merge_quast/merge_quast.py -i /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/input/quast_file -o /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/merge/merge_quast

gtdbtk.sh
source activate python3.6
gtdbtk classify_wf --genome_dir /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/fna --out_dir /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/gtdbtk --cpus 2 -x fa
source deactivate python3.6

cd /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/p
time srun -o gtdbtk.out -e gtdbtk.err -N 1 -c 20 -p slurm256 -w mnclient01 bash gtdbtk.sh &




export FGMP="/work/workspace/zhurj/software/FGMP-master"
export FGMPTMP="/work/workspace/zhurj/software/FGMP-master/tmp"
export PERL5LIB="$PERL5LIB:$FGMP/lib"

perl /work/workspace/zhurj/software/FGMP-master/src/fgmp.pl -g /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/fna/HN1.fa -o /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/fgmp/HN1.out

华农项目判断HN2 是否有细菌污染
checkm ssu_finder /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/fna/HN2.fa  /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/fna /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/ssu/HN2 -t 20 -x fa 

blastn -query /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/ssu/HN2/ssu_filter.fna  -db /work/database/ncbi/16S/current/16S_ribosomal_RNA -out /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/blastn/HN2/HN2.out -max_target_seqs 5 -num_threads 16 -outfmt 7
cat /work/database/ncbi/16S/current/16S_ribosomal_RNA.tab | grep -E "NR_028751.1|NR_170517.1|NR_170515.1|NR_170522.1|NR_170519.1"
fastANI -q /work/assembly/current/MNH/MNH260/MNH26050/genomic.fna -r /work/database/ncbi/genome/all/GCA/GCA_000/GCA_000020/GCA_000020225/genomic.fna.gz -t 16 -o /work/workspace/zhurj/project/2_swgs/taxcheck/fastANI/MNH26050_AKK

/work/workspace/zhurj/bin/miniconda3/bin/mash sketch -o /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/taxoMash/HN3 -l /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/input/mash_in

/work/workspace/zhurj/bin/miniconda3/bin/mash dist -t /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/taxoMash/HN3.msh /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/taxoMash/HN3.msh > /work/workspace/zhurj/project/2_swgs/huanong_zhuyongshan20210513/taxoMash/HN3_mash.dist

zhengjiao
8.3M    ./bin
178G    ./program
9.8T    ./project
939M    ./script
1.7T    ./download
3.2G    ./software
198G    ./develop
38M     ./deliver
12T    

seqtk subseq /work/database/ezbio/16S/current/Ezbio_16S_seqs.fa /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/ref.list > /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/ref.fasta
samtools faidx 


/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/43_extract_subseq_renameid_fa/subseq_renameid_fa.py -i /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/input/id_region -r /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/blast_ref/MNH05168/MNH05168.fasta -o /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/output/extract_4seq.fa --tmp /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/output

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/43_extract_subseq_renameid_fa/subseq_renameid_fa.py -i /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/input/test_region -r /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/blast_ref/MNH05168/MNH05168.fasta -o /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/output/test_4seq.fa --tmp /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/output

mlst
parallel -j 1 'cd /work/workspace/zhurj/project/2_swgs/MLST/db/pubmlst_sepidermidis_seqdef/loci/{} && wget https://rest.pubmlst.org/db/pubmlst_sepidermidis_seqdef/loci/{}/alleles_fasta ' ::: gtr mutS pyrR tpiA yqiL

parallel -j 1 'mkdir -p {1}/{2} && \
cd {1}/{2} && \
wget https://rest.pubmlst.org/db/pubmlst_shominis_seqdef/loci/{2}/alleles_fasta ' ::: /work/workspace/zhurj/project/2_swgs/MLST/db/pubmlst_shominis_seqdef/loci ::: arcC glpK gtr pta tpiA tuf

cat `find /work/workspace/zhurj/project/2_swgs/MLST/db/pubmlst_sepidermidis_seqdef | grep _fasta` > /work/workspace/zhurj/project/2_swgs/MLST/db/pubmlst_sepidermidis_seqdef/ref/merge.fasta
/work/program/current/ncbi-blast/bin/makeblastdb -in /work/workspace/zhurj/project/2_swgs/MLST/db/pubmlst_sepidermidis_seqdef/ref/merge.fasta -input_type fasta -dbtype 'nucl' -parse_seqids -out /work/workspace/zhurj/project/2_swgs/MLST/db/pubmlst_sepidermidis_seqdef/ref/merge -hash_index
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/MLST/data/input/test.fasta -task blastn -db /work/workspace/zhurj/project/2_swgs/MLST/db/pubmlst_sepidermidis_seqdef/ref/merge -out /work/workspace/zhurj/project/2_swgs/MLST/data/output/test.out -max_target_seqs 2 -num_threads 16 -outfmt 7

cat `find /work/workspace/zhurj/project/2_swgs/MLST/db/pubmlst_shominis_seqdef | grep _fasta` > /work/workspace/zhurj/project/2_swgs/MLST/db/pubmlst_shominis_seqdef/ref/merge.fasta
/work/program/current/ncbi-blast/bin/makeblastdb -in /work/workspace/zhurj/project/2_swgs/MLST/db/pubmlst_shominis_seqdef/ref/merge.fasta -input_type fasta -dbtype 'nucl' -parse_seqids -out /work/workspace/zhurj/project/2_swgs/MLST/db/pubmlst_shominis_seqdef/ref/merge -hash_index

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/mnh05168_4candidate.ffn -task blastn -db /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/blast_ref/MNH05168/MNH05168 -out /work/workspace/zhurj/project/2_swgs/MNH05168/orthovenn/blastn/mnh05168_4candidate.out -max_target_seqs 2 -num_threads 16


/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/22_SeqIO_format_convert/format_converter.py -i /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/need_reverse_complement.tab --in_format tab --out_format fasta -o /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/need_reverse_complement.fasta

source activate python3.6
parallel -j 1 '/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/22_SeqIO_format_convert/format_converter.py -i {1}/{2}.tab --in_format tab --out_format fasta -o {1}/{2}.fasta && \
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/muscle -in {1}/{2}.fasta -out {1}/{2}_align.fasta -maxiters 16  && \
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/fasttree -nt -gtr  {1}/{2}_align.fasta  >  {1}/{2}_fasttree.nwk ' ::: /work/workspace/zhurj/project/2_swgs/MLST/phylo/input ::: Sepidermidis_sam214 Shominis_sam41

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/fasttree -nt -gtr  /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/merge_muscle.fas   >  /work/workspace/zhurj/project/14_coworker/guozheng/phylum_newspecies_20210129/merge_fasttree.nwk


/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/MLST/data/input/MN228882.fasta -task blastn -db /work/workspace/zhurj/project/2_swgs/MLST/db/pubmlst_sepidermidis_seqdef/ref/merge -out /work/workspace/zhurj/project/2_swgs/MLST/data/output/MN228882_rmt7.out -max_target_seqs 5 -num_threads 16 -outfmt 7
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/MLST/data/input/MN228882.fasta -task blastn -db /work/workspace/zhurj/project/2_swgs/MLST/db/pubmlst_sepidermidis_seqdef/ref/merge -out /work/workspace/zhurj/project/2_swgs/MLST/data/output/MN228882.out -max_target_seqs 2 -num_threads 16

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/44_bacterial_MLST/mlst_typing.py -i /work/workspace/zhurj/project/2_swgs/MLST/data/input/test_in -o /work/workspace/zhurj/project/2_swgs/MLST/data/output/test_out --profile /work/workspace/zhurj/project/2_swgs/MLST/db/pubmlst_sepidermidis_seqdef/profiles_csv -r /work/workspace/zhurj/project/2_swgs/MLST/db/pubmlst_sepidermidis_seqdef/ref/merge --temp_folder /work/workspace/zhurj/project/2_swgs/MLST/data/tmp --blastn /work/program/current/ncbi-blast/bin/blastn
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/44_bacterial_MLST/mlst_typing.py -i /work/workspace/zhurj/project/2_swgs/MLST/data/input/Sepidermidis_sam214_in -o /work/workspace/zhurj/project/2_swgs/MLST/data/output/Sepidermidis_sam214_ST --profile /work/workspace/zhurj/project/2_swgs/MLST/db/pubmlst_sepidermidis_seqdef/profiles_csv -r /work/workspace/zhurj/project/2_swgs/MLST/db/pubmlst_sepidermidis_seqdef/ref/merge --temp_folder /work/workspace/zhurj/project/2_swgs/MLST/data/tmp --blastn /work/program/current/ncbi-blast/bin/blastn

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/44_bacterial_MLST/mlst_typing.py -i /work/workspace/zhurj/project/2_swgs/MLST/data/input/Shominis_sam41_in -o /work/workspace/zhurj/project/2_swgs/MLST/data/output/Shominis_sam41_ST --profile /work/workspace/zhurj/project/2_swgs/MLST/db/pubmlst_shominis_seqdef/profiles_csv -r /work/workspace/zhurj/project/2_swgs/MLST/db/pubmlst_shominis_seqdef/ref/merge --temp_folder /work/workspace/zhurj/project/2_swgs/MLST/data/tmp --blastn /work/program/current/ncbi-blast/bin/blastn


################################################################################################################################################################################################################
# MERCK 项目9个strain总结
2021.06.15 物种鉴定

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/21_multifile_defineLine_merge/merge_multifile_defineline.py --dir /work/classify/species/classify/MNH -o /work/workspace/zhurj/project/15_collaborator/merck_20210615/input --line 2 --file_suffix type.report

cat /work/workspace/zhurj/project/15_collaborator/merck_20210615/input/typereport_list | grep -P "MNH17714|MNH19893|MNH34106|MNH35904|MNH36333|MNH38402|MNH41089|MNH19716|MNH20865" > /work/workspace/zhurj/project/15_collaborator/merck_20210615/input/merge_mnlist
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/21_multifile_defineLine_merge/merge_multifile_defineline.py --batchfile /work/workspace/zhurj/project/15_collaborator/merck_20210615/input/merge_mnlist -o /work/workspace/zhurj/project/15_collaborator/merck_20210615/input --line 2 --file_suffix type.report

fastp: --poly_g_min_len 10 --poly_x_min_len 10 -q 15 -u 40 -n 5 -l 50
# QC 
# pipeline

mkdir -p /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC/rawdata
mkdir -p /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC/fastq
mkdir -p /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC/cleandata
mkdir -p /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC/assembly
mkdir -p /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC/checkm

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
perl /work/workspace/zhurj/script/5_wgspro/1_filter/mnid2rawdata.pl -i /work/workspace/zhurj/project/15_collaborator/merck_20210615/input/mnid9 -t 0 -r genome -o /work/workspace/zhurj/project/15_collaborator/merck_20210615/input -n rawdata9

# raw data in /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC/fastq
cat /work/workspace/zhurj/project/15_collaborator/merck_20210615/input/rawdata9_r1r2_1c.txt | awk -F "/" '{print "cp "$0" /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC/fastq"}' > /work/workspace/zhurj/project/15_collaborator/merck_20210615/p/cprawfastq.sh
echo -e "\n" >> /work/workspace/zhurj/project/15_collaborator/merck_20210615/p/cprawfastq.sh
bash /work/workspace/zhurj/project/15_collaborator/merck_20210615/p/cprawfastq.sh

find /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC/fastq | grep gz | sort > /work/workspace/zhurj/project/15_collaborator/merck_20210615/input/raw_fastq
# qcstat.pl doesn't work, need to find the causes
#############################################
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/25_merge_qc/merge_qc.py -i /work/workspace/zhurj/project/15_collaborator/merck_20210615/input/raw_fastq -o /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC/rawdata -t 16
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/perl /work/workspace/zhurj/script/1_tools/3_qcstat/qcstat.pl -i /work/workspace/zhurj/project/15_collaborator/merck_20210615/input/raw_fastq -o /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC/rawdata/qcstat --basecutoff 0.8 --cutoffq20 90 --cutoffq30 90

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/perl /work/workspace/zhurj/script/git/datafilter/datafilter.pl -i /work/workspace/zhurj/project/15_collaborator/merck_20210615/input/rawdata9_r1r2_2c.txt -o /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC
cd /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC/clean
rename .clean. . *.gz

find /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC/clean | grep gz | sort > /work/workspace/zhurj/project/15_collaborator/merck_20210615/input/clean_fastq

# clean data qc
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/25_merge_qc/merge_qc_clean.py -i /work/workspace/zhurj/project/15_collaborator/merck_20210615/input/clean_fastq -o /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC/cleandata -t 16 --rawqc_result /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC/rawdata/qcstat/qcmerge.txt



# assembled genome
# --print 可输出结果，-o为生成文件 --- 查原因
/work/workspace/zhurj/script/1_tools/10_mnid2assgenome/mnid2assgenome.py  \
-o /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC/input \
--batchfile /work/workspace/zhurj/project/15_collaborator/merck_20210615/input/mnid9 \
--type 1

/work/workspace/zhurj/script/1_tools/10_mnid2assgenome/mnid2assgenome.py  \
--print \
--batchfile /work/workspace/zhurj/project/15_collaborator/merck_20210615/input/mnid9 \
--type 1

cat /work/workspace/zhurj/project/15_collaborator/merck_20210615/input/mn2assemble.genome.txt| awk -F "/" '{print "cp "$0" /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC/assembly/"$(NF-1)".fna"}' > /work/workspace/zhurj/project/15_collaborator/merck_20210615/p/cpassembly.sh
bash /work/workspace/zhurj/project/15_collaborator/merck_20210615/p/cpassembly.sh

source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
# parallel -j 1 'mkdir -p {1}/quast/{2}/{3} ' \
::: /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC :::: /work/workspace/zhurj/project/15_collaborator/merck_20210615/input/mnid9

parallel -j 1 '/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/quast -l {2} -t 36 -o {1}/quast/{2} {1}/assembly/{2}.fna --silent' \
::: /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC :::: /work/workspace/zhurj/project/15_collaborator/merck_20210615/input/mnid9
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

cd /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC/p
time srun -o quast_7.out -e quast_7.err -N 1 -c 20 -p slurm256 -w mnclient01 bash quast_7.sh &

show_time = time.strftime("%Y-%m-%d %H:%M:%S",time.localtime())
find /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC/quast | grep transposed_report.tsv | awk -F "/" '{print $(NF-1)"\t"$0}' > /work/workspace/zhurj/project/15_collaborator/merck_20210615/input/quast_file

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/26_merge_quast/merge_quast.py -i /work/workspace/zhurj/project/15_collaborator/merck_20210615/input/quast_file -o /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC/merge/merge_quast

# checkm 
/work/workspace/zhurj/bin/miniconda3/bin/checkm lineage_wf /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC/assembly /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC/assembly -x fna -t 36 --tab_table -f /work/workspace/zhurj/project/15_collaborator/merck_20210615/QC/checkm/checkm.txt

/work/program/current/ncbi-blast/bin/makeblastdb -in /work/workspace/zhurj/project/2_swgs/MNH04863/assemble/MNH04863.fasta -input_type fasta -dbtype 'nucl' -parse_seqids -out /work/workspace/zhurj/project/2_swgs/MNH04863/assemble/MNH04863 -hash_index

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/MNH04863/gene/ClpB/ncbi_dataset/data/gene.fna  -task blastn -db /work/workspace/zhurj/project/2_swgs/MNH04863/assemble/MNH04863 -out /work/workspace/zhurj/project/2_swgs/MNH04863/gene/ClpB/result/clpB_MNH04863.out -max_target_seqs 20 -num_threads 16

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/2_swgs/MNH04863/gene/ClpB/ncbi_dataset/data/gene.fna  -task blastn -db /work/workspace/zhurj/project/2_swgs/MNH04863/assemble/MNH04863 -out /work/workspace/zhurj/project/2_swgs/MNH04863/gene/ClpB/result/clpB_MNH04863_outfmt7.out -max_target_seqs 20 -num_threads 16 -outfmt 7

/work/program/current/ncbi-blast/bin/makeblastdb -in /work/workspace/zhurj/project/2_swgs/MNH04863/protein/MNH04863E1L1S1.fasta -input_type fasta -dbtype 'prot' -parse_seqids -out /work/workspace/zhurj/project/2_swgs/MNH04863/protein/MNH04863E1L1S1 -hash_index

/work/program/current/ncbi-blast/bin/blastp -query /work/workspace/zhurj/project/2_swgs/MNH04863/gene/ClpB/ncbi_dataset/data/gene.fna  -task blastp -db /work/workspace/zhurj/project/2_swgs/MNH04863/protein/MNH04863E1L1S1 -out /work/workspace/zhurj/project/2_swgs/MNH04863/gene/ClpB/result/clpB_MNH04863_protein.out -max_target_seqs 5 -num_threads 16

/work/program/current/ncbi-blast/bin/blastp -query /work/workspace/zhurj/project/2_swgs/MNH04863/gene/ClpB/uniprot_ClpB_Hafnia_alvei.fasta  -task blastp -db /work/workspace/zhurj/project/2_swgs/MNH04863/protein/MNH04863E1L1S1 -out /work/workspace/zhurj/project/2_swgs/MNH04863/gene/ClpB/result/clpB_MNH04863_protein_outfmt7.out -max_target_seqs 5 -num_threads 16 -outfmt 7

'''
	
