2019.Nov07-13:13

去除环境变量里的重复项：
export PATH=$(echo $PATH | sed 's/:/\n/g' | sort | uniq | tr -s '\n' ':' | sed 's/:$//g')

vi /etc/hosts

vi /etc/group
database:x:1000:zhouyj,zhurj,liangyj,xiongz,liangzj,data,fastq,fasta
moon:x:1001:zhouyj,zhurj,liangyj,xiongz,liangzj,instal,data,fasta,fastq,moonis
bioinfor:x:1002:zhangdy,zhangcc
project:x:1003:zhouyj,zhurj,xiongz,liangzj,fasta,fastq,moonis
program:x:1004:zhouyj,zhurj,liangyj,xiongz,liangzj,instal


parallel -j 1 'usermod -a -G bioinfor {1}' ::: zhurj liangzj liangyj

test
useradd test
passwd test
# moon123456
parallel 'mkdir -p /work/workspace/test/{1} ' ::: bin develop download program project script software
usermod -a -G bioinfor test
chown -R test:bioinfor /work/workspace/test

parallel -j 1 'usermod -a -G {1} test' ::: database moon project program


root 权限删除 test
userdel -f -r test # -f 强制删除账户不管对方是否在登录，-r 删除账户下的目录
rm /work/workspace/test -rf

创建新用户
2021.05.13
# 查看已经安装的用户 cat /etc/passwd
IP: 192.168.2.20
useradd zhangdy
passwd zhangdy
# moon123456
parallel 'mkdir -p /work/workspace/zhangdy/{1} ' ::: bin develop download program project script software
usermod -a -G bioinfor zhangdy
parallel -j 1 'usermod -a -G {1} zhangdy' ::: database moon project program
chown -R zhangdy:bioinfor  /work/workspace/zhangdy
chown -R zhangdy:bioinfor  /home/zhangdy

2021.07.05
useradd zhangcc
passwd zhangcc
# moon123456
parallel 'mkdir -p /work/workspace/zhangcc/{1} ' ::: bin develop download program project script software
usermod -a -G bioinfor zhangcc
usermod -a -G database zhangcc
parallel -j 1 'usermod -a -G {1} zhangcc' ::: moon project program
chown -R zhangcc:bioinfor  /work/workspace/zhangcc
chown -R zhangcc:bioinfor  /home/zhangcc

2021.09.13
useradd lanzhou
passwd lanzhou
# moon123456
parallel 'mkdir -p /work/workspace/lanzhou/{1} ' ::: bin develop download program project script software
usermod -a -G bioinfor lanzhou
usermod -a -G database lanzhou
parallel -j 1 'usermod -a -G {1} lanzhou' ::: moon project program
chown -R lanzhou:bioinfor  /work/workspace/lanzhou
chown -R lanzhou:bioinfor  /home/lanzhou

服务器 IP: 192.168.2.20
账户： lanzhou
密码： moon123456
工作目录： /work/workspace/lanzhou


content in /home/zhurj/.condarc
```
channels:
  - defaults
show_channel_urls: true
channel_alias: https://mirrors.tuna.tsinghua.edu.cn/anaconda
default_channels:
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/pro
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/msys2
custom_channels:
  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  msys2: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  menpo: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  pytorch: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  simpleitk: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud


channels:
  - defaults
show_channel_urls: true
default_channels:
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
custom_channels:
  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  msys2: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  menpo: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  pytorch: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  simpleitk: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
```

清空缓存
# 查看内存使用情况
free -h
# root 操作
sync 
echo 1 > /proc/sys/vm/drop_caches
···
0：不释放（系统默认值）
 1：释放页缓存
 2：释放dentries和inodes
 3：释放所有缓存
···


source activate humann
humann2的reference： https://github.com/biobakery/humann/tree/3.0.0-alpha/humann2/data/misc
有用的论坛： https://groups.google.com/g/humann-users/c/V4SiiJfWw_I/m/10mvxRhrAgAJ
#  humann3 ： tend to stick with UniRef90 for anything human-associated and favor UniRef50 otherwise.
.libPaths()
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
R
install.packages("lattice")
install.packages("/work/workspace/zhurj/software/R/vegan_2.5-6.tar.gz",lib="/work/workspace/zhurj/bin/miniconda3/lib/R/library")

conda create -n R3.6 r=3.6.0
# R version 3.6.1 (2019-07-05)
# https://anaconda.org/conda-forge/r-vegan
conda install -c conda-forge r-vegan
# v2.5_6



**tcsh install
conda create -n tcsh
conda install -n tcsh -c conda-forge tcsh

login ruijuan@192.168.2.74
scp /mnt/d/home/ruijuan/workflows/software/kSNP3.1_Linux_package.zip zhurj@192.168.2.20:/work/workspace/zhurj/software
scp /mnt/d/home/ruijuan/workflows/software/FigTree_v1.4.4.tgz zhurj@192.168.2.20:/work/workspace/zhurj/software
scp /mnt/d/home/ruijuan/workflows/software/pyani-master.zip zhurj@192.168.2.20:/work/workspace/zhurj/software

login zhurj@192.168.2.20
cd /work/workspace/zhurj/software
unzip kSNP3.1_Linux_package.zip
unzip pyani-master.zip
tar zxvf FigTree_v1.4.4.tgz

# pyani installation
pip3 install pyani

# https://github.com/picrust/picrust2/wiki/PICRUSt2-Tutorial-(v2.1.4-beta)
# install picrust2
conda create -n picrust2 -c bioconda -c conda-forge picrust2=2.1.4_b

# https://bitbucket.org/biobakery/biobakery/wiki/lefse
# install lefse
conda create -n lefse 
conda install -c biobakery lefse

under qiime2-2019.10
qiime mmvec --help

# https://anaconda.org/bioconda/bwa
# under env base
conda install -c bioconda bwa

# CheckM
# https://github.com/Ecogenomics/CheckM/wiki (very good website!!!)
# installed under no conda environment
conda install -c bioconda checkm-genome 
# /work/workspace/zhurj/bin/miniconda3/bin/checkm
# example : https://linsalrob.github.io/ComputationalGenomicsManual/CheckM/
以下内容为python包，使用 pip3 freeze 查看是否安装，及安装的版本号
conda install -c conda-forge numpy
conda install -c conda-forge scipy
conda install -c conda-forge matplotlib
conda install -c bioconda pysam
conda install -c bioconda dendropy

# install conFindr
# https://olc-bioinformatics.github.io/ConFindr/install/#downloading-confindr-databases
conda install -c bioconda confindr

conda create -n python3.6
source activate python3.6
conda install python=3.6.5
conda install -c bioconda quast
conda install -c bioconda busco
conda install -c bioconda fastani
conda install -c bioconda augustus
conda install -c bioconda metaxa
conda install -c bioconda blast
conda install -c bioconda centrifuge
# https://ccb.jhu.edu/software/centrifuge/manual.shtml

conda create -n ipython
# https://anaconda.org/anaconda/ipython
conda install -c anaconda ipython
# https://www.php-master.com/post/355980.html
conda install -p /work/workspace/zhurj/bin/miniconda3/envs/ipython -c bioconda bioawk -y
conda install -c anaconda pandas -y

With a working anaconda installation, install the bioconda and conda-forge channels:
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install pyani

under env base
install infernal
download infernal-1.1.3-linux-intel-gcc.tar.gz to /work/workspace/zhurj/software
tar zxvf infernal-1.1.3-linux-intel-gcc.tar.gz
cd infernal-1.1.3-linux-intel-gcc/
configure:                 ./configure
build:                     make
automated tests:           make check
automated install:         make install


2020.08.05
ref: https://github.com/tseemann/prokka#installation
install package in linux(CentOS) without root user
# under no environment
# [zhurj@mnhead sh] /work/workspace/zhurj/software/prokka-master/bin/prokka
# python 3.6
# blastp 2.10.1+
#sudo yum install git perl-Time-Piece perl-XML-Simple perl-Digest-MD5 perl-App-cpanminus git java perl-CPAN perl-Module-Build
#sudo cpanm Bio::Perl
#git clone https://github.com/tseemann/prokka.git $HOME/prokka
#$HOME/prokka/bin/prokka --setupdb

mkdir /work/workspace/zhurj/instal/rpm
#yum list perl-Time-Piece
yumdownloader --destdir /work/workspace/zhurj/instal/rpm --resolve perl-Time-Piece
yumdownloader --destdir /work/workspace/zhurj/instal/rpm --resolve perl-XML-Simple
yumdownloader --destdir /work/workspace/zhurj/instal/rpm --resolve perl-Digest-MD5
yumdownloader --destdir /work/workspace/zhurj/instal/rpm --resolve perl-App-cpanminus
yumdownloader --destdir /work/workspace/zhurj/instal/rpm --resolve perl-CPAN
yumdownloader --destdir /work/workspace/zhurj/instal/rpm --resolve perl-Module-Build
mkdir /work/workspace/zhurj/instal/centos
cd /work/workspace/zhurj/instal/centos && rpm2cpio /work/workspace/zhurj/instal/rpm/*.rpm | cpio -id
# rpm2cpio outputs the .rpm file as a .cpio archive on stdout.
# cpio reads it from from stdin
# -i means extract (to the current directory)
# -d means create missing directory

export PATH="/work/workspace/zhurj/instal/centos/usr/sbin:/work/workspace/zhurj/instal/centos/usr/bin:/work/workspace/zhurj/instal/centos/bin:$PATH"

L='/lib:/lib64:/usr/lib:/usr/lib64'
export LD_LIBRARY_PATH="$L:/work/workspace/zhurj/instal/centos/usr/lib:/work/workspace/zhurj/instal/centos/usr/lib64"

# cpanm 的 mirrors （国内）
#http://mirrors.aliyun.com/CPAN/
# http://mirrors.163.com/cpan/ ftp://mirrors.sohu.com/CPAN/
export PERL_CPANM_OPT="--prompt -l /work/workspace/zhurj/lib/perl --mirror http://mirrors.aliyun.com/CPAN "
cpanm Bio::Perl -l /work/workspace/zhurj/lib/perl --mirror http://mirrors.aliyun.com/CPAN

export PERL5LIB=/work/workspace/zhurj/lib/perl/lib/perl5
cpanm XML::Simple -l /work/workspace/zhurj/lib/perl --mirror http://mirrors.aliyun.com/CPAN
cpanm Bio::Root::Version -l /work/workspace/zhurj/lib/perl --mirror http://mirrors.aliyun.com/CPAN

# export PERL5LIB=/home/foobar/code
# find /work/workspace/liangzj -type f -size +1G -exec du -lh {} \;

# blast 2.10
# ref: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

/work/workspace/zhurj/software/prokka-master/bin/prokka --setupdb

export PERL5LIB=/work/workspace/zhurj/lib/perl/lib/perl5
cpanm Bio::SearchIO::hmmer3 -l /work/workspace/zhurj/lib/perl --mirror http://mirrors.aliyun.com/CPAN

export PATH=$PATH:`pwd`/bin

# 2020.08.07
dRep: cluster genomes
ref: https://drep.readthedocs.io/en/latest/module_descriptions.html
# 1. To install dRep, simply run
pip install drep
# 2. To check which dependencies are installed on your system and accessible by dRep, run
dRep bonus testDir --check_dependencies
Checking dependencies
mash.................................... all good        (location = /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/mash)
nucmer.................................. !!! ERROR !!!   (location = None)
checkm.................................. all good        (location = /work/workspace/zhurj/bin/miniconda3/bin/checkm)
ANIcalculator........................... !!! ERROR !!!   (location = None)
prodigal................................ all good        (location = /work/workspace/zhurj/bin/prodigal)
centrifuge.............................. all good        (location = /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/centrifuge)
nsimscan................................ !!! ERROR !!!   (location = None)
fastANI................................. !!! ERROR !!!   (location = /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/fastANI)
'''
# https://www.ncbi.nlm.nih.gov/genome/annotation_prok/process/
# mash
RUN cd /root && wget https://github.com/marbl/Mash/releases/download/v1.1.1/mash-Linux64-v1.1.1.tar.gz  && \ 
	tar -xvzf mash-Linux64-v1.1.1.tar.gz && \
	mv mash-Linux64-v1.1.1/mash /usr/local/bin/
# Mummer
RUN cd /root && wget https://svwh.dl.sourceforge.net/project/mummer/mummer/3.23/MUMmer3.23.tar.gz  && \
	tar -xvzf MUMmer3.23.tar.gz && cd MUMmer3.23/ && \
	make install 
ENV PATH="$PATH:/root/MUMmer3.23/"
# gani 
RUN cd /root && wget https://ani.jgi-psf.org/download_files/ANIcalculator_v1.tgz && \
	tar -xvzf ANIcalculator_v1.tgz && \
	mv ANIcalculator_v1/ANIcalculator /usr/local/bin/ && \
	mv ANIcalculator_v1/nsimscan /usr/local/bin/ &&  \
	chmod 755 /usr/local/bin/ANIcalculator && \
	chmod 755 /usr/local/bin/nsimscan
# prodigal
RUN cd /root && wget https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux && \ 
	mv prodigal.linux /usr/local/bin/prodigal && \
	chmod 755 /usr/local/bin/prodigal
# centrifuge 
RUN cd /root && git clone https://github.com/infphilo/centrifuge && \
	cd centrifuge && make && \
	make install prefix=/usr/local && \
	cd .. && rm -r centrifuge/
# drep
RUN cd /root && git clone https://github.com/MrOlm/drep.git && \
	cd drep && pip3 install .
'''
cd /work/workspace/zhurj/software
tar -xvzf MUMmer3.23.tar.gz && cd MUMmer3.23/ && \
	make install 
#PATH="$PATH:/root/MUMmer3.23/"
export PATH=$PATH:/work/workspace/zhurj/software/MUMmer3.23

cd /work/workspace/zhurj/software
tar -xvzf ANIcalculator_v1.tgz
mv ANIcalculator /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin
mv nsimscan /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin
chmod 755 /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/ANIcalculator && \
chmod 755 /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/nsimscan

dRep compare -p 20 /work/workspace/zhurj/project/2_swgs/dRep/ChristensenellaStrain2_20200807/output -g /work/workspace/zhurj/project/2_swgs/dRep/ChristensenellaStrain2_20200807/indir/*.fna

# Roary pan genome analysis
ref: https://github.com/sanger-pathogens/Roary
install REF: https://www.jianshu.com/p/11b9d76babba

cd /work/workspace/zhurj/software
unzip Roary-master.zip

设置Roary的环境变量
export PERL5LIB=/work/workspace/zhurj/lib/perl/lib/perl5
export PERL5LIB=$PERL5LIB:/work/workspace/zhurj/software/Roary-master/bin/lib
cpanm Array::Utils Bio::Perl Exception::Class File::Basename File::Copy File::Find::Rule File::Grep File::Path File::Slurper File::Spec File::Temp File::Which FindBin Getopt::Long Graph Graph::Writer::Dot List::Util Log::Log4perl Moose Moose::Role Text::CSV PerlIO::utf8_strict Devel::OverloadInfo Digest::MD5::File  --mirror http://mirrors.aliyun.com/CPAN -l /work/workspace/zhurj/lib/perl/lib/perl5

# 安装必要的依赖项
# bedtools cd-hit blast mcl GNUparallel prank mafft fasttree
# python3.6 下 installed software： fasttree, mafft



CD-HIT is a very widely used program for clustering and comparing protein or nucleotide sequences.
ref: http://weizhongli-lab.org/cd-hit/
ref: https://github.com/weizhongli/cdhit


MMseqs2 (Many-against-Many sequence searching) is a software suite to search and cluster huge protein and nucleotide sequence sets.
ref: https://github.com/soedinglab/MMseqs2
ref: https://github.com/soedinglab/mmseqs2/wiki
cd /work/workspace/zhurj/software
# static build with AVX2 (fastest)
wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz; 
tar xvfz mmseqs-linux-avx2.tar.gz; 
#export PATH=$(pwd)/mmseqs/bin/:$PATH
ln -s /work/workspace/zhurj/software/mmseqs/bin/* /work/workspace/zhurj/bin/


#InterProScan 5
ref: https://github.com/ebi-pf-team/interproscan/wiki/InstallationRequirements
ref: https://github.com/ebi-pf-team/interproscan/wiki/HowToDownload
# ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.36-75.0/interproscan-5.36-75.0-64-bit.tar.gz， too slow to be downloaded
# download interproscan-5.36-75.0-64-bit.tar.gz from http://download.systemsbiology.nl/~jasperk/sapp/
# instruction of installation: https://www.jianshu.com/p/2a5c1e7866b1
# download panther-data-14.1.tar.gz ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/, 泽鑫帮忙下载的，下载速度太慢

cd /work/workspace/zhurj/software/my_interproscan
tar -pxvzf interproscan-5.36-75.0-64-bit.tar.gz
tar -pxvzf panther-data-14.1.tar.gz

OPTIONAL: Installing Panther in a different location
You may wish to install the Panther data files in a different location on the file system. If you do this, you will need to edit the [InterProScan5 home]/interproscan.properties file and set the following property to point to the correct location:

panther.temporary.file.directory=
panther.models.dir=PATH_TO/panther/14.1/
panther.hmm.path=PATH_TO/panther/14.1/panther.hmm
panther.names.tab=PATH_TO/panther/14.1/names.tab

# run interproscan5 
# ref： https://interproscan5-docs.readthedocs.io/en/latest/HowToRun.html ********************************************************
source activate py2
/work/workspace/zhurj/software/my_interproscan/interproscan-5.36-75.0/interproscan.sh -i /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/fna/PROKKA_08072020/interProScan_in/in.fna -f tsv -dp -b /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/fna/PROKKA_08072020/interProScan_out/test 

# ./interproscan.sh -appl CDD,COILS,Gene3D,HAMAP,MobiDBLite,PANTHER,Pfam,PIRSF,PRINTS,ProDom,PROSITEPATTERNS,PROSITEPROFILES,SFLD,SMART,SUPERFAMILY,TIGRFAM -i /path/to/sequences.fasta

# eggNOG mapper
# eggNOG注释的结果，包括了一些匹配和得分信息，以及GO，KEGG，BiGG，COG，KOG，NOG等功能注释结果。但不建议用它的GO和KEGG结果，因为这两个数据库是生信领域更新最快的，信息最全，eggNOG注释的结果可能会跟不上。可以采纳下它的COG、KOG、NOG的注释信息，因为COG/KOG几乎没有更新了，还停留在2003-2014年
ref: https://github.com/jhcepas/eggnog-mapper/releases.
ref: https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2#Installation
ref: http://eggnog5.embl.de/download/eggnog_5.0/ # annotation reference of hmmer ***********
ref: http://www.chenlianfu.com/?p=2804
instrocution: https://www.cnblogs.com/jessepeng/p/12753721.html
介绍
https://pianshen.com/article/5892362478/

referece download URL:
# http://eggnog5.embl.de/download/emapperdb-5.0.0/ # annotation reference of DIOMAND
# 时隔四年，NOG数据库更新啦！ https://zhuanlan.zhihu.com/p/64578206 # 介绍有用
# 数据下载信息： http://eggnog5.embl.de/download/eggnog_5.0/
gbff/                                              16-Apr-2019 08:02       -
id_mappings/                                       23-Apr-2020 09:14       -
per_tax_level/                                     19-Feb-2019 11:15       -
raw_data/                                          26-Mar-2020 19:02       -
e5.level_info.tar.gz                               19-Apr-2019 11:06      1M
e5.og_annotations.tsv                              15-Sep-2018 05:45    196M
e5.proteomes.faa                                   02-Mar-2019 05:58      9G
e5.sequence_aliases.tsv                            26-Jul-2018 11:15     14G
e5.taxid_info.tsv                                  30-Jul-2018 09:30      1M
e5.viruses.faa                                     15-Jul-2018 10:48     36M
其中e5.proteomes.faa为所有的蛋白组序列，e5.viruses.faa为所有的病毒蛋白序列，e5.taxid_info.tsv为Taxid对应的物种名称以及完整的谱系信息，e5.og_annotations.tsv为所有的NOG group信息，其第一列为Taxid，第二列为NOG groups，第三列为COG归属，第四列为Function。

2_members.tsv在per_tax_level/2 目录下
解压后的文件2_members.tsv有五列（如下所示），其中第一列为Taxid，因为我们下载的是细菌bacteria所以第一列均为2，第二列为NOG group，第三列为该NOG group所包含的蛋白序列数目，第四列为该NOG group所包含的物种数目，第五列为该NOG group所包含的蛋白序列id，第六列为该NOG group所包含的物种的Taxid
eggNOG diamond 建库

e5diamond_makedb.sh
```
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
diamond makedb --in /work/workspace/zhurj/database/eggnog/e5/e5_OG_fasta/eggnog5.proteins.all.faa -d /work/workspace/zhurj/database/eggnog/e5/diamond/eggnog5 -p 20
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

```
cd /work/workspace/zhurj/database/eggnog/e5/p
time srun -o e5diamond_makedb.out -e e5diamond_makedb.err -N 1 -c 20 -p slurm256 -w mnclient01 bash e5diamond_makedb.sh &

source activate py2
python setup.py build
python setup.py install
python download_eggnog_data.py

#cd /work/workspace/zhurj/software/eggnog-mapper-master/data && wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog.db.gz http://eggnogdb.embl.de/download/emapperdb-5.0.0/eggnog.db.gz && echo Decompressing... && gunzip eggnog.db.gz

#Downloading fasta files " at /work/workspace/zhurj/software/eggnog-mapper-master/data/eggnog_proteins.dmnd...
#cd /work/workspace/zhurj/software/eggnog-mapper-master/data && wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog_proteins.dmnd.gz  http://eggnogdb.embl.de/download/emapperdb-5.0.0/eggnog_proteins.dmnd.gz && echo Decompressing... && gunzip eggnog_proteins.dmnd.gz

# 下载数据库
# http://eggnog5.embl.de/download/emapperdb-4.5.1/
# http://eggnog5.embl.de/download/emapperdb-4.5.1/hmmdb_levels/
# http://eggnog5.embl.de/download/emapperdb-4.5.1/hmmdb_levels/bact_50/
cd eggnog-mapper ./download_eggnog_data.py bact  #euk真核，bact原核，arch古菌，viruses病毒

#注释
#python emapper.py -i test.fa --output ./ -d euk  #默认以HMMER搜索
#python emapper.py -m diamond -i test.fa --output ./ -d euk  #指定搜索数据库类型，可大类、小类
#python emapper.py -i test.fa --output ./ -d maNOG #哺乳动物NOG
#python emapper.py -i test.fa --output ./ -d maNOG --usemem --cpu 10  #内存和线程
worked command
source activate py2
python /work/workspace/zhurj/software/eggnog-mapper-master/emapper.py -m diamond -d bact -i /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/fna/PROKKA_08072020/interProScan_in/in.fna --output /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/fna/PROKKA_08072020/emapper/test --cpu 20

hmmer worked
export PATH=/work/workspace/zhurj/software/eggnog-mapper-master/bin/:$PATH
source activate py2
python /work/workspace/zhurj/software/eggnog-mapper-master/emapper.py -i /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/fna/PROKKA_08072020/interProScan_in/in.fna --output /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/fna/PROKKA_08072020/emapper/test_hmmer -d bact --usemem --cpu 20 --override
若有多组数据需要进行eggNOG注释，可以分两步运行。先将数据库放到内存中，再分别对多组数据进行运算，这样只读取一次数据库到内存中节约运算时间。
$ emapper.py -d euk --cpu 80 --servermode
待数据库加载完毕后，可以在其它终端中执行数据分析。
$ emapper.py -d euk:localhost:51400 -i proteins1.fasta -o eggNOG1 --no_file_comments
$ emapper.py -d euk:localhost:51400 -i proteins2.fasta -o eggNOG2 --no_file_comments

注意方法选择：diamond在序列少时相对较慢，但序列多时相对较快。HMMER方法对于亲源较远序列预测成功率更高，但数据量大时计算时间长，在线限制一次最多5000条序列。

mirror of R
https://mirrors.tuna.tsinghua.edu.cn/CRAN/  TUNA Team, Tsinghua University
https://mirrors.ustc.edu.cn/CRAN/   University of Science and Technology of China
https://mirror-hk.koddos.net/CRAN/  KoDDoS in Hong Kong
https://mirrors.e-ducation.cn/CRAN/ Elite Education
https://mirror.lzu.edu.cn/CRAN/ Lanzhou University Open Source Society
https://mirrors.tongji.edu.cn/CRAN/ Tongji University

install.packages("devtools", repo = "https://mirrors.ustc.edu.cn/CRAN/")
devtools::install_github('ggloor/CoDaSeq/CoDaSeq',host = "https://api.github.com")


# 2020.08.11
metaData
environment location: /work/workspace/zhurj/bin/miniconda3/envs/metaAna
bioBakery Workflows 
ref: https://bitbucket.org/biobakery/biobakery/wiki/biobakery_workflows
ref: https://bitbucket.org/biobakery/biobakery_workflows/wiki/Home#!requirements
bioBakery workflows is a collection of workflows and tasks for executing common microbial community analyses using standardized, validated tools and parameters. Quality control and statistical summary reports are automatically generated for most data types, which include 16S amplicons, metagenomes, and metatranscriptomes.

config conda
#$ conda config --add channels defaults
#$ conda config --add channels bioconda
#$ conda config --add channels conda-forge
#$ conda config --add channels biobakery
conda config --add channels biobakery
#conda config --remove channels https://mirrors.tuna.tsinghua.edu.cn/help/anaconda/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
conda install -c biobakery workflows
# Before installing the workflows with conda, make sure to configure your channels so biobakery is at the top of the list.

pip install -i https://pypi.tuna.tsinghua.edu.cn/simple biobakery_workflows

2020.08.12
https://huttenhower.sph.harvard.edu/humann

# install parallel
cd /work/workspace/zhurj/software
tar xjf parallel-latest.tar.bz2
cd /work/workspace/zhurj/software/parallel-20200722
/work/workspace/zhurj/software/parallel-20200722/configure --prefix=/work/workspace/zhurj/bin && make && make install


分泌蛋白预测
###### 满足有信号肽但无跨膜结构域的蛋白 → 经典分泌蛋白
方法介绍：https://zhuanlan.zhihu.com/p/52744204
ref: http://html.rhhz.net/SWJSTB/html/2016-8-129.htm
使用big-PI Predictor软件[23]对具有信号肽且含有0或1个跨膜结构域的氨基酸序列进行GPI锚定位点预测。最后，使用Protcomp软件完成对蛋白的亚细胞定位预测。
下载申请ref： https://services.healthtech.dtu.dk/software.php
信号肽预测
SignalP（Signal peptide and cleavage sites in gram+, gram- and eukaryotic amino acid sequences）
download ref： https://services.healthtech.dtu.dk/software.php
解释：

1. Sec/SPI: "standard" secretory signal peptides transported by the Sec translocon and cleaved by Signal Peptidase I (Lep)
   Sec/SPII: lipoprotein signal peptides transported by the Sec translocon and cleaved by Signal Peptidase II (Lsp)
   Tat/SPI: Tat signal peptides transported by the Tat translocon and cleaved by Signal Peptidase I (Lep)
2. SP(Sec/SPI) / LIPO(Sec/SPII) / TAT(Tat/SPI) (depending on what type of signal peptide is predicted), CS (the cleavage site) and OTHER (the probability that the sequence does not have any kind of signal peptide).
ref： http://www.cbs.dtu.dk/services/SignalP/（1.）
ref: http://www.cbs.dtu.dk/services/SignalP/instructions.php (2.)

command: signalp -org gram+ -fasta /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/prokka/MNH05168/MNH05168.faa -format short -prefix /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/signalp/MNH05168/MNH05168_short

膜结构域预测
TMHMM（Prediction of transmembrane helices in proteins）
cd /work/workspace/zhurj/software
tar zxvf /work/workspace/zhurj/software/tmhmm-2.0c.Linux.tar.gz && cd tmhmm-2.0c/bin
1. nsert the correct path(/work/workspace/zhurj/bin/miniconda3/bin/perl) for perl 5.x in the first line of the scripts
   bin/tmhmm and bin/tmhmmformat.pl (if not /usr/local/bin/perl)
2. ln -s /work/workspace/zhurj/software/tmhmm-2.0c/bin/tmhmm /work/workspace/zhurj/bin/tmhmm 

command: /work/workspace/zhurj/software/tmhmm-2.0c/bin/tmhmm /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/prokka/MNH05168/MNH05168.faa > /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/tmhmm/MNH05168_tmhmm.txt 

非经典分泌蛋白预测
SecretomeP（Prediction of non-classical protein secretion）
cd /work/workspace/zhurj/software
tar zxvf secretomep-1.0h.Linux.tar.gz


SecretomeP need prop-1.0c.Linux.tar.Z
cd /work/workspace/zhurj/software
tar xvf prop-1.0c.Linux.tar.Z
chmod 1777 tmp

vi prop
#! /work/workspace/zhurj/bin/miniconda3/envs/tcsh/bin/tcsh -f
```
else if ( $SYSTEM == "Linux" ) then             # typical Linux
   setenv AWK /usr/bin/gawk
   setenv ECHO "/bin/echo -e"
   setenv GNUPLOT /usr/bin/gnuplot              # /usr/local/bin/gnuplot-3.7
   setenv PPM2GIF /usr/bin/ppmtogif
   setenv SIGNALP /work/workspace/zhurj/bin/signalp
```



2020.08.13
under no ENV
export PERL5LIB=/home/foobar/code
cpanm Bio::DB::EUtilities -l /work/workspace/zhurj/lib/perl --mirror http://mirrors.aliyun.com/CPAN

cmake install
ref: https://blog.csdn.net/Anton8801/article/details/94282363
cd /work/workspace/zhurj/software

tar zxvf cmake-3.12.3.tar.gz
cd cmake-3.12.3
./bootstrap --prefix=/work/workspace/zhurj/bin
gmake
make
make install


pftools
ref： https://github.com/sib-swiss/pftools3/blob/master
cd /work/workspace/zhurj/software
unzip pftools3-master.zip
cd pftools3-master
mkdir build
cd build
cmake ..
cmake -DCMAKE_INSTALL_PREFIX:PATH=/work/workspace/zhurj/bin ..

# psort install
ref： https://psort.org/
install ref： 1. https://psort.org/downloads/index.html
              2. https://psort.org/downloads/INSTALL.html#installation

cd /work/workspace/zhurj/software
tar zxvf libpsortb-1.0.tar.gz
cd libpsortb-1.0
./configure --prefix=/work/workspace/zhurj/bin
make 
make install

cd /work/workspace/zhurj/software
tar zxvf bio-tools-psort-all.3.0.6.tar.gz
cd bio-tools-psort-all


2020.08.17
islandpath
ref:  https://github.com/brinkmanlab/islandpath
export PERL5LIB=/work/workspace/zhurj/lib/perl/lib/perl5
software required:
Data::Dumper
Log:Log4perl
Config::Simple
Moose
MooseX::Singleton
Bio::Perl

perl -e "use Data::Dumper;"
perl -e "use Log::Log4perl;"
perl -e "use Config::Simple;"
perl -e "use Bio::Perl;"
perl -e "use Moose;"
perl -e "use MooseX::Singleton;"


##
to check whether module installed
perl -e "use MODULE;"

Where is the MODULE
perldoc -l MODULE

find `perl -le 'print '@INC''` | grep Moose
##
/work/workspace/zhurj/lib/perl/lib/perl5/Log/Log4perl.pm
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/lib/5.26.2/x86_64-linux-thread-multi/Data/Dumper.pm
/work/workspace/zhurj/lib/perl/lib/perl5/Bio/Perl.pm
/work/workspace/zhurj/lib/perl/lib/perl5/Config/Simple.pm
/work/workspace/zhurj/lib/perl/lib/perl5/MooseX/Singleton.pm
/work/workspace/zhurj/lib/perl/lib/perl5/x86_64-linux-thread-multi/Moose.pm

cpanm Bio::Perl -l /work/workspace/zhurj/lib/perl --mirror http://mirrors.aliyun.com/CPAN
cpanm Config::Simple -l /work/workspace/zhurj/lib/perl --mirror http://mirrors.aliyun.com/CPAN
cpanm Moose -l /work/workspace/zhurj/lib/perl --mirror http://mirrors.aliyun.com/CPAN
cpanm MooseX::Singleton -l /work/workspace/zhurj/lib/perl --mirror http://mirrors.aliyun.com/CPAN


install MooseX::Singleton need TAP::Harness::Env 在 /work/workspace/zhurj/lib/perl 下
install MooseX::Singleton need Package::Stash，MRO::Compat 在 /work/workspace/zhurj/bin/miniconda3/lib/5.26.2
cpanm Package::Stash -l /work/workspace/zhurj/lib/perl --mirror http://mirrors.aliyun.com/CPAN
cpanm MRO::Compat -l /work/workspace/zhurj/lib/perl --mirror http://mirrors.aliyun.com/CPAN


## ref: https://mirrors.tuna.tsinghua.edu.cn/help/anaconda/
conda install -c bioconda islandpath


islandpath /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/prokka/MNH05168/MNH05168.gbk /work/workspace/zhurj/project/2_swgs/prokka/cm20200806/islandpath/gi_test.txt


构建进化树
http://blog.sciencenet.cn/home.php?mod=space&uid=656335&do=blog&id=1070129
基因组处理
http://blog.grimoire.cc/Bioinformatics-Phylogenetic-Tree/
常用细菌分析软件
https://www.microbialsystems.cn/zh/post/microbial_bioinformatics/
参考文献： A new genomic blueprint of the human gut microbiota（2019， Nature）
v.2.6.3 /work/workspace/zhurj/bin/miniconda3/bin/prodigal
muscle-3.8.1551 /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/muscle

进化树分析参考
https://www.jianshu.com/p/cbb36ac4f82a

RAxML 介绍
https://support.nesi.org.nz/hc/en-gb/articles/115001854444-RAxML
raxmlHPC-AVX or raxmlHPC-SSE3 with one task on only one CPU.
raxmlHPC-PTHREADS-AVX or raxmlHPC-PTHREADS-SSE3 with one task running on multiple CPUs.
raxmlHPC-MPI-AVX or raxmlHPC-MPI-SSE3 with multiple tasks, each running on one CPU.
raxmlHPC-HYBRID-AVX or raxmlHPC-HYBRID-SSE3 with multiple tasks, each of which runs on multiple CPUs.

Physift
https://phylosift.wordpress.com/tutorials/downloading-phylosift/
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/iqtree

IQ-TREE
参考： https://www.jianshu.com/p/fa057daadfbc

2020.08.21
MUMmer 3.23
ref: https://indexofire.github.io/pathongs/C06_Genome-Compare/01_mummer/
ref: https://www.jianshu.com/p/9701c9add82d
/work/assembly/current/MNH/MNH051/MNH05168/genomic.fna
/work/database/ncbi/genome/all/GCA/GCA_900/GCA_900447/GCA_900447735/genomic.fna.gz /

ref: http://blog.sciencenet.cn/blog-2970729-1174911.html
mummer /work/database/ncbi/genome/all/GCA/GCA_900/GCA_900447/GCA_900447735/genomic.fna.gz /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna
# step1: run nucmer for alignment
nucmer --prefix=ref_qry /work/workspace/zhurj/reference/NCBI/genome/GCA_900447735/genomic.fna /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna
nucmer --mum -p /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_NCTC7171/nucmer/MNH05168_NCTC7171 /work/workspace/zhurj/reference/NCBI/genome/GCA_900447735/genomic.fna /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna
delta-filter -i 85 -l 8000 -o 85 /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_NCTC7171/nucmer/MNH05168_NCTC7171.delta > /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_NCTC7171/nucmer/MNH05168_NCTC7171.best_delta
mummerplot -p /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_NCTC7171/nucmer/MNH05168_NCTC7171 /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_NCTC7171/nucmer/MNH05168_NCTC7171.best_delta -t postscript
ps2pdf /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_NCTC7171/nucmer/MNH05168_NCTC7171.ps /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_NCTC7171/nucmer/MNH05168_NCTC7171.pdf
convert -density 300 /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_NCTC7171/nucmer/MNH05168_NCTC7171.pdf /work/workspace/zhurj/project/2_swgs/prokka/MNH05168_NCTC7171/nucmer/MNH05168_NCTC7171.png

phiSpy (前噬菌体预测软件)
ref: https://www.jianshu.com/p/0bbbe1155950
ref： http://blog.sciencenet.cn/blog-3334560-1245695.html （PhiSpy：在细菌基因组中识别噬菌体）
整合在宿主基因组上的温和噬菌体的核酸称之为前噬菌体，能随宿主细菌DNA进行同步复制或分裂传代。前噬菌体序列的存在可能会使一些细菌获取抗生素抗性，增强对环境的适应性，提高粘附力或使细菌成为致病菌
PhiSpy.py -o prophages genome_prokka.gbk  # 默认
PhiSpy.py -o prophages genome_prokka.gbk --output_choice 255  # 输出所有结果


cpanm Switch --mirror http://mirrors.aliyun.com/CPAN

echo $PATH
/work/workspace/zhurj/bin/miniconda3/condabin:/work/workspace/zhurj/bin/miniconda3/bin:/home/zhurj/.local/bin/pip:/work/program/current:/work/program/bin:/work/program/instal/metaphlan2:/work/workspace/zhurj/bin:/work/workspace/zhurj/bin/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/home/zhurj/.local/bin:/home/zhurj/bin

/work/workspace/zhurj/bin/miniconda3/bin
/work/workspace/zhurj/bin/miniconda3/envs/roary
conda create -n roary
source activate roary
export PATH=/work/workspace/zhurj/bin/miniconda3/envs/roary/bin:/work/workspace/zhurj/bin/miniconda3/bin:/usr/bin:
conda install -c bioconda roary
cpan -f Bio::Roary --mirror http://mirrors.aliyun.com/CPAN

RepeatMasker is a program that screens DNA sequences for interspersed repeats and low complexity DNA sequences.
/work/workspace/zhurj/bin/miniconda3/envs/repeatmasker
#export /work/workspace/zhurj/bin/miniconda3/envs/repeatmasker/bin:/work/workspace/zhurj/bin/miniconda3/condabin:/work/workspace/zhurj/bin/miniconda3/envs/roary/bin:/work/workspace/zhurj/bin/miniconda3/bin:/usr/bin:.
export PATH=/work/workspace/zhurj/bin/miniconda3/envs/repeatmasker/bin:/usr/bin:
conda install -c bioconda repeatmasker

2020.08.28
circos
ref: ananconda.org/bioconda/circos
ref: http://www.chenlianfu.com/?p=2342
conda create -n circos
source activate circos
conda install -c bioconda circos

perl /work/workspace/zhurj/script/5_wgspro/1_filter/mnid2rawdata.pl
perl /work/workspace/zhurj/script/5_wgspro/1_filter/datafilter.pl -i /work/workspace/zhurj/project/1_metadata/metapipe/input/test_r1r2_2c.txt -o /work/workspace/zhurj/project/1_metadata/metapipe
perl /work/workspace/zhurj/script/5_wgspro/1_filter/mnid2rawdata.pl -i /work/workspace/zhurj/project/1_metadata/metapipe/input/ifile -t 0 -r metagenome -o /work/workspace/zhurj/project/1_metadata/metapipe/input -n test
/work/workspace/moonis/bin/genoassemble/soft/SPAdes-3.13.1-Linux/bin/spades.py

source activate python3.6
ref： https://anaconda.org/bioconda/cutadapt
v2.10
conda install -c bioconda cutadapt
cutadapt -h

source activate python3.6
ref: https://anaconda.org/bioconda/kneaddata
v0.6.1
conda install -c bioconda kneaddata
kneaddata -h

source activate python3.6
https://anaconda.org/bioconda/metabat2
v2.15
conda install -c bioconda metabat2
metabat2 -h

source activate python3.6
bwa
bowtie2

source activate python3.6
bwa
Version: 0.7.17-r1188

source activate python3.6
samtools 
Version: 1.3.1

source activate python3.6
checkM
CheckM v1.1.1

gtdbtk: version 1.2.0
GTDB-tk version 1.2.0

START2
https://pubmlst.org/software/analysis/start2/downloads.shtml

source activate python3.6
mlst
ref: https://anaconda.org/bioconda/mlst
ref: https://github.com/tseemann/mlst
conda install -c bioconda mlst
mlst
v2.16.1

mlst /work/workspace/zhurj/software/mlst-master/db 
mlst --label Anthrax  /work/database/ncbi/genome/all/GCA/GCA_001/GCA_001941/GCA_001941925/genomic.fna.gz

metaphlan3


source activate  python3.6
INFERNAL
REF: https://anaconda.org/bioconda/infernal
v1.1.2
conda install -c bioconda infernal

去除环境变量中重复项
export PATH=$(echo $PATH | tr : "\n"| sort | uniq | tr "\n" :)
'''
acount install
base                  *  /work/program/instal/Miniconda3-4.7.12
ETE-3                    /work/program/instal/Miniconda3-4.7.12/envs/ETE-3
GraPhlAn                 /work/program/instal/Miniconda3-4.7.12/envs/GraPhlAn
GraPhlAn-1.1.3           /work/program/instal/Miniconda3-4.7.12/envs/GraPhlAn-1.1.3
HUMAnN                   /work/program/instal/Miniconda3-4.7.12/envs/HUMAnN
HUMAnN-2                 /work/program/instal/Miniconda3-4.7.12/envs/HUMAnN-2
HUMAnN-2.8.1             /work/program/instal/Miniconda3-4.7.12/envs/HUMAnN-2.8.1
HUMAnN-3.0.0             /work/program/instal/Miniconda3-4.7.12/envs/HUMAnN-3.0.0
PhyloPhlAn               /work/program/instal/Miniconda3-4.7.12/envs/PhyloPhlAn
PhyloPhlAn-3.0           /work/program/instal/Miniconda3-4.7.12/envs/PhyloPhlAn-3.0
Prokka                   /work/program/instal/Miniconda3-4.7.12/envs/Prokka
Prokka-1.14              /work/program/instal/Miniconda3-4.7.12/envs/Prokka-1.14
antismash-5              /work/program/instal/Miniconda3-4.7.12/envs/antismash-5
antismash-5.1.1          /work/program/instal/Miniconda3-4.7.12/envs/antismash-5.1.1
circos                   /work/program/instal/Miniconda3-4.7.12/envs/circos
circos-0.69.8            /work/program/instal/Miniconda3-4.7.12/envs/circos-0.69.8
gnuplot                  /work/program/instal/Miniconda3-4.7.12/envs/gnuplot
gnuplot-5.2.7            /work/program/instal/Miniconda3-4.7.12/envs/gnuplot-5.2.7
                         /work/program/instal/miniconda3
                         /work/program/instal/miniconda3/envs/CheckM
                         /work/program/instal/miniconda3/envs/GNUPLOT
                         /work/program/instal/miniconda3/envs/MMvec
                         /work/program/instal/miniconda3/envs/PICRUSt2
                         /work/program/instal/miniconda3/envs/Prokka
                         /work/program/instal/miniconda3/envs/R.v360
                         /work/program/instal/miniconda3/envs/antiSMASH-4.2.0
                         /work/program/instal/miniconda3/envs/bioBakery
                         /work/program/instal/miniconda3/envs/circos
                         /work/program/instal/miniconda3/envs/qiime2-2019.7
'''
```
阿里云 http://mirrors.aliyun.com/pypi/simple/
中国科技大学 https://pypi.mirrors.ustc.edu.cn/simple/
豆瓣(douban) http://pypi.douban.com/simple/
清华大学 https://pypi.tuna.tsinghua.edu.cn/simple/
中国科学技术大学 http://pypi.mirrors.ustc.edu.cn/simple/
```
humann3
ref: HUMAnN 3.0 (alpha)安装及使用 https://zhuanlan.zhihu.com/p/240910229
ref: https://huttenhower.sph.harvard.edu/humann3/
humann3 version humann-3.0.0a4
numpy version:  numpy-1.19.1
diamond version 0.9.24

diamond ref: https://github.com/bbuchfink/diamond/releases/tag/v0.9.24

ref: https://blog.csdn.net/weixin_44223966/article/details/106875039
conda create -n humann python=3.7
conda activate humann
pip install humann -i https://pypi.mirrors.ustc.edu.cn/simple/
pip install numpy -i https://pypi.mirrors.ustc.edu.cn/simple/
pip install cython -i https://pypi.mirrors.ustc.edu.cn/simple/ # if exist, 
pip install biom-format -i https://pypi.mirrors.ustc.edu.cn/simple/
pip install h5py -i https://pypi.mirrors.ustc.edu.cn/simple/
pip install metaphlan -i https://pypi.tuna.tsinghua.edu.cn/simple/

https://github.com/bbuchfink/diamond/releases/tag/v0.8.38

```
我在主目录下面建立了一个humann文件夹，里面放了三个目录：chocophlan、uniref和misc三个目录，分别放置ChocoPhlAn、UniRef和Utility_Mapping三个数据库，而这三个数据库下载地址如下：
chocophlan : full = http://huttenhower.sph.harvard.edu/humann2_data/chocophlan/full_chocophlan.v296_201901.tar.gz
uniref : uniref50_diamond = http://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref_annotated/uniref50_annotated_v201901.tar.gz
uniref : uniref90_diamond = http://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref_annotated/uniref90_annotated_v201901.tar.gz
uniref : uniref50_ec_filtered_diamond = http://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref_ec_filtered/uniref50_ec_filtered_201901.tar.gz
uniref : uniref90_ec_filtered_diamond = http://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref_ec_filtered/uniref90_ec_filtered_201901.tar.gz
utility_mapping : full = http://huttenhower.sph.harvard.edu/humann2_data/full_mapping_v201901.tar.gz

上述文件下载解压后放到对应的三个目录下面即可。chocophlan下面是10669个物种的特征性序列文件，uniref下面是四个diamond的索引文件，utility_mapping是大17个mapping用的文件。随后在HUMAnN3的设置当中对上述目录进行相应的设置。
注意在前面diamond部分已经说过的，对uniref的uniref90_ec_filtered_201901数据索引文件使用0.8.38提取序列，再用0.9.24重建索引，否则不能顺利运行。所以需要从https://github.com/bbuchfink/diamond/releases/tag/v0.8.38这里下载该版本(diamond文件)，并使用上述的代码重新用0.9.24版建立索引。
/old/diamond getseq -d uniref90_201901.dmnd | /new/diamond makedb -d out_db
*************uniref90_ec_filtered_diamond 可以用 diamond 0.9.24 提取，不能用 0.3.38 提取

这些文件下载之后，按照代码本身的要求，需要放到script代码所在目录的metaphlan_databases目录下，具体位置应该是在~/anaconda3/envs/humann/lib/python3.7/site-packages/metaphlan/metaphlan_databases/这个目录里面。千万注意要把file_list.txt这个文件也一起放进去，因为我记得代码里一开始就是要识别这个文件在不在。
```

ln -s /work/workspace/zhurj/database/humann/201901/chocophlan/* /work/workspace/zhurj/database/humann/current/chocophlan/
source activate humann
#0.9.24
cd /work/workspace/zhurj/database/humann/201901/uniref/ 
/work/workspace/zhurj/software/diamondv0.8.38/diamond getseq --db uniref90_201901.dmnd > uniref90_201901.faa
```
#0.8.36
# diamond makedb --in uniref90_201901.faa --db uniref90_201901
```
cd /work/workspace/zhurj/database/humann/201901/uniref50/
# uniref50_201901.dmnd was uncompressed from uniref50_ec_filtered_201901.tar.gz
diamond getseq --db uniref50_201901.dmnd > uniref50_201901.faa

cd /work/workspace/zhurj/database/humann/201901/uniref50_all/
# uniref50_201901.dmnd was uncompressed from uniref50_annotated_v201901.tar.gz
diamond getseq --db uniref50_201901.dmnd > uniref50_201901.faa


ln -s /work/workspace/zhurj/database/humann/201901/uniref/* /work/workspace/zhurj/database/humann/current/uniref/
ln -s /work/workspace/zhurj/database/humann/201901/misc/* /work/workspace/zhurj/database/humann/current/misc/

# config file
# /work/workspace/zhurj/bin/miniconda3/envs/humann/lib/python3.7/site-packages/humann
```
nucleotide = data/chocophlan_DEMO
protein = data/uniref_DEMO
utility_mapping = data/misc

```
vi /work/workspace/zhurj/bin/miniconda3/envs/humann/lib/python3.7/site-packages/humann
nucleotide = /work/workspace/zhurj/database/humann/current/chocophlan
protein = /work/workspace/zhurj/database/humann/current/uniref
utility_mapping = /work/workspace/zhurj/database/humann/current/misc

```
humann_config --update database_folder nucleotide /work/workspace/zhurj/database/humann/current/chocophlan
humann_config --update database_folder protein /work/workspace/zhurj/database/humann/current/uniref
humann_config --update database_folder utility_mapping /work/workspace/zhurj/database/humann/current/misc
```


mkdir /work/workspace/zhurj/bin/miniconda3/envs/humann/lib/python3.7/site-packages/metaphlan/metaphlan_databases/
mv /work/workspace/zhurj/database/humann/201901/mpa/* /work/workspace/zhurj/bin/miniconda3/envs/humann/lib/python3.7/site-packages/metaphlan/metaphlan_databases -f

humann3 test
humann -i /work/workspace/zhurj/bin/miniconda3/envs/humann/lib/python3.7/site-packages/humann/tests/data/demo.fastq -o /work/workspace/zhurj/project/1_metadata/metapipe/humann3/test 
```
humann*
humann3*
humann3_databases*
humann_associate*
humann_barplot*
humann_benchmark*
humann_build_custom_database*
humann_config*
humann_databases*
humann_genefamilies_genus_level*
humann_humann1_kegg*
humann_infer_taxonomy*
humann_join_tables*
humann_reduce_table*
humann_regroup_table*
humann_rename_table*
humann_renorm_table*
humann_rna_dna_norm*
humann_split_stratified_table*
humann_split_table*
humann_strain_profiler*
humann_test*
humann_unpack_pathways*

```

CD-HIT
v4.8.1
conda install -c bioconda cd-hit
cd-hit

source activate /work/program/instal/miniconda/envs/Prokka
prokka 1.14.6

source activate python3.6
conda install -c bioconda rgi
v5.1.0
不能用

source activate py2
conda install -c bioconda rgi
v3.2.1
# worked

# 使用reference
https://github.com/arpcard/rgi
# 软件安装有先后顺序
conda create -n rgi python=3.6
source activate rgi
# Resistance Gene Identifier - 5.1.1
# RGI version5.1.1
conda install -c bioconda seaborn=0.9.0
conda install -c bioconda biopython=1.72 diamond=0.8.36 samtools=1.9
conda install -c bioconda bamtools=2.5.1 bedtools=2.27.1
conda install -c bioconda jellyfish=2.2.10 bowtie2=2.3.4.3 bwa=0.7.17 mock=2.0.0 prodigal=2.6.3

download rgi-master.zip from https://github.com/arpcard/rgi to /work/workspace/zhurj/software/
cd /work/workspace/zhurj/software && unzip rgi-master.zip
cd /work/workspace/zhurj/software/rgi-master
pip3 install pytest
pip3 install pyfaidx
pip3 install -i https://pypi.tuna.tsinghua.edu.cn/simple pyahocorasick
pip3 install -i https://pypi.tuna.tsinghua.edu.cn/simple filetype

cd /work/workspace/zhurj/reference/CARD/20200916
tar -xvf card-data.tar.bz2 ./card.json
rgi load --card_json /work/workspace/zhurj/reference/CARD/20200916/card.json --local
cp /work/workspace/zhurj/reference/CARD/20200916/card.json /work/workspace/zhurj/software/rgi-master/app/_data
pytest -v -rxs

rgi card_annotation -i /work/workspace/zhurj/reference/CARD/20200916/card.json > card_annotation.log 2>&1
rgi load -i /work/workspace/zhurj/reference/CARD/20200916/card.json --card_annotation /work/workspace/zhurj/reference/CARD/20200916/card_database_v3.1.0.fasta --local
#system
rgi load -i /work/workspace/zhurj/reference/CARD/20200916/card.json --card_annotation /work/workspace/zhurj/reference/CARD/20200916/card_database_v3.1.0.fasta

# wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants
mkdir /work/workspace/zhurj/reference/CARD/wildcard
#download wildcard_data.tar.bz2 to /work/workspace/zhurj/reference/CARD from https://card.mcmaster.ca/latest/variants
cd /work/workspace/zhurj/reference/CARD
tar -xjf card-prevalence.tar.bz2 -C wildcard
gunzip wildcard/*.gz

cd /work/workspace/zhurj/reference/CARD
rgi wildcard_annotation -i wildcard --card_json /work/workspace/zhurj/reference/CARD/20200916/card.json -v version_3.1.0 > wildcard_annotation.log 2>&1
rgi load --wildcard_annotation /work/workspace/zhurj/reference/CARD/wildcard_database_vversion_3.1.0.fasta --wildcard_index /work/workspace/zhurj/reference/CARD/wildcard/index-for-model-sequences.txt --card_annotation /work/workspace/zhurj/reference/CARD/20200916/card.json --local
#system
rgi load --wildcard_annotation /work/workspace/zhurj/reference/CARD/wildcard_database_vversion_3.1.0.fasta --wildcard_index /work/workspace/zhurj/reference/CARD/wildcard/index-for-model-sequences.txt --card_annotation /work/workspace/zhurj/reference/CARD/20200916/card.json 


# load the k-mer reference data:
cd /work/workspace/zhurj/reference/CARD
rgi load --kmer_database /work/workspace/zhurj/reference/CARD/wildcard/61_kmer_db.json --amr_kmers /work/workspace/zhurj/reference/CARD/wildcard/all_amr_61mers.txt --kmer_size 61 --local --debug > kmer_load.61.log 2>&1
# Local or working directory:
rgi database --version --local
# system
rgi load --kmer_database /work/workspace/zhurj/reference/CARD/wildcard/61_kmer_db.json --amr_kmers /work/workspace/zhurj/reference/CARD/wildcard/all_amr_61mers.txt --kmer_size 61 --debug > kmer_load.61.log 2>&1
# Local or working directory:
rgi clean --local

#run 
source activate rgi
rgi main --input_sequence /work/workspace/zhurj/project/2_swgs/MNH04863/assemble/MNH-4863.fna --output_file /work/workspace/zhurj/project/2_swgs/MNH04863/RGI/gri_result --input_type contig --clean

rgi main --input_sequence /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna --output_file /work/workspace/zhurj/project/2_swgs/MNH05168/RGI/gri_result --input_type contig --clean

2020.11.17
CARD: https://card.mcmaster.ca/download
1. Download CARD Data (see README), use https://card.mcmaster.ca/latest/data for automated downloads
2. Download Ontology Files (freely available, see README), use https://card.mcmaster.ca/latest/ontology for automated downloads
CRISPRCasFinder
download CRISPRCasFinder.zip to from https://github.com/dcouvin/CRISPRCasFinder

integron_finder 
INSTALLATION
ref: https://anaconda.org/bioconda/integron_finder
2.0rc6
conda install -c bioconda integron_finder
USAGE:
ref: https://github.com/gem-pasteur/Integron_Finder
integron_finder --local-max --func-annot --cpu 36 /work/assembly/current/MNH/MNH051/MNH05168/genomic.fna --outdir /work/workspace/zhurj/project/2_swgs/MNH05168/integron_finder
!!!!!!!!doesn't work

download rarlinux-x64-5.9.1.tar.gz from https://www.rarlab.com/download.html to /work/workspace/zhurj/software
cd /work/workspace/zhurj/software
tar zxvf rarlinux-x64-5.9.1.tar.gz
ln -s /work/workspace/zhurj/software/rar/rar /work/workspace/zhurj/bin
ln -s /work/workspace/zhurj/software/rar/unrar /work/workspace/zhurj/bin


source activate python3.6
# https://anaconda.org/bioconda/ecopy
# Species Diversity: https://ecopy.readthedocs.io/en/stable/diversity.html
conda install -c bioconda ecopy
pip install BeautifulSoup4


# phylophlan 3.0 
# PhyloPhlAn 3.0 微生物组系统发育分析: https://cloud.tencent.com/developer/article/1637229
# https://github.com/biobakery/phylophlan/wiki#configuration-file
conda create -n "phylophlan" -c bioconda phylophlan=3.0

conda activate phylophlan
phylophlan_write_default_configs.sh /work/workspace/zhurj/database/phylophlan/20201015
phylophlan_write_default_configs.sh  /work/workspace/zhurj/bin/miniconda3/envs/phylophlan/lib/python3.8/site-packages/phylophlan/phylophlan_configs
reference download: https://www.dropbox.com/s/x7cvma5bjzlllbt/phylophlan_databases.txt

cd /work/workspace/zhurj/bin/miniconda3/envs/phylophlan/lib/python3.8/site-packages/phylophlan/phylophlan_databases/phylophlan
phylophlan_setup_database -i phylophlan \
    -d phylophlan \
    -e faa.bz2 \
    -t a

phylophlan_setup_database -i amphora2 \
    -d amphora2 \
    -e faa.bz2 \
    -t a



Source activate R3.6
R
install.packages("Hmisc",repos="http://mirrors.ustc.edu.cn/CRAN/")
(rcorr work, 装了Hmisc以后)

install.packages("corrplot",repos="http://mirrors.ustc.edu.cn/CRAN/")
install.packages("PerformanceAnalytics",repos="http://mirrors.ustc.edu.cn/CRAN/")
install.packages("ggsci",repos="http://mirrors.ustc.edu.cn/CRAN/")
install.packages("ggpubr",repos="http://mirrors.ustc.edu.cn/CRAN/")
install.packages("psych",repos="http://mirrors.ustc.edu.cn/CRAN/")
install.packages("devtools",repos="http://mirrors.ustc.edu.cn/CRAN/")
install.packages("sparcc",repos="http://mirrors.ustc.edu.cn/CRAN/")
install.packages("BiocManager",repos="http://mirrors.ustc.edu.cn/CRAN/")
BiocManager::install("metagenomeSeq")
install.packages("heatmaply",repos="http://mirrors.ustc.edu.cn/CRAN/")
install.packages("optparse",repos="http://mirrors.ustc.edu.cn/CRAN/")
install.packages("tidyverse",repos="http://mirrors.ustc.edu.cn/CRAN/")
install.packages("ggthemes",repos="http://mirrors.ustc.edu.cn/CRAN/")
install.packages("randomcoloR",repos="http://mirrors.ustc.edu.cn/CRAN/")

install.packages("tidybayes",repos="http://mirrors.ustc.edu.cn/CRAN/")
install.packages("ggridges",repos="http://mirrors.ustc.edu.cn/CRAN/")
install.packages("rstan",repos="http://mirrors.ustc.edu.cn/CRAN/")  失败
install.packages("brms",repos="http://mirrors.ustc.edu.cn/CRAN/") 失败
install.packages("gganimate",repos="http://mirrors.ustc.edu.cn/CRAN/")
install.packages("heatmaply",repos="http://mirrors.ustc.edu.cn/CRAN/")
install.packages("pheatmap",repos="http://mirrors.ustc.edu.cn/CRAN/")
install.packages("hflights",repos="http://mirrors.ustc.edu.cn/CRAN/")
install.packages("MVN",repos="http://mirrors.ustc.edu.cn/CRAN/")
install.packages("diffEnrich",repos="http://mirrors.ustc.edu.cn/CRAN/") 
install.packages("plotrix",repos="http://mirrors.ustc.edu.cn/CRAN/")
install.packages("VennDiagram",repos="http://mirrors.ustc.edu.cn/CRAN/")
# 2021.03.15
# venn 参考  https://blog.csdn.net/elernino/article/details/42456871
# vennerable draw graph: https://www.jianshu.com/p/1a027cf12844
install.packages("Vennerable", repos="http://R-Forge.R-project.org")
install.packages("remotes",repos="http://mirrors.ustc.edu.cn/CRAN/")
remotes::install_github("js229/Vennerable")
options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))) 
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("RBGL")   
RBGL

# 2021.06.30
install.packages("randomForest",repos="http://mirrors.ustc.edu.cn/CRAN/")
# diffEnrich
# https://github.com/SabaLab/diffEnrich


library(devtools)
install_github("vqv/ggbiplot")

R3.6
conda install -c bioconda bioconductor-phyloseq
conda install -c bioconda bioconductor-microbiome

install.packages("microbiome",repos="http://mirrors.ustc.edu.cn/CRAN/")

library('ggpubr')
library('ggbiplot')
library(reshape2)
library(dplyr)

# Load data
data("mtcars")
my_data <- mtcars[, c(1,3,4,5,6,7)]
# print the first 6 rows
head(my_data, 6)
res <- cor(my_data)

symnum(res, abbr.colnames = FALSE)

corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

library(corrplot)
corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
Positive correlations are displayed in blue and negative correlations in red color. Color intensity and the size of the circle are proportional to the correlation coefficients. In the right side of the correlogram, the legend color shows the correlation coefficients and the corresponding colors.
The correlation matrix is reordered according to the correlation coefficient using “hclust” method.
tl.col (for text label color) and tl.srt (for text label string rotation) are used to change text colors and rotations.
Possible values for the argument type are : “upper”, “lower”, “full”

# Insignificant correlation are crossed
corrplot(res2$r, type="upper", order="hclust", 
         p.mat = res2$P, sig.level = 0.01, insig = "blank")
# Insignificant correlations are leaved blank
corrplot(res2$r, type="upper", order="hclust", 
         p.mat = res2$P, sig.level = 0.01, insig = "blank")

# https://github.com/slzhao/KEGGprofile
options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))) 
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("KEGGprofile")
BiocManager::install("KEGGREST")
BiocManager::install("clusterProfiler")

conda create -n R4.0
source activate R4.0
conda install -c conda-forge r-base
# R version 4.0.3 (2020-10-10)
安装DESeq2
# 1.判断是否有BiocManager包，若不存在则安装
options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))) #设置清华镜像，加速下载
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  BiocManager::install('annotate')
  BiocManager::install('GenomicRanges')
  BiocManager::install('genefilter')
  BiocManager::install('geneplotter')
  BiocManager::install('SummarizedExperiment')
if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install('DESeq2')  #通过BiocManager安装DESeq2 
BiocManager::install("pathview")
BiocManager::install("edgeR")
install.packages("vegan",repos="http://mirrors.ustc.edu.cn/CRAN/")
BiocManager::install("phyloseq")
BiocManager::install("biomformat")
install.packages("tidyverse",repos="http://mirrors.ustc.edu.cn/CRAN/")
install.packages("ggpubr",repos="http://mirrors.ustc.edu.cn/CRAN/")
install.packages("optparse",repos="http://mirrors.ustc.edu.cn/CRAN/") 
install.packages("ggbiplot",repos="http://mirrors.ustc.edu.cn/CRAN/") 

library(DESeq2) #加载library

pathview: https://learn.gencore.bio.nyu.edu/rna-seq-analysis/over-representation-analysis/



ref: https://anaconda.org/bioconda/kraken2
KrakenTools: https://github.com/jenniferlu717/KrakenTools
官网：https://ccb.jhu.edu/software/bracken/
需要的文件下载地址： https://benlangmead.github.io/aws-indexes/k2
Bracken: https://github.com/jenniferlu717/Bracken
Kraken2+Bracken: https://www.jianshu.com/p/dd8182e861cb
Kraken: https://indexofire.github.io/pathongs/C12_Metagenomics-Analysis/01_kraken/
conda create -n kraken2 -c bioconda kraken2
conda remove -n kraken2 bracken
source activate 
conda install -c bioconda braker2
v2.1.2
ln -s /work/workspace/zhurj/software/Bracken-master/bracken /work/workspace/zhurj/bin/miniconda3/envs/kraken2/bin
ln -s /work/workspace/zhurj/software/Bracken-master/bracken-build /work/workspace/zhurj/bin/miniconda3/envs/kraken2/bin
cd /work/workspace/zhurj/reference/kraken2/20200929
tar zxvf minikraken2_v1_8GB_201904_UPDATE.tgz
# Bracken数据库构建(2020.09.29问)
bracken-build -d /work/workspace/zhurj/reference/kraken2/20200929/minikraken2_v1_8GB -t 8 -k 35 -l 150 -x /work/workspace/zhurj/bin/miniconda3/envs/kraken2/bin
#结束后会在kraken数据库路径下生成database150mers.kmer_distrib
db: /work/workspace/zhurj/reference/kraken2/20200929/current

conda activate py2
# Microbial association network construction tutorial http://psbweb05.psb.ugent.be/conet/microbialnetworks/sparcc.php
conda install -c bioconda sparcc
SparCC.py -h
############################################# 
file: /work/workspace/zhurj/bin/miniconda3/envs/py2/bin/core_methods.py
replace line 293 with 295,296,297

SpiecEasi: https://github.com/zdk123/SpiecEasi
source activate R3.6
R
library(devtools)
install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(devtools)
install_github("hallucigenia-sparsa/seqtime") 
library(seqtime)

MEGAN 
MEGAN6 Download Page： https://software-ab.informatik.uni-tuebingen.de/download/megan6/welcome.html
宏基因组注释和可视化神器MEGAN入门： https://blog.csdn.net/woodcorpse/article/details/104381643
linux 安装位置： /work/workspace/zhurj/instal/megan/tools
nr: https://blog.csdn.net/woodcorpse/article/details/104381643
如何利用NR库快速进行物种鉴定: ttps://www.jianshu.com/p/45fdb5cf930a

download： ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz (82.3G)
time unpigz -k -p 16 nr.gz
time diamond makedb --in nr -d nr -p 32

# proteinortho_curves https://github.com/isabelschober/proteinortho_curves
/work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript /work/workspace/zhurj/software/proteinortho_curves-master/corePan_curve.r -p /work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/mergeAbundance -o /work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/corePan_curve

/work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript /work/workspace/zhurj/software/proteinortho_curves-master/corePan_curve.r -p /work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/mergeAbundance -i 20 -o /work/workspace/zhurj/project/1_metadata/mouse22/genecatalog/out/corePan_curve


# 什么！！！超70G的NT数据库文件一个小时搞定？: https://cloud.tencent.com/developer/article/1654480
# https://ftp.ncbi.nih.gov/blast/db/FASTA/
# https://ftp.ncbi.nih.gov/blast/db/v5/ # nr 数据库
conda create -n download
conda activate download
conda install -y -c hcc aspera-cli
# conda install -y -c bioconda sra-tools
# download nr https://ftp.ncbi.nih.gov/blast/db/FASTA/
# https://ftp.ncbi.nih.gov/blast/db/v5
# so so so cool
cd /work/workspace/zhurj/database/nr
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/FASTA/nr.gz ./
cd /work/workspace/zhurj/database/nr/db
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.00.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.01.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.02.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.03.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.04.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.05.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.06.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.07.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.08.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.09.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.10.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.11.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.12.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.13.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.14.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.15.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.16.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.17.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.18.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.19.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.20.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.21.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.22.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.23.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.24.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.25.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.26.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.27.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.28.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.29.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.30.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.31.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.32.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.33.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.34.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.35.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.36.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.37.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/v5/nr.38.tar.gz ./

ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/pub/taxonomy/taxdump.tar.gz ./
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/pub/taxonomy/accession2taxid/prot.accession2taxid.gz ./
gzip -dc prot.accession2taxid.gz > prot.accession2taxid
# mkdir ~/.taxonkit # 已存在，under env python3.6， taxonkit
tar zxf taxdump.tar.gz -C ~/.taxonkit
# 其主要有效文件有两个：
# names.dmp 记录物种名及其分类编号
# nodes.dmp 记录分类编号的节点信息
# 查看~/.taxonkit/names.dmp文件，使用关键词检索得到目标类的分类编号，例如：
# fungi 4751             # grep -P "\|\s+[fF]ungi\w*\s*\|" ~/.taxonkit/names.dmp
# plants 3193            # grep -P "\|\s+[pP]lant\w*\s*\|" ~/.taxonkit/names.dmp
# animals 33208          # grep -P "\|\s+[aA]nimal\w*\s*\|" ~/.taxonkit/names.dmp
# 提取古菌(2157)、细菌(2)和病毒(10239)这几个大类下的所有物种编号。
source deactivate download
source activate python3.6
taxonkit list -j 8 --ids 4751,2,2157,10239 > sub.meta.list
taxonkit list -j 8 -r --ids 4751,2,2157,10239 > sub.meta.rank.list
# 去除文件行前空格
sed -i 's/^[ ]*//g'  sub.meta.rank.list
cat sub.meta.rank.list | grep species | awk '{print $1}' > sub.meta.species.list
# blast db
# makeblastdb -in nr_meta.fasta -dbtype prot -title nr_meta -parse_seqids -out nr_meta_`date +%Y%m%d` -logfile nr_meta_`date +%Y%m%d`.log
#blastdbcmd -taxids 2696393 -db /work/workspace/zhurj/database/nr/db/nr  -dbtype prot -outfmt "%f" -out test.fasta
blastdbcmd -taxidlist /work/workspace/zhurj/database/nr/sub.meta.species.list -db /work/workspace/zhurj/database/nr/db/nr -dbtype prot -outfmt "%f" -out /work/workspace/zhurj/database/nr/microbe/nr_microbe.fasta
seqkit rmdup -s /work/workspace/zhurj/database/nr/microbe/nr_microbe.fasta > /work/workspace/zhurj/database/nr/microbe/nr_microbe_dup.fasta

cd /work/workspace/zhurj/database/p
time srun -o nr_micro_20201111.out -e nr_micro_20201111.err -N 1 -c 20 -p slurm256 -w mnclient02 bash nr_micro_20201111.sh &
nr_micro_20201111.sh
# 数据更新时间：2020.10.12
```
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
diamond makedb --in /work/workspace/zhurj/database/nr/microbe/nr_microbe.fasta -d /work/workspace/zhurj/database/nr/diamond/nr_microbe
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6


```

cazy
# download url: http://bcb.unl.edu/dbCAN2/download/
source activate python3.6
cd /work/workspace/zhurj/database/p
time srun -o cazy.diamond.07312020.out -e cazy.diamond.07312020.err -N 1 -c 20 -p slurm256 -w mnclient02 bash cazy.diamond.07312020.sh &

cazy.diamond.07312020.sh
```
source activate /work/workspace/zhurj/bin/miniconda3/envs/python3.6
diamond makedb --in /work/workspace/zhurj/database/cazy/20201111/CAZyDB.07312020.fa -d /work/workspace/zhurj/database/cazy/diamond/CAZyDB
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/python3.6

```

cog
https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/

/work/workspace/zhurj/instal/megan/jre/lib/libawt_xawt.so: libXtst.so.6

source activate metaAna
conda install -c bioconda krona
# krona reference
# https://ftp.ncbi.nih.gov/pub/taxonomy
Krona installed.  You still need to manually update the taxonomy
databases before Krona can generate taxonomic reports.  The update
script is ktUpdateTaxonomy.sh.  The default location for storing
taxonomic databases is /work/workspace/zhurj/bin/miniconda3/envs/metaAna/opt/krona/taxonomy

If you would like the taxonomic data stored elsewhere, simply replace
this directory with a symlink.  For example:

rm -rf /work/workspace/zhurj/bin/miniconda3/envs/metaAna/opt/krona/taxonomy
mkdir /path/on/big/disk/taxonomy
ln -s /path/on/big/disk/taxonomy /work/workspace/zhurj/bin/miniconda3/envs/metaAna/opt/krona/taxonomy

# bash /work/workspace/zhurj/bin/miniconda3/envs/metaAna/opt/krona/updateTaxonomy.sh

# krona test
source activate metaAna
ktImportText /work/workspace/zhurj/project/1_metadata/mouse22/krona/test/text.txt

source activate download
cd /work/workspace/zhurj/bin/miniconda3/envs/metaAna/opt/krona/taxonomy
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/pub/taxonomy/taxdump.tar.gz ./

source activate /work/workspace/zhurj/bin/miniconda3/envs/download
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/pub/taxonomy/accession2taxid/dead_nucl.accession2taxid.gz /work/workspace/zhurj/bin/miniconda3/envs/metaAna/opt/krona/taxonomy/accession2taxid/
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/pub/taxonomy/accession2taxid/dead_nucl.accession2taxid.gz.md5 /work/workspace/zhurj/bin/miniconda3/envs/metaAna/opt/krona/taxonomy/accession2taxid/
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/pub/taxonomy/accession2taxid/dead_prot.accession2taxid.gz /work/workspace/zhurj/bin/miniconda3/envs/metaAna/opt/krona/taxonomy/accession2taxid/
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/pub/taxonomy/accession2taxid/dead_prot.accession2taxid.gz.md5 /work/workspace/zhurj/bin/miniconda3/envs/metaAna/opt/krona/taxonomy/accession2taxid/
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/pub/taxonomy/accession2taxid/dead_wgs.accession2taxid.gz /work/workspace/zhurj/bin/miniconda3/envs/metaAna/opt/krona/taxonomy/accession2taxid/
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/pub/taxonomy/accession2taxid/dead_wgs.accession2taxid.gz.md5 /work/workspace/zhurj/bin/miniconda3/envs/metaAna/opt/krona/taxonomy/accession2taxid/
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz /work/workspace/zhurj/bin/miniconda3/envs/metaAna/opt/krona/taxonomy/accession2taxid/
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz.md5 /work/workspace/zhurj/bin/miniconda3/envs/metaAna/opt/krona/taxonomy/accession2taxid/
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz /work/workspace/zhurj/bin/miniconda3/envs/metaAna/opt/krona/taxonomy/accession2taxid/
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz.md5 /work/workspace/zhurj/bin/miniconda3/envs/metaAna/opt/krona/taxonomy/accession2taxid/
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/pub/taxonomy/accession2taxid/pdb.accession2taxid.gz /work/workspace/zhurj/bin/miniconda3/envs/metaAna/opt/krona/taxonomy/accession2taxid/
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/pub/taxonomy/accession2taxid/pdb.accession2taxid.gz.md5 /work/workspace/zhurj/bin/miniconda3/envs/metaAna/opt/krona/taxonomy/accession2taxid/
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/pub/taxonomy/accession2taxid/prot.accession2taxid.gz /work/workspace/zhurj/bin/miniconda3/envs/metaAna/opt/krona/taxonomy/accession2taxid/
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/pub/taxonomy/accession2taxid/prot.accession2taxid.gz.md5 /work/workspace/zhurj/bin/miniconda3/envs/metaAna/opt/krona/taxonomy/accession2taxid/
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/download

cd /work/workspace/zhurj/bin/miniconda3/envs/metaAna/opt/krona/p
srun -o rjdownloadTaxo.out -e rjdownloadTaxo.err -N 1 -c 20 -p slurm256 -w mnclient02 bash rjdownloadTaxo.sh &

source activate metaAna
# FastSpar | 用更快的 SparCC 进行微生物组相关性分析 https://cloud.tencent.com/developer/article/1727777
conda create -n SCNIC python=3 scnic
conda install -c bioconda -c conda-forge fastspar

uniprot2ko
ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/
KEGG FTP: https://www.kegg.jp/kegg/download/
# version 20201007
source activate download
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.uniprot.org:/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz /work/workspace/zhurj/database/kegg/acc/
ascp -v -k 1 -T -l 200m -i /work/workspace/zhurj/bin/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh anonftp@ftp.ebi.ac.uk:/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz /work/workspace/zhurj/database/kegg/acc/
source deactivate download

2020.11.30
VFDB 数据下载地址 http://www.mgc.ac.cn/VFs/download.htm
VFGxxxx 的注释信息 在 DNA sequence of full dataset中都有，文件在服务器上地址/work/workspace/zhurj/database/VFDB

2020.12.23
# 
conda install -c bioconda pyfastx

python3.6
pip install multiqc
# /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/multiqc

assembly software - 192.168.2.74
/mnt/d/home/x/data/M/coding/genomeAssemble

# -----------------------------------database update---------------------------------------------------------
# refseq 16S rRNA
# 2021.02.25
https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/

2021.05.13
# ref: https://anaconda.org/bioconda/fgmp
# FGMP: assessing fungal genome completeness
conda install -c bioconda fgmp
#cpanm -l /work/workspace/zhurj/software/FGMP-master/lib IO::All

# 系统命令
系统允许最大进程数：cat /proc/sys/kernel/pid_max
系统允许最大线程数：cat /proc/sys/kernel/threads-max
单用户允许最大线程数：ulimit -u
单进程允许最大线程数：cat /usr/include/bits/local_lim.h 的PTHREAD_THREADS_MAX
每个线程栈空间大小：ulimit -s
查看系统已用进程数：pstree -p|wc -l
查看一个进程中线程的数量：ps -Lf pid|wc -l
查看一个进程资源消耗情况：top -p pid

# 2021.06.17
blastdb 下载地址 https://ftp.ncbi.nlm.nih.gov/blast/db/

/work/workspace/zhurj/project/6_learning/python/Cookbook
https://onedrive.live.com/edit.aspx?action=editnew&resid=7A2E5A40CB933C63!827&ithint=file%2cxlsx&action=editnew&wdNewAndOpenCt=1624586271094&wdPreviousSession=6279bc17-8a71-40c4-aa74-e18c887e2ff8&wdOrigin=OFFICECOM-WEB.START.NEW

# 2021.06.29
GTDBTK v1.5.1
conda create -n GTDBr95
conda install -c bioconda gtdbtk
conda env remove -n GTDBr95
# Note, version gtdbtk v0.3.2

# 2021.07.08
# https://anaconda.org/conda-forge/falcon
# 使用FALCON对三代测序数据进行基因组组装： http://www.chenlianfu.com/?p=2755
# falcon组装及polish: https://www.jianshu.com/p/2c4a094c8675
conda create -n falcon python=2.7.9
source activate falcon
conda install -c bioconda pb-falcon
conda install -c bioconda dazz_db
conda install -c bioconda daligner
# deeptools 3.5.1
conda install -c bioconda deeptools
conda install -c bioconda bam2fastx
conda install -c bioconda circlator
conda install -c bioconda pbgcpp
conda install -c bioconda pbmm2  # doesn't work
conda install -c bioconda pb-assembly
/work/workspace/zhurj/bin/miniconda3/bin/samtools # samtools doesn't work


# canu: https://zhuanlan.zhihu.com/p/362068540

conda remove -n canu  --all
source activate canu
conda install -c bioconda canu
# bam2fastx : https://github.com/PacificBiosciences/bam2fastx
conda install -c bioconda bam2fastx
conda install -c bioconda pbmm2
# circlator 环化基因组
conda install -c bioconda circlator
conda install -c bioconda sambamba
# 2021.08.02
# 安装 https://anaconda.org/bioconda/repeatmasker
conda install -c bioconda repeatmasker
conda install -c bioconda trf
# 软件下载 从 http://www.repeatmasker.org/RepeatMasker/ 下载 RepeatMasker-4.1.2-p1.tar.gz
tar xf RepeatMasker-4.1.2-p1.tar.gz


2021.07.09
root 更新git
1. download “git-2.30.0.tar.gz” from https://cdn.kernel.org/pub/software/scm/git/
2. 
0）安装依赖软件
yum install -y curl-devel expat-devel gettext-devel openssl-devel zlib-devel gcc perl-ExtUtils-MakeMaker

1）卸载系统自带的底版本git（1.8.3.1）
[root@uatjenkins01 ~]# git --version
[root@uatjenkins01 ~]# yum remove git
 
2）编译安装最新的git版本
# put git-2.30.0.tar.gz under /usr/local/src/
[root@ ~]# cd /usr/local/src/
[root@ src]# tar -vxf git-2.30.0.tar.xz
[root@ src]# cd git-2.30.0
[root@ git-2.15.1]# make prefix=/usr/local/git all
[root@ git-2.15.1]# make prefix=/usr/local/git install
[root@ git-2.15.1]# echo "export PATH=$PATH:/usr/local/git/bin" >> /etc/profile
[root@ git-2.15.1]# source /etc/profile

[root@ ~]# git --version

2021.07.12
pacbio 用二代测序软件下载安装
# 下载文件存放路径： /work/workspace/zhurj/software
wget https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar
java -Xmx16G -jar pilon-1.24.jar

2021.07.21
网络设置查询： https://www.cnblogs.com/linhaifeng/articles/5937962.html#_label2
# windows 系统
打开cmd执行命令：ipconfig /all
全国通用DNS地址（国内用户推荐使用，速度较快！）
首先DNS服务器地址添：114.114.114.114  (位于北京人民英雄纪念碑）
备用DNS服务器地址添：114.114.115.115
全球通用DNS地址（此DNS地址为谷歌服务器的）
首选DNS服务器地址添：8.8.8.8
备用DNS服务器地址添：8.8.4.4

查看本地dns缓存命令：ipconfig /displaydns
清除本地dns缓存命令：ipconfig /flushdns


# 2021.07.21
download Xbin.tgz from http://ftp.xfree86.org/pub/XFree86/4.8.0/binaries/Linux-x86_64-glibc23/
cd /work/workspace/zhurj/software
# put Xbin.tgz under /work/workspace/zhurj/software
tar zxf Xbin.tgz -C /work/workspace/zhurj/bin/Xbin
cp /work/workspace/zhurj/bin/Xbin/lib64/libXtst.so.6 /work/workspace/zhurj/instal/megan/jre/lib
cp /work/workspace/zhurj/bin/Xbin/lib64/libXi.so.6 /work/workspace/zhurj/instal/megan/jre/lib
/work/workspace/zhurj/instal/megan/jre/lib

# 2021.07.23 megan reference 更新
# https://software-ab.informatik.uni-tuebingen.de/download/megan6/welcome.html
download megan-map-Jan2021.db.zip


env: python3.6
下载sqlite-tools-linux-x86-3360000.zip
ln -s /work/workspace/zhurj/software/sqlite-tools-linux-x86-3360000/sqlite3 /work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/sqlite3


download prot.accession2taxid.gz from ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz

2021.07.30
2.9T    ./zhurj
2.0T    ./liangyj
8.5T    ./liangzj
3.5T    ./zhangcc

# 2021.07.30
# python3.6
pip3 install pylint
pip3 install PasteScript
conda install -c bioconda trnascan-se
conda install -c bioconda 
# v1.1.4

# canu
conda install -c bioconda trnascan-se
conda install -c bioconda barrnap

# 查看服务器运行星狂
1. 进入root
2. ssh mnclient02 # 选择节点
3. top ， c （可以显示详细的命令）

# docker
126.com
shitou6
8docker

# 2021.09.24
# 9.9M gene download from https://db.cngb.org/microbiome/genecatalog/genecatalog_human/
IGC.fa.gz: Integrated non-redundant gene catalog (IGC, nucleotide sequences, fasta)
IGC.pep.gz: Integrated non-redundant gene catalog (IGC, amino acid sequences, fasta)
IGC.annotation_OF.summary.gz: IGC annotation and occurrence frequency summary table
linux dir： /work/workspace/zhurj/database/genecatalog

filezilla 下载  Unified Human Gastrointestinal Protein (UHGP) catalog
参考文献： 
Published: 20 July 2020
A unified catalog of 204,938 reference genomes from the human gut microbiome

下载地址：ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/uhgp_catalogue/
filezilla INPUT
HOST: ftp://ftp.ebi.ac.uk/
U: anonymous
W: 
P: 21
downloaded file: uhgp-90.tar.gz
linux dir： /work/workspace/zhurj/database/genecatalog

