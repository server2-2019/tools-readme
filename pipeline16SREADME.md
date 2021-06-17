# tools-readme
Record readme of often used software
2019.Jul02-14:05

# lefse : /mnt/d/home/ruijuan/.conda/envs/lefse

1. sign in github
account: zrjlyq@126.com
"New repository" => "shitou6-6/pipeline16S" => "Description: 16S rRNA gene sequencing analysis" =>  Private => "Initialize this repository with a README"

2. cd  /work/workspace/ruijuan/script/git
git clone git@github.com:shitou-6/pipeline16S.git

3. pipeline
(1) read filter -- 随机提取一定数据量的read进行后续分析，目的保证有足够数据的前提下，使每个样品read数相当
perl /work/workspace/ruijuan/script/git/pipeline16S/sampleRead.pl \
-i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/input/samplefile.txt  \
-o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/read30000 \
-n 30000

manifest creation
find /work/rawdata/run/guangzhou/2019/07/20190703/run00034/fastq -type l | grep gz |sort |awk -F "[/.]" 'BEGIN{print "sample-id,absolute-filepath,direction"}{if($(NF-2) == 1){print $(NF-3)","$0",forward"}else{print $(NF-3)","$0",reverse"}}' > /work/workspace/ruijuan/project/3_16S/anhuilung20190805/input/manifest10.txt

(2) qiime 
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/input/manifest \
  --output-path /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/reads.qza \
  --input-format PairedEndFastqManifestPhred33

dada2去燥,合并双端序列
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/reads.qza \
  --p-trim-left-f 19 \
  --p-trim-left-r 20 \
  --p-trunc-len-f 180 \
  --p-trunc-len-r 180 \
  --o-table /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/table.qza \
  --o-representative-sequences /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/rep-seqs.qza \
  --o-denoising-stats /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/denoising-stats.qza \
  --p-n-threads 20

  # 可视化denoising stats（qzv文件可在线查看）
  qiime metadata tabulate \
  --m-input-file /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/denoising-stats.qza \
  --o-visualization /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qzv/denoising-stats.qzv

  # filter FeatureTable
  qiime feature-table filter-features \
  --i-table /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/table.qza \
  --p-min-frequency 2 \
  --m-metadata-file /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/input/metadata \
  --p-exclude-ids \
  --o-filtered-table /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/table_filter.qza \
  --m-metadata-file /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/input/metadata

  # tree create for diversity analysis
  qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/rep-seqs.qza \
  --o-alignment /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/aligned-rep-seqs.qza \
  --o-masked-alignment /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/masked-aligned-rep-seqs.qza \
  --o-tree /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/unrooted-tree.qza \
  --o-rooted-tree /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/rooted-tree.qza \
  --p-n-threads 20

  ### ======================
  --output-dir /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza \
  ## output-dir 替代一下四行
  --o-alignment /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/aligned-rep-seqs.qza \
  --o-masked-alignment /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/masked-aligned-rep-seqs.qza \
  --o-tree /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/unrooted-tree.qza \
  --o-rooted-tree /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/rooted-tree.qza
  ###

  # alpha rarefaction 
  qiime diversity alpha-rarefaction \
  --i-table /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/table_filter.qza \
  --p-metrics  observed_otus --p-metrics shannon  --p-metrics faith_pd \
  --p-metrics goods_coverage --p-metrics chao1 --p-metrics simpson \
  --p-max-depth 20000 \
  --m-metadata-file /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/input/metadata \
  --i-phylogeny /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/rooted-tree.qza \
  --o-visualization /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qzv/alpha_rarefaction.qzv

  # diversity analysis
  qiime diversity core-metrics-phylogenetic \
  --i-phylogeny /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/rooted-tree.qza \
  --i-table /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/table_filter.qza \
  --p-sampling-depth 19508 \
  --m-metadata-file /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/input/metadata \
  --output-dir /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/core-metrics-results


  # training classifier
  qiime feature-classifier extract-reads \
  --i-sequences reference/99_otus.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \      #515 F引物
  --p-r-primer GGACTACHVGGGTWTCTAAT \  #806 R引物
  --o-reads /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/ref/ref-seqs.qza

  # Train the classifier（训练分类器）
  # 基于筛选的指定区段，生成实验特异的分类器
  qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/ref/ref-seqs.qza \
  --i-reference-taxonomy reference/ref-taxonomy.qza \
  --o-classifier /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/ref/classifier.qza

  # taxa classifier
  qiime feature-classifier classify-sklearn \
  --i-classifier /work/workspace/ruijuan/script/git/pipeline16S/reference/Greengenes_13_8_99%_OTUs_341F-805R_classifier.qza \
  --i-reads /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/rep-seqs.qza \
  --o-classification /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/taxonomy.qza

  # 结果可视化
  qiime metadata tabulate \
  --m-input-file /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/taxonomy.qza \
  --o-visualization /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qzv/taxonomy.qzv

  # 物种分类条形图
  qiime taxa barplot \
  --i-table /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/table_filter.qza \
  --i-taxonomy /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/taxonomy.qza \
  --m-metadata-file /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/input/metadata \
  --o-visualization /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qzv/taxa-bar-plots.qzv

  qiime tools export \
  --input-path /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qzv/taxa-bar-plots.qzv \
  --output-path /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qzv/taxa

  # acom analysis
  #按照 leve2 水平分析
  qiime taxa collapse \
  --i-table /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/table_filter.qza \
  --i-taxonomy /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/taxonomy.qza \
  --p-level 2 \
  --o-collapsed-table /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/acom/qza/taxa_table_l2.qza

  qiime composition add-pseudocount \
  --i-table /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/acom/qza/taxa_table_l2.qza \
  --o-composition-table /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/acom/qza/comp_table_l2.qza

  qiime composition ancom \
  --i-table /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/acom/qza/comp_table_l2.qza \
  --m-metadata-file /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/input/metadata \
  --m-metadata-column Condition \
  --p-transform-function clr \
  --o-visualization /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/acom/qzv/ancom_l2_clr_Condition.qzv

  # level 6
  qiime taxa collapse \
  --i-table /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/filter_table.qza \
  --i-taxonomy /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/acom/qza/taxa_table_l6.qza

  qiime composition add-pseudocount \
  --i-table /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/acom/qza/taxa_table_l6.qza \
  --o-composition-table /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/acom/qza/comp_table_l6.qza 

  qiime composition ancom \
  --i-table /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/acom/qza/comp_table_l6.qza \
  --m-metadata-file /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/input/metadata \
  --m-metadata-column Condition \
  --p-transform-function clr \
  --o-visualization /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/acom/qzv/ancom_l6_clr_Condition.qzv

  # level 7
  qiime taxa collapse \
  --i-table /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/filter_table.qza \
  --i-taxonomy /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/acom/qza/taxa_table_l7.qza

  qiime composition add-pseudocount \
  --i-table /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/acom/qza/taxa_table_l7.qza \
  --o-composition-table /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/acom/qza/comp_table_l7.qza

  qiime composition ancom \
  --i-table /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/acom/qza/comp_table_l7.qza \
  --m-metadata-file /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/input/metadata \
  --m-metadata-column Condition \
  --p-transform-function clr \
  --o-visualization /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/acom/qzv/ancom_l7_clr_Condition.qzv

  # picrust2
  qiime tools export \
  --input-path /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/rep-seqs.qza \
  --output-path /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/input/ref-seqs

  qiime tools export \
  --input-path /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qza/filter_table.qza \
  --output-path /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/input/table


  source deactivate 16S-2019.04
  source activate picrust2

  place_seqs.py \
  -s /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/input/ref-seqs/dna-sequences.fasta \
  -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/out.tre \
  -p 20
 
  hsp.py -i 16S \
  -t /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/out.tre \
  -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/marker_predicted_and_nsti.tsv.gz \
  -p 20

  hsp.py -i COG \
  -t /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/out.tre \
  -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/COG_predicted.tsv.gz \
  -p 20

  metagenome_pipeline.py -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/input/table/feature-table.biom \
  -m /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/marker_predicted_and_nsti.tsv.gz \
  -f /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/COG_predicted.tsv.gz \
  -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/COG \
  --strat_out

  add_descriptions.py -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/COG/pred_metagenome_unstrat.tsv \
  -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/COG/pred_metagenome_unstrat_descrip.tsv \
  -m COG

  hsp.py -i EC \
  -t /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/out.tre \
  -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/EC_predicted.tsv.gz \
  -p 20

  metagenome_pipeline.py -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/input/table/feature-table.biom \
  -m /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/marker_predicted_and_nsti.tsv.gz \
  -f /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/EC_predicted.tsv.gz \
  -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/EC \
  --strat_out

  add_descriptions.py -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/EC/pred_metagenome_unstrat.tsv \
  -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/EC/pred_metagenome_unstrat_descrip.tsv \
  -m EC

  hsp.py -i KO -t /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/out.tre \
  -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/KO_predicted.tsv.gz \
  -p 20

  metagenome_pipeline.py -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/input/table/feature-table.biom \
  -m /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/marker_predicted_and_nsti.tsv.gz \
  -f /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/KO_predicted.tsv.gz \
  -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/KO \
  --strat_out

  add_descriptions.py -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/KO/pred_metagenome_unstrat.tsv \
  -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/KO/pred_metagenome_unstrat_descrip.tsv \
  -m KO

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

 #--------------------------------------------------------------------------------------------------------------------------------------------------------------
  
source deactivate picrust2
source activate lefse

Rscript /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/lefsein_taxa.r /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qzv/taxa/level-2.csv Condition Name index 3 /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/lefse/taxa/level2/lefse.in
lefse-format_input.py /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/lefse/taxa/level2/lefse.in /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/lefse/taxa/level2/lefse_format.in -c 1 -s 2 -u 3 -f r -o 1000000
run_lefse.py /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/lefse/taxa/level2/lefse_format.in /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/lefse/taxa/level2/lefse_run.res -a 0.05 -w 0.05 -l 2 -y 0
lefse-plot_res.py /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/lefse/taxa/level2/lefse_run.res /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/lefse/taxa/level2/lefse_run.png --feature_font_size 8 --width 10 --format png
lefse-plot_cladogram.py /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/lefse/taxa/level2/lefse_run.res /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/lefse/taxa/level2/lefse_cladogram.png --format png
Rscript /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/lefsein_taxa.r /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qzv/taxa/level-6.csv Condition Name index 3 /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/lefse/taxa/level6/lefse.in
lefse-format_input.py /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/lefse/taxa/level6/lefse.in /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/lefse/taxa/level6/lefse_format.in -c 1 -s 2 -u 3 -f r -o 1000000
run_lefse.py /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/lefse/taxa/level6/lefse_format.in /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/lefse/taxa/level6/lefse_run.res -a 0.05 -w 0.05 -l 2 -y 0
lefse-plot_res.py /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/lefse/taxa/level6/lefse_run.res /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/lefse/taxa/level6/lefse_run.png --feature_font_size 8 --width 10 --format png
lefse-plot_cladogram.py /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/lefse/taxa/level6/lefse_run.res /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/lefse/taxa/level6/lefse_cladogram.png --format png
Rscript /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/lefsein_taxa.r /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/qzv/taxa/level-7.csv Condition Name index 3 /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/diff/lefse/taxa/level7/lefse.in





$lefse_avalue = ''; $lefse_wvalue = ''; $lefse_ldacutoff = ''; $lefse_strategy = '';
$lefse_graph_font_size --width $lefse_graph_width --format $lefse_oformat
kegg pathway LEfSe level3


lefse-format_input.py /work/workspace/ruijuan/project/3_16S/anhui20190527/all_30000_result/picrust2_out_pipeline/KEGG_pathways_l3/lefse_kegg_l3.tsv /work/workspace/ruijuan/project/3_16S/anhui20190527/all_30000_result/picrust2_out_pipeline/KEGG_pathways_l3/lefse_format.in -c 1 -s 2 -u 3 -f r -o 1000000
run_lefse.py /work/workspace/ruijuan/project/3_16S/anhui20190527/all_30000_result/picrust2_out_pipeline/KEGG_pathways_l3/lefse_format.in /work/workspace/ruijuan/project/3_16S/anhui20190527/all_30000_result/picrust2_out_pipeline/KEGG_pathways_l3/lefse_run.res -a 0.05 -w 0.05 -l 2.0 -y 0
lefse-plot_res.py /work/workspace/ruijuan/project/3_16S/anhui20190527/all_30000_result/picrust2_out_pipeline/KEGG_pathways_l3/lefse_run.res /work/workspace/ruijuan/project/3_16S/anhui20190527/all_30000_result/picrust2_out_pipeline/KEGG_pathways_l3/lefse_run.png --feature_font_size 8 --width 10 --format png 
lefse-plot_cladogram.py /work/workspace/ruijuan/project/3_16S/anhui20190527/all_30000_result/picrust2_out_pipeline/KEGG_pathways_l3/lefse_run.res /work/workspace/ruijuan/project/3_16S/anhui20190527/all_30000_result/picrust2_out_pipeline/KEGG_pathways_l3/lefse_cladogram.png --format png

our (@acom_levels, @lefse_items,@lefse_taxa_levels);
 66 our (%acomLevels, %lefseItems, %lefseTaxaLevels);
 67 our (%acom_levels_defines, %lefse_items_defines, %lefse_taxa_levels_defines);
 %acom_fun_defines




perl /work/workspace/ruijuan/script/git/pipeline16S/qiime2pro.pl -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000 --manifest /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/input/manifest --meta /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/input/metadata --trim-left-f 19 --trim-left-r 20 --trunc-len-f 250 --trunc-len-r 200

perl /work/workspace/ruijuan/script/git/pipeline16S/qiime2pro.pl -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000 --manifest /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/input/manifest --meta /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/input/metadata --trim-left-f 19 --trim-left-r 20 --trunc-len-f 250 --trunc-len-r 200 -ancom --ancom-column Condition --ancom-levels 2 --ancom-levels 6 --ancom-levels 7 -picrust2 -lefse --lefse-class Condition --lefse-subclass Name

Rscript
install.packages("/mnt/d/home/ruijuan/workflows/software/R/docopt_0.6.1.tar.gz") # USAGE
install.packages("/mnt/d/home/ruijuan/workflows/software/R/getopt_1.20.3.tar.gz") # 传参数

Rscript lefsein_taxa.r '/work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/tmp/taxa/level-7.csv' 'Condition' 'Name' 'index' 3 '/work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/tmp/lefse.in'
~/.conda/envs/16S-2019.04/bin/Rscript lefsein_taxa.r -i '/work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/tmp/taxa/level-7.csv' --class 'Condition' --subclass 'Name' --subject 'index' --cutoff 3 -o '/work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/tmp/lefse.in'
~/.conda/envs/16S-2019.04/bin/Rscript lefsein_taxa.r -i '/work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/tmp/taxa/level-7.csv' --class 'Condition' --subject 'index' --cutoff 3 -o '/work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/tmp/lefse.in'


library(ggplot2)
library(vegan)
library(ggpubr)
library('ggbiplot')
library(reshape2)
library('RColorBrewer')
library('ggsci')
library(dplyr)
library(docopt)
library(getopt)

LEfSe
# http://www.bioinfo-scrounger.com/archives/405 dplyr 数据处理
# http://www.bioinfo-scrounger.com/archives/405
df <- read.csv("/work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/tmp/taxa/level-7.csv",header = T,sep=",")
dff = df[,grep("index|Condition|Name|k__", names(df), value = TRUE)]
new_df <- dff %>% select("Condition","Name",everything())
myFun<-function(x){ 
return (length(which(x>0)))
}
count <- apply(new_df[,4:ncol(new_df)],2,myFun)
fcn <- names(which(count>3))
ndf <- new_df %>% select("Condition","Name","index",fcn)
rcount <- apply(ndf[,4:ncol(ndf)],1,sum)
rnames <- which(rcount>0)
ndff <- ndf[rnames,]
fdata <- as.data.frame(t(ndff))
rnames <- rownames(fdata)
rnames <- gsub('\\.','\\|',rnames)
rnames <- gsub('\\.','\\|',rnames)
rnames <- gsub('\\|f__\\|g__\\|s__$','\\_fgs',rnames)
rnames <- gsub('\\|g__\\|s__$','\\_gs',rnames)
rnames <- gsub('\\|s__$','\\_s',rnames)
rnames <- gsub('\\|__\\|__\\|__\\|__\\|__$','\\_cofgsUN',rnames)
rnames <- gsub('\\|__\\|__\\|__\\|__$','\\_ofgsUN',rnames)
rnames <- gsub('\\|__\\|__\\|__$','\\_fgsUN',rnames)
rnames <- gsub('\\|__\\|__$','\\_gsUN',rnames)
rnames <- gsub('\\|__$','\\_sUN',rnames)
rnames <- gsub('\\|\\|','\\|',rnames)
rnames <- gsub('p__\\|','p__',rnames)
write.table(fdata,"/work/workspace/ruijuan/script/git/pipeline16S/test/test.txt",quote=FALSE,sep="\t",col.names=FALSE,row.names=rnames)

~/.conda/envs/16S-2019.04/bin/Rscript lefsein_taxa.r -i /work/workspace/ruijuan/project/3_16S/anhui20190527/all_30000_result_190719/qzv/taxa/level-7.csv --class Group --subject index --cutoff 3 -o /work/workspace/ruijuan/project/3_16S/anhui20190527/all_30000_result_190719/qzv/taxa/diff/level7/lefse.in
lefse-format_input.py /work/workspace/ruijuan/project/3_16S/anhui20190527/all_30000_result_190719/qzv/taxa/diff/level7/lefse.in /work/workspace/ruijuan/project/3_16S/anhui20190527/all_30000_result_190719/qzv/taxa/diff/level7/lefse_format.in -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/ruijuan/project/3_16S/anhui20190527/all_30000_result_190719/qzv/taxa/diff/level7/lefse_format.in /work/workspace/ruijuan/project/3_16S/anhui20190527/all_30000_result_190719/qzv/taxa/diff/level7/lefse_run.res -a 0.05 -w 0.05 -l 2.0 -y 0
lefse-plot_res.py /work/workspace/ruijuan/project/3_16S/anhui20190527/all_30000_result_190719/qzv/taxa/diff/level7/lefse_run.res /work/workspace/ruijuan/project/3_16S/anhui20190527/all_30000_result_190719/qzv/taxa/diff/level7/lefse_run.png --feature_font_size 8 --width 10 --format png 
lefse-plot_cladogram.py /work/workspace/ruijuan/project/3_16S/anhui20190527/all_30000_result_190719/qzv/taxa/diff/level7/lefse_run.res /work/workspace/ruijuan/project/3_16S/anhui20190527/all_30000_result_190719/qzv/taxa/diff/level7/lefse_cladogram.png --format png

~/.conda/envs/16S-2019.04/bin/Rscript lefsein_taxa.r -i /work/workspace/ruijuan/project/3_16S/anhuilung20190805/result_30000/qzv/taxa/level-2.csv --class Condition --subject index --cutoff 3 -o /work/workspace/ruijuan/project/3_16S/anhuilung20190805/result_30000/diff/lefse/taxa/level2/lefse.in

~/.conda/envs/16S-2019.04/bin/Rscript lefsein_picrust.r -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/EC/pred_metagenome_unstrat.tsv -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/EC/lefse.in --metaf /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/input/metadata --class Condition --subject sampleid

rcount <- apply(ndf[,4:ncol(ndf)],1,sum)
data = read.table("/work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/EC/pred_metagenome_unstrat.tsv",sep="\t",header=T)
datat = t(data)
datat[1,1] = 'SampleId'
pdata = dataf[2:nrow(dataf),]
names(pdata) = t(as.vector(dataf[1,]))
meta = read.table("/work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/input/metaNCHFC1",header = T,sep="\t")
cla='Condition'
subject='SampleId'
pattern=paste(sutfbject,cla,sep = "|")
newmeta=meta[,grep(pattern,names(meta),value = TRUE)]
mdata=merge(pdata,newmeta,by='SampleId')
new_mdata <- mdata %>% select(cla,subject,everything())
f_mdata = as.data.frame(t(new_mdata))
write.table(f_mdata,"/work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/EC/lefse.in",quote=FALSE,sep="\t",col.names=FALSE,row.names=TRUE)

library(data.table)

myFun<-function(x){ 
return (length(which(x>0)))
}

#data = read.table("/work/workspace/ruijuan/script/git/pipeline16S/test/frequncy_in.txt",sep="\t",header = T)
# EC
~/.conda/envs/16S-2019.04/bin/Rscript lefsein_picrust.r -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/EC/pred_metagenome_unstrat.tsv -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/EC/lefse.in --metaf /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/input/metadata --class Condition --subject SampleId
# COG
~/.conda/envs/16S-2019.04/bin/Rscript lefsein_picrust.r -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/COG/pred_metagenome_unstrat.tsv -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/COG/lefse.in --metaf /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/input/metadata --class Condition --subject SampleId
# KO
~/.conda/envs/16S-2019.04/bin/Rscript lefsein_picrust.r -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/KO/pred_metagenome_unstrat.tsv -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/KO/lefse.in --metaf /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/input/metadata --class Condition --subject SampleId
# kegg module
# pathway
~/.conda/envs/16S-2019.04/bin/Rscript lefsein_picrust.r -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/kegg_module/path_abun_unstrat.tsv -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/kegg_module/lefse.in --metaf /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/input/metadata --class Condition --subject SampleId
# metacyc
# pathway
~/.conda/envs/16S-2019.04/bin/Rscript lefsein_picrust.r -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/metacyc/path_abun_unstrat.tsv -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/metacyc/lefse.in --metaf /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/input/metadata --class Condition --subject SampleId

# kegg level2
~/.conda/envs/16S-2019.04/bin/Rscript lefsein_picrust.r -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/kegg_path/level2/path_abun_unstrat.tsv -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/kegg_path/level2/lefse.in --metaf /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/input/metadata --class Condition --subject SampleId

# kegg level3
~/.conda/envs/16S-2019.04/bin/Rscript lefsein_picrust.r -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/kegg_path/level3/path_abun_unstrat.tsv -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/kegg_path/level3/lefse.in --metaf /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/input/metadata --class Condition --subject SampleId

# all 
~/.conda/envs/16S-2019.04/bin/Rscript lefsein_picrust.r -i /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/kegg_path/all/path_abun_unstrat.tsv -o /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/kegg_path/all/lefse.in --metaf /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/input/metadata --class Condition --subject SampleId


data = read.table("/work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/EC/pred_metagenome_unstrat.tsv",sep="\t",header = T)
fdata <- round(prop.table(data.matrix(data[1:nrow(data),2:ncol(data)]),2),digits=6)
fdata[is.nan(fdata)] = 0
rownames(fdata) <- data$function.
tdata <- t(fdata)
ftdata1 <- data.frame(tdata,stringsAsFactors = F)
ftdata <- as.data.frame(lapply(ftdata1,as.numeric))

count <- apply(ftdata,2,myFun)
fcn <- names(which(count>2))

ndf <- ftdata %>% select(fcn)
ndf$SampleId <- rownames(ftdata)
ndfo <- ndf %>% select('SampleId',fcn)
rcount <- apply(ndfo[,2:ncol(ndfo)],1,sum)
rnames <- which(rcount>=0)
ndff <- ndfo[rnames,]
#fdata <- as.data.frame(t(ndff))

meta = read.table("/work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/input/metaNCHFC1",header = T,sep="\t")
cla='Condition'
subject='SampleId'
pattern=paste(subject,cla,sep = "|")
newmeta=meta[,grep(pattern,names(meta),value = TRUE)]
mdata=merge(ndff,newmeta,by='SampleId')
new_mdata <- mdata %>% select(cla,subject,everything())
f_mdata = as.data.frame(t(new_mdata))
options(scipen = 200)
write.table(f_mdata,"test.txt",quote=FALSE,sep="\t",col.names=FALSE,row.names=TRUE)
write.table(f_mdata,"/work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/picrust2/EC/lefse.in",quote=FALSE,sep="\t",col.names=FALSE,row.names=TRUE)


perl /work/workspace/ruijuan/script/git/pipeline16S/qiime2pro.pl -o /work/workspace/ruijuan/project/3_16S/anhuiPR20190917/result --manifest /work/workspace/ruijuan/project/3_16S/anhuiPR20190917/input/manifest_30000.txt --meta /work/workspace/ruijuan/project/3_16S/anhuiPR20190917/input/metadata20190719.txt --trim-left-f 0 --trim-left-r 0 --trunc-len-f 151 --trunc-len-r 150 -ancom --ancom-column RIDCtag --ancom-levels 2 --ancom-levels 6 --ancom-levels 7 -picrust2 -lefse --lefse-class RIDCtag --lefse-subclass Name

perl /work/workspace/ruijuan/script/git/pipeline16S/qiime2pro.pl -o /work/workspace/ruijuan/project/3_16S/anhuiPR20190917/result --manifest /work/workspace/ruijuan/project/3_16S/anhuiPR20190917/input/manifest_30000.txt --meta /work/workspace/ruijuan/project/3_16S/anhuiPR20190917/input/metadata20190719.txt --trim-left-f 0 --trim-left-r 0 --trunc-len-f 151 --trunc-len-r 150 -ancom --ancom-column RIDCtag --ancom-levels 2 --ancom-levels 6 --ancom-levels 7 -picrust2 -lefse --lefse-class DIStag 


print PRO "$gb_con\n";

20191227
stage
perl /work/workspace/ruijuan/script/git/pipeline16S/qiime2pro.pl -o /work/workspace/ruijuan/project/3_16S/anhuiPRthreegroup20191227/resultstageg --manifest /work/workspace/ruijuan/project/3_16S/anhuiPRthreegroup20191227/input/manifest_30000.txt --meta /work/workspace/ruijuan/project/3_16S/anhuiPRthreegroup20191227/input/meta3group --trim-left-f 0 --trim-left-r 0 --trunc-len-f 151 --trunc-len-r 150 -ancom --ancom-column Stagetag --ancom-levels 2 --ancom-levels 6 --ancom-levels 7 -picrust2 -lefse --lefse-class Stagetag > /work/workspace/ruijuan/project/3_16S/anhuiPRthreegroup20191227/p/stagerun.sh

/work/workspace/ruijuan/project/3_16S/anhuiPRthreegroup20191227/resultstageg/qza/reads.qza

perl /work/workspace/ruijuan/script/git/pipeline16S/qiime2pro.pl -o /work/workspace/ruijuan/project/3_16S/anhuiPRthreegroup20191227/resultdistantg --manifest /work/workspace/ruijuan/project/3_16S/anhuiPRthreegroup20191227/input/manifest_30000.txt --meta /work/workspace/ruijuan/project/3_16S/anhuiPRthreegroup20191227/input/meta3group --trim-left-f 0 --trim-left-r 0 --trunc-len-f 151 --trunc-len-r 150 -ancom --ancom-column Distantag --ancom-levels 2 --ancom-levels 6 --ancom-levels 7 -picrust2 -lefse --lefse-class Distantag > /work/workspace/ruijuan/project/3_16S/anhuiPRthreegroup20191227/p/Distanrun.sh


perl /work/workspace/ruijuan/script/git/pipeline16S/qiime2pro.pl -o /work/workspace/ruijuan/project/3_16S/anhuiPRthreegroup20191227/resultil6g --manifest /work/workspace/ruijuan/project/3_16S/anhuiPRthreegroup20191227/input/manifest_30000.txt --meta /work/workspace/ruijuan/project/3_16S/anhuiPRthreegroup20191227/input/meta3group --trim-left-f 0 --trim-left-r 0 --trunc-len-f 151 --trunc-len-r 150 -ancom --ancom-column IL6tag --ancom-levels 2 --ancom-levels 6 --ancom-levels 7 -picrust2 -lefse --lefse-class IL6tag > /work/workspace/ruijuan/project/3_16S/anhuiPRthreegroup20191227/p/il6run.sh

full-length 16S cluster
# https://drive5.com/usearch/manual/cmd_sortbylength.html
usearch -sortbylength /work/workspace/zhurj/project/3_16S/clusterfle16S20200526/input/moonbioFL16S202004291127.fasta -fastaout /work/workspace/zhurj/project/3_16S/clusterfle16S20200526/output/moonbioFL16S202004291200_sorted.fasta -minseqlength 1200 -maxseqlength 1600
#https://drive5.com/usearch/manual/cmd_cluster_smallmem.html
usearch -cluster_smallmem /work/workspace/zhurj/project/3_16S/clusterfle16S20200526/output/moonbioFL16S202004291200_sorted.fasta -id 0.99 -centroids /work/workspace/zhurj/project/3_16S/clusterfle16S20200526/output/moonbioFL16S202004291200_nr.fasta -uc /work/workspace/zhurj/project/3_16S/clusterfle16S20200526/output/moonbioFL16S202004291200_clusters.uc
seqkit fx2tab -o /work/workspace/zhurj/project/3_16S/clusterfle16S20200526/output/moonbioFL16S202004291200_tab.txt /work/workspace/zhurj/project/3_16S/clusterfle16S20200526/output/moonbioFL16S202004291200_sorted.fasta


# EZBIO
usearch -sortbylength /work/database/ezbio/16S/current/Ezbio_16S_seqs.fa -fastaout /work/workspace/zhurj/project/3_16S/clusterfle16S20200526/ezbio/ezbio_sorted.fasta -minseqlength 1200 -maxseqlength 1600
usearch -cluster_smallmem /work/workspace/zhurj/project/3_16S/clusterfle16S20200526/ezbio/ezbio_sorted.fasta -id 0.99 -centroids /work/workspace/zhurj/project/3_16S/clusterfle16S20200526/ezbio/ezbio_nr.fasta -uc /work/workspace/zhurj/project/3_16S/clusterfle16S20200526/ezbio/ezbio_clusters.uc
seqkit fx2tab -o /work/workspace/zhurj/project/3_16S/clusterfle16S20200526/ezbio/Ezbio_16S_seqs_tab.txt /work/database/ezbio/16S/current/Ezbio_16S_seqs.fa 
Enterococcus faecium AJ301830

seqkit fx2tab -o /work/workspace/zhurj/project/3_16S/clusterfle16S20200526/ncbi/16S_ribosomal_RNA.txt /work/database/ncbi/16S/current/16S_ribosomal_RNA.fas
NR_112039.1     Enterococcus faecium strain JCM 5804 16S ribosomal RNA, partial sequence
NR_114742.1     Enterococcus faecium strain DSM 20477 16S ribosomal RNA, partial sequence
NR_115764.1     Enterococcus faecium strain ATCC 19434 16S ribosomal RNA, partial sequence
NR_042054.1     Enterococcus faecium strain LMG 11423 16S ribosomal RNA, partial sequence
NR_113903.1     Enterococcus faecium strain NBRC 100485 16S ribosomal RNA, partial sequence
NR_113904.1     Enterococcus faecium strain NBRC 100486 16S ribosomal RNA, partial sequence

python3.6
# http://sequenceconversion.bugaco.com/converter/biology/sequences/abi_to_fastq-illumina.php
seqret fastq::phred33 -sformat abi -osformat fastq -auto -stdout -sequence /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/ab1/A8-1-27F.ab1 > /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq/A8-1.1.fastq

seqret fastq::phred33 -sformat abi -osformat fastq -auto -stdout -sequence /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/ab1-20200601/A8-2-27F.ab1 > /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq-20200601/A8-2.1.fastq
seqret fastq::phred33 -sformat abi -osformat fastq -auto -stdout -sequence /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/ab1-20200601/A8-2-1492R_R.ab1 > /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq-20200601/A8-2.2.fastq


from Bio import SeqIO
records = SeqIO.parse("/work/workspace/zhurj/project/3_16S/sample3fl16S20200601/ab1-20200601/A8-2-27F.ab1", "abi")
count = SeqIO.write(records, "/work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq-20200601/A8-2.1.fastq", "fastq-solexa")
from Bio import SeqIO
records = SeqIO.parse("/work/workspace/zhurj/project/3_16S/sample3fl16S20200601/ab1-20200601/A8-2-1492R_R.ab1", "abi")
count = SeqIO.write(records, "/work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq-20200601/A8-2.2.fastq", "fastq-solexa")

seqret -sformat abi -osformat fastq -auto -stdout -sequence /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/ab1-20200601/A8-3-27F.ab1 > /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq-20200601/A8-3.1.fastq
seqret -sformat abi -osformat fastq -auto -stdout -sequence /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/ab1-20200601/A8-3-1492R.ab1 > /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq-20200601/A8-3.2.fastq

seqret -sformat abi -osformat fastq -auto -stdout -sequence /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/ab1-20200601/A8-4-27F.ab1 > /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq-20200601/A8-4.1.fastq
seqret -sformat abi -osformat fastq -auto -stdout -sequence /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/ab1-20200601/A8-4-1492R_R.ab1 > /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq-20200601/A8-4.2.fastq

seqret -sformat abi -osformat fastq -auto -stdout -sequence /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/ab1-20200601/A8-5-27F_R.ab1 > /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq-20200601/A8-5.1.fastq
seqret -sformat abi -osformat fastq -auto -stdout -sequence /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/ab1-20200601/A8-5-1492R.ab1 > /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq-20200601/A8-5.2.fastq

merger -asequence /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq-20200601/A8-2.1.fastq -bsequence /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq-20200601/A8-2.2.fastq -sbegin1 20  -send1 810 -sbegin2 20 630 -outfile /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/merger-20200601/A8-2.merger -outseq /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/merger-20200601/A8-2.fasta 

/work/workspace/zhurj/bin/miniconda3/bin/seqkit stats -a A8-2.1.fastq -j 16 -T
/work/workspace/zhurj/bin/miniconda3/bin/seqkit stats -a A8-2.2.fastq -j 16 -T

trimmomatic SE -threads 16 /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq-20200601/A8-2.1.fastq /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/clean-20200601/A8-2.clean.1.fastq LEADING:30 TRAILING:30 SLIDINGWINDOW:6:20
trimmomatic SE -threads 16 /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq-20200601/A8-2.2.fastq /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/clean-20200601/A8-2.clean.2.fastq LEADING:30 TRAILING:30 SLIDINGWINDOW:6:20

trimmomatic SE -threads 16 /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq-20200601/A8-2.1.fastq /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/clean-20200601/A8-2.clean.1.fastq CROP:810 HEADCROP:20
trimmomatic SE -threads 16 /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq-20200601/A8-2.2.fastq /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/clean-20200601/A8-2.clean.2.fastq CROP:710 HEADCROP:20

merger -asequence /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/clean-20200601/A8-2.clean.1.fastq -bsequence /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/clean-20200601/A8-2.clean.2.fastq -outfile /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/merger-20200601/A8-2.merger -outseq /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/merger-20200601/A8-2.fasta 
merger -asequence /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fasta-20200601/A8-2-F.fasta -bsequence /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fasta-20200601/A8-2-R.fasta -outfile /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/merger-20200601/A8-2-tmp.merger -outseq /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/merger-20200601/A8-2-tmp.fasta 
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/merger-20200601/A8-2-tmp.fasta -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/blast-20200601/A8-2-tmp.out -max_target_seqs 5 -num_threads 16 -outfmt 7 

flash /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/clean-20200601/A8-2.clean.1.fastq /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/clean-20200601/A8-2.clean.2.fastq -m 10 -x 0.2 -t 16 -o A8-2.merge.fastq -d /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/merger-20200601

from Bio import SeqIO

records = SeqIO.parse("/work/workspace/zhurj/project/3_16S/sample3fl16S20200601/merger-20200601/A8-2.merge.fastq", "fastq")
count = SeqIO.write(records, "/work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fasta-20200601/A8-2.fasta", "fasta")

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fasta-20200601/A8-2.fasta -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/blast-20200601/A8-2.out -max_target_seqs 5 -num_threads 16 -outfmt 7 

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fasta-20200601/A8-2-F.fasta -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/blast-20200601/A8-2-F.out -max_target_seqs 5 -num_threads 16 -outfmt 7
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fasta-20200601/A8-2-R.fasta -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/blast-20200601/A8-2-R.out -max_target_seqs 5 -num_threads 16 -outfmt 7 

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "BCVU01000117|CP001071|PJKF01000002|PJKB01000002"


# BLASTN 2.9.0+
# Query: A8-2-F
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
A8-2-F  GU983697  99.747  790 1 1 1 790 10  798 0.0 1413
A8-2-F  BCQB01000108  99.747  790 1 1 1 790 33  821 0.0 1413
A8-2-F  CP003504  99.620  790 2 1 1 790 33  821 0.0 1408
A8-2-F  KF621060  99.494  790 3 1 1 790 33  821 0.0 1404
A8-2-F  JXLE01000039  98.987  790 7 1 1 790 33  821 0.0 1386
# BLAST processed 1 queries
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "GU983697|KF621060|BCQB01000108|CP003504|JXLE01000039"
BCQB01000108    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus durans;NBRC 100479(T)
CP003504        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus hirae;ATCC 9790(T)
GU983697        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus lactis;BT159(T)
JXLE01000039    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus thailandicus;DSM 21767(T)
KF621060        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus xinjiangensis;48(T)


# BLASTN 2.9.0+
# Query: A8-2-R
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
A8-2-R  CP015992  99.420  690 4 0 1 690 740 1429  0.0 1227
A8-2-R  AP013068  99.420  690 4 0 1 690 740 1429  0.0 1227
A8-2-R  AJMR01000229  99.420  690 4 0 1 690 740 1429  0.0 1227
A8-2-R  CP032616  99.275  690 5 0 1 690 739 1428  0.0 1223
A8-2-R  SCOM01000009  99.130  690 6 0 1 690 740 1429  0.0 1218
# BLAST processed 1 queries
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "CP015992|AP013068|AJMR01000229|CP032616|SCOM01000009"
AJMR01000229    Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;Pseudomonas furukawaii;KF707(T)
AP013068        Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;NBRC 106553
CP015992        Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;TCU-HL1
CP032616        Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;DY-1
SCOM01000009    Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;SCOM_s;F(2018)


# BLASTN 2.9.0+
# Query: A8-2-merger
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
A8-2-F  GU983697        91.845  1398    110     3       1       1397    10      1404    0.0     1998
A8-2-F  KF621060        91.940  1402    101     8       1       1397    33      1427    0.0     1994
A8-2-F  BCQB01000108    91.702  1398    112     3       1       1397    33      1427    0.0     1989
A8-2-F  CP003504        91.631  1398    113     3       1       1397    33      1427    0.0     1984
A8-2-F  JXLB01000051    91.345  1398    117     3       1       1397    33      1427    0.0     1966
# BLAST processed 1 queries
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "GU983697|KF621060|BCQB01000108|CP003504|JXLB01000051"
BCQB01000108    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus durans;NBRC 100479(T)
CP003504        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus hirae;ATCC 9790(T)
GU983697        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus lactis;BT159(T)
JXLB01000051    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus ratti;DSM 15687(T)
KF621060        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus xinjiangensis;48(T)




from Bio import SeqIO
records = SeqIO.parse("/work/workspace/zhurj/project/3_16S/sample3fl16S20200601/ab1-20200601/A8-3-27F.ab1", "abi")
count = SeqIO.write(records, "/work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq-20200601/A8-3.1.fastq", "fastq-solexa")

records = SeqIO.parse("/work/workspace/zhurj/project/3_16S/sample3fl16S20200601/ab1-20200601/A8-3-1492R.ab1", "abi")
count = SeqIO.write(records, "/work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq-20200601/A8-3.2.fastq", "fastq-solexa")

trimmomatic SE -threads 16 /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq-20200601/A8-3.1.fastq /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/clean-20200601/A8-3.clean.1.fastq CROP:700 HEADCROP:30
trimmomatic SE -threads 16 /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq-20200601/A8-3.2.fastq /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/clean-20200601/A8-3.clean.2.fastq CROP:830 HEADCROP:30

merger -asequence /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fasta-20200601/A8-3-F.fasta -bsequence /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fasta-20200601/A8-3-R.fasta -outfile /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/merger-20200601/A8-3-merger.merger -outseq /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/merger-20200601/A8-3-merger.fasta 
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/merger-20200601/A8-3-merger.fasta -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/blast-20200601/A8-3-merger.out -max_target_seqs 5 -num_threads 16 -outfmt 7 
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fasta-20200601/A8-3-F.fasta -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/blast-20200601/A8-3-F.out -max_target_seqs 5 -num_threads 16 -outfmt 7
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fasta-20200601/A8-3-R.fasta -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/blast-20200601/A8-3-R.out -max_target_seqs 5 -num_threads 16 -outfmt 7 

# Query: A8-3-27F
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
A8-3-27F  GU983697  99.701  670 1 1 1 670 15  683 0.0 1196
A8-3-27F  BCQB01000108  99.701  670 1 1 1 670 38  706 0.0 1196
A8-3-27F  CP003504  99.701  669 1 1 2 670 39  706 0.0 1195
A8-3-27F  KF621060  99.403  670 3 1 1 670 38  706 0.0 1187
A8-3-27F  JXLE01000039  98.955  670 6 1 1 670 38  706 0.0 1174
# BLAST processed 1 queries

# Query: A8-3-R
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
A8-3-R  GU983697  99.750  800 2 0 1 800 618 1417  0.0 1434
A8-3-R  AJ301830  99.750  800 1 1 1 800 640 1438  0.0 1431
A8-3-R  JXKV01000056  99.500  800 4 0 1 800 641 1440  0.0 1425
A8-3-R  BCQB01000108  99.500  800 4 0 1 800 641 1440  0.0 1425
A8-3-R  LC473138  99.375  800 5 0 1 800 585 1384  0.0 1421
# BLAST processed 1 queries

# Query: A8-3-MERGER
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
A8-3-27F  GU983697  99.715  1404  3 1 1 1404  15  1417  0.0 2511
A8-3-27F  BCQB01000108  99.573  1404  5 1 1 1404  38  1440  0.0 2502
A8-3-27F  CP003504  99.501  1403  6 1 2 1404  39  1440  0.0 2496
A8-3-27F  AJ301830  99.359  1405  3 6 2 1404  38  1438  0.0 2472
A8-3-27F  JXLE01000039  99.074  1404  12  1 1 1404  38  1440  0.0 2470
# BLAST processed 1 queries

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "GU983697|BCQB01000108|CP003504|AJ301830|JXLE01000039|JXKV01000056|LC473138|KF621060"
AJ301830        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus faecium;LMG 11423(T)
BCQB01000108    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus durans;NBRC 100479(T)
CP003504        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus hirae;ATCC 9790(T)
GU983697        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus lactis;BT159(T)
JXKV01000056    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus mundtii;DSM 4838(T)
JXLE01000039    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus thailandicus;DSM 21767(T)
KF621060        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus xinjiangensis;48(T)
LC473138        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus hulanensis;190-7(T)


from Bio import SeqIO
records = SeqIO.parse("/work/workspace/zhurj/project/3_16S/sample3fl16S20200601/ab1-20200601/A8-4-27F.ab1", "abi")
count = SeqIO.write(records, "/work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq-20200601/A8-4.1.fastq", "fastq-solexa")

records = SeqIO.parse("/work/workspace/zhurj/project/3_16S/sample3fl16S20200601/ab1-20200601/A8-4-1492R_R.ab1", "abi")
count = SeqIO.write(records, "/work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq-20200601/A8-4.2.fastq", "fastq-solexa")

trimmomatic SE -threads 16 /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq-20200601/A8-4.1.fastq /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/clean-20200601/A8-4.clean.1.fastq CROP:800 HEADCROP:25
trimmomatic SE -threads 16 /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fastq-20200601/A8-4.2.fastq /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/clean-20200601/A8-4.clean.2.fastq CROP:800 HEADCROP:25

merger -asequence /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fasta-20200601/A8-4-F.fasta -bsequence /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fasta-20200601/A8-4-R.fasta -outfile /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/merger-20200601/A8-4-merger.merger -outseq /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/merger-20200601/A8-4-merger.fasta 
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/merger-20200601/A8-4-merger.fasta -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/blast-20200601/A8-4-merger.out -max_target_seqs 5 -num_threads 16 -outfmt 7 
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fasta-20200601/A8-4-F.fasta -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/blast-20200601/A8-4-F.out -max_target_seqs 5 -num_threads 16 -outfmt 7
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/fasta-20200601/A8-4-R.fasta -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/3_16S/sample3fl16S20200601/blast-20200601/A8-4-R.out -max_target_seqs 5 -num_threads 16 -outfmt 7 

# Query: A8-4-F
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
A8-4-F  GU983697  99.613  775 1 2 1 775 10  782 0.0 1378
A8-4-F  BCQB01000108  99.613  775 1 2 1 775 33  805 0.0 1378
A8-4-F  CP003504  99.484  775 2 2 1 775 33  805 0.0 1373
A8-4-F  KF621060  99.355  775 3 2 1 775 33  805 0.0 1369
A8-4-F  JXLE01000039  98.839  775 7 1 1 775 33  805 0.0 1355
# BLAST processed 1 queries
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "GU983697|BCQB01000108|CP003504|KF621060|JXLE01000039"
BCQB01000108    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus durans;NBRC 100479(T)
CP003504        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus hirae;ATCC 9790(T)
GU983697        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus lactis;BT159(T)
JXLE01000039    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus thailandicus;DSM 21767(T)
KF621060        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus xinjiangensis;48(T)

# Query: A8-4-R
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
A8-4-R  OLKI01000048  100.000 775 0 0 1 775 652 1426  0.0 1398
A8-4-R  NEIG01000032  100.000 775 0 0 1 775 652 1426  0.0 1398
A8-4-R  NC_021491 100.000 775 0 0 1 775 652 1426  0.0 1398
A8-4-R  MK680061  100.000 775 0 0 1 775 652 1426  0.0 1398
A8-4-R  MH517510  100.000 775 0 0 1 775 652 1426  0.0 1398
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "OLKI01000048|NEIG01000032|NC_021491|MK680061|MH517510"
MH517510        Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;Pseudomonas asiatica;RYU5(T)
MK680061        Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;Pseudomonas juntendi;BML3(T)
NC_021491       Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;CP005976_s;H8234
NEIG01000032    Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;NEIG_s;R17(2017)
OLKI01000048    Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;Pseudomonas shirazica;VM14(T)

# Query: A8-4-merge
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
A8-4-F  KF621060  92.516  1403  92  8 1 1398  33  1427  0.0 2033
A8-4-F  GU983697  91.761  1420  112 4 1 1419  10  1425  0.0 2021
A8-4-F  BCQB01000108  91.620  1420  114 4 1 1419  33  1448  0.0 2012
A8-4-F  CP003504  91.549  1420  115 4 1 1419  33  1448  0.0 2007
A8-4-F  JXLB01000051  91.708  1399  111 4 1 1398  33  1427  0.0 1987

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "KF621060|GU983697|BCQB01000108|CP003504|JXLB01000051"

2020.06.11
38个样品的分析
(2) qiime 
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/input/manifest \
  --output-path /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/reads.qza \
  --input-format PairedEndFastqManifestPhred33

dada2去燥,合并双端序列
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/reads.qza \
  --p-trim-left-f 25 \
  --p-trim-left-r 26 \
  --p-trunc-len-f 200 \
  --p-trunc-len-r 180 \
  --o-table /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/table.qza \
  --o-representative-sequences /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/rep-seqs.qza \
  --o-denoising-stats /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/denoising-stats.qza \
  --p-n-threads 20

  # 可视化denoising stats（qzv文件可在线查看）
  qiime metadata tabulate \
  --m-input-file /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/denoising-stats.qza \
  --o-visualization /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qzv/denoising-stats.qzv

# 先不filter feature -- 20200611
  # filter FeatureTable
  qiime feature-table filter-features \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/table.qza \
  --p-min-frequency 2 \
  --m-metadata-file /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/input/metadata \
  --p-exclude-ids \
  --o-filtered-table /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/table_filter.qza \

  # tree create for diversity analysis
  qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/rep-seqs.qza \
  --o-alignment /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/aligned-rep-seqs.qza \
  --o-masked-alignment /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/masked-aligned-rep-seqs.qza \
  --o-tree /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/unrooted-tree.qza \
  --o-rooted-tree /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/rooted-tree.qza \
  --p-n-threads 20

  ### ======================
  --output-dir /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza \
  ## output-dir 替代一下四行
  --o-alignment /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/aligned-rep-seqs.qza \
  --o-masked-alignment /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/masked-aligned-rep-seqs.qza \
  --o-tree /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/unrooted-tree.qza \
  --o-rooted-tree /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/rooted-tree.qza
  ###

  # alpha rarefaction 
  qiime diversity alpha-rarefaction \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/table.qza \
  --p-metrics  observed_otus --p-metrics shannon  --p-metrics faith_pd \
  --p-metrics goods_coverage --p-metrics chao1 --p-metrics simpson \
  --p-max-depth 50000 \
  --m-metadata-file /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/input/metadata \
  --i-phylogeny /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/rooted-tree.qza \
  --o-visualization /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qzv/alpha_rarefaction.qzv

  # diversity analysis
  qiime diversity core-metrics-phylogenetic \
  --i-phylogeny /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/rooted-tree.qza \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/table.qza \
  --p-sampling-depth 19508 \
  --m-metadata-file /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/input/metadata \
  --output-dir /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/core-metrics-results


  # training classifier
  # 已经train过， 本次38个样品无需再trainning
  qiime feature-classifier extract-reads \
  --i-sequences reference/99_otus.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \      #515 F引物
  --p-r-primer GGACTACHVGGGTWTCTAAT \  #806 R引物
  --o-reads /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/ref/ref-seqs.qza

  # Train the classifier（训练分类器）
  # 基于筛选的指定区段，生成实验特异的分类器
  qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/ref/ref-seqs.qza \
  --i-reference-taxonomy reference/ref-taxonomy.qza \
  --o-classifier /work/workspace/ruijuan/project/3_16S/fatmouse20190717/analy20190717/result_30000/ref/classifier.qza

  # taxa classifier
  qiime feature-classifier classify-sklearn \
  --i-classifier /work/workspace/ruijuan/script/git/pipeline16S/reference/Greengenes_13_8_99%_OTUs_341F-805R_classifier.qza \
  --i-reads /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/rep-seqs.qza \
  --o-classification /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/taxonomy.qza

  # 结果可视化
  qiime metadata tabulate \
  --m-input-file /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/taxonomy.qza \
  --o-visualization /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qzv/taxonomy.qzv

  # 物种分类条形图
  qiime taxa barplot \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/table.qza \
  --i-taxonomy /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/taxonomy.qza \
  --m-metadata-file /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/input/metadata \
  --o-visualization /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qzv/taxa-bar-plots.qzv

  qiime tools export \
  --input-path /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qzv/taxa-bar-plots.qzv \
  --output-path /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qzv/taxa

   qiime tools export \
  --input-path /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/rep-seqs.qza \
  --output-path /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/rep-seqs

  # acom analysis
  #按照 leve2 水平分析
  qiime taxa collapse \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/table.qza \
  --i-taxonomy /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/taxonomy.qza \
  --p-level 2 \
  --o-collapsed-table /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/acom/qza/taxa_table_l2.qza

  qiime composition add-pseudocount \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/acom/qza/taxa_table_l2.qza \
  --o-composition-table /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/acom/qza/comp_table_l2.qza

  qiime composition ancom \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/acom/qza/comp_table_l2.qza \
  --m-metadata-file /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/input/metadata \
  --m-metadata-column Time \
  --p-transform-function clr \
  --o-visualization /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/acom/qzv/ancom_l2_clr_Condition.qzv

  # level 6
  qiime taxa collapse \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/table.qza \
  --i-taxonomy /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/acom/qza/taxa_table_l6.qza

  qiime composition add-pseudocount \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/acom/qza/taxa_table_l6.qza \
  --o-composition-table /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/acom/qza/comp_table_l6.qza 

  qiime composition ancom \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/acom/qza/comp_table_l6.qza \
  --m-metadata-file /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/input/metadata \
  --m-metadata-column Time \
  --p-transform-function clr \
  --o-visualization /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/acom/qzv/ancom_l6_clr_Condition.qzv

  # level 7
  qiime taxa collapse \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/table.qza \
  --i-taxonomy /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/acom/qza/taxa_table_l7.qza

  qiime composition add-pseudocount \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/acom/qza/taxa_table_l7.qza \
  --o-composition-table /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/acom/qza/comp_table_l7.qza

  qiime composition ancom \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/acom/qza/comp_table_l7.qza \
  --m-metadata-file /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/input/metadata \
  --m-metadata-column Time \
  --p-transform-function clr \
  --o-visualization /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/acom/qzv/ancom_l7_clr_Condition.qzv

qiime tools export \
  --input-path /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/acom/qzv/ancom_l7_clr_Condition.qzv \
  --output-path /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/acom/qzv/ancom_l7

qiime tools export \
  --input-path /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/acom/qzv/ancom_l6_clr_Condition.qzv \
  --output-path /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/acom/qzv/ancom_l6
qiime tools export \
  --input-path /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/acom/qzv/ancom_l2_clr_Condition.qzv \
  --output-path /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/acom/qzv/ancom_l2
qiime tools export \
  --input-path /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/table.qza \
  --output-path /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/table


source deactivate 16S-2019.04
  source activate picrust2

  place_seqs.py \
  -s /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/rep-seqs/dna-sequences.fasta \
  -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/out.tre \
  -p 20
 
  hsp.py -i 16S \
  -t /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/out.tre \
  -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/marker_predicted_and_nsti.tsv.gz \
  -p 20

  hsp.py -i COG \
  -t /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/out.tre \
  -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/COG_predicted.tsv.gz \
  -p 20

  metagenome_pipeline.py -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/table/feature-table.biom \
  -m /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/marker_predicted_and_nsti.tsv.gz \
  -f /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/COG_predicted.tsv.gz \
  -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/COG \
  --strat_out

  add_descriptions.py -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/COG/pred_metagenome_unstrat.tsv \
  -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/COG/pred_metagenome_unstrat_descrip.tsv \
  -m COG

  hsp.py -i EC \
  -t /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/out.tre \
  -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/EC_predicted.tsv.gz \
  -p 20

  metagenome_pipeline.py -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/table/feature-table.biom \
  -m /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/marker_predicted_and_nsti.tsv.gz \
  -f /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/EC_predicted.tsv.gz \
  -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/EC \
  --strat_out

  add_descriptions.py -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/EC/pred_metagenome_unstrat.tsv \
  -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/EC/pred_metagenome_unstrat_descrip.tsv \
  -m EC

  hsp.py -i KO -t /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/out.tre \
  -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/KO_predicted.tsv.gz \
  -p 20

  metagenome_pipeline.py -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/table/feature-table.biom \
  -m /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/marker_predicted_and_nsti.tsv.gz \
  -f /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/KO_predicted.tsv.gz \
  -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/KO \
  --strat_out

  add_descriptions.py -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/KO/pred_metagenome_unstrat.tsv \
  -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/KO/pred_metagenome_unstrat_descrip.tsv \
  -m KO

  pathway_pipeline.py -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/EC/pred_metagenome_strat.tsv \
  -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/metacyc \
  -p 20

  add_descriptions.py -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/metacyc/path_abun_unstrat.tsv \
  -m METACYC \
  -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/metacyc/path_abun_unstrat_descrip.tsv
 
  pathway_pipeline.py -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/KO/pred_metagenome_strat.tsv \
  -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/kegg_module \
  --no_regroup \
  --map /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/ref/pathway_mapfiles/KEGG_modules_to_KO.tsv \
  -p 20

  add_descriptions.py -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/kegg_module/path_abun_unstrat.tsv \
  --custom_map_table /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/ref/description_mapfiles/KEGG_modules_info_adj.tsv \
  -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/kegg_module/path_abun_unstrat_descrip.tsv

  pathway_pipeline.py -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/KO/pred_metagenome_strat.tsv \
  -o //work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/kegg_path/all \
  --no_regroup \
  --map /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/ref/pathway_mapfiles/KEGG_pathways_to_KO.tsv \
  -p 20

  add_descriptions.py -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/kegg_path/all/path_abun_unstrat.tsv \
  --custom_map_table /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/ref/description_mapfiles/KEGG_pathways_info.tsv.gz \
  -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/kegg_path/all/path_abun_unstrat_descrip.tsv

  pathway_pipeline.py -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/KO/pred_metagenome_strat.tsv \
  -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/kegg_path/level2 \
  --no_regroup \
  --map /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/ref/kegg/kegg_path_ko_l2.txt \
  -p 20

  add_descriptions.py -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/kegg_path/level2/path_abun_unstrat.tsv \
  --custom_map_table /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/ref/kegg/kegg_path_description.txt \
  -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/kegg_path/level2/path_abun_unstrat_descrip.tsv

  pathway_pipeline.py -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/KO/pred_metagenome_strat.tsv \
  -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/kegg_path/level3 \
  --no_regroup \
  --map /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/ref/kegg/kegg_path_ko_l3.txt \
  -p 20

  add_descriptions.py -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/kegg_path/level3/path_abun_unstrat.tsv \
  --custom_map_table /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/ref/kegg/kegg_path_description.txt \
  -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/picrust2/kegg_path/level3/path_abun_unstrat_descrip.tsv

 #--------------------------------------------------------------------------------------------------------------------------------------------------------------
  
source deactivate picrust2
source activate lefse

mkdir -p /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level2
Rscript /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/lefsein_taxa.r -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qzv/taxa/level-2.csv -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level2/lefse.in --class Time -u index
lefse-format_input.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level2/lefse.in /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level2/lefse_format.in -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level2/lefse_format.in /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level2/lefse_run.res -a 0.05 -w 0.05 -l 2 -y 0
lefse-plot_res.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level2/lefse_run.res /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level2/lefse_run.png --feature_font_size 8 --width 10 --format png
lefse-plot_cladogram.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level2/lefse_run.res /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level2/lefse_cladogram.png --format png

mkdir -p /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level6
Rscript /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/lefsein_taxa.r -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qzv/taxa/level-6.csv -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level6/lefse.in --class Time -u index
lefse-format_input.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level6/lefse.in /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level6/lefse_format.in -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level6/lefse_format.in /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level6/lefse_run.res -a 0.05 -w 0.05 -l 2 -y 0
lefse-plot_res.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level6/lefse_run.res /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level6/lefse_run.png --feature_font_size 8 --width 10 --format png
lefse-plot_cladogram.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level6/lefse_run.res /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level6/lefse_cladogram.png --format png

mkdir -p /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level7
Rscript /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/lefsein_taxa.r -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qzv/taxa/level-7.csv -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level7/lefse.in --class Time -u index
lefse-format_input.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level7/lefse.in /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level7/lefse_format.in -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level7/lefse_format.in /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level7/lefse_run.res -a 0.05 -w 0.05 -l 2 -y 0
lefse-plot_res.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level7/lefse_run.res /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level7/lefse_run.png --feature_font_size 8 --width 10 --format png
lefse-plot_cladogram.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level7/lefse_run.res /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/level7/lefse_cladogram.png --format png

cd /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/
mkdir Time
mv level* Time

Type: Content vs ContentTissue
mkdir -p /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level2
Rscript /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/lefsein_taxa.r -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qzv/taxa/level-2.csv -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level2/lefse.in --class Type -u index
lefse-format_input.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level2/lefse.in /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level2/lefse_format.in -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level2/lefse_format.in /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level2/lefse_run.res -a 0.05 -w 0.05 -l 2 -y 0
lefse-plot_res.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level2/lefse_run.res /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level2/lefse_run.png --feature_font_size 8 --width 10 --format png
lefse-plot_cladogram.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level2/lefse_run.res /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level2/lefse_cladogram.png --format png

mkdir -p /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level6
Rscript /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/lefsein_taxa.r -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qzv/taxa/level-6.csv -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level6/lefse.in --class Type -u index
lefse-format_input.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level6/lefse.in /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level6/lefse_format.in -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level6/lefse_format.in /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level6/lefse_run.res -a 0.05 -w 0.05 -l 2 -y 0
lefse-plot_res.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level6/lefse_run.res /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level6/lefse_run.png --feature_font_size 8 --width 10 --format png
lefse-plot_cladogram.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level6/lefse_run.res /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level6/lefse_cladogram.png --format png

mkdir -p /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level7
Rscript /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/lefsein_taxa.r -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qzv/taxa/level-7.csv -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level7/lefse.in --class Type -u index
lefse-format_input.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level7/lefse.in /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level7/lefse_format.in -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level7/lefse_format.in /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level7/lefse_run.res -a 0.05 -w 0.05 -l 2 -y 0
lefse-plot_res.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level7/lefse_run.res /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level7/lefse_run.png --feature_font_size 8 --width 10 --format png
lefse-plot_cladogram.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level7/lefse_run.res /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Type/level7/lefse_cladogram.png --format png

Location: cecum vs colon
mkdir -p /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level2
Rscript /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/lefsein_taxa.r -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qzv/taxa/level-2.csv -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level2/lefse.in --class Location -u index
lefse-format_input.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level2/lefse.in /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level2/lefse_format.in -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level2/lefse_format.in /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level2/lefse_run.res -a 0.05 -w 0.05 -l 2 -y 0
lefse-plot_res.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level2/lefse_run.res /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level2/lefse_run.png --feature_font_size 8 --width 10 --format png
lefse-plot_cladogram.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level2/lefse_run.res /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level2/lefse_cladogram.png --format png

mkdir -p /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level6
Rscript /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/lefsein_taxa.r -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qzv/taxa/level-6.csv -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level6/lefse.in --class Location -u index
lefse-format_input.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level6/lefse.in /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level6/lefse_format.in -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level6/lefse_format.in /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level6/lefse_run.res -a 0.05 -w 0.05 -l 2 -y 0
lefse-plot_res.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level6/lefse_run.res /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level6/lefse_run.png --feature_font_size 8 --width 10 --format png
lefse-plot_cladogram.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level6/lefse_run.res /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level6/lefse_cladogram.png --format png

mkdir -p /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level7
Rscript /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/lefsein_taxa.r -i /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qzv/taxa/level-7.csv -o /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level7/lefse.in --class Location -u index
lefse-format_input.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level7/lefse.in /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level7/lefse_format.in -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level7/lefse_format.in /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level7/lefse_run.res -a 0.05 -w 0.05 -l 2 -y 0
lefse-plot_res.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level7/lefse_run.res /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level7/lefse_run.png --feature_font_size 8 --width 10 --format png
lefse-plot_cladogram.py /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level7/lefse_run.res /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/diff/lefse/taxa/Location/level7/lefse_cladogram.png --format png




echo "Enterococcus faecium" | taxonkit name2taxid | taxonkit lineage --taxid-field 2
echo "Enterococcus faecium" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t
echo "Megasphaera massiliensis" | taxonkit name2taxid | taxonkit lineage --taxid-field 2 |cut -f 2,3 | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species | csvtk pretty -t

/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/rep-seqs/dna-sequences.fasta

from Bio import SeqIO
records = SeqIO.parse("/work/rawdata/ab1/16S/MNH/MNH279/MNH27992/E1/L1/S1/MNH27992-27F.ab1", "abi")
count = SeqIO.write(records, "/work/workspace/zhurj/project/3_16S/sam38baoxiaMM20200608/MNH27992/fastq/MNH27992.1.fastq", "fastq-solexa")

records = SeqIO.parse("/work/rawdata/ab1/16S/MNH/MNH279/MNH27992/E1/L1/S1/MNH27992-1492R.ab1", "abi")
count = SeqIO.write(records, "/work/workspace/zhurj/project/3_16S/sam38baoxiaMM20200608/MNH27992/fastq/MNH27992.2.fastq", "fastq-solexa")

trimmomatic SE -threads 16 /work/workspace/zhurj/project/3_16S/sam38baoxiaMM20200608/MNH27992/fastq/MNH27992.1.fastq /work/workspace/zhurj/project/3_16S/sam38baoxiaMM20200608/MNH27992/fastq/MNH27992.clean.1.fastq CROP:800 HEADCROP:30
trimmomatic SE -threads 16 /work/workspace/zhurj/project/3_16S/sam38baoxiaMM20200608/MNH27992/fastq/MNH27992.2.fastq /work/workspace/zhurj/project/3_16S/sam38baoxiaMM20200608/MNH27992/fastq/MNH27992.clean.2.fastq CROP:800 HEADCROP:31

merger -asequence /work/workspace/zhurj/project/3_16S/sam38baoxiaMM20200608/MNH27992/fasta/MNH27992-27F.fasta -bsequence /work/workspace/zhurj/project/3_16S/sam38baoxiaMM20200608/MNH27992/fasta/MNH27992-1492R.fasta -outfile /work/workspace/zhurj/project/3_16S/sam38baoxiaMM20200608/MNH27992/fasta/MNH27992-merger.merger -outseq /work/workspace/zhurj/project/3_16S/sam38baoxiaMM20200608/MNH27992/fasta/MNH27992-merger.fasta 
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/3_16S/sam38baoxiaMM20200608/MNH27992/fasta/MNH27992-merger.fasta -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/3_16S/sam38baoxiaMM20200608/MNH27992/blastn/MNH27992-merger.out -max_target_seqs 5 -num_threads 16 -outfmt 7 

# BLASTN 2.9.0+
# Query: MNH27992-27F
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
MNH27992-27F    JX424772        99.074  1404    13      0       1       1404    39      1442    0.0     2474
MNH27992-27F    HE576794        97.013  1406    40      2       1       1404    39      1444    0.0     2340
MNH27992-27F    HM990965        96.167  1409    49      4       1       1404    39      1447    0.0     2285
MNH27992-27F    GQ480004        94.554  1414    64      4       1       1404    39      1449    0.0     2197
MNH27992-27F    HQ801062        94.160  1404    80      2       3       1404    27      1430    0.0     2156

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "JX424772|HE576794"
HE576794        Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Megasphaera;Megasphaera elsdenii;DSM 20460(T)
JX424772        Bacteria;Firmicutes;Negativicutes;Veillonellales;Veillonellaceae;Megasphaera;Megasphaera massiliensis;NP3(T)

/work/program/current/ncbi-blast/bin/makeblastdb -in /work/workspace/zhurj/project/3_16S/sam38baoxiaMM20200608/resultall_20200609/fasta/features.fasta -input_type fasta -dbtype 'nucl' -parse_seqids -out features -hash_index
/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/3_16S/sam38baoxiaMM20200608/MNH27992/fasta/v4noprimer.fasta -task blastn -db /work/workspace/zhurj/project/3_16S/sam38baoxiaMM20200608/resultall_20200609/fasta/features -out /work/workspace/zhurj/project/3_16S/sam38baoxiaMM20200608/resultall_20200609/blastn/MNH27992-v4.out -max_target_seqs 5 -num_threads 16 -outfmt 7 

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/3_16S/sam38baoxiaMM20200608/MNH27992/fasta/v4noprimer.fasta -task blastn -db /work/workspace/zhurj/project/3_16S/sam38baoxiaMM20200608/resultall_20200609/fasta/features -out /work/workspace/zhurj/project/3_16S/sam38baoxiaMM20200608/resultall_20200609/blastn/MNH27992-v4-format.out -max_target_seqs 5 -num_threads 16 

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/3_16S/sam38baoxiaMM20200608/resultall_20200609/MMsimilar/mmsimilar.fasta -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/3_16S/sam38baoxiaMM20200608/resultall_20200609/blastn/mmsimilar.out -max_target_seqs 5 -num_threads 16 -outfmt 7 

# BLASTN 2.9.0+
# Query: 2747bd94a19fbfa82f1e90e326d2f100
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
2747bd94a19fbfa82f1e90e326d2f100  PAC001560 100.000 253 0 0 1 253 500 752 3.81e-128 457
2747bd94a19fbfa82f1e90e326d2f100  EU469921  98.814  253 3 0 1 253 500 752 8.39e-124 444
2747bd94a19fbfa82f1e90e326d2f100  GQ451293  98.024  253 5 0 1 253 500 752 4.34e-121 434
2747bd94a19fbfa82f1e90e326d2f100  AF371828  96.838  253 8 0 1 253 500 752 2.74e-117 421
2747bd94a19fbfa82f1e90e326d2f100  HM124281  96.443  253 9 0 1 253 499 751 1.17e-115 416
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "PAC001560|EU469921"
EU469921        Bacteria;Firmicutes;Clostridia;Clostridiales;Oscillospiraceae;PAC000672_g;GOR_aag74d05
PAC001560       Bacteria;Firmicutes;Clostridia;Clostridiales;Oscillospiraceae;PAC000672_g;None

# BLASTN 2.9.0+
# Query: aff516c121bd4574170ec019010c5d86
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
aff516c121bd4574170ec019010c5d86  HM124015  98.024  253 5 0 1 253 499 751 4.34e-121 434
aff516c121bd4574170ec019010c5d86  AM405167  97.628  253 6 0 1 253 497 749 5.29e-120 430
aff516c121bd4574170ec019010c5d86  EU794153  96.838  253 8 0 1 253 499 751 2.74e-117 421
aff516c121bd4574170ec019010c5d86  PAC001245 96.047  253 10  0 1 253 497 749 1.42e-114 412
aff516c121bd4574170ec019010c5d86  EF445183  96.047  253 10  0 1 253 496 748 1.42e-114 412
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "HM124015"
HM124015        Bacteria;Firmicutes;Clostridia;Clostridiales;Oscillospiraceae;HM124015_g;G47

# BLASTN 2.9.0+
# Query: 41b4cb52f742b40791edeaa8c577c716
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
41b4cb52f742b40791edeaa8c577c716  PAC002348 99.605  253 1 0 1 253 498 750 1.62e-126 453
41b4cb52f742b40791edeaa8c577c716  PAC000684 96.443  253 9 0 1 253 498 750 1.17e-115 416
41b4cb52f742b40791edeaa8c577c716  PAC001776 96.047  253 10  0 1 253 498 750 1.42e-114 412
41b4cb52f742b40791edeaa8c577c716  PAC001705 96.047  253 10  0 1 253 498 750 1.42e-114 412
41b4cb52f742b40791edeaa8c577c716  EU344160  96.047  253 10  0 1 253 491 743 1.42e-114 412
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "PAC002348"
PAC002348       Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;KE159571_g;None

# BLASTN 2.9.0+
# Query: 7963654625b74816cd1eb1373a31352b
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
7963654625b74816cd1eb1373a31352b  PAC001245 100.000 253 0 0 1 253 497 749 3.81e-128 457
7963654625b74816cd1eb1373a31352b  EF445183  98.419  253 4 0 1 253 496 748 1.02e-122 439
7963654625b74816cd1eb1373a31352b  PAC002888 97.233  253 7 0 1 253 496 748 2.25e-118 425
7963654625b74816cd1eb1373a31352b  DQ793380  96.838  253 8 0 1 253 497 749 2.74e-117 421
7963654625b74816cd1eb1373a31352b  KY978733  96.443  253 9 0 1 253 497 749 1.17e-115 416
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "PAC001245"
PAC001245       Bacteria;Firmicutes;Clostridia;Clostridiales;Oscillospiraceae;Monoglobus;None

# BLASTN 2.9.0+
# Query: cbdc0f4fcc7b1e56aacc302c2b9b870c
# Database: /work/database/ezbio/16S/current/Ezbio_16S_seqs
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
cbdc0f4fcc7b1e56aacc302c2b9b870c  PAC001776 100.000 253 0 0 1 253 498 750 3.81e-128 457
cbdc0f4fcc7b1e56aacc302c2b9b870c  EU455375  98.419  253 4 0 1 253 496 748 1.02e-122 439
cbdc0f4fcc7b1e56aacc302c2b9b870c  PAC001705 97.233  253 7 0 1 253 498 750 2.25e-118 425
cbdc0f4fcc7b1e56aacc302c2b9b870c  EU462940  97.233  253 7 0 1 253 496 748 2.25e-118 425
cbdc0f4fcc7b1e56aacc302c2b9b870c  PAC002348 96.443  253 9 0 1 253 498 750 1.17e-115 416
# BLAST processed 5 queries
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "PAC001776"
PAC001776       Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;KE159571_g;None

/mnt/d/home/ruijuan/.conda/envs/R/bin/Rscript
/mnt/d/home/ruijuan/.conda/envs/16S-2019.04/bin/Rscript

python3 /work/workspace/ruijuan/script/4_16Sseq/1_rplot/lib/betadiv.py \
-m /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/input/metadata \
-g Time \
/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/core-metrics-results \
/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/R/beta/Time


python3 /work/workspace/ruijuan/script/4_16Sseq/1_rplot/lib/betadiv.py \
-m /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/input/metadata \
-g Location \
/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/core-metrics-results \
/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/R/beta/Location
/mnt/d/home/ruijuan/.conda/envs/16S-2019.04/bin/Rscript /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/R/beta/Location/p/bboxplot.r

python3 /work/workspace/ruijuan/script/4_16Sseq/1_rplot/lib/betadiv.py \
-m /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/input/metadata \
-g Type \
/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/core-metrics-results \
/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/R/beta/Type
/mnt/d/home/ruijuan/.conda/envs/16S-2019.04/bin/Rscript /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/R/beta/Type/p/bboxplot.r

PCOA
 python3 /work/workspace/ruijuan/script/4_16Sseq/1_rplot/lib/pcoa.py \
-m /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/input/metadata \
-g Time \
/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/core-metrics-results \
/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/R/pcoa/Time
/mnt/d/home/ruijuan/.conda/envs/16S-2019.04/bin/Rscript /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/R/pcoa/Time/p/pcoaplot.r

 python3 /work/workspace/ruijuan/script/4_16Sseq/1_rplot/lib/pcoa.py \
-m /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/input/metadata \
-g Location \
/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/core-metrics-results \
/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/R/pcoa/Location

/mnt/d/home/ruijuan/.conda/envs/16S-2019.04/bin/Rscript /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/R/pcoa/Location/p/pcoaplot.r

python3 /work/workspace/ruijuan/script/4_16Sseq/1_rplot/lib/pcoa.py \
-m /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/input/metadata \
-g Type \
/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/core-metrics-results \
/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/R/pcoa/Type

/mnt/d/home/ruijuan/.conda/envs/16S-2019.04/bin/Rscript /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/R/pcoa/Type/p/pcoaplot.r

qiime tools export --input-path bray_curtis_distance_matrix.qza --output-path pcoa_bray_curtis

import numpy as np
a = np.loadtxt('speciesmatrix.txt')
b= a/a.sum(axis=1)[:,np.newaxis]
np.savetxt('speciesfre_matrix.txt',b,fmt='%f',delimiter="\t",newline="\n")

a = np.loadtxt('genusmatrix.txt')
b= a/a.sum(axis=1)[:,np.newaxis]
np.savetxt('genusfre_matrix.txt',b,fmt='%f',delimiter="\t",newline="\n")

library(ggplot2)
library(vegan)
library(ggpubr)
library('ggbiplot')
library(reshape2)
library('RColorBrewer')
library('ggsci')
library(dplyr)
library(docopt)
library(getopt)
data = read.table("/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/R/taxa/input/genusType.txt",header = T)
compare_means(log2Abundance ~ Group, data=data, group.by = "Genus")
d <- ggboxplot(data, x="Group", y="log2Abundance", fill = "Group",palette = c("#00AFBB", "#FC4E07"), facet.by = "Genus",ylim = c(0, 16), nrow=2) + 
scale_x_discrete(name="", breaks = c("Content", "ContentTissue"),labels=c("Content" = "Content", "ContentTissue" = "ContentTissue")) +
stat_compare_means(aes(label = ..p.signif..),label.y = 15.5) +
guides(fill=FALSE)
dif <- ggarrange(d, ncol=1,nrow=1,labels=c("A"))
ggsave("/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/R/taxa/output/genusType_20200623.jpg", dif, width = 16, height = 8)

data = read.table("/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/R/taxa/input/genusLocation.txt",header = T,sep="\t")
compare_means(log2Abundance ~ Group, data=data, group.by = "Genus")
d <- ggboxplot(data, x="Group", y="log2Abundance", fill = "Group",palette = c("#00AFBB", "#FC4E07"), facet.by = "Genus",ylim = c(0, 16), nrow=2) + 
scale_x_discrete(name="", breaks = c("colon", "cecum"),labels=c("colon" = "colon", "cecum" = "cecum")) +
stat_compare_means(aes(label = ..p.signif..),label.y = 15.5) +
guides(fill=FALSE)
dif <- ggarrange(d, ncol=1,nrow=1,labels=c("A"))
ggsave("/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/R/taxa/output/genusLocation_20200623.jpg", dif, width = 16, height = 8)

data = read.table("/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/R/taxa/input/genusTime.txt",header = T,sep="\t")
compare_means(log2Abundance ~ Group, data=data, group.by = "Genus")
d <- ggboxplot(data, x="Group", y="log2Abundance", fill = "Group",palette = c("#00AFBB", "#FC4E07"), facet.by = "Genus",ylim = c(0, 16), nrow=1) + 
scale_x_discrete(name="", breaks = c("after", "before"),labels=c("after" = "after", "before" = "before")) +
stat_compare_means(aes(label = ..p.signif..),label.y = 15.5) +
guides(fill=FALSE)
dif <- ggarrange(d, ncol=1,nrow=1,labels=c("A"))
ggsave("/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/R/taxa/output/genusTime_20200623.jpg", dif, width = 4, height = 4)

data = read.table("/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/R/taxa/input/speciesLocation.txt",header = T,sep="\t")
compare_means(log2Abundance ~ Group, data=data, group.by = "Genus")
d <- ggboxplot(data, x="Group", y="log2Abundance", fill = "Group",palette = c("#00AFBB", "#FC4E07"), facet.by = "Genus",ylim = c(0, 16), nrow=2) + 
scale_x_discrete(name="", breaks = c("colon", "cecum"),labels=c("colon" = "colon", "cecum" = "cecum")) +
stat_compare_means(aes(label = ..p.signif..),label.y = 15.5) +
guides(fill=FALSE)
dif <- ggarrange(d, ncol=1,nrow=1,labels=c("A"))
ggsave("/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/R/taxa/output/speciesLocation_20200623.jpg", dif, width = 8, height = 8)

data = read.table("/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/R/taxa/input/speciesType.txt",header = T)
compare_means(log2Abundance ~ Group, data=data, group.by = "Genus")
d <- ggboxplot(data, x="Group", y="log2Abundance", fill = "Group",palette = c("#00AFBB", "#FC4E07"), facet.by = "Genus",ylim = c(0, 16), nrow=2) + 
scale_x_discrete(name="", breaks = c("Content", "ContentTissue"),labels=c("Content" = "Content", "ContentTissue" = "ContentTissue")) +
stat_compare_means(aes(label = ..p.signif..),label.y = 15.5) +
guides(fill=FALSE)
dif <- ggarrange(d, ncol=1,nrow=1,labels=c("A"))
ggsave("/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/R/taxa/output/speciesType_20200623.jpg", dif, width = 10, height = 8)

data = read.table("/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/R/taxa/input/speciesTime.txt",header = T,sep="\t")
compare_means(log2Abundance ~ Group, data=data, group.by = "Genus")
d <- ggboxplot(data, x="Group", y="log2Abundance", fill = "Group",palette = c("#00AFBB", "#FC4E07"), facet.by = "Genus",ylim = c(0, 16), nrow=1) + 
scale_x_discrete(name="", breaks = c("after", "before"),labels=c("after" = "after", "before" = "before")) +
stat_compare_means(aes(label = ..p.signif..),label.y = 15.5) +
guides(fill=FALSE)
dif <- ggarrange(d, ncol=1,nrow=1,labels=c("A"))
ggsave("/work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/R/taxa/output/speciesTime_20200623.jpg", dif, width = 2, height = 4)


from Bio import SeqIO

records = SeqIO.parse("/work/workspace/zhurj/project/3_16S/sample3fl16S20200601/merger-20200601/A8-4-merger.fasta", "fasta")
count = SeqIO.write(records, "/work/workspace/zhurj/project/3_16S/sample3fl16S20200601/merger-20200601/A8-4-merger.txt", "tab")
print("Converted %i records" % count)

from Bio import SeqIO

records = SeqIO.parse("/work/workspace/zhurj/project/3_16S/sample3fl16S20200601/merger-20200601/A8-3-merger.fasta", "fasta")
count = SeqIO.write(records, "/work/workspace/zhurj/project/3_16S/sample3fl16S20200601/merger-20200601/A8-3-merger.txt", "tab")

records = SeqIO.parse("/work/workspace/zhurj/project/3_16S/sample3fl16S20200601/merger-20200601/A8-2-merger.fasta", "fasta")
count = SeqIO.write(records, "/work/workspace/zhurj/project/3_16S/sample3fl16S20200601/merger-20200601/A8-2-merger.txt", "tab")


2020.07.29
35个样品
source activate 16S-2019.04

qiime feature-classifier extract-reads \
  --i-sequences reference/99_otus.qza \
  --p-f-primer CCTAYGGGRBGCASCAG \      #341 F引物
  --p-r-primer GGACTACHVGGGTWTCTAAT \  #806 R引物
  --o-reads /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/ref/ref-seqs.qza

qiime feature-classifier extract-reads --i-sequences /work/workspace/ruijuan/reference/16S/gg_13_8_otus/qza/99_otus.qza --p-f-primer CCTAYGGGRBGCASCAG  --p-r-primer GGACTACHVGGGTWTCTAAT  --o-reads /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/ref/ref-seqs.qza

  # Train the classifier（训练分类器）
  # 基于筛选的指定区段，生成实验特异的分类器
  qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/ref/ref-seqs.qza --i-reference-taxonomy /work/workspace/ruijuan/reference/16S/gg_13_8_otus/qza/ref-taxonomy.qza --o-classifier /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/ref/classifier.qza


(2) qiime 
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/input/manifest \
  --output-path /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/reads.qza \
  --input-format PairedEndFastqManifestPhred33

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/reads.qza \
  --p-trim-left-f 25 \
  --p-trim-left-r 26 \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --o-table /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/table.qza \
  --o-representative-sequences /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/rep-seqs.qza \
  --o-denoising-stats /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/denoising-stats.qza \
  --p-n-threads 20

  # 可视化denoising stats（qzv文件可在线查看）
  qiime metadata tabulate \
  --m-input-file /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/denoising-stats.qza \
  --o-visualization /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qzv/denoising-stats.qzv

# 先不filter feature -- 20200611
  # filter FeatureTable
  qiime feature-table filter-features \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/table.qza \
  --p-min-frequency 2 \
  --m-metadata-file /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/input/metadata \
  --p-exclude-ids \
  --o-filtered-table /work/workspace/ruijuan/project/3_16S/baoxia38sam20200609/result_all/qza/table_filter.qza \

  # tree create for diversity analysis
  qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/rep-seqs.qza \
  --o-alignment /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/aligned-rep-seqs.qza \
  --o-masked-alignment /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/masked-aligned-rep-seqs.qza \
  --o-tree /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/unrooted-tree.qza \
  --o-rooted-tree /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/rooted-tree.qza \
  --p-n-threads 20

  # alpha rarefaction 
  qiime diversity alpha-rarefaction \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/table.qza \
  --p-metrics  observed_otus --p-metrics shannon  --p-metrics faith_pd \
  --p-metrics goods_coverage --p-metrics chao1 --p-metrics simpson \
  --p-max-depth 50000 \
  --m-metadata-file /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/input/meta35 \
  --i-phylogeny /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/rooted-tree.qza \
  --o-visualization /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qzv/alpha_rarefaction.qzv

# diversity analysis
  qiime diversity core-metrics-phylogenetic \
  --i-phylogeny /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/rooted-tree.qza \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/table.qza \
  --p-sampling-depth 19508 \
  --m-metadata-file /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/input/meta35 \
  --output-dir /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/core-metrics-results
 
# taxa classifier
  qiime feature-classifier classify-sklearn \
  --i-classifier /work/workspace/ruijuan/script/git/pipeline16S/reference/Greengenes_13_8_99%_OTUs_341F-805R_classifier.qza \
  --i-reads /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/rep-seqs.qza \
  --o-classification /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/taxonomy.qza

# 结果可视化
  qiime metadata tabulate \
  --m-input-file /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/taxonomy.qza \
  --o-visualization /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qzv/taxonomy.qzv

qiime taxa barplot --i-table /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/table.qza --i-taxonomy /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/taxonomy.qza --m-metadata-file /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/input/meta35 --o-visualization /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qzv/taxa-bar-plots.qzv

qiime tools export \
  --input-path /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qzv/taxa-bar-plots.qzv \
  --output-path /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qzv/taxa

qiime tools export \
  --input-path /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/rep-seqs.qza \
  --output-path /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/rep-seqs


# acom analysis
  #按照 leve2 水平分析
  qiime taxa collapse \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/table.qza \
  --i-taxonomy /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/taxonomy.qza \
  --p-level 2 \
  --o-collapsed-table /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/acom/qza/taxa_table_l2.qza

  qiime composition add-pseudocount \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/acom/qza/taxa_table_l2.qza \
  --o-composition-table /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/acom/qza/comp_table_l2.qza

  qiime composition ancom \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/acom/qza/comp_table_l2.qza \
  --m-metadata-file /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/input/meta35 \
  --m-metadata-column Condition \
  --p-transform-function clr \
  --o-visualization /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/acom/qzv/ancom_l2_clr_Condition.qzv

  # level 6
  qiime taxa collapse \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/table.qza \
  --i-taxonomy /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/acom/qza/taxa_table_l6.qza

  qiime composition add-pseudocount \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/acom/qza/taxa_table_l6.qza \
  --o-composition-table /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/acom/qza/comp_table_l6.qza 

  qiime composition ancom \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/acom/qza/comp_table_l6.qza \
  --m-metadata-file /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/input/meta35 \
  --m-metadata-column Condition \
  --p-transform-function clr \
  --o-visualization /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/acom/qzv/ancom_l6_clr_Condition.qzv

  # level 7
  qiime taxa collapse \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/table.qza \
  --i-taxonomy /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/acom/qza/taxa_table_l7.qza

  qiime composition add-pseudocount \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/acom/qza/taxa_table_l7.qza \
  --o-composition-table /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/acom/qza/comp_table_l7.qza

  qiime composition ancom \
  --i-table /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/acom/qza/comp_table_l7.qza \
  --m-metadata-file /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/input/meta35 \
  --m-metadata-column Condition \
  --p-transform-function clr \
  --o-visualization /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/acom/qzv/ancom_l7_clr_Condition.qzv

qiime tools export \
  --input-path /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/acom/qzv/ancom_l7_clr_Condition.qzv \
  --output-path /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/acom/qzv/ancom_l7

qiime tools export \
  --input-path /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/acom/qzv/ancom_l6_clr_Condition.qzv \
  --output-path /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/acom/qzv/ancom_l6
qiime tools export \
  --input-path /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/acom/qzv/ancom_l2_clr_Condition.qzv \
  --output-path /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/acom/qzv/ancom_l2
qiime tools export \
  --input-path /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/table.qza \
  --output-path /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/table

 source activate picrust2
mkdir -p /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2
  place_seqs.py \
  -s /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/rep-seqs/dna-sequences.fasta \
  -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/out.tre \
  -p 20
 
  hsp.py -i 16S \
  -t /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/out.tre \
  -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/marker_predicted_and_nsti.tsv.gz \
  -p 20

  hsp.py -i COG \
  -t /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/out.tre \
  -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/COG_predicted.tsv.gz \
  -p 20

  metagenome_pipeline.py -i /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/table/feature-table.biom \
  -m /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/marker_predicted_and_nsti.tsv.gz \
  -f /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/COG_predicted.tsv.gz \
  -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/COG \
  --strat_out

  add_descriptions.py -i /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/COG/pred_metagenome_unstrat.tsv \
  -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/COG/pred_metagenome_unstrat_descrip.tsv \
  -m COG

  hsp.py -i EC \
  -t /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/out.tre \
  -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/EC_predicted.tsv.gz \
  -p 20

  metagenome_pipeline.py -i /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/table/feature-table.biom \
  -m /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/marker_predicted_and_nsti.tsv.gz \
  -f /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/EC_predicted.tsv.gz \
  -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/EC \
  --strat_out

  add_descriptions.py -i /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/EC/pred_metagenome_unstrat.tsv \
  -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/EC/pred_metagenome_unstrat_descrip.tsv \
  -m EC

  hsp.py -i KO -t /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/out.tre \
  -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/KO_predicted.tsv.gz \
  -p 20

  metagenome_pipeline.py -i /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/table/feature-table.biom \
  -m /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/marker_predicted_and_nsti.tsv.gz \
  -f /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/KO_predicted.tsv.gz \
  -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/KO \
  --strat_out

  add_descriptions.py -i /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/KO/pred_metagenome_unstrat.tsv \
  -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/KO/pred_metagenome_unstrat_descrip.tsv \
  -m KO

  pathway_pipeline.py -i /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/EC/pred_metagenome_strat.tsv \
  -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/metacyc \
  -p 20

  add_descriptions.py -i /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/metacyc/path_abun_unstrat.tsv \
  -m METACYC \
  -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/metacyc/path_abun_unstrat_descrip.tsv
 
  pathway_pipeline.py -i /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/KO/pred_metagenome_strat.tsv \
  -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/kegg_module \
  --no_regroup \
  --map /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/ref/pathway_mapfiles/KEGG_modules_to_KO.tsv \
  -p 20

  add_descriptions.py -i /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/kegg_module/path_abun_unstrat.tsv \
  --custom_map_table /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/ref/description_mapfiles/KEGG_modules_info_adj.tsv \
  -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/kegg_module/path_abun_unstrat_descrip.tsv

  pathway_pipeline.py -i /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/KO/pred_metagenome_strat.tsv \
  -o //work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/kegg_path/all \
  --no_regroup \
  --map /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/ref/pathway_mapfiles/KEGG_pathways_to_KO.tsv \
  -p 20

  add_descriptions.py -i /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/kegg_path/all/path_abun_unstrat.tsv \
  --custom_map_table /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/ref/description_mapfiles/KEGG_pathways_info.tsv.gz \
  -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/kegg_path/all/path_abun_unstrat_descrip.tsv

  pathway_pipeline.py -i /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/KO/pred_metagenome_strat.tsv \
  -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/kegg_path/level2 \
  --no_regroup \
  --map /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/ref/kegg/kegg_path_ko_l2.txt \
  -p 20

  add_descriptions.py -i /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/kegg_path/level2/path_abun_unstrat.tsv \
  --custom_map_table /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/ref/kegg/kegg_path_description.txt \
  -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/kegg_path/level2/path_abun_unstrat_descrip.tsv

  pathway_pipeline.py -i /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/KO/pred_metagenome_strat.tsv \
  -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/kegg_path/level3 \
  --no_regroup \
  --map /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/ref/kegg/kegg_path_ko_l3.txt \
  -p 20

  add_descriptions.py -i /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/kegg_path/level3/path_abun_unstrat.tsv \
  --custom_map_table /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/ref/kegg/kegg_path_description.txt \
  -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/picrust2/kegg_path/level3/path_abun_unstrat_descrip.tsv

#--------------------------------------------------------------------------------------------------------------------------------------------------------------
source activate lefse
install.packages("dplyr",lib="/mnt/d/home/ruijuan/.conda/envs/lefse/lib/R/library")
.libPaths(c("/mnt/d/home/ruijuan/.conda/envs/lefse/lib/R/library"))
install.packages("",lib="",repos= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
export R_LIBS=/mnt/d/home/ruijuan/.conda/envs/lefse/lib/R/library
install.packages("qiimer", repos= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"), lib="/mnt/d/home/ruijuan/.conda/envs/lefse/lib/R/library")

mkdir -p /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level2
~/.conda/envs/16S-2019.04/bin/Rscript /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/lefsein_taxa.r -i /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qzv/taxa/level-2.csv -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level2/lefse.in --class Condition -u index
lefse-format_input.py /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level2/lefse.in /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level2/lefse_format.in -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level2/lefse_format.in /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level2/lefse_run.res -a 0.05 -w 0.05 -l 2 -y 0
lefse-plot_res.py /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level2/lefse_run.res /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level2/lefse_run.png --feature_font_size 8 --width 10 --format png
lefse-plot_cladogram.py /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level2/lefse_run.res /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level2/lefse_cladogram.png --format png

mkdir -p /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level6
~/.conda/envs/16S-2019.04/bin/Rscript /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/lefsein_taxa.r -i /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qzv/taxa/level-6.csv -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level6/lefse.in --class Condition -u index
lefse-format_input.py /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level6/lefse.in /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level6/lefse_format.in -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level6/lefse_format.in /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level6/lefse_run.res -a 0.05 -w 0.05 -l 2 -y 0
lefse-plot_res.py /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level6/lefse_run.res /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level6/lefse_run.png --feature_font_size 8 --width 10 --format png
lefse-plot_cladogram.py /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level6/lefse_run.res /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level6/lefse_cladogram.png --format png

mkdir -p /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level7
~/.conda/envs/16S-2019.04/bin/Rscript /mnt/d/work/workspace/ruijuan/script/git/pipeline16S/lefsein_taxa.r -i /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qzv/taxa/level-7.csv -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level7/lefse.in --class Condition -u index
lefse-format_input.py /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level7/lefse.in /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level7/lefse_format.in -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level7/lefse_format.in /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level7/lefse_run.res -a 0.05 -w 0.05 -l 2 -y 0
lefse-plot_res.py /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level7/lefse_run.res /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level7/lefse_run.png --feature_font_size 8 --width 10 --format png
lefse-plot_cladogram.py /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level7/lefse_run.res /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/diff/lefse/taxa/level7/lefse_cladogram.png --format png



cd /work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/blastn/db
/work/program/current/ncbi-blast/bin/makeblastdb -in /work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/blastn/db/feature.fasta -input_type fasta -dbtype 'nucl' -parse_seqids -out feature -hash_index

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/blastn/targetEF/efnoprimer.fasta -task blastn -db /work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/blastn/db/feature -out /work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/blastn/targetEF/output/MNH05168-v3v4-format.out -max_target_seqs 10  -num_threads 16 

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/blastn/targetEF/efnoprimer.fasta -task blastn -db /work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/blastn/db/feature -out /work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/blastn/targetEF/output/MNH05168-v3v4-format7.out -max_target_seqs 10 -num_threads 16 -outfmt 7 

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/blastn/targetEF/input/efeature3.fasta -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/blastn/targetEF/output/feature3-format7.out -max_target_seqs 20 -num_threads 16 -outfmt 7 

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "NGMM01000009|MIEK01000071|JXLB01000051|CP003504BCQ|B01000108|AJAN01000023"
AJAN01000023    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus villorum;ATCC 700913(T)
BCQB01000108    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus durans;NBRC 100479(T)
JXLB01000051    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus ratti;DSM 15687(T)
MIEK01000071    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus rivorum;LMG 25899(T)
NGMM01000009    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;NGMM_s;9E7_DIV0242
cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "KF621060|JXLE01000039|GU983697|AJAT01000017|AY321376|AJ301830|LC127059|JXKG01000044|AHYR01000005"
AHYR01000005    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus dispar;ATCC 51266(T)
AJ301830        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus faecium;LMG 11423(T)
AJAT01000017    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus phoeniculicola;ATCC BAA-412(T)
AY321376        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus sanguinicola;SS-1729(T)
GU983697        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus lactis;BT159(T)
JXKG01000044    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus canintestini;DSM 21207(T)
JXLE01000039    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus thailandicus;DSM 21767(T)
KF621060        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus xinjiangensis;48(T)
LC127059        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus saigonensis;VE80(T)

误: package or namespace load failed for ‘ggplot2’:
 package ‘scales’ was installed by an R version with different internals; it needs to be reinstalled for use with this R version


~/.conda/envs/16S-2019.04/bin/Rscript
source deactivate 16S-2019.04
[ruijuan@localhost level3]$
mkdir -p /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/R/beta/Condition
python3 /work/workspace/ruijuan/script/4_16Sseq/1_rplot/lib/betadiv.py \
-m /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/input/meta35 \
-g Condition \
/work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/core-metrics-results \
/work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/R/beta/Condition


python3 /work/workspace/ruijuan/script/4_16Sseq/1_rplot/lib/betadiv.py \
-m /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/input/metadata \
-g Location \
/work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/core-metrics-results \
/work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/R/beta/Location
/mnt/d/home/ruijuan/.conda/envs/16S-2019.04/bin/Rscript /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/R/beta/Location/p/bboxplot.r

python3 /work/workspace/ruijuan/script/4_16Sseq/1_rplot/lib/betadiv.py \
-m /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/input/metadata \
-g Type \
/work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/core-metrics-results \
/work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/R/beta/Type
/mnt/d/home/ruijuan/.conda/envs/16S-2019.04/bin/Rscript /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/R/beta/Type/p/bboxplot.r

PCOA
mkdir -p /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/R/pcoa/Condition
 python3 /work/workspace/ruijuan/script/4_16Sseq/1_rplot/lib/pcoa.py \
-m /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/input/meta35 \
-g Condition \
/work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/core-metrics-results \
/work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/R/pcoa/Condition
#/mnt/d/home/ruijuan/.conda/envs/16S-2019.04/bin/Rscript /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/R/pcoa/Time/p/pcoaplot.r

 python3 /work/workspace/ruijuan/script/4_16Sseq/1_rplot/lib/pcoa.py \
-m /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/input/metadata \
-g Location \
/work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/core-metrics-results \
/work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/R/pcoa/Location

/mnt/d/home/ruijuan/.conda/envs/16S-2019.04/bin/Rscript /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/R/pcoa/Location/p/pcoaplot.r

python3 /work/workspace/ruijuan/script/4_16Sseq/1_rplot/lib/pcoa.py \
-m /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/input/metadata \
-g Type \
/work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/qza/core-metrics-results \
/work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/R/pcoa/Type

/mnt/d/home/ruijuan/.conda/envs/16S-2019.04/bin/Rscript /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/R/pcoa/Type/p/pcoaplot.r

qiime tools export --input-path bray_curtis_distance_matrix.qza --output-path pcoa_bray_curtis

192.168.2.20
import numpy as np
a = np.loadtxt('/work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/taxonomyBoxplot/input/speciesmatrix.txt')
b= a/a.sum(axis=1)[:,np.newaxis]
np.savetxt('/work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/taxonomyBoxplot/input/speciesfre_matrix.txt',b,fmt='%f',delimiter="\t",newline="\n")

a = np.loadtxt('/work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/taxonomyBoxplot/input/genusmatrix.txt')
b= a/a.sum(axis=1)[:,np.newaxis]
np.savetxt('/work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/taxonomyBoxplot/input/genusfre_matrix.txt',b,fmt='%f',delimiter="\t",newline="\n")

library(ggplot2)
library(vegan)
library(ggpubr)
library('ggbiplot')
library(reshape2)
library('RColorBrewer')
library('ggsci')
library(dplyr)
library(docopt)
library(getopt)
data = read.table("/work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/R/taxa/input/genusCondition.txt",header = T,sep="\t")
compare_means(log2Abundance ~ Group, data=data, group.by = "Genus")
d <- ggboxplot(data, x="Group", y="log2Abundance", fill = "Group",palette = c("#00AFBB", "#FC4E07"), facet.by = "Genus",ylim = c(0, 16), nrow=2) + 
scale_x_discrete(name="", breaks = c("Control", "Test"),labels=c("Control" = "Control", "Test" = "Test")) +
stat_compare_means(aes(label = ..p.signif..),label.y = 15.5) +
guides(fill=FALSE)
dif <- ggarrange(d, ncol=1,nrow=2,labels=c("A"))
ggsave("/work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/R/taxa/output/genusCondition_20200730.jpg", dif, width = 10, height = 12)


data = read.table("/work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/R/taxa/input/speciesCondition.txt",header = T,sep="\t")
compare_means(log2Abundance ~ Group, data=data, group.by = "Genus")
d <- ggboxplot(data, x="Group", y="log2Abundance", fill = "Group",palette = c("#00AFBB", "#FC4E07"), facet.by = "Genus",ylim = c(0, 16), nrow=2) + 
scale_x_discrete(name="", breaks = c("Control", "Test"),labels=c("Control" = "Control", "Test" = "Test")) +
stat_compare_means(aes(label = ..p.signif..),label.y = 15.5) +
guides(fill=FALSE)
dif <- ggarrange(d, ncol=1,nrow=1,labels=c("A"))
ggsave("/work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/R/taxa/output/speciesCondition_20200730.jpg", dif, width = 12, height = 6)


heatmap
install.packages("caTools",repos= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
install.packages("gplots",repos= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
install.packages("ggplot2",repos= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"),lib="/mnt/d/home/ruijuan/.conda/envs/16S-2019.04/lib/R/library")
install.packages("reshape",repos= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
rm(list = ls())
export R_LIBS=/mnt/d/home/ruijuan/R/x86_64-redhat-linux-gnu-library/3.5
python /work/workspace/zhengjiao/script/profileTest/corr_heatmap.py -i /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/corr/input/l6fre.txt -p /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/corr/input/clinic -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/corr/output/table1only -x 15 -y 15 -f1 0.35

python /work/workspace/zhengjiao/script/profileTest/corr_heatmap.py -i /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/corr/input/l6fre.txt -c /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/corr/input/clinic -a cross -o /work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/corr/output/cross -x 15 -y 15 -f1 0.35

data = read.table("/work/workspace/ruijuan/project/3_16S/baoxia35sam20200729/result_all/corr/input/spegroup.txt",header=T,row.names=1,sep="\t",check.names=F,comment.char="")   
method = 'spearman'

dat = as.matrix(data)
var_names = colnames(dat)
cor_result = as.data.frame(matrix(data=NA,nrow=length(var_names)-1,ncol=1))
p_result = cor_result

coli = 'Group'
for (rowi in var_names){
  if (rowi==coli){
    cor_result[rowi,coli] = 0
    p_result[rowi,coli] = 1
  }else{
    tem_cor = cor.test(as.vector(dat[rowi]),as.vector(dat[coli]),method=method)
    cor_result[rowi,coli] = tem_cor$estimate
    p_result[rowi,coli] = tem_cor$p.value
  }
}


/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/blastn/targetEF/input/efeature3.fasta -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/blastn/targetEF/output/feature3-format7.out -max_target_seqs 20 -num_threads 16 -outfmt 7 


/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/rename_levels_name.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm/genus_sam5_abun0.005 -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm


```
source activate /work/workspace/zhurj/bin/miniconda3/envs/py2
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm/name_filter_genus_sam5_abun0.005 /work/workspace/zhurj/project/3_16S/anhuilung20201229/input/meta /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm/lefse_in
lefse-format_input.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm/lefse_in /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm/lefse_format.in -c 1 -u 2 -f r -o 1000000
run_lefse.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm/lefse_format.in /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm/lefse_run.res -a 0.05 -w 0.05 -l 2 -y 0
lefse-plot_res.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm/lefse_run.res  /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm/lefse_run.pdf --feature_font_size 8 --width 15 --format pdf
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse_fitler.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm/lefse_run.res /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm/lefse_run_filter.res
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/toheadmapin.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm/name_filter_genus_sam5_abun0.005 /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm/lefse_run_filter.res /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm/dif_levels_all_data
source deactivate /work/workspace/zhurj/bin/miniconda3/envs/py2


import os
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import subprocess

ifile = '/work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm/name_filter_genus_sam5_abun0.005'
metaf = '/work/workspace/zhurj/project/3_16S/anhuilung20201229/input/meta'
odir = '/work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm/test'
pvalue_Anova_test = 0.05
pvalue_Wilcoxon_test = 0.05
LDA_score_cutoff = 2.0
feature_font_size = 8
width = 8
# graph_format options {png,svg,pdf}
graph_format = 'pdf'
programf = os.path.join(odir,'lefse_pro.py')


con = program = ''

ofp = open(programf,'w')
con = 'source activate /work/workspace/zhurj/bin/miniconda3/envs/py2\n'
ofp.write(con)
# lefse_in
infile = ifile
program = '/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/lefse_in.py'
ofile = os.path.join(odir,'lefse_in')
con = "{} {} {} {} \n".format(program,infile,metaf,ofile)
ofp.write(con)

# lefse_format.in
program = 'lefse-format_input.py'
ifile = ofile
ofile = os.path.join(odir,'lefse_format.in')
con = "{} {} {} -c 1 -u 2 -f r -o 1000000 \n".format(program,ifile,ofile)
ofp.write(con)

# lefse_run.res
program = 'run_lefse.py'
ifile = ofile
ofile = os.path.join(odir,'lefse_run.res')
con = "{} {} {} -a {} -w {} -l {} -y 0 \n".format(program,ifile,ofile,pvalue_Anova_test,pvalue_Wilcoxon_test,LDA_score_cutoff)
ofp.write(con)

# lefse_run.pdf
program = 'lefse-plot_res.py'
ifile = ofile
ofile = os.path.join(odir,'lefse-plot_res.pdf')
con = "{} {} {}  --feature_font_size {} --width {} --format {} \n".format(program,ifile,ofile,feature_font_size,width,graph_format)
ofp.write(con)

# lefse_run_filter.res
program = '/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse_fitler.py'
ifile = os.path.join(odir,'lefse_run.res')
ofile = os.path.join(odir,'lefse_run_filter.res')
con = "{} {} {} \n".format(program,ifile,ofile)
ofp.write(con)

# lefse_run_filter.res
program = '/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/toheadmapin.py'
ifile = ofile
ofile = os.path.join(odir,'dif_levels_all_data')
con = "{} {} {} {}\n".format(program,infile,ifile,ofile)
ofp.write(con)

con = 'source deactivate /work/workspace/zhurj/bin/miniconda3/envs/py2\n'
ofp.write(con)
ofp.close()

con = 'bash {}\n'.format(programf)
subprocess.run(con,shell=True)

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm/name_filter_genus_sam5_abun0.005 --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/input/meta -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm/test -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf



cd /work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/taxoreana20210111
# genus
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/qiime2_levels_to_filtered_frequency.py -i /work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/taxoreana20210111/input/level6 -o /work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/taxoreana20210111/genus --sample_num 35
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/rename_levels_name.py -i /work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/taxoreana20210111/genus/genus_sam3_abun0.005 -o /work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/taxoreana20210111/genus

# family
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/qiime2_levels_to_filtered_frequency.py -i /work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/taxoreana20210111/input/level5 -o /work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/taxoreana20210111/family --sample_num 35 --remove_num_col 0 --level family

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/rename_levels_name.py -i /work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/taxoreana20210111/family/family_sam3_abun0.005 -o /work/workspace/zhurj/project/3_16S/sam35baoxiaEF20200730/taxoreana20210111/family


from Bio import SeqIO

records = SeqIO.parse("seq.txt", "tab")
count = SeqIO.write(records, "seq.fasta", "fasta")
print("Converted %i records" % count)

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/22_SeqIO_format_convert/format_converter.py -i /work/workspace/zhurj/project/14_coworker/baojia/blastseq4_20210113/seq.txt -o /work/workspace/zhurj/project/14_coworker/baojia/blastseq4_20210113/seq.fasta
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/22_SeqIO_format_convert/format_converter.py -i /work/workspace/zhurj/project/14_coworker/baojia/blastseq4_20210113/out --in_format fasta -o /work/workspace/zhurj/project/14_coworker/baojia/blastseq4_20210113/out.txt --out_format tab
cd /work/workspace/zhurj/project/14_coworker/baojia/blastseq4_20210113/ref
makeblastdb -dbtype nucl -in /work/workspace/zhurj/project/14_coworker/baojia/blastseq4_20210113/ref/seq.fasta -out /work/workspace/zhurj/project/14_coworker/baojia/blastseq4_20210113/ref/seq

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/14_coworker/baojia/blastseq4_20210113/seq.fasta -task blastn -db /work/workspace/zhurj/project/14_coworker/baojia/blastseq4_20210113/ref/seq -out /work/workspace/zhurj/project/14_coworker/baojia/blastseq4_20210113/seq.out -num_threads 16 -outfmt 7 

/work/program/current/ncbi-blast/bin/blastn -query /work/workspace/zhurj/project/14_coworker/baojia/blastseq4_20210113/seq.fasta -task blastn -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/14_coworker/baojia/blastseq4_20210113/in_16S.out -max_target_seqs 5 -num_threads 16 -outfmt 7

cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "GU983697|AJ301830|BCQB01000108|CP003504|JXLE01000039|AY321376"

AJ301830        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus faecium;LMG 11423(T)
AY321376        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus sanguinicola;SS-1729(T)
BCQB01000108    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus durans;NBRC 100479(T)
CP003504        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus hirae;ATCC 9790(T)
GU983697        Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus lactis;BT159(T)
JXLE01000039    Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus thailandicus;DSM 21767(T)

cd /work/workspace/zhurj/project/14_coworker/baojia/blastseq4_20210113
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/muscle -in seq.fasta -clwstrict -out musleclw.clw
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/clustalo -i musleclw.clw --percent-id --distmat-out=pim.txt --full --force

infile = '/work/workspace/zhurj/project/14_coworker/wanglin/genus_num'
ofile = '/work/workspace/zhurj/project/14_coworker/wanglin/unique_genus_count'
df = pd.read_csv(infile,sep='\t',header=None,index_col=None)
df.columns = ['genus','count']
new_df = df.groupby(['genus']).sum()
new_df.to_csv(ofile,sep='\t',header=True,index=True)


/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/rename_levels_name.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm/genus_sam5_abun0.005 -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_norm
# 先不用重命名的规则，直接用filter后的结果，先手动修改名称
#/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/rename_levels_name.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/genus_sam3_abun0.0001 -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/27_frequency2count/frequency2count.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/used_genus_frequency_sam3_abun0.005 -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/used_genus_count50000
/work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/used_genus_count50000 header 改为 “#OTU ID”


co-abundance
genus
#/work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_count_100000
#'fastspar':'/work/workspace/liangzj/program/Miniconda3/envs/local/bin/fastspar',
#'fastspar_bootstrap':'/work/workspace/liangzj/program/Miniconda3/envs/local/bin/fastspar_bootstrap',
#'fastspar_pvalues':'/work/workspace/liangzj/program/Miniconda3/envs/local/bin/fastspar_pvalues'

mkdir /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/bootstrap_counts
mkdir /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/bootstrap_correlation
/work/workspace/liangzj/program/Miniconda3/envs/local/bin/fastspar -c /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/used_genus_count50000 -r /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/cor_sparcc.fastspar.txt -a /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/cov_sparcc.fastspar.txt -t 20 

/work/workspace/liangzj/program/Miniconda3/envs/local/bin/fastspar_bootstrap --otu_table /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/used_genus_count50000 --n 100 --prefix /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/bootstrap_counts/bootstrap -t 20


parallel /work/workspace/liangzj/program/Miniconda3/envs/local/bin/fastspar --otu_table {} --correlation /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/bootstrap_correlation/cor_{/} --covariance /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/bootstrap_correlation/cov_{/} -i 50 -t 20 ::: /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/bootstrap_counts/*

/work/workspace/liangzj/program/Miniconda3/envs/local/bin/fastspar_pvalues --otu_table /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/used_genus_count50000 --correlation /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/cor_sparcc.fastspar.txt --prefix /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/bootstrap_correlation/cor_bootstrap --permutations 100 --outfile /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/pvalues.tsv -t 20

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/28_count2fastspar/count2fastspar.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/used_genus_count50000 -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/cor_pvalue_to_edge_node_list.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/cor_sparcc.fastspar.txt -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226 --pfile /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/pvalues.tsv --p_cutoff 0.05 --cor_cutoff 0.6

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/used_genus_frequency_sam3_abun0.005 --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/input/meta -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/lefse -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf



R
# REFERENCE
# https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
-------------------------------------Example---------------------------------------
library(corrplot)
## corrplot 0.84 loaded
M <- cor(mtcars)
corrplot(M, method = "circle")

res1 <- cor.mtest(mtcars, conf.level = .95)
corrplot(M, p.mat = res1$p, insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "white")

-------------------------------------Example---------------------------------------
library(corrplot)
infile <- '/work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/correlation/genus_clinical'
ofile <- '/work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/correlation/old_dif_genus_clinical_corPvalue.pdf'
data <- read.table(infile,header = TRUE, sep = "\t",row.names= 1)
M <- cor(data, use="complete.obs",method = "spearman")
diag(M) = NA
# corrplot(M, method = "circle")
res1 <- cor.mtest(data, conf.level = .95)
pdf(ofile,width=8,height=8)
corrplot(M, p.mat = res1$p, insig = "label_sig",  type = "upper", sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "white",na.label = "-")
dev.off()

library(corrplot)
infile <- '/work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/correlation/genus_abun_clinic_20210301'
ofile <- '/work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/correlation/genus_clinical_corPvalue_20210301.pdf'
data <- read.table(infile,header = TRUE, sep = "\t",row.names= 1)
M <- cor(data, use="complete.obs",method = "spearman")
diag(M) = NA
# corrplot(M, method = "circle")
res1 <- cor.mtest(data, conf.level = .95)
pdf(ofile,width=8,height=8)
corrplot(M, p.mat = res1$p, insig = "label_sig",  type = "upper", sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "white",na.label = "-")
dev.off()


library(Hmisc)
library(ggplot2)
library(dplyr)

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}

plot_width <- 10
plot_height <- 4
text_size <- 12
infile <- '/work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/correlation/genus_clinical'
ofile <- '/work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/tmp/corr/test.pdf'
rowf <- '/work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/tmp/corr/row'
columnf <- '/work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/tmp/corr/column'
data <- read.table(infile,header = TRUE, sep = "\t",row.names= 1)
row_df <- read.table(rowf)
column_df <- read.table(columnf)

res2 <- rcorr(as.matrix(data),type="spearman")
format_df <- flattenCorrMatrix(res2$r, res2$P)
format_df$stars <- cut(format_df$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
plot.data <- subset(format_df, row %in% row_df$V1 & column %in% column_df$V1)
#data <- subset(pre_data, group %in% unique_group)


# Plot everything
my_corrploat <- ggplot(aes(x=column, y=row, fill=cor), data=plot.data) +
    geom_tile() + 
    scale_fill_gradient2(low="#D7191C", mid="white", high="#2C7BB6",breaks=c(-1,-0.5,0,0.5,1), limits=c(-1,1)) + 
    #   geom_text(aes(label=stars, color=cor), size=8) + scale_colour_gradient(low="grey30", high="white", guide="none") +
    geom_text(aes(label=stars), color="black", size=5) + 
    labs(y=NULL, x=NULL, fill="Spearman Correlation") + 
    #geom_vline(xintercept=1.5, size=1.5, color="grey50") + 
    theme_bw() + 
    theme(axis.text.x=element_text(angle = -45, hjust = 0),text=element_text(size=12))


ggsave(ofile,my_corrploat,width = plot_width,height = plot_height)


/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/33_corr_heatmap/corr_heatmap_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/correlation/genus_clinical -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/tmp/corr/test.pdf -r /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/tmp/corr/row -c /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/tmp/corr/column --method spearman --plot_width 10 --plot_height 5 --text_size 12

lefse
phylum
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/qiime2_levels_to_filtered_frequency.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/phylum_tab -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/phylum/input --cutoff_sample 3 --abundance_cutoff 0.005 --level phylum

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/phylum_CAG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy_meta_CAG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy/lefse/phylum/CAG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
# no significant phylum
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/phylum_CBG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy_meta_CBG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy/lefse/phylum/CBG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
# no significant phylum
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/phylum_CAG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/side_effect_meta_CAG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/phylum/CAG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
# 0 significant phylum
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/phylum_CBG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/side_effect_meta_CBG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/phylum/CBG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
# 0 significant phylum

genus
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/genus_CAG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy_meta_CAG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy/lefse/genus/CAG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
# 6 genera are significant different

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/genus_CBG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy_meta_CBG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy/lefse/genus/CBG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
# NO sinificant genus -- CBG

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/genus_CAG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/side_effect_meta_CAG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/genus/CAG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
# 6 significant phylum

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/genus_CBG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/side_effect_meta_CBG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/genus/CBG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
# 2 significant phylum

species
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/qiime2_levels_to_filtered_frequency.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/species_tab -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/species/input --cutoff_sample 3 --abundance_cutoff 0.005 --level species

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/species_CAG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy_meta_CAG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy/lefse/species/CAG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
# 4 species are significant different

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/species_CBG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy_meta_CBG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy/lefse/species/CBG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
# 3 sinificant species -- CBG

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/species_CAG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/side_effect_meta_CAG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/species/CAG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
# 6 significant phylum

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/species_CBG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/side_effect_meta_CBG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/species/CBG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
# 2 significant phylum

kegg
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/kegg_CAG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy_meta_CAG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy/lefse/kegg/CAG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
# 5 species are significant different

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/kegg_CBG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy_meta_CBG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy/lefse/kegg/CBG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
# 2 sinificant species -- CBG

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/kegg_CAG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/side_effect_meta_CAG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/kegg/CAG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
# 0 significant phylum

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/kegg_CBG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/side_effect_meta_CBG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/kegg/CBG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
# 0 significant phylum


metacyc
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/metacyc_CAG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy_meta_CAG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy/lefse/metacyc/CAG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
# 4 species are significant different

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/metacyc_CBG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy_meta_CBG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy/lefse/metacyc/CBG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
# 3 sinificant species -- CBG

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/metacyc_CAG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/side_effect_meta_CAG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/metacyc/CAG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
# 5 significant phylum

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/metacyc_CBG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/side_effect_meta_CBG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/metacyc/CBG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
# 5 significant phylum


beta diversity
---------side_effect
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/vedist_to_dist.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/genus_CAG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/beta_diversity/genus/CAG/genus_braycurtis_dist -m bray
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/distformat.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/beta_diversity/genus/CAG/genus_braycurtis_dist /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/side_effect_meta_CAG  /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/beta_diversity/genus/CAG/genus_beta_diversity_in
/work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/beta_diversity_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/beta_diversity/genus/CAG/genus_beta_diversity_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/beta_diversity/genus/CAG -n genus_beta -t Bray_Curtis_Distance --plot_width 4 --plot_height 4 --y_max 1.5
/work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript /work/workspace/zhurj/script/1_tools/29_beta_diversity/betaDiversity_defineComparegroup_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/beta_diversity/genus/CAG/genus_beta_diversity_in -g /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/side_effect_compare -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/beta_diversity/genus/CAG -n genus_beta_two -t Bray_Curtis_Distance --plot_width 3 --plot_height 3 --y_max 1.2

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/vedist_to_dist.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/genus_CBG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/beta_diversity/genus/CBG/genus_braycurtis_dist -m bray
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/distformat.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/beta_diversity/genus/CBG/genus_braycurtis_dist /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/side_effect_meta_CBG  /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/beta_diversity/genus/CBG/genus_beta_diversity_in
/work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/beta_diversity_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/beta_diversity/genus/CBG/genus_beta_diversity_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/beta_diversity/genus/CBG -n genus_beta -t Bray_Curtis_Distance --plot_width 4 --plot_height 4 --y_max 1.5

/work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript /work/workspace/zhurj/script/1_tools/29_beta_diversity/betaDiversity_defineComparegroup_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/beta_diversity/genus/CBG/genus_beta_diversity_in -g /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/side_effect_compare -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/beta_diversity/genus/CBG -n genus_beta_two -t Bray_Curtis_Distance --plot_width 3 --plot_height 3 --y_max 1.2


---------efficacy
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/vedist_to_dist.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/genus_CAG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy/beta_diversity/genus/CAG/genus_braycurtis_dist -m bray
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/distformat.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy/beta_diversity/genus/CAG/genus_braycurtis_dist /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy_meta_CAG  /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy/beta_diversity/genus/CAG/genus_beta_diversity_in
/work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript /work/workspace/zhurj/script/1_tools/29_beta_diversity/betaDiversity_defineComparegroup_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy/beta_diversity/genus/CAG/genus_beta_diversity_in -g /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy_compare -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy/beta_diversity/genus/CAG -n genus_beta -t Bray_Curtis_Distance --plot_width 4 --plot_height 4 --y_max 1.3

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/vedist_to_dist.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/genus_CBG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy/beta_diversity/genus/CBG/genus_braycurtis_dist -m bray
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/distformat.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy/beta_diversity/genus/CBG/genus_braycurtis_dist /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy_meta_CBG  /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy/beta_diversity/genus/CBG/genus_beta_diversity_in
/work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript /work/workspace/zhurj/script/1_tools/29_beta_diversity/betaDiversity_defineComparegroup_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy/beta_diversity/genus/CBG/genus_beta_diversity_in -g /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy_compare -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy/beta_diversity/genus/CBG -n genus_beta -t Bray_Curtis_Distance --plot_width 4 --plot_height 4 --y_max 1.4

/work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript /work/workspace/zhurj/script/1_tools/29_beta_diversity/betaDiversity_defineComparegroup_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/beta_diversity/genus/CBG/genus_beta_diversity_in -g /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/side_effect_compare -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/beta_diversity/genus/CBG -n genus_beta_two -t Bray_Curtis_Distance --plot_width 3 --plot_height 3 --y_max 1.2


****************************************************************************************************************
测试
ref <- '/work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/side_effect_compare'
my_comparisons <- strsplit(readLines(ref), "[[:space:]]+")
allgroup_df <- read.table(ref,header=F,sep='\t')
unique_group <- unique(sort(as.vector(as.matrix(allgroup_df))))

infile <- '/work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/beta_diversity/genus/CBG/genus_beta_diversity_in'
pre_data <- read.table(infile,sep='\t',header=T)
data <- subset(pre_data, group %in% unique_group)
****************************************************************************************************************

diff bar graph creation
library(ggplot2)
library(dplyr)

infile <- '/work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/side_effect_metacyc_CAG'
ofile <- '/work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/side_effect_metacyc_CAG.jpg'
label_size <- 2.5
x_axis_min = -5
x_axis_max = 5
plot_width = 8
plot_height = 2
legend_x_axis = 0.9
legend_y_axis = 1


df <- read.table(infile,header = T,sep="\t")
mycolor = c("#00AFBB", "#FC4E07","#E7B800","#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
df <- df %>%
  mutate(
    item = factor(item, levels = item[order(value, decreasing = FALSE)]),
    label_y = ifelse(value < 0, 0.2, -0.2),
    label_hjust = ifelse(value < 0, 0, 1)
  )
unique_group = as.vector(unique(df$class))
color_num = length(unique_group)
select_color = mycolor[1:color_num]
names(select_color) = unique_group

my_plot <- ggplot(df, aes(x = item, y = value, fill = class)) +
  geom_bar(stat = "identity") +
  #geom_bar(stat = "identity", col = "black") +
  geom_text(aes(y = label_y, label = item, hjust = label_hjust),size=label_size) +
  coord_flip() +
  scale_fill_manual(values = select_color) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = c(legend_x_axis,legend_y_axis),
        legend.justification = c(legend_x_axis,legend_y_axis),
        legend.direction='vertical',
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey80", linetype = "dashed"),
        panel.grid.minor.x = element_blank()) +
  scale_y_continuous(expression(log[10](italic("LDA score"))),
                     breaks = x_axis_min:x_axis_max, limits = c(x_axis_min, x_axis_max))

ggsave(ofile,my_plot,width = plot_width,height = plot_height)

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/31_diff_barplot/dif_barplot_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/side_effect_metacyc_CAG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/side_effect_metacyc_CAG.jpg --x_axis_min -5 --x_axis_max 5 --label_size 2.5 --legend_x 0.9 --legend_y 1 --plot_width 8 --plot_height 2
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/31_diff_barplot/dif_barplot_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/side_effect_metacyc_CBG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/side_effect_metacyc_CBG.jpg --x_axis_min -6 --x_axis_max 5 --label_size 2.4 --legend_x 0.9 --legend_y 1 --plot_width 8 --plot_height 2

boxplot 
#############################################################################################################

library('ggpubr')

infile <- "/work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/dif_genus_graph_log2_in"
ofile <- "/work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/side_effect_dif_genus.jpg"
comparefile <- "/work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/side_effect_compare"
title <- "Boxplot of genus"
plot_width <- 10
plot_height <- 10
y_min <- 0
y_step <- 2
y_max <- 14
num_col <- 3

my_comparisons <- strsplit(readLines(comparefile), "[[:space:]]+")
allgroup_df <- read.table(comparefile,header=F,sep='\t')
unique_group <- unique(sort(as.vector(as.matrix(allgroup_df))))
mycolor = c("#00AFBB", "#FC4E07","#E7B800","#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
color_num = length(unique_group)
select_color = mycolor[1:color_num]

Data = read.table(infile,header=T,sep='\t')
p <- ggboxplot(Data, x = "Group", y = "Abundance", fill = "Group", facet.by = "Species", ylim = c(y_min, y_max),ncol=num_col) + 
ggtitle(title) +
theme(legend.position="none",plot.title = element_text(face="bold.italic", hjust=0.5),)+
scale_y_continuous(name="Relative Abundance", breaks=seq(y_min,y_max,y_step)) +
scale_x_discrete(name="") +
stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", label = "p.signif") +
guides(fill=FALSE)
ggsave(ofile, p, width = plot_width, height =plot_height)

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/toheadmapin.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/genus_CAG /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/genus/CAG/lefse_run_filter.res /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/genus/CAG/dif_genus_in
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse2boxplot.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/genus/CAG/dif_genus_in  /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/side_effect_meta_CAG /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/genus/CAG/dif_genus_graph_in
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/32_relativeAbun2log2/relative_abundance2log2multiple.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/genus/CAG/dif_genus_graph_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/dif_genus_CAG_log2_in --column 3 --with_header --time 10000
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/30_diff_boxplot/dif_boxplot_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/dif_genus_CAG_log2_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/side_effect_dif_genus_CAG.jpg -g /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/side_effect_compare -t Different_genus --plot_width 10 --plot_height 10 --y_min 0 --y_step 2 --y_max 14 --ncol 3


/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/toheadmapin.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/genus_CBG /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/genus/CBG/lefse_run_filter.res /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/genus/CBG/dif_genus_in
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse2boxplot.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/genus/CBG/dif_genus_in  /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/side_effect_meta_CBG /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/genus/CBG/dif_genus_graph_in
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/32_relativeAbun2log2/relative_abundance2log2multiple.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/genus/CBG/dif_genus_graph_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/dif_genus_CBG_log2_in --column 3 --with_header --time 10000
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/30_diff_boxplot/dif_boxplot_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/dif_genus_CBG_log2_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/side_effect_dif_genus_CBG.jpg -g /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/side_effect_compare -t Different_genus --plot_width 5 --plot_height 5 --y_min 0 --y_step 2 --y_max 9 --ncol 2


/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/toheadmapin.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/species_CAG /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/species/CAG/lefse_run_filter.res /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/species/CAG/dif_species_in
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse2boxplot.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/species/CAG/dif_species_in  /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/side_effect_meta_CAG /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/species/CAG/dif_species_graph_in
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/32_relativeAbun2log2/relative_abundance2log2multiple.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/species/CAG/dif_species_graph_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/dif_species_CAG_log2_in --column 3 --with_header --time 10000
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/30_diff_boxplot/dif_boxplot_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/dif_species_CAG_log2_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/side_effect_dif_species_CAG.jpg -g /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/side_effect_compare -t Different_genus --plot_width 10 --plot_height 10 --y_min 0 --y_step 2 --y_max 14 --ncol 3


/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/toheadmapin.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/species_CBG /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/species/CBG/lefse_run_filter.res /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/species/CBG/dif_species_in
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse2boxplot.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/species/CBG/dif_species_in  /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/side_effect_meta_CBG /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/species/CBG/dif_species_graph_in
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/32_relativeAbun2log2/relative_abundance2log2multiple.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/species/CBG/dif_species_graph_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/dif_species_CBG_log2_in --column 3 --with_header --time 10000
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/30_diff_boxplot/dif_boxplot_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/dif_species_CBG_log2_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/side_effect_dif_species_CBG.jpg -g /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/side_effect_compare -t Different_species --plot_width 5 --plot_height 5 --y_min 0 --y_step 2 --y_max 9 --ncol 2


/work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/side_effect/lefse/phylum/CBG/lefse_run_filter.res


efficacy38
lefse
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/phylum_38 --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/meta38 -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/lefse/phylum -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/genus_38 --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/meta38 -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/lefse/genus -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/species_38 --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/meta38 -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/lefse/species -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/kegg_38 --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/meta38 -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/lefse/kegg -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/metacyc_38 --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/meta38 -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/lefse/metacyc -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/metacyc_38 --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/meta38 -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/lefse/metacyc001 -a 0.01 -w 0.01 -l 2.0 --feature_font_size 8 --width 8 -f pdf

beta diversity
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/vedist_to_dist.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/genus_38 -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/beta_diversity/genus_braycurtis_dist -m bray
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/distformat.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/beta_diversity/genus_braycurtis_dist /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/meta38  /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/beta_diversity/genus_beta_diversity_in
/work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/beta_diversity_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/beta_diversity/genus_beta_diversity_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/beta_diversity -n genus_beta -t Bray_Curtis_Distance --plot_width 4 --plot_height 4 --y_max 1.4


boxplot
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/toheadmapin.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/phylum_38 /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/lefse/phylum/lefse_run_filter.res /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/lefse/phylum/dif_phylum_in
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse2boxplot.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/lefse/phylum/dif_phylum_in  /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/meta38 /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/lefse/phylum/dif_phylum_graph_in
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/32_relativeAbun2log2/relative_abundance2log2multiple.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/lefse/phylum/dif_phylum_graph_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/dif_efficacy38_phylum_log2_in --column 3 --with_header --time 10000
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/30_diff_boxplot/dif_boxplot_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/dif_efficacy38_phylum_log2_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/efficacy38_dif_phylum.jpg -g /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy38_compare -t Different_phylum --plot_width 5 --plot_height 5 --y_min 0 --y_step 2 --y_max 16 --ncol 2

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/toheadmapin.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/genus_38 /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/lefse/genus/lefse_run_filter.res /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/lefse/genus/dif_genus_in
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse2boxplot.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/lefse/genus/dif_genus_in  /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/meta38 /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/lefse/genus/dif_genus_graph_in
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/32_relativeAbun2log2/relative_abundance2log2multiple.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/lefse/genus/dif_genus_graph_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/dif_efficacy38_genus_log2_in --column 3 --with_header --time 10000
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/30_diff_boxplot/dif_boxplot_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/dif_efficacy38_genus_log2_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/efficacy38_dif_genus.jpg -g /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy38_compare -t Different_genus --plot_width 10 --plot_height 5 --y_min 0 --y_step 2 --y_max 16 --ncol 4

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/toheadmapin.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/species_38 /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/lefse/species/lefse_run_filter.res /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/lefse/species/dif_species_in
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse2boxplot.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/lefse/species/dif_species_in  /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/meta38 /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/lefse/species/dif_species_graph_in
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/32_relativeAbun2log2/relative_abundance2log2multiple.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy38/lefse/species/dif_species_graph_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/dif_efficacy38_species_log2_in --column 3 --with_header --time 10000
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/30_diff_boxplot/dif_boxplot_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/dif_efficacy38_species_log2_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/efficacy38_dif_species.jpg -g /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy38_compare -t Different_species --plot_width 10 --plot_height 5 --y_min 0 --y_step 2 --y_max 16 --ncol 4


barplot
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/31_diff_barplot/dif_barplot_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy38_kegg -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/efficacy38_kegg.jpg --x_axis_min -5 --x_axis_max 5 --label_size 3.5 --legend_x 0.9 --legend_y 1 --plot_width 8 --plot_height 8 --bar_width 0.2

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/31_diff_barplot/dif_barplot_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy38_metacyc001 -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/efficacy38_metacyc001.jpg --x_axis_min -5 --x_axis_max 5 --label_size 2.5 --legend_x 0.9 --legend_y 1 --plot_width 8 --plot_height 8 --bar_width 0.1

# pvalue 设置 0.05,0.05 差异pathway 92 个
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/31_diff_barplot/dif_barplot_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy38_metacyc -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/efficacy38_metacyc.jpg --x_axis_min -7 --x_axis_max 7 --label_size 3.5 --legend_x 0.9 --legend_y 1 --plot_width 8 --plot_height 25 --bar_width 0.1



efficacy19
lefse
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_phylum_CAG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_meta_CAG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/lefse/phylum/CAG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_phylum_CBG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_meta_CBG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/lefse/phylum/CBG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_genus_CAG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_meta_CAG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/lefse/genus/CAG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_genus_CBG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_meta_CBG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/lefse/genus/CBG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_species_CAG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_meta_CAG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/lefse/species/CAG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_species_CBG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_meta_CBG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/lefse/species/CBG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_metacyc_CAG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_meta_CAG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/lefse/metacyc/CAG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_metacyc_CBG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_meta_CBG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/lefse/metacyc/CBG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_kegg_CAG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_meta_CAG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/lefse/kegg/CAG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/abun_to_lefse.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_kegg_CBG --metaf /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_meta_CBG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/lefse/kegg/CBG -a 0.05 -w 0.05 -l 2.0 --feature_font_size 8 --width 8 -f pdf


beta diversity
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/vedist_to_dist.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_genus_CAG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/beta_diversity/CAG/genus_braycurtis_dist -m bray
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/distformat.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/beta_diversity/CAG/genus_braycurtis_dist /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_meta_CAG  /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/beta_diversity/CAG/genus_beta_diversity_in
/work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/beta_diversity_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/beta_diversity/CAG/genus_beta_diversity_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/beta_diversity/CAG -n genus_beta -t Bray_Curtis_Distance --plot_width 4 --plot_height 4 --y_max 1.3

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/vedist_to_dist.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_genus_CBG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/beta_diversity/CBG/genus_braycurtis_dist -m bray
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/distformat.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/beta_diversity/CBG/genus_braycurtis_dist /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_meta_CBG  /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/beta_diversity/CBG/genus_beta_diversity_in
/work/workspace/zhurj/bin/miniconda3/envs/R4.0/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/beta_diversity_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/beta_diversity/CBG/genus_beta_diversity_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/beta_diversity/CBG -n genus_beta -t Bray_Curtis_Distance --plot_width 4 --plot_height 4 --y_max 1.3


boxplot
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/toheadmapin.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_genus_CAG /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/lefse/genus/CAG/lefse_run_filter.res /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/lefse/genus/CAG/dif_genus_in
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse2boxplot.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/lefse/genus/CAG/dif_genus_in  /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_meta_CAG /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/lefse/genus/CAG/dif_genus_graph_in
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/32_relativeAbun2log2/relative_abundance2log2multiple.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/lefse/genus/CAG/dif_genus_graph_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/dif_efficacy19_CAG_genus_log2_in --column 3 --with_header --time 10000
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/30_diff_boxplot/dif_boxplot_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/dif_efficacy19_CAG_genus_log2_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/efficacy19_CAG_dif_genus.jpg -g /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy19_compare -t Different_genus --plot_width 12 --plot_height 10 --y_min 0 --y_step 2 --y_max 12 --ncol 4



/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/2_metapro/metapipe/toheadmapin.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_species_CAG /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/lefse/species/CAG/lefse_run_filter.res /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/lefse/species/CAG/dif_species_in
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/lib/python/script/lefse2boxplot.py /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/lefse/species/CAG/dif_species_in  /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/efficacy19_meta_CAG /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/lefse/species/CAG/dif_species_graph_in
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/32_relativeAbun2log2/relative_abundance2log2multiple.py -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/efficacy19/lefse/species/CAG/dif_species_graph_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/dif_efficacy19_CAG_species_log2_in --column 3 --with_header --time 10000
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/30_diff_boxplot/dif_boxplot_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/dif_efficacy19_CAG_species_log2_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/efficacy19_CAG_dif_species.jpg -g /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy19_compare -t Different_species --plot_width 12 --plot_height 10 --y_min 0 --y_step 2 --y_max 14 --ncol 4


barplot
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/31_diff_barplot/dif_barplot_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy19_kegg_CAG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/efficacy19_kegg_CAG.jpg --x_axis_min -4 --x_axis_max 4 --label_size 3.5 --legend_x 0.95 --legend_y 1 --plot_width 8 --plot_height 2 --bar_width 0.2
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/31_diff_barplot/dif_barplot_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy19_kegg_CBG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/efficacy19_kegg_CBG.jpg --x_axis_min -4 --x_axis_max 4 --label_size 3.3 --legend_x 0.95 --legend_y 1 --plot_width 8 --plot_height 2 --bar_width 0.1
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/31_diff_barplot/dif_barplot_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy19_metacyc_CAG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/efficacy19_metacyc_CAG.jpg --x_axis_min -4 --x_axis_max 6 --label_size 3 --legend_x 0.95 --legend_y 1 --plot_width 8 --plot_height 3 --bar_width 0.2
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/31_diff_barplot/dif_barplot_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy19_metacyc_CBG -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/efficacy19_metacyc_CBG.jpg --x_axis_min -7 --x_axis_max 6 --label_size 3 --legend_x 0.95 --legend_y 1 --plot_width 8 --plot_height 3 --bar_width 0.2

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/33_corr_heatmap/corr_heatmap_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy38_corr_taxa_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/efficacy38_corr_taxa.pdf -r /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy38_corr_taxa_row -c /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy38_corr_column --method spearman --plot_width 10 --plot_height 5 --text_size 12


/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/33_corr_heatmap/corr_heatmap_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy38_corr_pathway_in2 -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/efficacy38_corr_pathway.pdf -r /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy38_corr_pathway_row2 -c /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy38_corr_column --method spearman --plot_width 8 --plot_height 10 --text_size 12

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/33_corr_heatmap/corr_heatmap_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy19_corr_CAG_taxa_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/efficacy19_corr_CAG_taxa.pdf -r /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy19_corr_CAG_taxa_row -c /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy38_corr_column --method spearman --plot_width 10 --plot_height 6 --text_size 12

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/33_corr_heatmap/corr_heatmap_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy19_corr_CAG_pathway_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/efficacy19_corr_CAG_pathway.pdf -r /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy19_corr_CAG_pathway_row -c /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy38_corr_column --method spearman --plot_width 10 --plot_height 5 --text_size 12

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/33_corr_heatmap/corr_heatmap_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy19_corr_CBG_pathway_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/efficacy19_corr_CBG_pathway.pdf -r /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy19_corr_CBG_pathway_row -c /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy38_corr_column --method spearman --plot_width 10 --plot_height 5 --text_size 12

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/33_corr_heatmap/corr_heatmap_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/side_effect_corr_CAG_taxa_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/side_effect_corr_CAG_taxa.pdf -r /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/side_effect_corr_CAG_taxa_row -c /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy38_corr_column --method spearman --plot_width 10 --plot_height 5 --text_size 12

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/33_corr_heatmap/corr_heatmap_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/side_effect_corr_CAG_pathway_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/side_effect_corr_CAG_pathway.pdf -r /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/side_effect_corr_CAG_pathway_row -c /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy38_corr_column --method spearman --plot_width 10 --plot_height 3 --text_size 12

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/33_corr_heatmap/corr_heatmap_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/side_effect_corr_CBG_taxa_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/side_effect_corr_CBG_taxa.pdf -r /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/side_effect_corr_CBG_taxa_row -c /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy38_corr_column --method spearman --plot_width 10 --plot_height 3 --text_size 12

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/33_corr_heatmap/corr_heatmap_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/side_effect_corr_CBG_pathway_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/side_effect_corr_CBG_pathway.pdf -r /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/side_effect_corr_CBG_pathway_row -c /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy38_corr_column --method spearman --plot_width 10 --plot_height 3 --text_size 12

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/33_corr_heatmap/corr_heatmap_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/side_effect_corr_CBG_taxa_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/side_effect_corr_CBG_taxa.pdf -r /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/side_effect_corr_CBG_taxa_row -c /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy38_corr_column --method spearman --plot_width 10 --plot_height 3 --text_size 12

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/33_corr_heatmap/corr_heatmap_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/sam52_taxa_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/sam52_taxa.pdf -r /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/sam52_taxa_row -c /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/efficacy38_corr_column --method spearman --plot_width 10 --plot_height 5 --text_size 12

2021.05.18
分析临床指标与差异phyla和genera的关联分析
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/33_corr_heatmap/corr_heatmap_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/sam52_taxa_CBG_clinic_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/sam52_taxa_CBG_clinic.pdf -r /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/sam52_taxa_clinic_row -c /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/sam52_taxa_clinic_column --method spearman --plot_width 10 --plot_height 4 --text_size 12

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/33_corr_heatmap/corr_heatmap_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/sam52_taxa_CAG_clinic_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/sam52_taxa_CAG_clinic.pdf -r /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/sam52_taxa_clinic_row -c /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/sam52_taxa_clinic_column --method spearman --plot_width 10 --plot_height 4 --text_size 12

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/30_diff_boxplot/dif_boxplot_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/dif_genus_CBG_gender_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/gender_CBG_dif_genus.jpg -g /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/gender_compare -t Different_genus --plot_width 5 --plot_height 10 --y_min 0 --y_step 2 --y_max 12 --ncol 1

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/30_diff_boxplot/dif_boxplot_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/dif_genus_CAG_gender_in -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/graph/gender_CAG_dif_genus.jpg -g /work/workspace/zhurj/project/3_16S/anhuilung20201229/taxa_samcutoff3_20210226/input/graph_in/gender_compare -t Different_genus --plot_width 5 --plot_height 10 --y_min 4 --y_step 2 --y_max 16 --ncol 1

2021.03.12
CAG_CBG genus
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/30_diff_boxplot/dif_boxplot_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/CAG_CBG52_20210312/graph_in/cagcbg52_genus -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/CAG_CBG52_20210312/graph/cagcbg52_dif_genus.jpg -g /work/workspace/zhurj/project/3_16S/anhuilung20201229/CAG_CBG52_20210312/input/cag_cbg_compare -t Different_genus --plot_width 14 --plot_height 5 --y_min 0 --y_step 2 --y_max 14 --ncol 6

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/30_diff_boxplot/dif_boxplot_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/CAG_CBG52_20210312/graph_in/cagcbg52_phylum -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/CAG_CBG52_20210312/graph/cagcbg52_dif_phylum.jpg -g /work/workspace/zhurj/project/3_16S/anhuilung20201229/CAG_CBG52_20210312/input/cag_cbg_compare -t Different_phylum --plot_width 4 --plot_height 5 --y_min 4 --y_step 2 --y_max 16 --ncol 2
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/30_diff_boxplot/dif_boxplot_creation.r -i /work/workspace/zhurj/project/3_16S/anhuilung20201229/CAG_CBG52_20210312/graph_in/cagcbg52_species -o /work/workspace/zhurj/project/3_16S/anhuilung20201229/CAG_CBG52_20210312/graph/cagcbg52_dif_species.jpg -g /work/workspace/zhurj/project/3_16S/anhuilung20201229/CAG_CBG52_20210312/input/cag_cbg_compare -t Different_species --plot_width 16 --plot_height 5 --y_min 0 --y_step 2 --y_max 16 --ncol 10

陈博合作项目的数据分析：

  qiime taxa barplot \
  --i-table /work/workspace/liangzj/project/16S/2021/CPACPB/newCPACPB0223/2.qiime_data/table.qza \
  --i-taxonomy /work/workspace/liangzj/project/16S/2021/CPACPB/newCPACPB0223/2.qiime_data/taxonomy.qza \
  --m-metadata-file /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/input/metadata \
  --o-visualization /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/taxa-bar-plots.qzv

  qiime tools export \
  --input-path /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/taxa-bar-plots.qzv \
  --output-path /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/qzv


/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/qiime2_levels_to_filtered_frequency.py -i /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/data/phylum_tab -o /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/data/ --cutoff_sample 3 --abundance_cutoff 0.005 --level phylum --remove_num_col 1
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/7_16S/qiime2_levels_to_filtered_frequency.py -i /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/data/genus_tab -o /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/data/ --cutoff_sample 3 --abundance_cutoff 0.001 --level genus --remove_num_col 1

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/27_frequency2count/frequency2count.py -i /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/data/used_phylum_sam3_abun0.005 -o /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/data/used_phylum_count50000
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/27_frequency2count/frequency2count.py -i /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/data/used_genus_sam3_abun0.001 -o /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/data/used_genus_count50000

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/34_count2log2/count2log.py -i /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/data/used_phylum_count50000 -o /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/data/phylum_log2
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/34_count2log2/count2log.py -i /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/data/used_genus_count50000 -o /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/data/genus_log2

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/heatmap_taxo_topn_addrowanno.r -i /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/species_geneabun_heatmapin_top150 -o /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species -n 'dif genus abundance' -r /work/workspace/zhurj/project/1_metadata/mouse22/input/sgb_group --plot_width 12  --plot_height 18 --rowf /work/workspace/zhurj/project/1_metadata/mouse22/cdhit/taxo/taxo_lefse_dg/species/status

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/2_metapro/metapipe/heatmap_taxo_topn_addrowanno.r -i /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/data/dif_genus_log2 -o /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/graph/genus -n 'dif genus abundance' -r /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/input/group --plot_width 18  --plot_height 8 --rowf /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/input/rowmeta

/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/35_heatmap_anno/heatmap_taxo_addrowanno.r -i /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/data/dif_genus_log2 -o /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/graph/genus -n 'dif genus abundance' -r /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/input/group --plot_width 18  --plot_height 8 --rowf /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/input/rowmeta --treerow_n 1 

metacyc
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/34_count2log2/count2log.py -i /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/data/metacyc_tab -o /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/data/metacyc_log2

heatmap
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/35_heatmap_anno/heatmap_taxo_addrowanno.r -i /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/data/dif_metacyc_log2 -o /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/graph/metacyc -n 'dif metacyc' -r /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/input/group --plot_width 18  --plot_height 8 --rowf /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/input/rowmeta --treerow_n 1 

KO
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/34_count2log2/count2log.py -i /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/data/KO_tab -o /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/data/KO_log2

heatmap 
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/35_heatmap_anno/heatmap_taxo_addrowanno.r -i /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/data/dif_KO_log2 -o /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/graph/KO -n 'dif KO gene' -r /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/input/group --plot_width 15  --plot_height 25 --rowf /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/input/rowmeta --treerow_n 1 

phylum
/work/workspace/zhurj/bin/miniconda3/envs/R3.6/bin/Rscript /work/workspace/zhurj/script/1_tools/35_heatmap_anno/heatmap_taxo_addrowanno.r -i /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/data/dif_phylum_log2 -o /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/taxa/graph/phylum -n 'dif phylum' -r /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/input/group --plot_width 15  --plot_height 2 --rowf /work/workspace/zhurj/project/3_16S/chenjuanCPACPB20210318/input/rowmeta --treerow_n 1 

# 2021.04.28
# MNC-168 专利申请， 与其他11个菌株进化树制作
# 除MNC-168 外， 其他11个菌株的16S 序列
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/seqtk subseq /work/database/ezbio/16S/current/Ezbio_16S_seqs.fa /work/workspace/zhurj/project/16_patent/MNC168/entercoccusPhylo_20210428/namelist > /work/workspace/zhurj/project/16_patent/MNC168/entercoccusPhylo_20210428/strain11.fa

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/22_SeqIO_format_convert/format_converter.py -i /work/workspace/zhurj/project/16_patent/MNC168/entercoccusPhylo_20210428/strain11.fa --in_format fasta -o/work/workspace/zhurj/project/16_patent/MNC168/entercoccusPhylo_20210428/strain11.tab --out_format tab

# MNC-168 16S 序列 fasta 转 tab
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/22_SeqIO_format_convert/format_converter.py -i /work/workspace/zhurj/project/16_patent/MNC168/entercoccusPhylo_20210428/mnc168.fa --in_format fasta -o/work/workspace/zhurj/project/16_patent/MNC168/entercoccusPhylo_20210428/mnc168.tab --out_format tab

# 12个菌株16S 序列（包括MNC-168的序列） tab 转 fasta
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/22_SeqIO_format_convert/format_converter.py -i /work/workspace/zhurj/project/16_patent/MNC168/entercoccusPhylo_20210428/strain12.tab --in_format tab -o/work/workspace/zhurj/project/16_patent/MNC168/entercoccusPhylo_20210428/strain12.fa --out_format fasta

2021.05.07 
赵博画进化树
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/22_SeqIO_format_convert/format_converter.py -i /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/input_tab --in_format tab -o /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/input.fa --out_format fasta
muscle -in /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/input.fa -out /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/alignment.fas
fasttree -nt -gtr  /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/alignment.fas   >  /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/alignment_fasttree.nwk

/usr/bin/bash /work/workspace/zhurj/script/1_tools/38_species2level/genus_to_level.sh /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/test_genus /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/test_level
/usr/bin/bash /work/workspace/zhurj/script/1_tools/38_species2level/genus_to_level.sh /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/genus /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/genus_level

mycolors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")

R3.6
library(ggsci)
library(vegan)
mycolors <- unique(c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                      pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
                      pal_locuszoom("default")(7),pal_igv("default")(51),
                      pal_uchicago("default")(9),pal_startrek("uniform")(7),
                      pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
                      pal_simpsons("springfield")(16),pal_gsea("default")(12)))

mycolors[1:15]

mycolors[16:49]

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/seqtk subseq input.fa  /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/incorrect_strain_id > /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/needre_classify_strain.fna 
blastn -query /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/needre_classify_strain.fna  -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/needre_classify_strain.out -max_target_seqs 2 -num_threads 16 -outfmt 7
blastn -query /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/input.fa  -db /work/database/ezbio/16S/current/Ezbio_16S_seqs -out /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/all_1.out -max_target_seqs 1 -num_threads 16 -outfmt 7

blastn -query /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/input.fa  -db /work/database/ncbi/16S/current/16S_ribosomal_RNA -out /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/all_ribosomalRNA_1.out -max_target_seqs 1 -num_threads 16 -outfmt 7

# 将反向的16S rRNA序列调整方向
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/22_SeqIO_format_convert/format_converter.py -i /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/new/reverse_seq.tab --in_format tab -o /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/new/reverse.fa --out_format fasta
seqtk seq -r /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/new/reverse.fa > /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/new/reverse_to_forward.fa
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/22_SeqIO_format_convert/format_converter.py -i /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/new/forward_seq.tab --in_format tab -o /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/new/forward.fa --out_format fasta


cat /work/database/ezbio/16S/current/Ezbio_16S_seqs.tab | grep -E "AB696397|JF719608|AY442267|EU861881|X77437"

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/seqtk subseq /work/database/ncbi/16S/current/../2021.Apr28/16S_ribosomal_RNA.fas /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/new/ref_strain_id > /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/new/ref_strain.fa

/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/seqtk subseq /work/database/ncbi/16S/current/../2021.Apr28/16S_ribosomal_RNA.fas /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/new/ref_dedup_id > /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/new/ref_dedup_strain.fa
/work/workspace/zhurj/bin/miniconda3/envs/python3.6/bin/python /work/workspace/zhurj/script/1_tools/39_fasta_id_rename/id_rename_fasta.py -i /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/new/ref_dedup_strain.fa -o /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/new/ref_rename_id.fa --delimiter space
cat forward.fa reverse_to_forward.fa ref_rename_id.fa > merge_forward.fa

muscle -in /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/new/merge_forward.fa -out /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/new/merge_forward_alignment.fas
#fasttree -nt -gtr -notop -bionj -slow /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/new/merge_forward_alignment.fas   >  /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/new/merge_forward_alignment_fasttree.nwk
fasttree -nt -gtr -notop /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/new/merge_forward_alignment.fas   >  /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/new/merge_forward_alignment_fasttree.nwk

/usr/bin/bash /work/workspace/zhurj/script/1_tools/38_species2level/genus_to_level.sh /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/ref_genus /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/ref_genus_level

cat /work/workspace/zhurj/project/3_16S/phylo/zhaobo_20210507/ref_genus_level | grep -v "Eukaryota" > ref_bacteria_level
