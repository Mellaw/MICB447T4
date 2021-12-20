#Data quality control and diversity analysis in QIIME2
#Date last modified: Dec 19, 2021

#Load data 
#link relevant files to the team folder using symbolic links
ln -s /mnt/datasets/project_2/infant

#using manifest important data into qiime
cd ~/data/t4_project2/infant

qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path ./infant_manifest.txt \
  --output-path /data/t4_project2/script_outputs/demux_seqs.qza
 
#Conversion of qza to qzv
cd ~/data/t4_project2/script_outputs

qiime demux summarize \
  --i-data ./demux_seqs.qza \
  --o-visualization /data/t4_project2/script_outputs/demux_seqs.qzv

#Quality control 
#At truncation length 129, Phred for Q3 & Q4 = 37, 38
qiime dada2 denoise-single \
  --i-demultiplexed-seqs ./demux_seqs.qza \
  --p-trunc-len 129 \
  --o-table ./dada2_table.qza \
  --o-representative-sequences ./dada2_rep_set.qza \
  --o-denoising-stats ./dada2_stats.qza


#FILTER FOR DESIRED METADATA CATEGORIES (check w Karen from her code)
qiime feature-table filter-samples \
  --i-table ./dada2_table.qza \
  --m-metadata-file /data/t4_project2/infant/infant_metadata.txt \
  --p-where "[life_stage]='Infant'" \
  --o-filtered-table infant-filtered-table.qza

#Review quality control step 
qiime metadata tabulate \
  --m-input-file ./dada2_stats.qza  \
  --o-visualization ./dada2_stats.qzv
  
#Generate a features table 
qiime feature-table summarize \
  --i-table ./dada2_table.qza \
  --m-sample-metadata-file /data/t4_project2/infant/infant_metadata.txt \
  --o-visualization ./dada2_table.qzv

#Download SILVA database and upload to server
# URL; https://docs.qiime2.org/2021.8/data-resources/#sepp-reference-databases
wget \
  -O "sepp-refs-silva-128.qza" \
  "https://docs.qiime2.org/2021.8/data-resources/#sepp-reference-databases"

#Create session to run phylogenetic tree building in background 
screen -S phylotree

#Generate phylogentic trees
qiime fragment-insertion sepp \
  --i-representative-sequences ./dada2_rep_set.qza \
  --i-reference-database sepp-refs-silva-128.qza\
  --o-tree ./tree.qza \
  --o-placements ./tree_placements.qza \
  --p-threads 8  # update to a higher number if you can



#Diversity Analysis 
# Generate alpha rarefaction curve
qiime diversity alpha-rarefaction \  
--i-table /data/t4_project2/script_outputs/infant-filtered-table.qza \
--m-metadata-file /data/t4_project2/infant/infant_metadata.txt \  
--o-visualization /data/t4_project2/script_outputs/alpha_rarefaction_curves.qzv \   
--p-min-depth 1 \  
--p-max-depth 59159

# Generate phylogenetic tree
qiime fragment-insertion sepp \
  --i-representative-sequences ./dada2_rep_set.qza \
  --i-reference-database /data/mouse_tutorial/sepp-refs-silva-128.qza\
  --o-tree ./tree.qza \
  --o-placements ./tree_placements.qza \
  --p-threads 8

# Diversity analysis
qiime diversity core-metrics-phylogenetic \
  --i-table /data/t4_project2/script_outputs/infant-filtered-table.qza \
  --i-phylogeny ./tree.qza \
  --m-metadata-file /data/t4_project2/infant/infant_metadata.txt \
  --p-sampling-depth 26745 \
  --output-dir ./core-metrics-results

# Generate faith_pd statistics
qiime diversity alpha-group-significance \
  --i-alpha-diversity ./core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file /data/t4_project2/infant/infant_metadata.txt \
  --o-visualization ./core-metrics-results/faiths_pd_statistics.qzv

# Generate Pielou’s evenness statistics
qiime diversity alpha-group-significance \
 --i-alpha-diversity ./core-metrics-results/evenness_vector.qza \
 --m-metadata-file /data/t4_project2/infant/infant_metadata.txt \
 --o-visualization ./core-metrics-results/evenness_statistics.qzv

# Beta Diversity
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/t4_project2/infant/infant_metadata.txt \
  --m-metadata-column feed \
  --o-visualization core-metrics-results/unweighted-unifrac-feed-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file /data/t4_project2/infant/infant_metadata.txt \
  --m-metadata-column feed \
  --o-visualization core-metrics-results/weighted-unifrac-feed-significance.qzv

#Taxonomic classification and taxonomy barchart
qiime feature-classifier classify-sklearn \
  --i-reads ./dada2_rep_set.qza \
  --i-classifier /mnt/datasets/classifiers/silva-138-99-515-806-nb-classifier.qza  \
  --o-classification ./taxonomy.qza

#Tabulate taxonomy associated with the sequences with metadata
qiime metadata tabulate \
  --m-input-file ./taxonomy.qza \
  --o-visualization ./taxonomy.qzv

#Tabulate the representative sequence
qiime feature-table tabulate-seqs \
  --i-data ./dada2_rep_set.qza \
  --o-visualization ./dada2_rep_set.qzv

#Filter out the sample with sequence fewer than the rarefaction depth set in the previous step (26745)
qiime feature-table filter-samples \
  --i-table /data/t4_project2/script_outputs/infant-filtered-table.qza \
  --p-min-frequency 26745 \
  --o-filtered-table ./table_2k.qza

#Use rarefied table to build an interactive bar plot of taxonomy for each sample
qiime taxa barplot \
  --i-table ./table_2k.qza \
  --i-taxonomy ./taxonomy.qza \
  --m-metadata-file /data/t4_project2/infant/infant_metadata.txt \
  --o-visualization ./taxa_barplot.qzv


#Export infant filtered feature tables as .biom file
qiime tools export \
  --input-path ./infant-filtered-table.qza \
  --output-path /root/data/t4_project2/script_outputs/qiime_out

#Export taxonmic classification for each ASV as a tsv file
qiime tools export \
  --input-path ./taxonomy.qza \
  --output-path /root/data/t4_project2/script_outputs/qiime_out

#Export the rooted phylogenetic tree as a .nwk file
qiime tools export \
  --input-path ./tree.qza \
  --output-path /root/data/t4_project2/script_outputs/qiime_out
