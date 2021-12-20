#Additional commands to stratify samples by AD status and calculate alpha and beta diversity in QIIME2
#Date last modified: Dec 19 2021

# To filter for infant observations 
qiime feature-table filter-samples \
  --i-table ./dada2_table.qza \
  --m-metadata-file /data/t4_project2/infant/infant_metadata.txt \
  --p-where "[life_stage]='Infant'" \
  --o-filtered-table infant-filtered-table.qza

#To filter for AD samples and visualize feature table
qiime feature-table filter-samples \
  --i-table ./infant-filtered-table.qza \
  --m-metadata-file /data/t4_project2/infant/infant_metadata.txt \
  --p-where "[eczema_inf]='yes'" \
  --o-filtered-table /data/t4_project2/AD_outputs/AD-inf-filtered-table.qza

qiime feature-table summarize \
  --i-table ./AD-inf-filtered-table.qza \
  --m-sample-metadata-file /data/t4_project2/infant/infant_metadata.txt \
  --o-visualization ./AD-inf-filtered-table.qzv

#to filter for non AD samples and visualize feature table
qiime feature-table filter-samples \
  --i-table ./infant-filtered-table.qza \
  --m-metadata-file /data/t4_project2/infant/infant_metadata.txt \
  --p-where "[eczema_inf]='no'" \
  --o-filtered-table /data/t4_project2/Non_AD_outputs/Non-AD-inf-filtered-table.qza

qiime feature-table summarize \
  --i-table ./Non-AD-inf-filtered-table.qza \
  --m-sample-metadata-file /data/t4_project2/infant/infant_metadata.txt \
  --o-visualization ./Non-AD-inf-filtered-table.qzv

# Generate alpha rarefaction curve for AD infants 
qiime diversity alpha-rarefaction \  
--i-table ./AD-inf-filtered-table.qza \
--m-metadata-file /data/t4_project2/infant/infant_metadata.txt \  
--o-visualization /data/t4_project2/AD_outputs/alpha_rarefaction_curves.qzv \   
--p-min-depth 1 \  
--p-max-depth 49972

# Generate alpha rarefaction curve for Non- AD infants 
qiime diversity alpha-rarefaction \  
--i-table /data/t4_project2/Non_AD_outputs/Non-AD-inf-filtered-table.qza \
--m-metadata-file /data/t4_project2/infant/infant_metadata.txt \  
--o-visualization /data/t4_project2/Non-AD_outputs/alpha_rarefaction_curves.qzv \   
--p-min-depth 1 \  
--p-max-depth 49972

# Generate phylogenetic tree
qiime fragment-insertion sepp \
  --i-representative-sequences ./dada2_rep_set.qza \
  --i-reference-database /data/mouse_tutorial/sepp-refs-silva-128.qza\
  --o-tree ./tree.qza \
  --o-placements ./tree_placements.qza \
  --p-threads 8

# Diversity analysis (core metrics) for AD samples
qiime diversity core-metrics-phylogenetic \
  --i-table ./AD-inf-filtered-table.qza \
  --i-phylogeny /data/t4_project2/script_outputs/tree.qza \
  --m-metadata-file /data/t4_project2/infant/infant_metadata.txt \
  --p-sampling-depth 26750 \
  --output-dir ./core-metrics-results

# Diversity analysis (core metrics) for non AD samples
qiime diversity core-metrics-phylogenetic \
  --i-table ./Non-AD-inf-filtered-table.qza \
  --i-phylogeny /data/t4_project2/script_outputs/tree.qza \
  --m-metadata-file /data/t4_project2/infant/infant_metadata.txt \
  --p-sampling-depth 26750 \
  --output-dir ./core-metrics-results

# Generate faith_pd statistics for AD samples
qiime diversity alpha-group-significance \
  --i-alpha-diversity ./core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file /data/t4_project2/infant/infant_metadata.txt \
  --o-visualization ./core-metrics-results/faiths_pd_statistics.qzv

# Generate faith_pd statistics for non AD samples
qiime diversity alpha-group-significance \
  --i-alpha-diversity ./core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file /data/t4_project2/infant/infant_metadata.txt \
  --o-visualization ./core-metrics-results/faiths_pd_statistics.qzv

# Generate Pielou’s evenness statistics for AD samples
qiime diversity alpha-group-significance \
 --i-alpha-diversity ./core-metrics-results/evenness_vector.qza \
 --m-metadata-file /data/t4_project2/infant/infant_metadata.txt \
 --o-visualization ./core-metrics-results/evenness_statistics.qzv

# Generate Pielou’s evenness statistics for non AD samples
qiime diversity alpha-group-significance \
 --i-alpha-diversity ./core-metrics-results/evenness_vector.qza \
 --m-metadata-file /data/t4_project2/infant/infant_metadata.txt \
 --o-visualization ./core-metrics-results/evenness_statistics.qzv

# Beta Diversity for AD samples
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

# Beta Diversity for non AD samples
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

#Taxonomic classification for AD
qiime feature-classifier classify-sklearn \
  --i-reads /data/t4_project2/script_outputs/dada2_rep_set.qza \
  --i-classifier /mnt/datasets/classifiers/silva-138-99-515-806-nb-classifier.qza  \
  --o-classification ./taxonomy.qza

#Taxonomic classification for non AD
qiime feature-classifier classify-sklearn \
  --i-reads /data/t4_project2/script_outputs/dada2_rep_set.qza \
  --i-classifier /mnt/datasets/classifiers/silva-138-99-515-806-nb-classifier.qza  \
  --o-classification ./taxonomy.qza

#Tabulate taxonomy associated with the sequences with metadata for AD
qiime metadata tabulate \
  --m-input-file ./taxonomy.qza \
  --o-visualization ./taxonomy.qzv

#Tabulate taxonomy associated with the sequences with metadata for non AD
qiime metadata tabulate \
  --m-input-file ./taxonomy.qza \
  --o-visualization ./taxonomy.qzv

#Tabulate the representative sequence for AD
qiime feature-table tabulate-seqs \
  --i-data /data/t4_project2/script_outputs/dada2_rep_set.qza \
  --o-visualization ./dada2_rep_set.qzv

#Tabulate the representative sequence for non AD
qiime feature-table tabulate-seqs \
  --i-data /data/t4_project2/script_outputs/dada2_rep_set.qza \
  --o-visualization ./dada2_rep_set.qzv

#Filter out the sample with sequence fewer than the rarefaction depth set in the previous step (26750) for AD
qiime feature-table filter-samples \
  --i-table /data/t4_project2/script_outputs/infant-filtered-table.qza \
  --p-min-frequency 26750 \
  --o-filtered-table ./table_2k.qza

#Filter out the sample with sequence fewer than the rarefaction depth set in the previous step (26750) for non AD
qiime feature-table filter-samples \
  --i-table /data/t4_project2/script_outputs/infant-filtered-table.qza \
  --p-min-frequency 26750 \
  --o-filtered-table ./table_2k.qza

#Use rarefied table to build an interactive bar plot of taxonomy for AD samples
qiime taxa barplot \
  --i-table ./table_2k.qza \
  --i-taxonomy ./taxonomy.qza \
  --m-metadata-file /data/t4_project2/infant/infant_metadata.txt \
  --o-visualization ./taxa_barplot.qzv

#Use rarefied table to build an interactive bar plot of taxonomy for non AD samples
qiime taxa barplot \
  --i-table ./table_2k.qza \
  --i-taxonomy ./taxonomy.qza \
  --m-metadata-file /data/t4_project2/infant/infant_metadata.txt \
  --o-visualization ./taxa_barplot.qzv
