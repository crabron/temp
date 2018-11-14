#!/bin/bash

set -e

# while getops "ibm" opt ;
# do
#     case $opt in    
#             b) BASE=$OPTARG;
        
#         esac 
#     done

qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path $PWD\gz_fold
    --output-path demux.qza
    --source-format PairedEndFastqManifestPhred33 

qiime dada2 denoise-single \
    --i-demultiplexed-seqs demux.qza \
    --p-trim-left 20 \
    --p-trunc-len 180 \
    --o-representative-sequences rep-seqs.qza \
    --o-table table.qza \
    --o-denoising-stats stats.qza

qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences rep-seqs.qza \
    --o-alignment aligned-rep-seqs.qza \
    --o-masked-alignment masked-aligned-rep-seqs.qza \
    --o-tree unrooted-tree.qza \
    --o-rooted-tree rooted-tree.qza

qiime diversity alpha-phylogenetic \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-metric faith_pd \
  --o-alpha-diversity faith_pd_vector.qza

qiime diversity alpha \
  --i-table table.qza \
  --p-metric observed_otus \
  --o-alpha-diversity observed_otus_vector.qza

