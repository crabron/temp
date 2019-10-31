path = "in/" 
path_trein_set = "data/silva_nr_v132_train_set.fa"
path_trein_set_species = "data/silva_species_assignment_v132.fa"
maxEE = "2,5" 
mult = TRUE
mlt = NULL
require(dada2)
require(Biostrings)
require(DECIPHER)
require(phyloseq)
require(seqinr)
require(data.table)
require(metagMisc)
require(tibble)

fn <- sort(list.files(path, pattern=".fastq", full.names = TRUE))
sample.names <-sub('\\.fastq$', '', basename(fn))
on.exit()
filtFn <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq.gz"))
out <- filterAndTrim(fn, filtFn, truncLen=250, maxEE=5, multithread=TRUE)
err<- learnErrors(filtFn, multithread=mult)
derepFn <- derepFastq(filtFn, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFn) <- sample.names
dadaFn <- dada(derepFn, err=err, multithread=TRUE, pool="pseudo")
seqtab <- makeSequenceTable(dadaFn)
dim(seqtab)
table(nchar(getSequences(seqtab)))
getN <- function(x) sum(getUniques(x))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=mult, verbose=TRUE)

briefToSeq <- colnames(seqtab.nochim)
names(briefToSeq) <- paste0("Seq", seq(ncol(seqtab.nochim))) 
st.brief <- seqtab.nochim
colnames(st.brief) <- names(briefToSeq) 

write.table(briefToSeq, file = "454.brief.Cd.tsv", sep = "\t")

      #dna <- DNAStringSet(briefToSeq)
      #alignment <- AlignSeqs(DNAStringSet(dna), anchor=NA,verbose=FALSE, processors = mlt)
      #writeXStringSet(alignment, file="align.fasta")
      
taxa.dada2 <- assignTaxonomy(briefToSeq,path_trein_set , multithread=mult)
taxa.dada2.species <- assignSpecies(briefToSeq, path_trein_set_species) # maybe use not seqs but brief 
rownames(taxa.dada2.species) <- rownames(briefToSeq)
briefToSeq.df <- data.frame(briefToSeq)
rownames(taxa.dada2.species) <- rownames(briefToSeq.df)
rownames(taxa.dada2) <- rownames(taxa.dada2.species)
taxa <- cbind2(taxa.dada2, taxa.dada2.species[,2])
colnames(taxa)[7] <- "Species"

st.brief.t.df <- data.frame(t(st.brief))
write.table(st.brief.t.df, file = "otu_table.txt", sep= "\t", col.names = NA, quote=FALSE)

briefToSeq.ls <- as.list(briefToSeq.df[,c("briefToSeq")])
briefToSeq.names <- as.list(rownames(briefToSeq.df))
write.fasta( briefToSeq.ls, briefToSeq.names , "rep_seq.fasta", as.string = FALSE)
      
      system("conda activate qiime2-2019.1
      sed -i '1s/^/#OTU_ID/' otu_table.txt 
      biom convert -i otu_table.txt -o otu_table.biom --to-hdf5
      qiime tools import \
      --input-path rep_seq.fasta \
      --output-path rep_seq.qza \
      --type 'FeatureData[Sequence]'
      qiime fragment-insertion sepp --i-representative-sequences rep_seq.qza  --o-tree data/insertion-tree.qza  --o-placements data/insertion-placements.qza --p-threads 20
      qiime tools import \
      --input-path otu_table.biom \
      --type 'FeatureTable[Frequency]' \
      --input-format BIOMV210Format \
      --output-path feature-table.qza
      qiime fragment-insertion filter-features \
      --i-table feature-table.qza \
      --i-tree data/insertion-tree.qza \
      --o-filtered-table filtered_table.qza \
      --o-removed-table removed_table.qza
      unzip -p filtered_table.qza */data/* > filtered_table.biom
      biom convert -i  filtered_table.biom -o filtered_table.txt --to-tsv
      sed -i '1d' filtered_table.txt
      unzip -p data/insertion-tree.qza */data/* > tree.nwk
             ")
      
filt.otu <-t(as.data.frame(fread("filtered_table.txt")))
first <-  filt.otu[1,]
filt.otu <- filt.otu[-c(1),]
colnames(filt.otu) <-  first
class(filt.otu) <- "numeric"
filt.otu.matrix <- as.matrix(filt.otu)
#head.col <- scan("head.txt", character(), quote = "")
#rownames(filt.otu.matrix) <- head.col
tree <- read_tree(treefile="tree.nwk")

mapp <- read.csv("map.csv" , header=TRUE, sep="\t")
map <- data.frame(row.names="ID", mapp)

ps <- phyloseq(otu_table(filt.otu.matrix, taxa_are_rows=FALSE), 
               sample_data(map), 
               tax_table(taxa),
               phy_tree(tree))



pop_taxa <- function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

delete_mit_chl <- function(ps){
  badTaxa <- taxa_names(subset_taxa(ps, Order=="Chloroplast"))
  ps <- pop_taxa(ps, badTaxa)
  badTaxa <- taxa_names(subset_taxa(ps, Family=="Mitochondria"))
  ps <- pop_taxa(ps, badTaxa)
  return(ps)
}

ps.f <- delete_mit_chl(ps)

# data postprocessing. Delete and combine samples
ps.f.i <- prune_samples(sample_data(ps.f)$Inconceivability != c("pos"), ps.f)
ps.f.i <- prune_taxa(taxa_sums(ps.f.i) > 0, ps.f.i)
ps.f.i.Cdt <- prune_samples(sample_data(ps.f.i)$Genotype %in% c("Cdt"), ps.f.i)
ps.f.i.Cdt <- prune_taxa(taxa_sums(ps.f.i.Cdt) > 0, ps.f.i.Cdt)
ps.f.i.wt <- prune_samples(sample_data(ps.f.i)$Genotype %in% c("wt"), ps.f.i)
ps.f.i.wt <- prune_taxa(taxa_sums(ps.f.i.wt) > 0, ps.f.i.wt)

ps.f.i.wt.noCd <- prune_samples(sample_data(ps.f.i.wt)$Cd %in% c("neg"), ps.f.i.wt)
ps.f.i.wt.noCd <- prune_taxa(taxa_sums(ps.f.i.wt.noCd) > 0, ps.f.i.wt.noCd) 
ps.f.i.wt.Cd <- prune_samples(sample_data(ps.f.i.wt)$Cd %in% c("pos"), ps.f.i.wt)
ps.f.i.wt.Cd <- prune_taxa(taxa_sums(ps.f.i.wt.Cd) > 0, ps.f.i.wt.Cd) 
ps.f.i.Cdt.noCd <- prune_samples(sample_data(ps.f.i.Cdt)$Cd %in% c("neg"), ps.f.i.Cdt)
ps.f.i.Cdt.noCd <- prune_taxa(taxa_sums(ps.f.i.Cdt.noCd) > 0, ps.f.i.Cdt.noCd) 
ps.f.i.Cdt.Cd <- prune_samples(sample_data(ps.f.i.Cdt)$Cd %in% c("pos"), ps.f.i.Cdt)
ps.f.i.Cdt.Cd <- prune_taxa(taxa_sums(ps.f.i.Cdt.Cd) > 0, ps.f.i.Cdt.Cd)


amp_boxplot(amp.f,
            group_by = "Drought",
            tax_show = 10,
            tax_aggregate = "Phylum",
            tax_class = "Proteobacteria"
)

#plot nice beta

beta_custom_norm_NMDS <- function(ps, seed = 6788, normtype="vst", color="Cd", shape="Drought"){
  require(phyloseq)
  require(ggplot2)
  require(ggpubr)
  require(DESeq2)
  
  # beta_NMDS <- function(){
  #normalisation. unifrac - rarefaction; wunifrac,bray - varstab

  diagdds = phyloseq_to_deseq2(ps, ~ Repeats)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  if (normtype =="vst")
    pst <- varianceStabilizingTransformation(diagdds)
  if (normtype =="log") 
    pst <- rlogTransformation(diagdds)
  
  pst.dimmed <- t(as.matrix(assay(pst))) 
  pst.dimmed[pst.dimmed < 0.0] <- 0.0
  ps.varstab <- ps
  otu_table(ps.varstab) <- otu_table(pst.dimmed, taxa_are_rows = FALSE) 
  
  ps.rand <- rarefy_even_depth(ps, rngseed = seed)
  
  #beta and ordination
  
  ordination.b <- ordinate(ps.varstab, "NMDS", "bray")
  ordination.u <- ordinate(ps.rand, "NMDS", "unifrac")
  ordination.w <- ordinate(ps.varstab, "NMDS", "wunifrac")
  
  #plotting
  p = plot_ordination(ps, ordination.b, type="sample", color, shape, title="NMDS - Bray", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 10)) + geom_point(size = 3) 
  p1.1 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.u, type="sample", color, shape, title="NMDS - unifrac", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 10)) + geom_point(size = 3) 
  p1.2 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.w, type="sample", color, shape, title="NMDS - wunifrac", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 10)) + geom_point(size = 3) 
  p1.3 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  #merge by ggpubr
  
  p.all <- ggarrange(p1.1, p1.2, p1.3, ncol = 3 , nrow = 1, common.legend = TRUE, legend = "right", font.label = list(size = 12, face = "bold", color ="black"))
  
  return(p.all)
}

p1 <- beta_custom_norm_NMDS(ps.f.i.wt, seed =5544, color = "Drought", shape = "Cd", normtype = "vst")
p2 <- beta_custom_norm_NMDS(ps.f.i.Cdt, seed =5544, color = "Drought", shape = "Cd", normtype = "vst")
ggarrange(p1, p2, ncol = 1 , nrow = 2, label.y = 0.75, label.x = 0.935,  labels = c("SGE.wt","SGE.Cdt"), font.label = list(size = 14, face = "bold", color ="black"))

#alpha pd - richness - shannon

library(picante)
ps.merged.wt <- merge_samples(ps.f.i.wt, "Repeats")
pd(ps.merged.wt@otu_table@.Data, ps.merged.wt@phy_tree)

ps.merged.Cdt <- merge_samples(ps.f.i.Cdt, "Repeats")
pd(ps.merged.Cdt@otu_table@.Data, ps.merged.Cdt@phy_tree)

alpha.custom <- function(ps.f.i, arrga = "PD"){
  pd <- pd(ps.f.i@otu_table@.Data, ps.f.i@phy_tree)
  
  pd1 <- estimate_richness(ps.f.i, measures=c("Observed", "InvSimpson", "Shannon"))
  pd <- cbind(pd,ps.f.i@sam_data, pd1 )
  
  bmi <- levels(pd$Вариант)
  bmi.pairs <- combn(seq_along(bmi), 2, simplify = FALSE, FUN = function(i)bmi[i])
  p1 <- ggviolin(pd, x = "Вариант", y = arrga,
                 add = "boxplot", fill = "Вариант") + stat_compare_means(comparisons = bmi.pairs, method = "wilcox.test") +
  scale_x_discrete(limits = c("конт полив", "конт засуха","Сd полив", "Cd засуха")) + theme(axis.title.x = element_blank(), legend.title = element_blank())
  return(p1)
}

p.alpha.Cdt.oo <- alpha.custom(ps.f.i.Cdt, arrga = "Observed")
p.alpha.Cdt.sh <- alpha.custom(ps.f.i.Cdt, arrga = "Shannon")
p.alpha.Cdt.is <- alpha.custom(ps.f.i.Cdt, arrga = "InvSimpson")
p.alpha.Cdt.pd <- alpha.custom(ps.f.i.Cdt, arrga = "PD")

p.alpha.wt.oo <- alpha.custom(ps.f.i.wt, arrga = "Observed")
p.alpha.wt.sh <- alpha.custom(ps.f.i.wt, arrga = "Shannon")
p.alpha.wt.is <- alpha.custom(ps.f.i.wt, arrga = "InvSimpson")
p.alpha.wt.pd <- alpha.custom(ps.f.i.wt, arrga = "PD")

p1 <- ggarrange(p.alpha.Cdt.oo, p.alpha.Cdt.sh,p.alpha.Cdt.is, p.alpha.Cdt.pd , ncol = 4 ,label.x = 0.105, nrow = 1, common.legend = TRUE)
p2 <- ggarrange(p.alpha.wt.oo, p.alpha.wt.sh,p.alpha.wt.is, p.alpha.wt.pd , ncol = 4 ,label.x = 0.105, nrow = 1, common.legend = TRUE)
ggarrange(p1, p2, ncol = 1 , nrow = 2, label.y = 0.9855, label.x = 0.32,  labels = c("SGE.wt","SGE.Cdt"), font.label = list(size = 14, face = "bold", color ="black"))  


# ampvis heatmap on phylum level
amp_heatmap_grouped <- function(ps){
  require(ampvis2)
  require(phyloseq)
  require(ggpubr)
  
  # split the original dataset into two then convert the goddamnned phyloseq object into ampvis2 asshanded class
  ps.Cd <- prune_samples(sample_data(ps)$Cd %in% c("pos"), ps)
  ps.Cd <- prune_taxa(taxa_sums(ps.Cd) > 0, ps.Cd)  
  amp.Cd <- phyloseq_to_amp(ps.Cd)
  
  ps.wCd <- prune_samples(sample_data(ps)$Cd %in% c("neg"), ps)
  ps.wCd <- prune_taxa(taxa_sums(ps.wCd) > 0, ps.wCd)  
  amp.wCd <- phyloseq_to_amp(ps.wCd)
  
  # plot some ampvised heatmaps
  p.heat.Cd <- amp_boxplot(amp.wCd,
              group_by = "Drought",
              tax_show = 8,
              tax_aggregate = "Phylum",
              tax_class = "Proteobacteria") + labs(title = "Cd-", color = "drought") + theme_bw() + theme(text = element_text(size=14))

  p.heat.wCd <- amp_boxplot(amp.Cd,
              group_by = "Drought",
              tax_show = 8,
              tax_aggregate = "Phylum",
              tax_class = "Proteobacteria") + labs(title = "Cd+", color = "drought") + theme_bw()  + theme(text = element_text(size=14))
  
  p <- ggarrange(p.heat.Cd, p.heat.wCd, ncol = 2 ,label.x = 0.105, nrow = 1, common.legend = TRUE)
  return(p)
  }

p1 <- amp_heatmap_grouped(ps.f.i.wt)
p2 <- amp_heatmap_grouped(ps.f.i.Cdt)
ggarrange(p1, p2, ncol = 1 , nrow = 2, label.y = 0.985, label.x = 0.35,  labels = c("SGE.wt","SGE.Cdt"), font.label = list(size = 14, face = "bold", color ="black"))                   

permanova.454 <- function(ps, dist = "bray"){
  require(phyloseq)
  require(vegan)
  dist <- phyloseq::distance(ps, dist)
  metadata <- as(sample_data(ps), "data.frame")
  ad <- adonis2(dist ~ Drought, data = metadata, permutations = 10000)
  return(ad)
}



library(phyloseq)             
library(vegan)

# функция преобразует otu-table из phyloseq объекта в otu-table необходимый для пакета vegan
  veganifyOTU <- function(physeq){
       require(phyloseq)
       if(taxa_are_rows(physeq)){physeq <- t(physeq)}
       return(as(otu_table(physeq), "matrix"))
    }

#вместо ps необходимый объект класса phyloseq
otus.ps.vegan <- veganifyOTU(ps.f.i.meta)
metadata <- as(sample_data(ps.f.i.meta), "data.frame")
#собственно cca(у rda аналогичный синтаксис), указано какие факторы из метаданных учитывать
vare.cca <- cca(otus.ps.vegan ~ B+Ca+Fe+K+Mg+M+Na+P+S+Zn+X15N+N.tot+X15N.nakop+Ntot.nakop, data=metadata)

kate.ggcca.sites <- function(ps){
  require(ggvegan)
  require(vegan)
  require(ggplot2)
  require(dplyr)
  require(ggrepel)
  require(phyloseq)
  require(DESeq2)
  
  diagdds = phyloseq_to_deseq2(ps, ~ Repeats)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  pst <- varianceStabilizingTransformation(diagdds)
  pst.dimmed <- t(as.matrix(assay(pst))) 
  pst.dimmed[pst.dimmed < 0.0] <- 0.0
  ps.varstab <- ps
  otu_table(ps.varstab) <- otu_table(pst.dimmed, taxa_are_rows = FALSE)

  # функция преобразует otu-table из phyloseq объекта в otu-table необходимый для пакета vegan
  veganifyOTU <- function(physeq){
    require(phyloseq)
    if(taxa_are_rows(physeq)){physeq <- t(physeq)}
    return(as(otu_table(physeq), "matrix"))
  }
  
  #вместо ps необходимый объект класса phyloseq

  otus.ps.vegan <- veganifyOTU(ps.varstab)
  metadata <- as(sample_data(ps.varstab), "data.frame")
  #собственно cca(у rda аналогичный синтаксис), указано какие факторы из метаданных учитывать
  vare.cca <- cca(otus.ps.vegan ~ B+Ca+Fe+K+Mg+M+Na+P+S+Zn+X15N+N.tot+X15N.nakop+Ntot.nakop, data=metadata)
  
  otu_table(ps.varstab) <- otu_table(pst.dimmed, taxa_are_rows = FALSE)
  fdat <- fortify(vare.cca)
  p.sites <- ggplot(fdat %>% filter(Score %in% c("sites","biplot"))) + geom_point(data = fdat %>% dplyr::filter(Score == "sites"), mapping = aes(x=CCA1, y=CCA2, colour = factor(Score))) + geom_segment(data = fdat %>% dplyr::filter(Score == "biplot"), aes(x = 0, xend = CCA1, y = 0, yend = CCA2), alpha=0.8,
                                                                                                                                                                                                         color = "red",arrow = arrow(angle = 3))  + 
    geom_text_repel(aes(x=CCA1, y=CCA2, label= Label),size=4) 
  p <- p.sites + theme(legend.position = "none", panel.background = element_rect(fill = "white", colour = "grey50"))
  return(p)
}

kate.ggcca.species <- function(ps){
  require(ggvegan)
  require(vegan)
  require(ggplot2)
  require(dplyr)
  require(ggrepel)
  require(phyloseq)
  require(DESeq2)
  require(microbiomeSeq)
  
  diagdds = phyloseq_to_deseq2(ps, ~ Repeats)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  pst <- varianceStabilizingTransformation(diagdds)
  pst.dimmed <- t(as.matrix(assay(pst))) 
  pst.dimmed[pst.dimmed < 0.0] <- 0.0
  ps.varstab <- ps
  otu_table(ps.varstab) <- otu_table(pst.dimmed, taxa_are_rows = FALSE)
  
  # функция преобразует otu-table из phyloseq объекта в otu-table необходимый для пакета vegan
  veganifyOTU <- function(physeq){
    require(phyloseq)
    if(taxa_are_rows(physeq)){physeq <- t(physeq)}
    return(as(otu_table(physeq), "matrix"))
  }
  
  #вместо ps необходимый объект класса phyloseq
  ps.varstab <- taxa_level(ps.varstab, "Phylum")
  otus.ps.vegan <- veganifyOTU(ps.varstab)
  metadata <- as(sample_data(ps.varstab), "data.frame")
  #собственно cca(у rda аналогичный синтаксис), указано какие факторы из метаданных учитывать
  vare.cca <- cca(otus.ps.vegan ~ B+Ca+Fe+K+Mg+M+Na+P+S+Zn+X15N+N.tot+X15N.nakop+Ntot.nakop, data=metadata)
  
  otu_table(ps.varstab) <- otu_table(pst.dimmed, taxa_are_rows = FALSE)
  fdat <- fortify(vare.cca)
  p.sites <- ggplot(fdat %>% filter(Score %in% c("species","biplot"))) + geom_point(data = fdat %>% dplyr::filter(Score == "species"), mapping = aes(x=CCA1, y=CCA2, colour = factor(Score))) + geom_segment(data = fdat %>% dplyr::filter(Score == "biplot"), aes(x = 0, xend = CCA1, y = 0, yend = CCA2), alpha=0.8,
                                                                                                                                                                                                             color = "red",arrow = arrow(angle = 3))  + 
    geom_text_repel(aes(x=CCA1, y=CCA2, label= Label),size=4) 
  p <- p.sites + theme(legend.position = "none", panel.background = element_rect(fill = "white", colour = "grey50"))
  return(p)
}

# mantel. vegan wrapper. oop elements.

mantel.all <- setClass("mantel.all", slots = list())

mantel.all <- function(ps, dist = "bray"){
  require(purrr)
  require(phyloseq)
  require(vegan)
  require(tibble)
  require(dplyr)
  
  mantel.local <- function(ps, element){
    require(vegan)
    require(phyloseq)
    map <- as.data.frame(ps@sam_data)
    # map.num <- mapply(map, FUN=as.numeric)
    dist <- phyloseq::distance(ps, dist)
    env.dist <- vegdist(scale(map[[element]]) , "euclid")
    mant <- mantel(dist, env.dist, method="pearson", permutations = 9999)
    return(list(element,mant$statistic, mant$signif))
  }
  
  l <- colnames(ps@sam_data)
  f <- function(x) {mantel.local(ps, x)}
  l.mantel <- lapply(l, purrr::possibly(f, "NaN"))
  l.mantel <- lapply(l.mantel, function(x) x[x != "NaN"])
  d.mantel <- do.call(rbind.data.frame, l.mantel)
  colnames(d.mantel) <- c("ooo","statistics", "significance")
  d.mantel <- remove_rownames(d.mantel) %>% column_to_rownames(var = "ooo")
  return(d.mantel)
}

