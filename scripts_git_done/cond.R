func.more <- function(x){ 
  ress <- x[1] > x[2] | x[1] > x[3] | x[3] > x[4]
  return(ress)
}
otus.lapp <- as.data.frame(apply(otus.t@.Data, 1, func.more))

colnames(otus.lapp) <- c("Cond")
View(otus.lapp)
View(otus.t@.Data)
head(otus.t)

library(phyloseq)             
library(vegan)

# функция преобразует otu-table из phyloseq объекта в otu-table необходимый для пакета vegan
veganifyOTU <- function(physeq){
  require(phyloseq)
  if(taxa_are_rows(physeq)){physeq <- t(physeq)}
  return(as(otu_table(physeq), "matrix"))
}

#вместо ps необходимый объект класса phyloseq
otus.ps.vegan <- veganifyOTU(ps.s.merged)
metadata <- as(sample_data(ps.s.merged), "data.frame")
#собственно cca(у rda аналогичный синтаксис), указано какие факторы из метаданных учитывать
vare.cca <- cca(otus.ps.vegan ~ Mean_depth + N.total + Р2О5 + К2О , data=metadata)

colnames(merged.metadata)

svir.ggcca.species <- function(vare.cca){
  require(vegan)
  require(ggplot2)
  require(dplyr)
  require(ggrepel)
  require(ggvegan)
  fdat <- fortify(vare.cca)
  p.sites <- ggplot(fdat %>% filter(Score %in% c("species","biplot"))) + geom_point(data = fdat %>% dplyr::filter(Score == "species"), mapping = aes(x=CCA1, y=CCA2, colour = factor(Score))) + geom_segment(data = fdat %>% dplyr::filter(Score == "biplot"), aes(x = 0, xend = CCA1, y = 0, yend = CCA2), alpha=0.8,
                                                                                                                                                                                                             color = "red",arrow = arrow(angle = 3))  + 
    geom_text_repel(aes(x=CCA1, y=CCA2, label= Label),size=4) 
  p <- p.sites + theme(legend.position = "none", panel.background = element_rect(fill = "white", colour = "grey50"))
  return(p)
}


ps.s.5 <- prune_samples(sample_data(ps.s)$Site %in% c("5"), ps.s)
ps.s.bf <- prune_samples(sample_data(ps.s)$Horizont %in% c("BF"), ps.s)
ps.s.merged.5 <- prune_samples(sample_data(ps.s.merged)$Run %in% c("5"), ps.s.merged)
sample_data(ps.s.merged)

Des.ls <- function(ps){
  require(DESeq2)
  require(phyloseq)
  diagdds = phyloseq_to_deseq2(ps, ~  Site)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  diagdds = DESeq(diagdds)
  samp <-sample_data(ps)
  dds.counts <- diagdds@assays@.xData$data$counts
  dds.counts.df <- as.data.frame(dds.counts)
  aggdata <- t(aggregate.data.frame(t(dds.counts.df), by=list(samp$Site), median))
  colnames(aggdata) <- aggdata[1,]
  aggdata <- aggdata[-1,]
  res = results(diagdds)
  res.df <- as.data.frame(res)
  nice <- cbind(res.df,as.data.frame(tax_table(ps)[rownames(res.df),]), as.data.frame(aggdata)[rownames(res.df),])               
  return(nice)
}  

kate.ggcca.sites <- function(vare.cca){
  require(vegan)
  require(ggplot2)
  require(dplyr)
  require(ggrepel)
  fdat <- fortify(vare.cca)
  p.sites <- ggplot(fdat %>% filter(Score %in% c("sites","biplot"))) + geom_point(data = fdat %>% dplyr::filter(Score == "sites"), mapping = aes(x=CCA1, y=CCA2, colour = factor(Score))) + geom_segment(data = fdat %>% dplyr::filter(Score == "biplot"), aes(x = 0, xend = CCA1, y = 0, yend = CCA2), alpha=0.4,
                                                                                                                                                                                                         color = "red",arrow = arrow(angle = 3))  + 
    geom_text_repel(aes(x=CCA1, y=CCA2, label= Label),size=4) 
  p <- p.sites + theme(legend.position = "none", panel.background = element_rect(colour = "grey50"))
  return(p)
}

ps.s.5.pr <- prune_taxa(taxa_sums(ps.s.5) > 10, ps.s.5)

amp_heatmap(amp, group_by = "Site", facet_by = "Horizont", tax_show = 40,tax_aggregate = "Genus", tax_add = "Family", color_vector = c("white","azure4"))



beta_NMDS <- function(ps){
  require(phyloseq)
  require(ggplot2)
  require(DESeq2)

  diagdds <- phyloseq_to_deseq2(ps, ~ Horizont)                  
  diagdds = estimateSizeFactors(diagdds, type="poscounts")
  diagdds = estimateDispersions(diagdds, fitType = "local") 
  pst <- varianceStabilizingTransformation(diagdds)
  pst.dimmed <- t(as.matrix(assay(pst))) 
  pst.dimmed[pst.dimmed < 0.0] <- 0.0 
  ps.varstab <- ps            
  otu_table(ps.varstab) <- otu_table(pst.dimmed, taxa_are_rows = FALSE)
  ps	<- ps.varstab
  
  ordination.b <- ordinate(ps, "NMDS", "bray")
  p = plot_ordination(ps, ordination.b, type="sample", color="Horizont", shape = "Site", title="NMDS - bray", axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.1 <- p + stat_ellipse( type="norm", alpha=0.7)
  return(p1.1)
}

ps.s.5.oe.pr <- prune_samples(sample_data(ps.s.5.pr)$Horizont %in% c("O", "E"), ps.s.5.pr)
des.horisont.nh4 <- Des.ls(ps.s.5.oe.pr)
length(which(des.horisont$padj > 0.05))
length(which(des.horisont$padj < 0.05))
length(which(des.horisont.nh4$padj > 0.05))
length(which(des.horisont.nh4$padj < 0.05))

ps.s.5.pr <- prune_taxa(taxa_sums(ps.s.5) > 100, ps.s.5)

library(DESeq2)

diagdds <- phyloseq_to_deseq2(ps.s.5.pr, ~ Horizont )                  
diagdds = estimateSizeFactors(diagdds, type="poscounts")
diagdds = estimateDispersions(diagdds, fitType = "local") 
pst <- varianceStabilizingTransformation(diagdds)
ps.varstab <- t(as.matrix(assay(pst))) 
ps.varstab <- ps.s.5.pr
otu_table(ps.varstab) <- otu_table(pst.dimmed, taxa_are_rows = FALSE) 

cbind.tables <- function(des.horisont){
  temp.df.more <- des.horisont[des.horisont$padj > 0.05,]
  table.df.more <- table(temp.df.more$Phylum)
  sum.more <- sum(table.df.more)
  table.df.more <- table.df.more/sum.more
 temp.df.less <- des.horisont[des.horisont$padj < 0.05,]
  table.df.less <- table(temp.df.less$Phylum)
  sum.less <- sum(table.df.less)
  table.df.less <- table.df.less/sum.less
  table.df <- as.data.frame(cbind(table.df.less, table.df.more))
  table.df <- table.df[rowSums(table.df) > 0,] 
  return(table.df)
}

View(table.df)
table.df <- cbind.tables(des.oe)

library(tibble)

table.fam.all <- rownames_to_column(.data = table.df, var = "Phylum")
dfm <- melt(table.df, id=c(table.fam.all$table.df.less, table.fam.all$table.df.more))
dfm <- cbind(dfm, table.fam.all$Phylum)
p <- ggplot(data=dfm, aes(x=dfm$`table.fam.all$Phylum`, y=value)) + geom_col(aes(fill = variable), position = "dodge") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p <- ggplot(dfm,aes(x = dfm$`table.fam.all$Family`,y = value)) +   geom_col(aes(fill = variable),stat = "identity",position = "dodge") 

#make DeSeq2 table from bf for 5-6-7
#ps.s.bf.pr.100 <- prune_taxa(taxa_sums(ps.s.bf) > 100, ps.s.bf)
Des.tr <- function(ps.s.bf){
  des.bf.all <- Des.ls(ps.s.bf)
  des.bf.all.2 <- des.bf.all
  des.bf.all.2 <- des.bf.all.2[des.bf.all.2$padj < 0.05,] 
  des.bf.all.2 <- des.bf.all.2[des.bf.all.2$baseMean > 30,]
  des.bf.all.2 <- des.bf.all.2[!is.na(des.bf.all.2$padj),]
  return(des.bf.all.2)
}

Des.ch <- function(ps){
  require(qgraph)
  des.all <- Des.tr(ps)
  ps.meangfl <- prune_taxa(rownames(des.all), ps)
  ps.5.6 <- prune_samples(sample_data(ps.meangfl)$Site %in% c("5", "6"), ps.meangfl)
  ps.6.7 <- prune_samples(sample_data(ps.meangfl)$Site %in% c("6", "7"), ps.meangfl)
  des.5.6 <- Des.ls(ps.5.6)
  des.6.7 <- Des.ls(ps.6.7)
  des.res <- cbind(des.5.6$log2FoldChange, des.6.7$log2FoldChange)
  rownames(des.res) <- rownames(des.all)
  des.res[is.na(des.res) == TRUE] <- 0
  Q <- qgraph(cor(t(des.res)), layout = "spring")
  return(des.res)
}



filter.des.tables <- function(des.horisont){
  temp.df.more <- des.horisont[des.horisont$padj > 0.05,]
  table.df.more <- table(temp.df.more$Phylum)
  sum.more <- sum(table.df.more)
  table.df.more <- table.df.more/sum.more
  temp.df.less <- des.horisont[des.horisont$padj < 0.05,]
  table.df.less <- table(temp.df.less$Phylum)
  sum.less <- sum(table.df.less)
  table.df.less <- table.df.less/sum.less
  table.df <- as.data.frame(cbind(table.df.less, table.df.more))
  table.df <- table.df[rowSums(table.df) > 0,] 
  return(table.df)
}


require(DESeq2)
require(phyloseq)
diagdds = phyloseq_to_deseq2(ps.s.bf, ~  Site)                  
diagdds = estimateSizeFactors(diagdds, type="poscounts")
diagdds = estimateDispersions(diagdds, fitType = "local") 
diagdds = DESeq(diagdds)

write.table(des.bf.all, file = "~/storage/svir/LS/not.nice.reads/some.results/des.bf.all.tsv", sep= "\t", col.names = NA, quote=FALSE)
