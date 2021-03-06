#for fires article the lastest endlesses code

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
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.1 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.u, type="sample", color, shape, title="NMDS - unifrac", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.2 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  p = plot_ordination(ps, ordination.w, type="sample", color, shape, title="NMDS - wunifrac", 
                      axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
  p1.3 <- p + stat_ellipse( type="norm", alpha=0.7)
  
  #merge by ggpubr
  
  p.all <- ggarrange(p1.1, p1.2, p1.3, ncol = 3 , nrow = 1, common.legend = TRUE, legend = "right")
  
  return(ordination.b)
}
library(ggforce)
p = plot_ordination(ps.f, ord.test.b, type="sample", shape = "Site", title="NMDS - Bray-Curtis", 
                    axes = c(1,2) ) + theme_bw() + theme(text = element_text(size = 14)) + geom_point(size = 3) 
p
p + geom_mark_ellipse(aes(group = Repeats, label = Repeats))



alpha.custom <- function(ps.f.i, arrga = "PD"){
  pd <- pd(ps.f.i@otu_table@.Data, ps.f.i@phy_tree)
  
  pd1 <- estimate_richness(ps.f.i, measures=c("Observed", "InvSimpson", "Shannon"))
  pd <- cbind(pd,ps.f.i@sam_data, pd1 )
  
  bmi <- levels(pd$Вариант)
  bmi.pairs <- combn(seq_along(bmi), 2, simplify = FALSE, FUN = function(i)bmi[i])
  p1 <- ggviolin(pd, x = "Horizont", y = arrga,
                 add = "boxplot", fill = "Horizont") + stat_compare_means(comparisons = bmi.pairs, method = "wilcox.test") +
    scale_x_discrete(limits = c("AY", "AC","C")) + theme(axis.title.x = element_blank(), legend.title = element_blank())
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


x = rnorm(50)
hist(ps.f@otu_table, breaks = seq(min(x), max(x), length.out = 11))


ps.group <- merge_samples(ps.f, "Repeats", fun=median)
otus.td.group <- as.data.frame(t(ps.group@otu_table@.Data))
colnames(otus.td.group)

part.hist.fire <- function(otus, column, nameee){
  x <- otus.td.group[column]/sum(otus.td.group[column])*100
  x <- x[x > 0]
  cut.vals <- cut(x, breaks = c(0,0.005, 0.01,0.05, 0.1, 0.5, 1, 10, Inf), right = FALSE)
  xy <- data.frame(x, cut = cut.vals)
  
  p <- ggplot(xy, aes(x = cut)) +
    theme_bw() + ggtitle(nameee) + 
    geom_bar() +
    scale_x_discrete(drop = FALSE) + theme(axis.text.x = element_text(angle = 90), axis.title = element_blank(), panel.grid = element_blank()) 
  return(p)
}
library(ggpubr)
p1 <- part.hist.fire(otus.td.group, "AY.Control", "AY control")
p2 <- part.hist.fire(otus.td.group, "AY.High_fire", "AY crown fire")
p3 <- part.hist.fire(otus.td.group, "AY.Low_fire","AY surface fire")
p4 <- part.hist.fire(otus.td.group, "AC.Control","AC control")
p5 <- part.hist.fire(otus.td.group, "AC.High_fire","AC crown fire")
p6 <- part.hist.fire(otus.td.group, "AC.Low_fire", "AC surface fire")

p <- ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2)
p

real <- read.csv("Real_eng.csv" , header=TRUE, sep="\t")
real <- data.frame(real)
real$wildfire <- factor(real$wildfire,levels = c("control", "crown", "surface"))
real$horizon <- factor(real$horizon,levels = c("AY", "AC", "C"))
ggplot(real, aes(x=horizon, y=Amount)) + 
  geom_point() + geom_errorbar(real, mapping=aes( ymin = Ошибка_низ, ymax=Ошибка_верх),width=0.1, size=0.5, )+ facet_grid(wildfire~Domain) +
  theme_bw()+ theme(panel.grid = element_blank()) 
library(ggplot2)
