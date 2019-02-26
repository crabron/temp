''' use script CRT_Functions_v1.1.R from https://github.com/ShadeLab/ConditionallyRareTaxa.git 
'''

SimpleRareToPrev.f(otu_fp="av_tp_r_otu.txt", abund_thresh=0.01, abund_thresh_ALL=FALSE, b_thresh=0.90, rdp_lastcol=FALSE)


write.table(av_tp_r_dat, "av_tp_r_otu.txt", quote = FALSE, sep = "\t")

av_tp_r_dat <- data.frame(t(otu_table(AY.ps)))
write.table(av_tp_r_dat, "AY.ps.txt", quote = FALSE, sep = "\t")
SimpleRareToPrev.f(otu_fp="AY.ps.txt", abund_thresh=0.001, abund_thresh_ALL=FALSE, b_thresh=0.90, rdp_lastcol=FALSE)





mantel(dist.AY.ps.with.crt, dist.AY.ps.without.crt, method="pearson")
dist.AY.ps.with.crt <- distance(AY.ps.with.crt, "bray")
dist.AY.ps.without.crt <- distance(AY.ps.without.crt, "bray")
AY.ps.with.crt.pcoa <- pcoa(dist.AY.ps.with.crt, correction="none")
AY.ps.without.crt.pcoa <- pcoa(dist.AY.ps.without.crt, correction="none")
p1 <- plot_ordination(AY.ps, AY.ps.without.crt.pcoa, type = "samples", axes = 1:2,color="Site")
p2 <- plot_ordination(AY.ps, AY.ps.with.crt.pcoa, type = "samples", axes = 1:2,color="Site")
mantel(dist.AY.ps.with.crt, dist.AY.ps.without.crt, method="pearson")
p <- ggarrange(p2, p1,  labels= c("CRT","все остальное"), ncol = 2, nrow = 1,label.y = 0.08)