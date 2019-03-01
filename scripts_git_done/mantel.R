ps.gr.Arhaea <- subset_taxa(ps.gr, Kingdom=="Archaea")

mantel.res <- bioenv(otu_table(ps.gr), gr.map.d, method="pearson")
summary(mantel.res)

dist.ps.gr <- distance(ps.gr, "bray")
env.dist <- vegdist(scale(gr.map.num[,9]), "euclid")
mantel(dist.ps.gr, env.dist, method="pearson", permutations = 9999)