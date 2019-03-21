diagdds_more = function(ps){
    diagdds = phyloseq_to_deseq2(ps, ~ Repeats)
    diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
    res = results(diagdds)
    res = res[order(res$padj, na.last=NA), ]
    sigtab = res[(res$padj < 0.1), ]
    sigtab = sigtab[(sigtab$log2FoldChange > 1),]
    sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
    return(sigtab)
}


diagdds_less = function(ps){
    diagdds = phyloseq_to_deseq2(ps, ~ Repeats)
    diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
    res = results(diagdds)
    res = res[order(res$padj, na.last=NA), ]
    sigtab = res[(res$padj < 0.1), ]
    sigtab <- sigtab[(sigtab$log2FoldChange < 1),]
    sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
    return(sigtab)
}

change.prop <-prop.table(table(sigtab$Phylum))

boxplot(log10(assays(diagdds)[["cooks"]]), range=0, las=2)

cts <- counts(dds)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds <- estimateSizeFactors(dds, geoMeans=geoMeans)
    
    
    
cook <- dds@assays[["cooks"]]
                  
For Lise:
                  

diagdds_Lise = function(ps, name){
    diagdds <- phyloseq_to_deseq2(ps, ~ Repeats)
    diagdds <- DESeq(diagdds, test="Wald", fitType="parametric")
    samp <-sample_data(ps)
    dds.counts <- diagdds@assays@.xData$data$counts
    dds.counts.df <- as.data.frame(dds.counts)
    aggdata <- t(aggregate.data.frame(t(dds.counts.df), by=list(samp$Repeats), median))
    colnames(aggdata) <- aggdata[1,]
    aggdata <- aggdata[-1,]
    res = results(diagdds)
    res.df <- as.data.frame(res)
    nice <- cbind(res.df, as.data.frame(tax_table(ps)[rownames(res.df),]), as.data.frame(aggdata)[rownames(res.df),])
    
    return(nice)
}  

d1 <- diagdds_Lise(ps.1, "site1_diagdds.csv")                  
d2 <- diagdds_Lise(ps.2, "site2_diagdds.csv")                   
d3 <- diagdds_Lise(ps.3, "site3_diagdds.csv")                  
d4 <-diagdds_Lise(ps.4, "site4_diagdds.csv")                  
d5 <-diagdds_Lise(ps.15 "site5_diagdds.csv")                   

                  
nice <- cbind(as.data.frame(sigtab), as.data.frame(tax_table(ps.1)[rownames(sigtab),]), as.data.frame(aggdata[rownames(sigtab),])) 
                  
                  EF517956.1.1666
	

                  
                  
                     