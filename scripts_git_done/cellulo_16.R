path = "in/"
path_trein_set = "data/silva_nr_v132_train_set.fa"
path_trein_set_species = "data/silva_species_assignment_v132.fa"
name_Brief = "cellulo.nameBrief.txt"
truncLen = "220,180"
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
  
  fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
  fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
  sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
  on.exit()
  filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
  
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,180), trimLeft=20, maxN=0, maxEE=c(2,5), rm.phix=TRUE, compress=TRUE, multithread=mult)
  errF <- learnErrors(filtFs, multithread=mult)
  errR <- learnErrors(filtRs, multithread=mult)
  derepFs <- derepFastq(filtFs, verbose=TRUE)
  derepRs <- derepFastq(filtRs, verbose=TRUE)
  # Name the derep-class objects by the sample names
  names(derepFs) <- sample.names
  names(derepRs) <- sample.names
  dadaFs <- dada(derepFs, err=errF, multithread=mult, pool=TRUE)
  dadaRs <- dada(derepRs, err=errR, multithread=mult, pool=TRUE)
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
  seqtab <- makeSequenceTable(mergers)
  dim(seqtab)
  table(nchar(getSequences(seqtab)))
  getN <- function(x) sum(getUniques(x))
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=mult, verbose=TRUE)
  
  
  briefToSeq <- colnames(seqtab.nochim)
  names(briefToSeq) <- paste0("Seq", seq(ncol(seqtab.nochim))) 
  st.brief <- seqtab.nochim
  colnames(st.brief) <- names(briefToSeq) 
  
  write.table(briefToSeq, file = name_Brief, sep = "\t")
  
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
  
  write.table(taxa, file = "taxa.txt", sep= "\t", col.names = NA, quote=FALSE)
  
  #some SEPP bash tree magic
  
  filt.otu <-t(as.data.frame(fread("filtered_table.txt")))
  first <-  filt.otu[1,]
  filt.otu <- filt.otu[-c(1),]
  colnames(filt.otu) <-  first
  class(filt.otu) <- "numeric"
  filt.otu.matrix <- as.matrix(filt.otu)
  #head.col <- scan("head.txt", character(), quote = "")
  #rownames(filt.otu.matrix) <- head.col
  tree <- read_tree(treefile="tree.nwk")
  
  mapp <- read.csv("" , header=TRUE, sep="\t")
  map <- data.frame(row.names="ID", mapp)
  
  ps <- phyloseq(otu_table(filt.otu.matrix, taxa_are_rows=FALSE), 
                 sample_data(map), 
                 tax_table(taxa),
                 phy_tree(tree))
  
  