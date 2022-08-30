### Analyse 16S data with DADA2

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.14")

library(dada2); packageVersion("dada2")
library(tidyverse)

setwd("~/My Drive/16S/16S.2/trimmed.truseq")
path <- "~/My Drive/16S/16S.2/trimmed.truseq"

# Forward and reverse fastq filenames have format:
fnFs <- sort(list.files(path, pattern="_1.trimmed.fq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.trimmed.fq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:10]) # cut off 20 bp
plotQualityProfile(fnRs[1:10]) # cut off 70 bp

## Filter and Trim

# Place filtered files in filtered/subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,230), trimLeft=c(20,20), 
                     maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

#write.csv(out, 'filter.trimm.summary.csv')

path2 <- "~/My Drive/16S/16S.1/trimmed.truseq/filtered"
# Forward and reverse fastq filenames have format:
fnFs <- sort(list.files(path2, pattern="_F_filt.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path2, pattern="_R_filt.fastq.gz", full.names = TRUE))

plotQualityProfile(fnFs[1:10])
plotQualityProfile(fnRs[1:10])

## Learn (machine learn) the error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

## Sample inference

dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool = 'pseudo')
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool = 'pseudo')
# inspecting returned dada-class object
dadaFs[[7]]

## Merge paired reads

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

## Construct sequence table

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 300:443]
table(nchar(getSequences(seqtab2)))

## Remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

seqtab2.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab2.nochim)
sum(seqtab2.nochim)/sum(seqtab2)

#write.csv(seqtab2.nochim, "dada2.ASV.sample.csv")

## Track reads through the pipeline

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab2.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#write.csv(track, 'dada2.track.summary.csv')

## Assign taxonomy

taxa <- assignTaxonomy(seqtab2.nochim, "~/My Drive/16S/16S.1/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "~/My Drive/16S/16S.1/silva_species_assignment_v138.1.fa.gz")
taxa.print <- taxa # Removing sequence row names for display only
rownames(taxa.print) <- NULL
head(taxa.print)

## upload genomes
setwd('~/My Drive/WGS/ACassemblies')
ac83 <- getSequences('AC83_scaffolds.fasta')
match83 <- sum(sapply(names(seqtab2.nochim), function(x) any(grepl(x, ac83))))

#write.csv(taxa, 'dada2.16S2.taxa.csv')

tseq <- t(seqtab2.nochim)
trn <- rownames(tseq)
tseq <- as_tibble(tseq)
tseq$ASV <- trn

taxrn <- rownames(taxa)
taxa <- as_tibble(taxa)
taxa$ASV <- taxrn

dada2out <- left_join(taxa, tseq, by = 'ASV')

#write_csv(dada2out, 'dada2.16S2.output.csv')
### Try out the IDTAXA algorithm sometime!!

## Hand off to phyloseq

if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")

library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

samples.out <- rownames(seqtab2.nochim)
environment <- rep(c('endo', 'epi'), each=5)
plant <- rep(1:5, times=2)
samdf <- data.frame(environment=environment, plant=plant)
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(seqtab2.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

plot_richness(ps, x="plant", measures=c("Shannon", "Simpson"), color="environment")+
  theme(axis.text.x = element_text(angle = 0))+
  scale_color_manual(values = c('#BF812D', '#35978F'))
  
plot_richness(ps, x="environment", measures=c("Shannon", "Simpson"))+
  theme(axis.text.x = element_text(angle = 0))

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="environment", title="Bray NMDS")

# Bar plot
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="plant", fill="Order") + facet_wrap(~environment, scales="free_x")

plot_bar(ps, x="plant", fill="Order") + facet_wrap(~environment, scales="free_x")

ps.no <- ps %>% subset_taxa( Family!= "Mitochondria" | is.na(Family) & Order!="Chloroplast" | is.na(Order))
  

plot_bar(ps.no, x="plant", fill="Order") + facet_wrap(~environment, scales="free_x")+
theme(axis.text.x = element_text(angle = 0))

top200 <- names(sort(taxa_sums(ps.no), decreasing=TRUE))[1:200]
ps.top200 <- transform_sample_counts(ps.no, function(OTU) OTU/sum(OTU))
ps.top200 <- prune_taxa(top200, ps.top200)
plot_bar(ps.top200, x="plant", fill="Order") + facet_wrap(~environment, scales="free_x")+
  theme(axis.text.x = element_text(angle = 0))

ps.no@tax_table
