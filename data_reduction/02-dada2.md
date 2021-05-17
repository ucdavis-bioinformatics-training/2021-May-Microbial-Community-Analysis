# Preprocessed reads to phyloseq object using dada2

Some of you should also launch an interactive session on the cluster before performing the below opperation.

```
cd /share/workshop/mca_workshop/$USER
module load R
R
```

**From here on out the commands will be within R**

## R Environment Setup

We are going to use a library directory for this workshop that is located inside the workshop folder. All packages were installed in the last few days.

<div class="r_output">getwd()

.libPaths("/share/workshop/mca_htstream/r_lib")

## CODE to install the needed packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager:::install(c("igraph"))
# BiocManager:::install(c("knitr", "kableExtra", "tidyr")
# BiocManager::install(c("dada2", "biomformat", "phyloseq", "pheatmap"))
# BiocManager::install(c("DECIPHER", "phangorn"))
##

library(dada2); packageVersion("dada2")
library(biomformat); packageVersion("biomformat")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(DECIPHER); packageVersion("DECIPHER")
library(phangorn); packageVersion("phangorn")
library(pheatmap); packageVersion("pheatmap")
library(tidyr); packageVersion("pheatmap")
library(kableExtra); packageVersion("pheatmap")
library(Biostrings); packageVersion("pheatmap")

options(stringsAsFactors=F)
</div>

## First task is to get the samples we will be processing

Here we plan to process all the samples that were preprocessed.

    path <- "01-HTS_Preproc"
    list.files(path)

Read in "SE" files names

    fnSEs = sort(list.files(path, pattern="*_SE.fastq", full.names=T))
    sample.names <-  sapply(strsplit(basename(fnSEs), "_SE"), `[`, 1)

Examine quality profile plots:

    png("fnSEs_quality.png", width=2000, height=1000)
    plotQualityProfile(fnSEs)
    dev.off()


***Note:*** that filtering and trimming should not be necessary if reads are preprocessed with htsteam and filtered for any reads containing an "N".

## DADA2 Denoising

We get right to denoising. DADA2 will generate an error model based on \~100 million basepair (by default the first n reads across N samples that sum to 100M), this error model is to estimate specific error-signature of our run. If your dataset spans multiple runs, this process should be ran for each run separately.


    errU <- learnErrors(fnSEs, multithread=TRUE)

    png("errors-profiles.png", width=2500, height=1500)
        plotErrors(errU, nominalQ=TRUE)
    dev.off()

Then DADA2 denoises the dataset using the error model generated above. First, dereplicating the sequences, then clustering and aligning, then denoiseing. It does this by incorporating the consensus quality profiles and abundances of each sequence, and then figuring out if each sequence is more likely to be of biological origin or a technical artifact.

This step can be run on individual samples, all samples pooled or pseudo. Running it on individual samples is the least computationally intensive, followed by pseudo and the pool. Pooling the samples allow the algorithm to identify a lowly abundant sequence in 1 sample that is more abundant in another.


    dadaUs = dada(fnSEs, err=errU, multithread=TRUE, pool="pseudo")
    gc()


## Generate an ASV count table

Now we aggregate all the denoised sequence data into a single table.

# Construct sequence table:
<div class="r_output">seqtab <- makeSequenceTable(dadaUs)
rownames(seqtab) = sample.names
save(seqtab, file="seqtab.RData")

dim(seqtab)
</div>

## Chimera removal

DADA2 recommends a chimera removal step. This step takes ASVs with very high read count support and then looks for low read count support ASVs that have perfect matches to the high-support ASVs over part of their length. A mixture of 2 other ASVs.

<img src="https://drive5.com/usearch/manual/chimera.gif" alt="chimera" width="250px"/>


<div class="r_output">seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
save(seqtab.nochim, file="seqtab.nochim.RData")

dim(seqtab.nochim)

# What fraction of total reads are retained?
sum(seqtab.nochim)/sum(seqtab)
# [1] 0.9503154


# Sequence lengths before nochim:
table(nchar(getSequences(seqtab)))


# Look at distribution of fragment lengths discarded by the chimera removal step.
table(nchar(setdiff(getSequences(seqtab),getSequences(seqtab.nochim))))
</div>

### Track reads through the pipeline

<div class="r_output">getN <- function(x) sum(getUniques(x))

track = data.frame(input = sapply(dadaUs, getN),
                   denoised = rowSums(seqtab),
                   nochim = rowSums(seqtab.nochim),
                   ASVs = rowSums(seqtab>0),
                   ASVs.nochim = rowSums(seqtab.nochim>0))

rownames(track) <- sample.names[1:3]
track

write.table(data.frame(sample=rownames(track), track), file="read_tracking.tsv", sep='\t', row.names=F)
</div>

### Assign taxonomy to nochim data:


https://zenodo.org/record/4587955#.YKKouZNKgq0

<div class="r_output">dir.create("SILVA")
download.file("https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1", "SILVA/silva_nr99_v138.1_train_set.fa.gz")
download.file("https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz?download=1", "SILVA/silva_species_assignment_v138.1.fa.gz")
</div>

<div class="r_output">taxa.silva.nochim = assignTaxonomy(seqtab.nochim, "SILVA/silva_nr99_v138.1_train_set.fa.gz", multithread = T, minBoot=0)

taxa.silva.nochim = addSpecies(taxa.silva.nochim, "SILVA/silva_species_assignment_v138.1.fa.gz")

save(taxa.silva.nochim, file="taxa.silva.nochim.RData")
gc()
colSums(!is.na(taxa.silva.nochim))
</div>

## Create Amplicon Sequence Variant (ASV) table with counts:
<div class="r_output">
CountTable.silva = data.frame(taxa.silva.nochim[colnames(seqtab.nochim), ], t(seqtab.nochim))
CountTable.silva = CountTable.silva[order(rowSums(t(seqtab.nochim)), decreasing=T), ]
write.table(data.frame(Sequence=rownames(CountTable.silva), CountTable.silva), row.names=F, file="AllSamples_AmpliconSequenceVariant_count_SILVA.tsv", sep='\t')
</div>



<div class="r_output">split_names <- strsplit(sample.names,split="_")
grp <- sapply(split_names, "[[",1L)
temp <- sapply(split_names,"[[",2L)
replicate <- sapply(split_names,"[[",3L)

mdata <- data.frame("SampleID"=sample.names, "Group"=grp, "Temp"=temp, "Replicate"=replicate)
mdata2 = mdata[match(rownames(seqtab.nochim), mdata$SampleID), ]
rownames(mdata2) = mdata2$SampleID
</div>


## Make a phylogenetic tree from all of the ASV sequences

  * optim.pml optimizes the different model parameters, this essentially searches for a better tree using a bunch of stochastic rearrangement strategies.

<div class="r_output">library("DECIPHER")
library("phangorn")

ASVs.nochim = DNAStringSet(colnames(seqtab.nochim))
names(ASVs.nochim) = paste0("ASV", 1:ncol(seqtab.nochim))

alignment = AlignSeqs(DNAStringSet(ASVs.nochim), anchor=NA, processors=30)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA") # turn into phyDat format
dm <- dist.ml(phang.align) #
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)  # computes the likelihood of a phylogenetic tree given a sequence alignment and a model

fitGTR <- update(fit, k=4, inv=0.2)

fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))
save(fitGTR, file="fitGTR-3sample.RData")
</div>

<div class="r_output"># Make Phyloseq objects
library(phyloseq)
library(Biostrings)

load("fitGTR.RData")
load("seqtab.nochim.RData")
load('taxa.silva.nochim.RData')

# create a Phyloseq object with Silva annotation:
ASVs.nochim = DNAStringSet(colnames(seqtab.nochim))
names(ASVs.nochim) = paste0("ASV", 1:ncol(seqtab.nochim))

tmp.seqtab = seqtab.nochim
colnames(tmp.seqtab) = names(ASVs.nochim)
tmp.taxa = taxa.silva.nochim
rownames(tmp.taxa) = names(ASVs.nochim)

ps.silva.nochim = phyloseq(
             otu_table(tmp.seqtab, taxa_are_rows=FALSE),
             sample_data(mdata2),
             tax_table(tmp.taxa),
             refseq(ASVs.nochim),
             phy_tree(fitGTR$tree))
save(ps.silva.nochim, file="phyloseq_nochim_silva.RData")
</div>
