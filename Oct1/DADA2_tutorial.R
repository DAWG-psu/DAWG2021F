### R studio basic tutorial ###
list.files() #This command shows files in your working directory
getwd()
setwd("~/storage/work/tuc289/DAWG2021F")
list.files()

### You can also set up the working directory by clicking 
###'Session -> Set Working Directory -> Choose directory ..'

##Import data##
table <- read.table("DAWG_tutorial_table.csv", sep = ",", header=T, row.names =1)
rownames(table) #check row names of the table
dim(table) #check dimension of the table
colnames(table) #check column names of the table

##check the data structure and data type
str(table) #check structure of the table

##functions
table <- read.table("DAWG_tutorial_table.csv", sep = ",", header=T, row.names =1)
?read.table

##Output file and export your data
save.image("filename.Rdata")
load("filename.Rdata")

## Save a single object to a file
saveRDS(table, "table.rds")
# Restore it under a different name
table <- readRDS("table.rds")



##DADA2 tutorial pipelines##
##Check DADA2 documents - https://benjjneb.github.io/dada2/index.html

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", update = TRUE)

library(dada2); packageVersion("dada2")

##Let's get our tutorial dataset##

##Click Terminal -> type 'cd /storage/work/[your PSU id]/'
##follow below command lines line by line (copy and paste to Terminal)

##mkdir DAWG2021F
##cd DAWG2021F
##wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1uEfDchKzojpXuFXvhcsJ4wNmVlrhH-iH' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1uEfDchKzojpXuFXvhcsJ4wNmVlrhH-iH" -O fastq.zip && rm -rf /tmp/cookies.txt
##unqip fastq.zip
##rm fastq.zip

##type 'pwd' and you will see '/storage/work/[your PSU ID]/DAWG2021F
##copy that!

##Go back to 'Console'

path <- "/storage/work/tuc289/DAWG2021F"
#change directory with your [PSU ID]
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_1.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
# IF your file has different pattern, you need to change 'pattern'
# For example, if your files looks like 'XXX_XXX_L001_R1_001.fastq' and 'XXX_XXX_L001_R2_001.fastq'
# You should chnage above lines as follow
# fnFs <- sort(list.files(path, pattern= "R1_001.fastq", full.names = TRUE))
# fnRs <- sort(list.files(path, pattern= "R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
# Same here, if your fastq file name has different pattern, you need to edit the line above

#Visualizing the quality profiles of the resd (first two samples, but you can check all of them too)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
#Quality score (phred score) - https://en.wikipedia.org/wiki/Phred_quality_score
# Q = -10 logP or P = 10^(-Q/10)
# Quality score is logarithmically related to the base-calling error probability P
# For example, if Q-score = 10 , meaning that prabability of incorrect base call is 1/10 (90% accuracy)
# If Q-score = 30, meaning that error rate is 1 in 1000 (99.9% accuracy)

#From this plot, we need to decide the position for truncation (truncLen)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
# Understaning 'paste0' function
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,240),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
#maxN, maxEE, trunQ, rm.phix, truncLen .. what are they?
head(out)

#Alternatives: Zymo Research has recently developed a tool called Figaro that can help you choose DADA2 truncation length parameters: 
#https://github.com/Zymo-Research/figaro#figaro


#It takes several minutes
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaFs[[1]]

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaRs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

##Let's download taxonomy database (SILVA v138.1 (updated at March 2021))
##wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
##wget https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# Making phyloseq object for downstream analysis
library(phyloseq); packageVersion("phyloseq")
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               tax_table(taxa))
ps
