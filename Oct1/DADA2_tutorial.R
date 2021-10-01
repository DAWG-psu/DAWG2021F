### R studio basic tutorial ###
list.files() #This command shows files in your working directory
getwd()
setwd("/storage/work/tuc289/DAWG2021F/")
list.files()

### You can also set up the working directory by clicking 
###'Session -> Set Working Directory -> Choose directory ..'

##Import data##
##Go to 'terminal'
##wget https://github.com/DAWG-psu/DAWG2021F/raw/main/Oct1/DAWG_tutorial_table.csv
table <- read.table("DAWG_tutorial_table.csv", sep = ",", header=T, row.names =1)
rownames(table) #check row names of the table
dim(table) #check dimension of the table
colnames(table) #check column names of the table

##check the data structure and data type
str(table) #check structure of the table
table$SAMPLE_ID <- as.factor(table$SAMPLE_ID)
table$Variable.B..group. <- as.factor(table$Variable.B..group.)
table$Variable.A <- as.factor(table$Variable.A)
str(table)

##Examine data frame
table[1,] ##Only showing first row
table[1:2,] ##Showing first two row
table[,1] ##First column
table[,2:4] ##Second to fourth columns
##IF your dataframe has 'header'
table$SAMPLE_ID

##functions
?read.table

##Output file and export your data
save.image("filename.Rdata")
load("filename.Rdata")

## Save a single object to a file
saveRDS(table, "table.rds")
# Restore it under a different name
table <- readRDS("table.rds")

##Export your plot and table##
write.csv(table, "table_output.csv")
plot(table[,4:5])
##plot -> export -> save as ...
##There are other ways to export high resolution image files




##DADA2 tutorial pipelines##
##Check DADA2 documents - https://benjjneb.github.io/dada2/index.html

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", update = TRUE)
BiocManager::install("ShortRead")
library(dada2); packageVersion("dada2")
library(ShortRead)
##Let's get our tutorial dataset##

##Click Terminal -> type 'cd /storage/work/[your PSU id]/'
##follow below command lines line by line (copy and paste to Terminal)

##mkdir DAWG2021F
##cd DAWG2021F
##wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1uEfDchKzojpXuFXvhcsJ4wNmVlrhH-iH' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1uEfDchKzojpXuFXvhcsJ4wNmVlrhH-iH" -O fastq.zip && rm -rf /tmp/cookies.txt
##unzip fastq.zip
##rm fastq.zip


##type 'pwd' and you will see '/storage/work/[your PSU ID]/DAWG2021F
##copy that!
## ans paste here! : /storage/work/tuc289/DAWG2021F
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

#In gray-scale is a heat map of the frequency of each quality score at each base position. 
#The mean quality score at each position is shown by the green line, 
#and the quartiles of the quality score distribution by the orange lines. 
#The red line shows the scaled proportion of reads that extend to at least that position 

#From this plot, we need to decide the position for truncation (truncLen)

#Check for primers ####
FWD <- "GTGYCAGCMGCCGCGGTAA"
REV <- "GGACTACNVGGGTWTCTAAT"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
# Understaning 'paste0' function
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,240),
                     maxN=0, maxEE=c(2,2), truncQ=2,
                     compress=TRUE, multithread=TRUE,
                     trimLeft = c(19, 20)) # On Windows set multithread=FALSE
#maxN, maxEE, trunQ, rm.phix, truncLen .. what are they?
head(out)

#Alternatives: Zymo Research has recently developed a tool called Figaro that can help you choose DADA2 truncation length parameters: 
#https://github.com/Zymo-Research/figaro#figaro

#check if primers are removed
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = filtFs[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = filtRs[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = filtFs[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = filtRs[[1]]))

#Learning error rates
#It takes several minutes
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

##Main function, 
#The dada function takes as input dereplicated amplicon sequencing reads and returns the inferred composition of the sample (or samples). 
#Put another way, dada removes all sequencing errors to reveal the members of the sequenced community. 

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaFs[[1]]
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaRs[[1]]

#Merge forward and reverse reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])

#Making sequencetable
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Remove chimeric reads (A bimera is a two-parent chimera, in which the left side is made up of one parent sequence, and the right-side made up of a second parent sequence.)
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

