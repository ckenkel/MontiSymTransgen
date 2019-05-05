#####################################

# installing packages: execute just once when first using a script:

source("https://bioconductor.org/biocLite.R")
biocLite("ShortRead")
#OR
#
install.packages("~/Downloads/dada2-1.4", #to install from source, just indicate pkg download location
                 repos = NULL,
                 type = "source",
                 dependencies = c("Depends", "Suggests","Imports"))

#####################################

library(dada2); packageVersion("dada2"); citation("dada2")
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")
library(phyloseq); packageVersion("phyloseq")
library(ShortRead)

#Prior to continuing with R analysis, complete the following steps in terminal
#Prepend year in filename
#Open Terminal window, navigate to folder with fastq.gz files 
#Type:
#for f in *.gz; do mv "$f" "2016-$f"; done

#Trim filenames

#Repeat for 2016 folder

#Argh. Files also differentially named either R1/R2 OR I1/I2
#Find and replace I's with R's
#Type:
#find . -depth -name '*I1*' -execdir bash -c 'mv -i "$1" "${1//I1/R1}"' bash {} \;
#find . -depth -name '*I2*' -execdir bash -c 'mv -i "$1" "${1//I2/R2}"' bash {} \;
#now all files use R1/R2 convention

#Also - some pre-trimming. Retain only PE reads that match amplicon primer. Remove reads containing Illumina sequencing adapters

#Still in terminal, type the following:
#ls *R1_001.fastq | cut -d '_' -f 1 > samples.list
#for file in $(cat samples.list); do  mv ${file}_*R1*.fastq ${file}_R1.fastq; mv ${file}_*R2*.fastq ${file}_R2.fastq; done 
#for file in $(cat samples.list); do bbduk.sh in1=${file}_R1.fastq in2=${file}_R2.fastq ref=adaptors.fasta k=12 out1=${file}_R1_NoIll.fastq out2=${file}_R2_NoIll.fastq ; done &>bbduk_NoIll.log
#for file in $(cat samples.list); do bbduk.sh in1=${file}_R1_NoIll.fastq in2=${file}_R2_NoIll.fastq restrictleft=21 k=10 literal=GTGAATTGCAGAACTCCGTG,CCTCCGCTTACTTATATGCTT outm1=${file}_R1_NoIll_NoITS.fastq outu1=${file}_R1_check.fastq outm2=${file}_R2_NoIll_NoITS.fastq outu2=${file}_R2_check.fastq; done &>bbduk_NoITS.log

#seq-tk for random subsampling of read data...to verify that 2015 run is not different from 2016
##################For random subsampling of 2016 dataset to bring it down to 2015 size

#TooHighSams.list is list of 2016 sample names that need to be subsampled to 7000 reads; note s100 sets same random seed so that the same F and reverse reads are selected. In case below, re-setting random seed each time, but retaining it for both R1 and R2 processing
#for file in $(cat TooHighSams.list); do n=$RANDOM; echo $n; seqtk sample -s $n ${file}_R1_NoIll_NoITS.fastq 7000 >${file}_R1_sub.fastq; seqtk sample -s $n ${file}_R2_NoIll_NoITS.fastq 7000 >${file}_R2_sub.fastq; done 
#Replace original samples with these _sub files and run the analysis below as before


############Once above preprocessing steps completed, begin DADA2 anaylsis here

#Set path to trimmed, renamed fastq files
path <- "~/Dropbox/AIMSpostdoc/KateMontiSpawnExpt/MontiAllFastq" # CHANGE ME to the directory containing the fastq files after unzipping.
fns <- list.files(path)
fns

################################
##### Trimming/Filtering #######
################################

fastqs <- fns[grepl(".fastq$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
fnRs <- fastqs[grepl("_R2", fastqs)] # Just the reverse read files


# Get sample names, assuming files named as so: SAMPLENAME_XXX.fastq; OTHERWISE MODIFY
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1) #the last number will select the field for renaming
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)


#########Visualize Raw data

#First, lets look at quality profile of R1 reads
 
plotQualityProfile(fnFs[c(1,2,3,4)])
plotQualityProfile(fnFs[c(111,112,113,114)])

#Then look at quality profile of R2 reads

quartz()
plotQualityProfile(fnRs[c(1,2,3,4)])
plotQualityProfile(fnRs[c(111,112,113,114)])

#The reverse reads are worse quality for the 2015 data, especially at the end, which is common in Illumina sequencing. 2016 looks really good. This isn’t too worrisome, DADA2 incorporates quality information into its error model which makes the algorithm more robust, but trimming as the average qualities crash is still a good idea as long as our reads will still overlap. 
#For Pochon ITS2 primers, have 160 bp overlap. Can trim quite a bit

#The distribution of quality scores at each position is shown as a grey-scale heat map, with dark colors corresponding to higher frequency. 
#green is the mean, orange is the median, and the dashed orange lines are the 25th and 75th quantiles.
#Recommend trimming where quality profile crashes - in this case, forward reads mostly fine up to 210; for reverse >160 bases it gets below 30; still should leave ~30-50bp overlap

#If using this workflow on your own data: Your reads must still overlap after truncation in order to merge them later! If you are using a less-overlapping primer set, your truncLen must be large enough to maintain 20 + biological.length.variation nucleotides of overlap between them.
#BUT: For common ITS amplicon strategies, it is undesirable to truncate reads to a fixed length due to the large amount of length variation at that locus. That is OK, just leave out truncLen. Make sure you removed the forward and reverse primers from both the forward and reverse reads though! 
 #Not too many indels in ITS2; at least not super long ones

# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "trimmedSubTest")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Filter
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, #truncLen=c(210,160), #leaves ~30bp overlap
              maxN=0, #DADA does not allow Ns
              maxEE=c(1,1), #allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
              truncQ=2, 
              trimLeft=c(20,21), #N nucleotides to remove from the start of each read: ITS2 primers = F 20bp; R 21bp 
              rm.phix=TRUE, #remove reads matching phiX genome
              matchIDs=TRUE, #enforce matching between id-line sequence identifiers of F and R reads
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

head(out)
tail(out)


#A word on Expected Errors vs a blanket quality threshold
#Take a simple example: a read of length two with quality scores Q3 and Q40, corresponding to error probabilities P=0.5 and P=0.0001. The base with Q3 is much more likely to have an error than the base with Q40 (0.5/0.0001 = 5,000 times more likely), so we can ignore the Q40 base to a good approximation. Consider a large sample of reads with (Q3, Q40), then approximately half of them will have an error (because of the P=0.5 from the Q2 base). We express this by saying that the expected number of errors in a read with quality scores (Q3, Q40) is 0.5.
#As this example shows, low Q scores (high error probabilities) dominate expected errors, but this information is lost by averaging if low Qs appear in a read with mostly high Q scores. This explains why expected errors is a much better indicator of read accuracy than average Q.  
################################
##### Learn Error Rates #######
################################
#DADA2 learns its error model from the data itself by alternating estimation of the error rates and the composition of the sample until they converge on a jointly consistent solution (this is similar to the E-M algorithm)
#As in many optimization problems, the algorithm must begin with an initial guess, for which the maximum possible error rates in this data are used (the error rates if only the most abundant sequence is correct and all the rest are errors).


setDadaOpt(MAX_CONSIST=30) #increase number of cycles to allow convergence
errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)

#sanity check: visualize estimated error rates
#error rates should decline with increasing qual score
#red line is based on definition of quality score alone
#black line is estimated error rate after convergence
#dots are observed error rate for each quality score

plotErrors(errF, nominalQ=TRUE) #some issues with C2G and G2C variants being underestimated, but not terrible
plotErrors(errR, nominalQ=TRUE) #again, worse with G2C and C2G; a little T2G, but rest err on the side of being conservative (above red line)

#why do values increase at Q40 in some plots?
#this artefact exists b/c in many sequencing runs there are almost no Q=40 bases. The loess smoothing hits the edge and a lack of observations, causing weird behavior. BUT as there essentially aren't (almost) any Q=40 bases to correct anyway and at the worst, the error rates are overestimated, so it's actually conservative for calling new variants

################################
##### Dereplicate reads #######
################################
#Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance”: the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#DADA2 retains a summary of the quality information associated with each unique sequence. The consensus quality profile of a unique sequence is the average of the positional qualities from the dereplicated reads. These quality profiles inform the error model of the subsequent denoising step, significantly increasing DADA2’s accuracy.
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


################################
##### Infer Sequence Variants #######
################################

#Must change some of the DADA options b/c original program optomized for ribosomal data, not ITS - from github, "We currently recommend BAND_SIZE=32 for ITS data." leave as default for 16S/18S

setDadaOpt(BAND_SIZE=32)

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#now, look at teh dada class objects by sample
#will tell how many 'real' variants in unique input seqs
#By default, the dada function processes each sample independently, but pooled processing is available with pool=TRUE and that may give better results for low sampling depths at the cost of increased computation time. See our discussion about pooling samples for sample inference. 
dadaFs[[109]]
dadaRs[[109]]

################################
##### Merge paired reads #######
################################

#To further cull spurious sequence variants
#Merge the denoised forward and reverse reads
#Paired reads that do not exactly overlap are removed

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[119]])

summary((mergers[[119]]))

#We now have a data.frame for each sample with the merged $sequence, its $abundance, and the indices of the merged $forward and $reverse denoised sequences. Paired reads that did not exactly overlap were removed by mergePairs.

################################
##### Construct sequence table #######
################################
#a higher-resolution version of the “OTU table” produced by classical methods

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

plot(table(nchar(getSequences(seqtab)))) #real variants appear to be right in that 294-304 window

#The sequence table is a matrix with rows corresponding to (and named by) the samples, and 
#columns corresponding to (and named by) the sequence variants. 
#Do merged sequences all fall in the expected range for amplicons? ITS2 Pochon ~340bp-41bp primers; accept 294-304
#Sequences that are much longer or shorter than expected may be the result of non-specific priming, and may be worth removing

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(294,304)] #again, being fairly conservative wrt length

table(nchar(getSequences(seqtab2)))
dim(seqtab2)

################################
##### Remove chimeras #######
################################
#The core dada method removes substitution and indel errors, but chimeras remain. 
#Fortunately, the accuracy of the sequences after denoising makes identifying chimeras easier 
#than it is when dealing with fuzzy OTUs: all sequences which can be exactly reconstructed as 
#a bimera (two-parent chimera) from more abundant sequences.

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab2)
#The fraction of chimeras varies based on factors including experimental procedures and sample complexity, 
#but can be substantial. Here chimeras make up about 36% of the inferred sequence variants (138-89 = 49 => 49/138), 
#BUT those variants account for only about 0.5% of the total sequence reads
#Most of your reads should remain after chimera removal (it is not uncommon for a majority of sequence variants to be removed though)

################################
##### Track Read Stats #######
################################

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab2), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
tail(track) #again, 2016 better than 2015

write.csv(track,file="ReadFilterStats_SubsampledData_21Nov17.csv",row.names=TRUE,quote=FALSE)


################################
##### Assign Taxonomy #######
################################

#It is common at this point, especially in 16S/18S/ITS amplicon sequencing, to classify sequence variants taxonomically. 
#DADA2 provides a native implementation of the RDP's naive Bayesian classifier. The assignTaxonomy function takes a set of sequences and a training set of taxonomically classified sequences, and outputs the taxonomic assignments with at least minBoot bootstrap confidence.
#Here, I have supplied a modified version of the GeoSymbio ITS2 database (Franklin et al. 2012)
#
taxa <- assignTaxonomy(seqtab.nochim, "~/Dropbox/AIMSpostdoc/KateMontiSpawnExpt/MontiAllFastq/Training/GeoSymbio_ITS2_LocalDatabase_verForPhyloseq.fasta", minBoot=5,multithread=TRUE,tryRC=TRUE,outputBootstraps=FALSE)
unname(head(taxa, 30))
unname(taxa)

#Lowered bootstrap threshold from 50 to 5. Was not returning hits for many sequences. But reducing to 5 improved sequence return and identities largely match separate blastn search against the same database

#Now, save outputs so can come back to the analysis stage at a later point if desired
saveRDS(seqtab.nochim, file="21Nov_SubSam_seqtab_nochim.rds")
saveRDS(taxa, file="21_Nov_SubSam_taxa_blastCorrected.rds")

#If you need to read in previously saved datafiles
seqtab.nochim <- readRDS("21Sep_seqtab_nochim.rds")
taxa <- readRDS("21_Sep_taxa_blastCorrected.rds")

#back in terminal, run blast against Symbio database to compare
#makeblastdb -in GeoSymbio_ITS2_LocalDatabase.fasta -dbtype nucl
#blastn -query Sep21_OTUs_All.fasta -db /Users/drcarl/Dropbox/AIMSpostdoc/KateMontiSpawnExpt/MontiAllFastq/Training/GeoSymbio_ITS2_LocalDatabase_verForPhyloseq.fasta -num_descriptions 5 -num_alignments 5 -out NODES_All.br
#grep -A 12 'Query=' NODES_All.br

##########################
####using LULU to collapse intragenomic variants if present
####Problem: co-occurrence of ASVs across samples is expected for this dataset! We're  
####looking at heritability of symbionts in parents and eggs.
####therefore, just use co-occurrence in independent adult coral samples across years  
####to remove highly sequence similar, co-occurring ITS2 types, if present
############

#######IN TERMINAL#########

#First produce a blastdatabase with the ASVs
#makeblastdb -in Sep21_OTUs_All.fasta -parse_seqids -dbtype nucl

# Then blast the ASVs against the database to produce the match list
# which provides information about sequence similarity for collapsing
# ITS2 diversity 
#blastn -db Sep21_OTUs_All.fasta -outfmt '6 qseqid sseqid pident' -out match_list.txt -qcov_hsp_perc 90 -perc_identity 84 -query Sep21_OTUs_All.fasta

#####Back in R###### 


#first, read in ASV table
alldat<-read.csv("Sep21_OutputDADA_AllOTUs.csv")
head(alldat)

#And match list
matchList<-read.table("match_list.txt")
head(matchList)

#Reformat ASV table to desired LULU format
rownames(alldat)<-alldat$X
ASVs<-data.frame(t(alldat[,2:90])) 
head(ASVs)

#Use Parent corals only to ID intragenomic types 
sub<-ASVs[,c(1,13,26,39,52,65,71,83,96,109,110,111,124,125,126,139,140,141,154,155,156)]

#are zeros messing up lulu algorithm? Yes. 
#Must remove zero count ASVs (those which only appear in eggs). 
#these are all low low confidence types anyway (don't appear in many samples period)
sub$zeros = apply(sub[,1:21],1,function(x){sum(x<1)}) #number of ASVs with less than 1 read per sample (i.e. the number of zero count samples)
sub2<-sub[-which(sub$zeros>20),] #if all 21 samples are zero count, toss these ASVs

#Now, run the LULU curation

curated_result <- lulu(sub2, matchList,minimum_relative_cooccurence=0.7) #using just 70% co-occurrence across samples and 84% sequence similarity. => Samples are highly similar sequence-wise, problem is co-occurrence. these types just don't co-occur frequently enough. 
summary(curated_result)#does not discard any ASVs! 

#nothing to collapse. Proceed with full dataset

################################
##### handoff 2 phyloseq #######
################################

#import dataframe holding sample information
samdf<-read.csv("monti2yrs_VARIABLESTable_forphyloseq.csv")
head(samdf)
rownames(samdf) <- samdf$Samplename

# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

ps


#Visualize alpha-diversity - ***Should be done on raw, untrimmed dataset***
#total species diversity in a landscape (gamma diversity) is determined by two different things, the mean species diversity in sites or habitats at a more local scale (alpha diversity) and the differentiation among those habitats (beta diversity)
#Shannon:Shannon entropy quantifies the uncertainty (entropy or degree of surprise) associated with correctly predicting which letter will be the next in a diverse string. Based on the weighted geometric mean of the proportional abundances of the types, and equals the logarithm of true diversity. When all types in the dataset of interest are equally common, the Shannon index hence takes the value ln(actual # of types). The more unequal the abundances of the types, the smaller the corresponding Shannon entropy. If practically all abundance is concentrated to one type, and the other types are very rare (even if there are many of them), Shannon entropy approaches zero. When there is only one type in the dataset, Shannon entropy exactly equals zero (there is no uncertainty in predicting the type of the next randomly chosen entity).
#Simpson:equals the probability that two entities taken at random from the dataset of interest represent the same type. equal to the weighted arithmetic mean of the proportional abundances pi of the types of interest, with the proportional abundances themselves being used as the weights. Since mean proportional abundance of the types increases with decreasing number of types and increasing abundance of the most abundant type, λ obtains small values in datasets of high diversity and large values in datasets of low diversity. This is counterintuitive behavior for a diversity index, so often such transformations of λ that increase with increasing diversity have been used instead. The most popular of such indices have been the inverse Simpson index (1/λ) and the Gini–Simpson index (1 − λ).

plot_richness(ps, x="ColonyID", measures=c("Shannon", "Simpson"), color="Type") + theme_bw()

#Ordinate Samples
ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray",k=20)


#Exploratory Bar-plots

top30 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:30]
ps.top30 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top30 <- prune_taxa(top30, ps.top30)
plot_bar(ps.top30, x="Type", fill="Class") + facet_wrap(~ColonyID+Timepoint, scales="free_x")

btm30 <- names(sort(taxa_sums(ps), decreasing=FALSE))[1:30]
ps.btm30 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.btm30 <- prune_taxa(btm30, ps.btm30)
plot_bar(ps.btm30, x="Type", fill="Class") + facet_wrap(~ColonyID+Timepoint, scales="free_x")



################################
##### output 'OTU' table #######
################################

#seqtab.nochim is the 'OTU' table...but is a little unwieldy
#For Symbiodinium, sequence classification is not so great...
#want fasta file of 'OTUs' and table designated by 'OTU'

#First, output fasta file for 'OTUs'
path='~/Dropbox/AIMSpostdoc/KateMontiSpawnExpt/Sep21_OTUs_All.fasta'
uniquesToFasta(seqtab.nochim, path, ids = NULL, mode = "w", width = 20000)

#then, rename output table and write it out
ids <- paste0("sq", seq(1, length(colnames(seqtab.nochim))))
colnames(seqtab.nochim)<-ids

write.csv(seqtab.nochim,file="Sep21_OutputDADA_AllOTUs.csv",quote=F)

str(seqtab.nochim)

#For our purposes, we also want to focus on the corals used in both 2015 and 2016
#subset data
focus = subset_samples(ps, Focus== "Yes")

seqtab<-otu_table(focus)
ids <- paste0("sq", seq(1, length(colnames(seqtab))))
colnames(seqtab)<-ids
head(seqtab)
write.csv(seqtab,file="Nov21_SubSam_OutputDADA_AllOTUs_FocusYesOnly.csv",quote=F)

#############################################

# Principal coordinate analysis 

library(vegan)
library(MCMC.OTU)
library(ggfortify)
library(cluster)
library(labdsv)
#if necessary, rerun lines 273-296 above to regenerate ps

alldat<-read.csv("Sep21_OutputDADA_AllOTUs_FocusYesOnly_forMCMC.csv") #Sep21_OutputDADA_AllOTUs_FocusYesOnly_forMCMC.csv
names(alldat)
str(alldat)
head(alldat) #Note: 2015-M11r11 has no read data left. Must remove from dataset or zero.cut eval will fail below
dat<-alldat[c(1:3,5:110),] #remove r11

#props=apply(dat[,c(5:93)],2,function(x){sum(x)/sum(dat[,c(5:93)])})
#tdat<-as.data.frame(t(dat[5:93]))
#tdat$zeros = apply(tdat,1,function(x){sum(x<1)}) #making new column counting number of samples with counts <=1 within each isogroup (host)
#tdat$low= tdat$zeros/110

# purging under-sequenced samples; and ASVs represented less than 3 unique times
goods2=purgeOutliers(dat,count.columns=5:93,otu.cut=0,zero.cut=0.027) #cols=5:93; zero cut=0.027; also 2 recruit samples removed for having too little read data

taxa_names(ps)<-names(dat)[5:93]

#subset to fully replicated samples, with sufficient read depth (based on purging step above)

ps<-prune_samples(sample_data(ps)$Samplename %in% goods2$sample, ps)

# creating a log-transfromed normalized dataset for PCoA:
goods.log=logLin(data=goods2,count.columns=5:length(names(goods2)))

rownames(goods.log)<-goods2$sample

#replace phyloseq ASV table with logLin transform
OTU<-otu_table(goods.log,taxa_are_rows=FALSE)
otu_table(ps)<-OTU


#plot PCoA

ord<-ordinate(ps,"PCoA","manhattan")


p1=plot_ordination(ps,ord,type="taxa",color="Class")
print(p1)

p2=plot_ordination(ps,ord,type="samples",color="Time",shape="Time")+geom_point(size=3,alpha=0.75)+scale_colour_brewer(type="qual",palette="Set2")
p2

p3=plot_ordination(ps,ord,type="samples",color="Mom",shape="Mom")+geom_point(size=3,alpha=0.75)+scale_colour_brewer(type="qual",palette="Set1")
p3

p4=plot_ordination(ps,ord,type="samples",color="Type",shape="Type")+geom_point(size=3,alpha=0.75)+scale_colour_brewer(type="qual",palette="Set1")
p4

require(gridExtra)
grid.arrange(p2,p3,p4)

unname(taxa)



################################
###### DESeq for Stats #######
################################


#load deseq

library("DESeq"); packageVersion("DESeq")
library(genefilter)


#must swap rows/columns to match DESeq format

head(goods2)
rownames(goods2)<-goods2$sample
names(goods2)

counts<-data.frame(t(goods2[,5:38])) #5:38 only the columns with count data

########################Creating table of conditions for your experiment, 

Year=Type=Family=c(1:length(names(counts)))
Year[grep("2015",names(counts))]="yr2015"
Year[grep("2016",names(counts))]="yr2016" 
Type[grep("r",names(counts))]="larvae"
Type[grep("r",names(counts),invert=TRUE)]="adult"
m7<-c("M7","6.7")
m9<-c("M9","6.9")
m11<-c("M11","6.11")
Family[grep(paste(m7,collapse="|"),names(counts))]="M7"
Family[grep(paste(m9,collapse="|"),names(counts))]="M9"
Family[grep(paste(m11,collapse="|"),names(counts))]="M11"
Family[grep("24",names(counts))]="M24"

conditions=data.frame(cbind(Year,Type,Family))
head(conditions)

real=newCountDataSet(counts,conditions) 
real=estimateSizeFactors(real)

sizeFactors(real)


#build the DESeq object
real=estimateDispersions(real,sharingMode="gene-est-only",method="pooled",fitType="parametric")  #using a pooled dispersion estimate, parametric fit; NOTE: must use local fit for subsampled dataset
#can use gene-est-only for ASV data, should have >7 reps per factor level

goods=t(counts(real,normalized=TRUE))

# what is the proportion of samples with data for these ASVs?
withData=apply(goods,2,function(x){sum(x>0)/length(x)})

hist(withData) #Most ASVs not well represented across samples; there are 5 with counts in more than 60% of samples

# what percentage of total counts does each ASV represent?
props=apply(goods,2,function(x){sum(x)/sum(goods)})

barplot(props)

props#: sq1=66%; sq2=17%; sq3=5.2%; sq4=5.0%; sq5=1.4%; rest<1% NOTE: for subsampled dataset, recovered ASV proportions are almost identical



 vsd=getVarianceStabilizedData(real) #variance stabilized counts => good for heatmaps

 normal=counts(real,normalized=TRUE) #just extracting normalized counts for abundance bargraphs


######################## Now for the real Model Testing - account for family; but interested in type*year effects; allow 30 iterations

fit0=fitNbinomGLMs(real, count ~ Family, glmControl=list(maxit=30)) 
fit1=fitNbinomGLMs(real, count ~ Family+Type, glmControl=list(maxit=30))
fit2=fitNbinomGLMs(real, count ~ Family+Year, glmControl=list(maxit=30))

fit3=fitNbinomGLMs(real, count ~ Family+Type+Year, glmControl=list(maxit=30))
fit4=fitNbinomGLMs(real, count ~ Family+Type*Year, glmControl=list(maxit=30))

fit5=fitNbinomGLMs(real, count ~ Year, glmControl=list(maxit=30))

# testing section

pvals.y<-nbinomGLMTest(fit2,fit0) #testing significance of year term
pvals.t<-nbinomGLMTest(fit1,fit0) #testing significance of type term
pvals.i<-nbinomGLMTest(fit4,fit3) #testing significance of interaction term


#adjusting for multiple testing AND
#making non-convergent model p-values NA's -NOTE: both models must have converged 
adjp.t<-p.adjust(pvals.t,method="BH")
adjp.t=data.frame(adjp.t)
pvals.t=data.frame(pvals.t)
converged<-as.data.frame(cbind(fit0$converged,fit1$converged))
converged$test<-converged$V1+converged$V2
rownames(fit1)->rownames(pvals.t); rownames(fit1)->rownames(converged);rownames(fit1)->rownames(adjp.t);
converged$test<-apply(converged,1,all)
badmods<-subset(converged,(!converged$test))
for ( i in rownames(badmods)){pvals.t[i,1]<-NA}
for ( i in rownames(badmods)){adjp.t[i,1]<-NA}

adjp.y<-p.adjust(pvals.y,method="BH")
adjp.y=data.frame(adjp.y)
pvals.y=data.frame(pvals.y)
converged<-as.data.frame(cbind(fit0$converged,fit2$converged))
converged$test<-converged$V1+converged$V2
rownames(fit2)->rownames(pvals.y); rownames(fit2)->rownames(converged);rownames(fit2)->rownames(adjp.y);
converged$test<-apply(converged,1,all)
badmods<-subset(converged,(!converged$test))
for ( i in rownames(badmods)){pvals.y[i,1]<-NA}
for ( i in rownames(badmods)){adjp.y[i,1]<-NA}


adjp.i<-p.adjust(pvals.i,method="BH")
adjp.i=data.frame(adjp.i)
pvals.i=data.frame(pvals.i)
converged<-as.data.frame(cbind(fit4$converged,fit3$converged))
rownames(fit4)->rownames(pvals.i); rownames(fit4)->rownames(converged);rownames(fit4)->rownames(adjp.i);
converged$test<-apply(converged,1,all)
badmods<-subset(converged,(!converged$test))
for ( i in rownames(badmods)){pvals.i[i,1]<-NA}
for ( i in rownames(badmods)){adjp.i[i,1]<-NA}

summary(adjp.t)
summary(adjp.y)
summary(adjp.i)


#creating table of all multiple test corrected p-values with variance stabilized count data 
PPV_VSD<-cbind(vsd, "adjp.t" = adjp.t$adjp.t, "adjp.y" = adjp.y$adjp.y,"adjp.i" = adjp.i$adjp.i,"pval.t" = pvals.t$pvals.t, "pval.y" = pvals.y$pvals.y, "pval.i" = pvals.i$pvals.i)  

PPV_Norm<-cbind(normal, "adjp.t" = adjp.t$adjp.t, "adjp.y" = adjp.y$adjp.y,"adjp.i" = adjp.i$adjp.i,"pval.t" = pvals.t$pvals.t, "pval.y" = pvals.y$pvals.y, "pval.i" = pvals.i$pvals.i)  


write.csv(PPV_VSD, file="VSDandPVALS_FocusOnly_SubSam_nov21.csv", quote=F) #writing an output file of vsd plus p-values
write.csv(PPV_Norm, file="NormCtsandPVALS_FocusOnly_SubSam_nov21.csv", quote=F) #

########################## counting, venn diagram:

p<-data.frame(PPV_VSD) #OR
p<-read.csv("VSDandPVALS_FocusOnly_oct19.csv"); rownames(p)<-p$X

inter=row.names(p[p$adjp.i<=0.05 & !is.na(p$adjp.i),])#all rows where adjusted pvalue is less than or equal to 0.05 and is not an NA
year=row.names(p[p$adjp.y<=0.05 & !is.na(p$adjp.y),])
type=row.names(p[p$adjp.t<=0.05 & !is.na(p$adjp.t),])
none=c("sq1","sq6","sq7","sq9","sq13","sq14","sq15","sq16","sq17","sq19","sq22","sq23","sq24","sq25","sq26","sq27","sq28","sq29","sq31","sq37","sq50","sq54")


candidates=list("Age"=type,"Year"=year,"Age*Year"=inter,"None"=none)
library(gplots)
quartz()
venn(candidates)
inter
year
type

#

#############################################
# exploring correlations between sig ASVs 

# creating a log-transformed normalized dataset ignoring zero counts:
nl=startedLog(data=goods2,count.columns=c(6:9,12,14:16,22,24,25,35),logstart=0)


panel.cor.pval2<-function (x, y, digits = 2, cex.cor, p.cut = 1) 
{
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor.test(x, y, use = "na.or.complete")$p.value
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    if (missing(cex.cor)) 
        cex.cor <- 0.8/strwidth(txt)
    if (r > p.cut) {
        txt = ""
    }
    text(0.5, 0.5, txt)
}


# displaying a matrix of scatterplots and p-values of ASV correlations
# (onlu p-values better than 0.1 are displayed)
pairs(nl,lower.panel=panel.smooth,upper.panel=panel.cor.pval2)

pvals=c(2.1e-06,0.005, 0.4, 1.2e-05,1.9e-07,0.00055, 0.222, 0.25, 0.0015, 0.12, 0.0039, 4.4e-08,0.0014, 5.5e-10,2.4e-05,0.035,0.0091,0.14,0.0019,0.43,0.089,0.49,0.0057,0.54,0.11,0.051,0.13,0.11,0.70,0.72,0.0041,0.00061,0.0084,0.07,0.90,0.024,0.35,0.00056,0.00097,0.01,0.15,0.38,0.00013,0.60,0.11,0.0087,0.98,0.086,0.00005,0.61,0.78,0.64,0.012,0.00054,0.034,0.21,0.0047,0.23,0.16,0.64,0.12,0.0011,0.15,0.21,0.31,0.77)
p.adjust(pvals,method="BH")

pairs(nl,lower.panel=panel.smooth,upper.panel=panel.cor) 
#need to do multiple test correction on pvals

#----------------------------------------


###############################################
##### plotting  #######
###############################################

#rerun lines 273-296 and 365-383 above if need be to regenerate dat and ps
#rename taxa (going back to phyloseq object)


#Try phyloseq PCoA step here...or nmds?

#subset to significant ASVs (see DESeq analysis below with MCMC pre-filter step, Venn diagram)

all<-c("sq11","sq4"  ,"sq8",   "sq5" ,"sq10", "sq21" ,"sq3", "sq12", "sq2" ,"sq20","sq18" ,"sq33")

SigYear<-prune_taxa(taxa_names(ps) %in% all, ps)

#establish sample order for plotting

SamSortY<-rownames(sample_data(SigYear)[order(sample_data(SigYear)$Time,sample_data(SigYear)$Mom,sample_data(SigYear)$Type),])
SamSortA<-rownames(sample_data(SigYear)[order(sample_data(SigYear)$Type,sample_data(SigYear)$Mom),])

#plot a sample heatmap
theme_set(theme_bw())
plot_heatmap(SigYear,method=NULL,sample.order=SamSortA,taxa.order=rev(all))


###################Plot a barplot

p<-read.csv("NormCtsandPVALS_FocusOnly_oct19.csv"); rownames(p)<-p$X
normDat<-as.data.frame(t(p[,2:108]))
head(normDat)
rownames(normDat)<-rownames(otu_table(ps))

#replace phyloseq ASV table with normalized counts
OTU<-otu_table(normDat,taxa_are_rows=FALSE)
otu_table(ps)<-OTU

#subset to significant ASVs (see DESeq analysis below with MCMC pre-filter step, Venn diagram)

all<-c("sq11","sq4"  ,"sq8",   "sq5" ,"sq10", "sq21" ,"sq3", "sq12", "sq2" ,"sq20","sq18" ,"sq33")

SigYear<-prune_taxa(taxa_names(ps) %in% all, ps)

otus<-rownames(tax_table(SigYear))

# over-write the former kingdom slot in the tax table with the functional group
tax_table(SigYear)[,"Kingdom"] <- otus


plot_bar(SigYear,fill="Kingdom")+labs(fill='SeqVars')


###################Plot a CarlPlot: VSD values

p<-read.csv("VSDandPVALS_FocusOnly_oct19.csv"); rownames(p)<-p$X #approximates log2 transform
normDat<-as.data.frame(t(p[,2:108]))
head(normDat)
rownames(normDat)<-rownames(otu_table(ps))

#to regenerate ps, go to lines 273-293 above

#replace phyloseq ASV table with normalized counts
OTU<-otu_table(normDat,taxa_are_rows=FALSE)
otu_table(ps)<-OTU

#subset to significant ASVs (see DESeq analysis below with MCMC pre-filter step, Venn diagram)

all<-c("sq11","sq4"  ,"sq8",   "sq5" ,"sq10", "sq21" ,"sq3", "sq12", "sq2" ,"sq20","sq18" ,"sq33")

SigYear<-prune_taxa(taxa_names(ps) %in% all, ps)

otus<-rownames(tax_table(SigYear))

# over-write the former kingdom slot in the tax table with the functional group
tax_table(SigYear)[,"Kingdom"] <- otus

df<-psmelt(SigYear)
df$Type<-factor(df$Type,levels=c("egg","adult"))
head(df)


source("summarySE.R")

# int<-subset(df,Kingdom=="sq4"|Kingdom=="sq5"|Kingdom=="sq8")
# yr<-subset(df,Kingdom=="sq3"|Kingdom=="sq10"|Kingdom=="sq11"|Kingdom=="sq12"|Kingdom=="sq18"|Kingdom=="sq20"|Kingdom=="sq21"|Kingdom=="sq33")
# stage<-subset(df,Kingdom=="sq2"|Kingdom=="sq4")


all.summ=summarySE(df,measurevar="Abundance",groupvars=c("OTU","Time","Type"))

pd <- position_dodge(.3)
ggplot(all.summ,aes(x=Time,y=Abundance,fill=Type))+
	geom_point(aes(group=Type,pch=Type),position=pd,size=2.5)+
	geom_errorbar(aes(ymin=Abundance-se,ymax=Abundance+se),lwd=0.4,width=0.3,position=pd)+
	geom_line(aes(group=Type,linetype=Type),position=pd)+
	scale_color_manual(values=c("lightblue","deepskyblue"))+
	facet_wrap(~OTU,scales="free_y",ncol=6)+
	theme_bw()

all.summ=summarySE(yr,measurevar="Abundance",groupvars=c("OTU","Time"))

quartz()
ggplot(all.summ,aes(x=Time,y=Abundance))+
	geom_point(aes(),position=pd,size=2.5)+
	geom_errorbar(aes(ymin=Abundance-se,ymax=Abundance+se),lwd=0.4,width=0.3,position=pd)+
	geom_line(aes(),position=pd)+
	#scale_fill_manual(values=c("lightblue","deepskyblue"))+
	facet_wrap(~OTU,scales="free_y",ncol=8)+
	theme_bw()

all.summ=summarySE(stage,measurevar="Abundance",groupvars=c("OTU","Type"))
quartz()
pd <- position_dodge(.3)
ggplot(all.summ,aes(x=Type,y=Abundance))+
	geom_point(aes(),position=pd,size=2.5)+
	geom_errorbar(aes(ymin=Abundance-se,ymax=Abundance+se),lwd=0.4,width=0.3,position=pd)+
	geom_line(aes(),position=pd)+
	#scale_fill_manual(values=c("lightblue","deepskyblue"))+
	facet_wrap(~OTU,scales="free_y",ncol=6)+
	theme_bw()



###########################
#Relationships among sig ASVs?
##########################
library(pegas)
library(phangorn)
library(seqinr)
input <- "Cvars_ClustalAlignment_WithGeoSym_WithGenBankRefs.fasta"
d <- ape::read.dna(input, format='fasta')
e <- dist.dna(d)
h <- pegas::haplotype(d,labels=c("D1","sq33","D1a")) #,labels=c("D1","sq33","D1a")
h <- sort(h, what = "label")
(net <- pegas::haploNet(h))
ind.hap<-with(
stack(setNames(attr(h, "index"), rownames(h))),
table(hap=ind, pop=rownames(d)[values])
)
plot(net, size=c(0.5,0.5,1), scale.ratio=2,cex=0.8, pie=ind.hap,fast=TRUE)
legend(-8, 0, colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=19, ncol=2)

input <- "Cvars_C15only_ClustalAlignment_WithGeoSym_WithGenBankRefs.fasta"
d <- ape::read.dna(input, format='fasta')
e <- dist.dna(d)
h <- pegas::haplotype(d,labels=c("sq20","sq3","sq21","C15","sq18","sq11","sq10","sq8","sq5","sq4","sq2","sq12")) #,labels=c("sq20","sq3","sq21","C15","sq18","sq11","sq10","sq8","sq5","sq4","sq2","sq12")
#h <- sort(h, what = "label")
(net <- pegas::haploNet(h))
ind.hap<-with(
stack(setNames(attr(h, "index"), rownames(h))),
table(hap=ind, pop=rownames(d)[values])
)
plot(net,  scale.ratio=1, pie=ind.hap)
legend(-8, 0, colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=19, ncol=2)

#plot phylogenetic tree
input <- "Cvars_ClustalAlignment_WithGeoSym_WithGenBankRefs.fasta"
d <- ape::read.dna(input, format='fasta')
d_phyDat<-phyDat(d,type="DNA",levels=NULL) #convert sequence data to phyDat object

dna_dist<-dist.ml(d_phyDat) #build tree with distance based methods
seqs_UPGMA<-upgma(dna_dist)
seqs_NJ<-NJ(dna_dist)
plot(seqs_UPGMA)
plot(seqs_NJ)

mt<-modelTest(d_phyDat)
print(mt)

fit<-pml(seqs_UPGMA,d_phyDat)
fitHKY<-optim.pml(fit,model="HKY",optInv=FALSE,optGamma=FALSE,rearrangement="NNI",control=pml.control(trace=0)) 


bs <- bootstrap.pml(fitHKY, bs=100,optNni=TRUE,control=pml.control(trace=0))

plotBS(fitHKY$tree, bs, p=50, type="p")


###Plotting bar graphs for family 7 as example of new ASVs by year and intra colony spatial var
#first sum each row in ASV columns and obtain mean sum, then multiply each value in matrix by mean sum

head(PPV_Norm)
names(PPV_Norm)

goods<-data.frame(t(PPV_Norm[,c(2:108)]))

fam7_15<-goods[c(24:35),]
fam7_15$type<-c("adult",rep("larvae",11))
fam7_15$sample<-c("A","L1","L10","L11","L2","L3","L4","L5","L6","L7","L8","L9") #fix names with strsplit to get rid of year
fam7_16<-goods[c(79,78,80:92),] #re-order so adult samples are Left, Center, Right
fam7_16$type<-c(rep("adult",3),rep("larvae",12))
fam7_16$sample<-c("A1L","A2C","A3R","L1","L10","L11","L12","L2","L3","L4","L5","L6","L7","L8","L9")

######Now plot normalized read proportions by year for family 7
gss_7_15=otuStack(fam7_15,count.columns=c(1:34),condition.columns=c(35:36))[1:408,] #must remove 'summ' 
head(gss_7_15)

gss_7_16=otuStack(fam7_16,count.columns=c(1:34),condition.columns=c(35:36))[1:510,] #must remove 'summ' 

ggplot(gss_7_16, aes(x=sample, y=count, fill = otu))+geom_bar(position = "stack",stat="identity")+theme_bw()

