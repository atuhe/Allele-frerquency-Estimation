#Determining the minor/major allele and their frequencies

# Load required packages
#install.packages("genetics")
library(genetics)

fms <- read.table("FMS_data.txt")
attach(fms) #We use the attach() function so that we can call each variable by its name without having to indicate the corresponding dataframe.
#A dataframe must be re- attached at the start of a new R session for the corresponding variable names to be recognized
View(fms)
head(fms)
fms[1:5,1:10]
names(fms)

##Identifying the minor allele and its frequency
#Suppose we are interested in determining the minor allele for the SNP 
#labeled actn3 rs540874 in the FAMuSS data.

#To do this, we need to calculate corresponding allele frequencies

# First we determine the number of observations with each genotype
#attach(fms)
GenoCount <- summary(actn3_rs540874)
GenoCount
#Or
table(fms$actn3_rs540874)
table(fms$actn3_rs540874, exclude = NULL)
#we see n = 226 individuals have the AA genotype, 
#n = 595 individuals have the GA genotype and n = 395 
#individuals have the GG genotype.
#An additional n = 181 individuals are missing the genotype
#For simplicity, we assume that this missingness is 
#non- informative. 
#That is, we make the strong assumption that our estimates 
#of the allele frequencies would be the same had we observed 
#the genotypes for these individuals.

#To calculate the allele frequencies, we begin by determining
#our reduced sample size (that is, the number of individuals 
#with complete data):
sum(is.na(actn3_rs540874))
nrow(fms)
nrow(fms)-sum(is.na(actn3_rs540874))
sum(complete.cases(fms$actn3_rs540874))

NumbObs <- sum(!is.na(actn3_rs540874))
NumbObs

#The genotype frequencies for AA, GA and GG are then given 
#respectively by
GenoFreq <- as.vector(GenoCount/NumbObs)
GenoFreq

GenoFreq <- as.vector(table(fms$actn3_rs540874))/NumbObs
GenoFreq

#The frequencies of the A and G alleles are calculated as 
#follows:
#Recall: For diploid organisms, if f(AA), f(AB) and f(BB)
#are the frequencies  of the 3 genotypes at a locus with
#2 alleles, then the frequency, p of the A-allele and the
#frequency q of the B-allele in the population are given by:
#p=f(AA)+1/2f(AB)=frequency of A
#q=f(BB)+1/2f(AB)=frequency of B
#Note: AA organisms do have 2 A alleles; 
#AB orbanisms do have 1 A allele and 
#BB do have 0 Z alleles. 
FreqA <- GenoFreq[1] + GenoFreq[2]/2
FreqA
(2*GenoFreq[1] + GenoFreq[2])/2 #Divide by total number of alles = 2

FreqG <- GenoFreq[3] + GenoFreq[2]/2
FreqG

#check: p + q = FreqA + FreqG = 1
FreqA + FreqG

#Thus, we report A is the minor allele at this SNP locus, 
#with a frequency of 0.43.

#Alternatively, we can achieve the same result using 
#the genotype() and summary() functions within the genetics 
#package

#install.packages("genetics")
#library(genetics)

#We then create a genotype object and summarize the 
#corresponding genotype and allele frequencies:
Geno <- genotype(actn3_rs540874,sep="")
summary(Geno)

obj <- summary(Geno)
obj$allele.freq
class(obj$allele.freq) 

#Here we again see that A corresponds to the minor allele 
#at this SNP locus, with a frequency of 0.43, while G is 
#the major allele, with a greater frequency of 0.57


# Do it for all SNPs in the actn3 gene
library(stringr)
regx <- str_detect(names(fms), "actn3")
regx

snpsdf <- fms[,regx]
names(snpsdf)

#Read in the SNPs in the act3 gene SNPs as genotype variables and create a corresponding dataframe
Actn3Snp1 <- genotype(actn3_r577x,sep="")
Actn3Snp2 <- genotype(actn3_rs540874,sep="")
Actn3Snp3 <- genotype(actn3_rs1815739,sep="")
Actn3Snp4 <- genotype(actn3_1671064,sep="")
Actn3AllSnps <- data.frame(Actn3Snp1,Actn3Snp2,Actn3Snp3,Actn3Snp4)

snps.list <- as.list(snpsdf)

GenoCount.list <- lapply(snps.list, genotype, sep="")
GenoCount.list
GenoCount.list[[2]]


obj <- lapply(GenoCount.list, summary)
obj

lapply(obj, function(x) x$allele.freq)
class(obj[[1]]) 






