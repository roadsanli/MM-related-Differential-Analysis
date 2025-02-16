library(tidyr)
library(stringr)
library(dplyr)
library(qqman)
library(qqman)
library(lattice)
source("/share/home/xxwu/Xiaohui/reference_1000g/prepareSS.R")
source("/share/home/xxwu/UKB_data/qqplot.r")
bim <- read.table("/share/home/xxwu/Xiaohui/reference_1000g/g1000_eur_msu.bim")
#### hail outcome prepare(other gwas not) ֻ
trait <- "female_Multiple_myeloma"   
ss <- NULL
i <- "X"
#gzfile <- paste0("gunzip /share/home/xxwu/Xiaohui/GWAS_Hail/#res/",trait,"_chr",i,".tsv.gz")
#system(gzfile)
s <- read.table(paste0("/share/home/xxwu/Xiaohui/GWAS_Hail/res/",trait,"_chr",i,".tsv"),header=T)
ss <- rbind(ss,s)
RM <- paste0(paste0("rm /share/home/xxwu/Xiaohui/GWAS_Hail/res/",trait,"_chr",i,".tsv") )
#system(RM)

b <- separate(data = ss, col = locus, into = c("chr", "BP"), sep = ":")
c <- separate(data = b, col = alleles, into = c("1","A1","2","A2","4"), sep = "\"")
##A1  REF
##A2  ALT
c <- c[,-c(3,5,7)]
gwas <- c

######### 01  merge gwas snp with 1000g by chr bp
##gwas file :gwas file include snp      chr /bp /A1(minor) /A2 /beta /se /pvalue  
##gwas file format
#bp chr 
colnames(gwas)[2] <- "bp"
colnames(gwas)[9] <- "se"
colnames(gwas)[11] <- "pvalue"

Info <- 0 
MAF  <- 0 
N <- 0  

###
#gwas$chr <- as.numeric(gwas$chr)
gwas$bp <- as.numeric(gwas$bp)

##output:chr bp rs and so on

########## 02 gwas cleaning

##remove error length allele snp
ss <- gwas
#ss <- ss[nchar(ss$A1)==1,]
#ss <- ss[nchar(ss$A2)==1,]

reference <- readRDS("reference.rds")
ss$rs <- reference[match(ss$bp, reference$BP),5]
##remove duplicate snp
snp <- unique(ss$rs)
ss <- ss[match(snp,ss$rs),]


##############summary
C <- data.frame(SNP  = ss$rs,
                chr  = ss$chr,
                BP   = ss$bp,
                A1   = ss$A1,
                A2   = ss$A2,
                beta = ss$beta,
                se   = ss$se, #se
                pvalue=ss$pvalue)##p
C <- C[!is.na(C$SNP),]
C <- C[!is.na(C$chr),]
C <- C[!is.na(C$BP),]
C <- C[!is.na(C$A1),]
C <- C[!is.na(C$A2),]
C <- C[!is.na(C$beta),]
C <- C[!is.na(C$se),]
C <- C[!is.na(C$pvalue),]
write.table(C,paste0("/share/home/xxwu/Xiaohui/GWAS_clean/",trait,"_clean_gwas.txt"),row.names=FALSE,quote=FALSE)

