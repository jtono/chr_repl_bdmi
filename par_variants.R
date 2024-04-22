library(dplyr)
#####cer#####
#read in vcf file - remove one comment line first so get header
cer9 <- read.table("filtered_par/cerW_diSNP9_filtered.vcf", header=FALSE, comment.char="#")
cer10 <- read.table("filtered_par/cerW_diSNP10_filtered.vcf", header=FALSE, comment.char="#")
cer15 <- read.table("filtered_par/cerW_diSNP15_filtered.vcf", header=FALSE, comment.char="#")

#check that it worked
head(cer15)

#give names to cols
names(cer9) <- c("chr", "pos", "id","ref","alt","qual","filter","info","format","cerw")
names(cer10) <- c("chr", "pos", "id","ref","alt","qual","filter","info","format","cerw")
names(cer15) <- c("chr", "pos", "id","ref","alt","qual","filter","info","format","cerw")

####cer - region A####
#range for A - YIL166C to YIL156W is region, did transformations YIL151C to YIL169C - chrIX
#- go a bit outside (1kb)
#for W303 - 24534	to	62140
cer9A <- subset(cer9, pos>23534&pos<63140)

#only get columns want
#info want:
#chr, pos, ref, alt, qual, cerw
#last one GT:AD:DP:GQ:PL
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
cer9A <- select(cer9A, c(chr, pos, ref, alt, qual, cerw))

#split last column
cer9A$cerw <- as.character(cer9A$cerw)

##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (high-quality bases)"> - ref, alt
cer9A$cerwref <- 0
cer9A$cerwalt <- 0
for (i in 1:length(cer9A$chr)){
  AD <- strsplit(strsplit(cer9A$cerw, ":")[[i]][2],",")
  cer9A$cerwref[i] <- AD[[1]][1]
  cer9A$cerwalt[i] <- AD[[1]][2]
}

#calculate frequencies
cer9A$altf <- as.numeric(cer9A$cerwalt)/(as.numeric(cer9A$cerwref)+as.numeric(cer9A$cerwalt))
#none over 90%
#only 4 actually varying, so could just do all

####cer - region C####
#range for C - YJL218W to YJL165C is region, did transformations YJL164C to YJL219W - chrX
#- go a bit outside (1kb)
#for W303 - 23187	24890	to	113654	114847
cer10C <- subset(cer10, pos>22187&pos<115847)
#none

####cer - region E####
#range for E rep2 - YOL160W to YOL157C is region, did transformations pau20 (YOL161C) to zps1 (YOL154W)\ - chrXV
#for W303 - 13641	14003 to 36750	37499\
#Go a bit outside (1kb on either end).
cer15E <- subset(cer15, pos>12641&pos<38499)

#only get columns want
#info want:
#chr, pos, ref, alt, qual, cerw
#last one GT:AD:DP:GQ:PL
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
cer15E <- select(cer15E, c(chr, pos, ref, alt, qual, cerw))

#split last column
cer15E$cerw <- as.character(cer15E$cerw)

##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (high-quality bases)"> - ref, alt
cer15E$cerwref <- 0
cer15E$cerwalt <- 0
for (i in 1:length(cer15E$chr)){
  AD <- strsplit(strsplit(cer15E$cerw, ":")[[i]][2],",")
  cer15E$cerwref[i] <- AD[[1]][1]
  cer15E$cerwalt[i] <- AD[[1]][2]
}

#calculate frequencies
cer15E$altf <- as.numeric(cer15E$cerwalt)/(as.numeric(cer15E$cerwref)+as.numeric(cer15E$cerwalt))
#none over 90%
#only 5 actually varying, so could just do all

####cer - region F####
#range for F - YOL138C to YOL024W is region, did transformations YOL020W to YOL141W - chrXV
#for W303 - 288189	289967	to	58476	60563
#Go a bit outside (1kb on either end).
cer15F <- subset(cer15, pos>57476&pos<290967)

#only get columns want
#info want:
#chr, pos, ref, alt, qual, cerw
#last one GT:AD:DP:GQ:PL
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
cer15F <- select(cer15F, c(chr, pos, ref, alt, qual, cerw))

#split last column
cer15F$cerw <- as.character(cer15F$cerw)

##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (high-quality bases)"> - ref, alt
cer15F$cerwref <- 0
cer15F$cerwalt <- 0
for (i in 1:length(cer15F$chr)){
  AD <- strsplit(strsplit(cer15F$cerw, ":")[[i]][2],",")
  cer15F$cerwref[i] <- AD[[1]][1]
  cer15F$cerwalt[i] <- AD[[1]][2]
}

#calculate frequencies
cer15F$altf <- as.numeric(cer15F$cerwalt)/(as.numeric(cer15F$cerwref)+as.numeric(cer15F$cerwalt))
#all fixed but only 5 of them


####par - region A####
#read in vcf file - remove one comment line first so get header
par9 <- read.table("filtered/IXparN_diSNP9_filtered.vcf", header=FALSE, comment.char="#")

#check that it worked
head(par9)

#give names to cols
names(par9) <- c("chr", "pos", "id","ref","alt","qual","filter","info","format","IXparn")

#range for A - YIL166C to YIL156W is region, did transformations YIL151C to YIL169C - chrIX
#- go a bit outside (1kb)
#For N17 - (??) end to 41322\
par9A <- subset(par9, pos<42322)

#only get columns want
#info want:
#chr, pos, ref, alt, qual, cerw
#last one GT:AD:DP:GQ:PL
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
par9A <- select(par9A, c(chr, pos, ref, alt, qual, IXparn))

#split last column
par9A$IXparn <- as.character(par9A$IXparn)

##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (high-quality bases)"> - ref, alt
par9A$parnref <- 0
par9A$parnalt <- 0
for (i in 1:length(par9A$chr)){
  AD <- strsplit(strsplit(par9A$IXparn, ":")[[i]][2],",")
  par9A$parnref[i] <- AD[[1]][1]
  par9A$parnalt[i] <- AD[[1]][2]
}

#calculate frequencies
par9A$altf <- as.numeric(par9A$parnalt)/(as.numeric(par9A$parnref)+as.numeric(par9A$parnalt))
#3 over 90%, only those 3 over 50%
#only 28 in total


####par - region C####
#read in vcf file - remove one comment line first so get header
par10 <- read.table("filtered/XparN_diSNP10_filtered.vcf", header=FALSE, comment.char="#")

#check that it worked
head(par10)

#give names to cols
names(par10) <- c("chr", "pos", "id","ref","alt","qual","filter","info","format","Xparn")

#range for C - YJL218W to YJL165C is region, did transformations YJL164C to YJL219W - chrX
#- go a bit outside (1kb)
#for N17 - 1422 to 102835\
par10C <- subset(par10, pos>422&pos<103835)

#only get columns want
#info want:
#chr, pos, ref, alt, qual, cerw
#last one GT:AD:DP:GQ:PL
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
par10C <- select(par10C, c(chr, pos, ref, alt, qual, Xparn))

######left off here#####

#split last column
par9A$IXparn <- as.character(par9A$IXparn)

##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (high-quality bases)"> - ref, alt
par9A$parnref <- 0
par9A$parnalt <- 0
for (i in 1:length(par9A$chr)){
  AD <- strsplit(strsplit(par9A$IXparn, ":")[[i]][2],",")
  par9A$parnref[i] <- AD[[1]][1]
  par9A$parnalt[i] <- AD[[1]][2]
}

#calculate frequencies
par9A$altf <- as.numeric(par9A$parnalt)/(as.numeric(par9A$parnref)+as.numeric(par9A$parnalt))
#3 over 90%, only those 3 over 50%
#only 28 in total


