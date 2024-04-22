library(dplyr)
#####XVpar mapped to W303#####
#read in vcf file - remove one comment line first so get header
xvparw <- read.table("filtered_par/XVparW_diSNP15_filtered.vcf", header=FALSE, comment.char="#")

#give names to cols
names(xvparw) <- c("chr", "pos", "id","ref","alt","qual","filter","info","format","xvpar")


####xvpar - region E####
#range for E - YOL160W to YOL157C is region, did transformations pau20 (YOL161C) to zps1 (YOL154W)\ - chrXV
#for W303 - 13641	14003 to 36750	37499\
#Go a bit outside (1kb on either end).
xvparwE <- subset(xvparw, pos>12641&pos<38499)

#only get columns want
#info want:
#chr, pos, ref, alt, qual, cerw
#last one GT:AD:DP:GQ:PL
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
xvparwE <- select(xvparwE, c(chr, pos, ref, alt, qual, xvpar))

#split last column
xvparwE$xvpar <- as.character(xvparwE$xvpar)
##second col - FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (high-quality bases)"> - ref, alt
##second last col - ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##third col - ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
xvparwE$nref <- 0
xvparwE$nalt <- 0
xvparwE$gq <- 0
xvparwE$dp <- 0
for (i in 1:length(xvparwE$chr)){
  AD <- strsplit(strsplit(xvparwE$xvpar, ":")[[i]][2],",")
  xvparwE$nref[i] <- as.numeric(AD[[1]][1])
  xvparwE$nalt[i] <- as.numeric(AD[[1]][2])
  xvparwE$gq[i] <- as.numeric(strsplit(strsplit(xvparwE$xvpar, ":")[[i]][4],",")[[1]])
  xvparwE$dp[i] <- as.numeric(strsplit(strsplit(xvparwE$xvpar, ":")[[i]][3],",")[[1]])
}

#calculate frequencies
xvparwE$falt <- xvparwE$nalt/(xvparwE$nref+xvparwE$nalt)
#calculate qual/dp
xvparwE$qd <- xvparwE$qual/xvparwE$dp

############decide which snps to use##########
length(xvparwE$chr)
#858

#use only with falt ==1
Esnps <- xvparwE[xvparwE$falt==1,]
length(Esnps$chr)
#774

#use only with QUAL > 50
Esnps <- Esnps[Esnps$qual>50,]
length(Esnps$chr)
#769

#use only with GQ > 20
Esnps <- Esnps[Esnps$gq>20,]
length(Esnps$chr)
#706

#use only with QUAL/DP > 20
Esnps <- Esnps[Esnps$qd>10,]
length(Esnps$chr)
#706

#use only with DP > 20
Esnps <- Esnps[Esnps$dp>20,]
length(Esnps$chr)
#608

write.csv(Esnps, "Ereg_snps.csv")

####xvpar - region F####
#range for F - YOL138C to YOL024W is region, did transformations YOL020W to YOL141W - chrXV
#for W303 - 288189	289967	to	58476	60563
#Go a bit outside (1kb on either end).
xvparwF <- subset(xvparw, pos>57476&pos<290967)

#only get columns want
#info want:
#chr, pos, ref, alt, qual, cerw
#last one GT:AD:DP:GQ:PL
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
xvparwF <- select(xvparwF, c(chr, pos, ref, alt, qual, xvpar))

#split last column
xvparwF$xvpar <- as.character(xvparwF$xvpar)
##second col - FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (high-quality bases)"> - ref, alt
##second last col - ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##third col - ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
xvparwF$nref <- 0
xvparwF$nalt <- 0
xvparwF$gq <- 0
xvparwF$dp <- 0
for (i in 1:length(xvparwF$chr)){
  AD <- strsplit(strsplit(xvparwF$xvpar, ":")[[i]][2],",")
  xvparwF$nref[i] <- as.numeric(AD[[1]][1])
  xvparwF$nalt[i] <- as.numeric(AD[[1]][2])
  xvparwF$gq[i] <- as.numeric(strsplit(strsplit(xvparwF$xvpar, ":")[[i]][4],",")[[1]])
  xvparwF$dp[i] <- as.numeric(strsplit(strsplit(xvparwF$xvpar, ":")[[i]][3],",")[[1]])
}


#calculate frequencies
xvparwF$falt <- xvparwF$nalt/(xvparwF$nref+xvparwF$nalt)
#calculate qual/dp
xvparwF$qd <- xvparwF$qual/xvparwF$dp

############decide which snps to use##########
length(xvparwF$chr)
#16839

#use only with falt ==1
Fsnps <- xvparwF[xvparwF$falt==1,]
length(Fsnps$chr)
#16833

#use only with QUAL > 50
Fsnps <- Fsnps[Fsnps$qual>50,]
length(Fsnps$chr)
#16757

#use only with GQ > 20
Fsnps <- Fsnps[Fsnps$gq>20,]
length(Fsnps$chr)
#12604

#use only with QUAL/DP > 20
Fsnps <- Fsnps[Fsnps$qd>10,]
length(Fsnps$chr)
#12604

#use only with DP > 20
Fsnps <- Fsnps[Fsnps$dp>20,]
length(Fsnps$chr)
#4971

write.csv(Fsnps, "Freg_snps.csv")






