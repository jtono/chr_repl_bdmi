library(reshape)
library(dplyr)
is.sorted = Negate(is.unsorted)

########functions#########
#get the coverage of all genes (for first three bases) on a certain chromosome from a certain coverage file
get_cov_gene <- function(genedatfr, colname.c, colname.p, coverage, chr){
  cov.c <- subset(coverage, Chr==paste("W303",chr, sep="."))
  cov.p <- subset(coverage, Chr==paste("N_17",chr, sep="."))
  genedatfr <- subset(genedatfr, seqid==chr)

  genedatfr[[colname.c]] <- 0
  genedatfr[[colname.p]] <- 0
  for (i in 1:length(genedatfr$ID)){
    if(genedatfr$strand.c[i]=="-"){
      genedatfr[[colname.c]][i] <- mean(cov.c$depth[(genedatfr$end.c[i]-2):(genedatfr$end.c[i])])
    }else{
      genedatfr[[colname.c]][i] <- mean(cov.c$depth[(genedatfr$start.c[i]):(genedatfr$start.c[i]+2)])
    }
    if(genedatfr$strand.p[i]=="-"){
      genedatfr[[colname.p]][i] <- mean(cov.p$depth[(genedatfr$end.p[i]-2):(genedatfr$end.p[i])])
    }else{
      genedatfr[[colname.p]][i] <- mean(cov.p$depth[(genedatfr$start.p[i]):(genedatfr$start.p[i]+2)])
    }
  }

  return(genedatfr)
}
######load in files######
#load in cer coverage
cer <- read.table("coverage/cer.coverage",header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
# renames the header
cer <- rename(cer, c("Chr"=V1, "locus"=V2, "depth"=V3))

#load in IXpar coverage
IXpar <- read.table("coverage/IXpar.coverage",header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
# renames the header
IXpar <- rename(IXpar, c("Chr"=V1, "locus"=V2, "depth"=V3))

#load in Xpar coverage
Xpar <- read.table("coverage/Xpar.coverage",header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
# renames the header
Xpar <- rename(Xpar, c("Chr"=V1, "locus"=V2, "depth"=V3))

#load in XVpar coverage
XVpar <- read.table("coverage/XVpar.coverage",header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
# renames the header
XVpar <- rename(XVpar, c("Chr"=V1, "locus"=V2, "depth"=V3))

#load in gene files
genes.c <- read.csv("W303.genes.csv", head=TRUE)
genes.p <- read.csv("N17.genes.csv", head=TRUE)

#get rid of weird W303.scplasm1
cer <- cer[-which(cer$Chr=="W303.scplasm1"),]
IXpar <- IXpar[-which(IXpar$Chr=="W303.scplasm1"),]
Xpar <- Xpar[-which(Xpar$Chr=="W303.scplasm1"),]
XVpar <- XVpar[-which(XVpar$Chr=="W303.scplasm1"),]

#####A samples#######
A1 <- read.table("coverage/A1.coverage",header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
# renames the header
A1 <- rename(A1, c("Chr"=V1, "locus"=V2, "depth"=V3))

A2 <- read.table("coverage/A2.coverage",header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
# renames the header
A2 <- rename(A2, c("Chr"=V1, "locus"=V2, "depth"=V3))

#########find genes for region A - chr9########
#limit to only genes on chr interested in
genes9.c <- subset(genes.c, seqid=="chr09")
genes9.p <- subset(genes.p, seqid=="chr09")


#find which genes in both and in same order
c.geneind9 <-c()
for (i in 1:length(genes9.c$ID)){
  c.geneind9 <- c(c.geneind9, which(genes9.c$ID%in%genes9.p$ID[i]))
}

#in order
is.sorted(c.geneind9)

#make n17 (par) equivalent
p.geneind9 <- c()
for (i in 1:length(c.geneind9)){
  id <- as.character(genes9.c$ID[c.geneind9[i]])
  p.geneind9 <- c(p.geneind9, which(genes9.p$ID==id))
}
#in order
is.sorted(p.geneind9)

#pull out only genes that are in both
genes.c.both9 <- genes9.c[c.geneind9,]
genes.p.both9 <- genes9.p[p.geneind9,]

#make data frame with genes and between genes regions
ID <- 1
start.c <- 1
end.c <- c()
start.p <- 1
end.p <- c()
strand.c <- "X"
strand.p <- "X"
for (i in 1:length(genes.c.both9$X)){
  ID <- c(ID, as.character(genes.c.both9$ID[i]), i+1)
  start.c <- c(start.c, genes.c.both9$start[i], (genes.c.both9$end[i]+1))
  end.c <- c(end.c, (genes.c.both9$start[i]-1), genes.c.both9$end[i])
  start.p <- c(start.p, genes.p.both9$start[i], (genes.p.both9$end[i]+1))
  end.p <- c(end.p, (genes.p.both9$start[i]-1), genes.p.both9$end[i])
  strand.c <- c(strand.c, genes.c.both9$strand[i], "X")
  strand.p <- c(strand.p, genes.p.both9$strand[i], "X")
}
cer9.c <- subset(cer, Chr=="W303.chr09")
cer9.p <- subset(cer, Chr=="N_17.chr09")
end.c <- c(end.c, max(cer9.c$locus))
end.p <- c(end.p, max(cer9.p$locus))
cov9.cp <- data.frame(start.c, end.c, strand.c, start.p, end.p, strand.p)

#########figure out coverage of each gene for each sample - region A######
#pull out only chr of interest
#'for A - chr9
#'do by species
cer9.c <- subset(cer, Chr=="W303.chr09")
cer9.p <- subset(cer, Chr=="N_17.chr09")
IXpar9.c <- subset(IXpar, Chr=="W303.chr09")
IXpar9.p <- subset(IXpar, Chr=="N_17.chr09")
Xpar9.c <- subset(Xpar, Chr=="W303.chr09")
Xpar9.p <- subset(Xpar, Chr=="N_17.chr09")
XVpar9.c <- subset(XVpar, Chr=="W303.chr09")
XVpar9.p <- subset(XVpar, Chr=="N_17.chr09")


#get depth of A, T, G (or first 3 bases)
cov9.cp$depth.cer.c <- 0
cov9.cp$depth.cer.p <- 0
cov9.cp$depth.IXpar.c <- 0
cov9.cp$depth.IXpar.p <- 0
cov9.cp$depth.Xpar.c <- 0
cov9.cp$depth.Xpar.p <- 0
cov9.cp$depth.XVpar.c <- 0
cov9.cp$depth.XVpar.p <- 0
for (i in 1:length(cov9.cp$ID)){
  if(cov9.cp$strand.c[i]=="-"){
    cov9.cp$depth.cer.c[i] <- mean(cer9.c$depth[(cov9.cp$end.c[i]-2):(cov9.cp$end.c[i])])
    cov9.cp$depth.IXpar.c[i] <- mean(IXpar9.c$depth[(cov9.cp$end.c[i]-2):(cov9.cp$end.c[i])])
    cov9.cp$depth.Xpar.c[i] <- mean(Xpar9.c$depth[(cov9.cp$end.c[i]-2):(cov9.cp$end.c[i])])
    cov9.cp$depth.XVpar.c[i] <- mean(XVpar9.c$depth[(cov9.cp$end.c[i]-2):(cov9.cp$end.c[i])])
  }else{
    cov9.cp$depth.cer.c[i] <- mean(cer9.c$depth[(cov9.cp$start.c[i]):(cov9.cp$start.c[i]+2)])
    cov9.cp$depth.IXpar.c[i] <- mean(IXpar9.c$depth[(cov9.cp$start.c[i]):(cov9.cp$start.c[i]+2)])
    cov9.cp$depth.Xpar.c[i] <- mean(Xpar9.c$depth[(cov9.cp$start.c[i]):(cov9.cp$start.c[i]+2)])
    cov9.cp$depth.XVpar.c[i] <- mean(XVpar9.c$depth[(cov9.cp$start.c[i]):(cov9.cp$start.c[i]+2)])
  }
  if(cov9.cp$strand.p[i]=="-"){
    cov9.cp$depth.cer.p[i] <- mean(cer9.p$depth[(cov9.cp$end.p[i]-2):(cov9.cp$end.p[i])])
    cov9.cp$depth.IXpar.p[i] <- mean(IXpar9.p$depth[(cov9.cp$end.p[i]-2):(cov9.cp$end.p[i])])
    cov9.cp$depth.Xpar.p[i] <- mean(Xpar9.p$depth[(cov9.cp$end.p[i]-2):(cov9.cp$end.p[i])])
    cov9.cp$depth.XVpar.p[i] <- mean(XVpar9.p$depth[(cov9.cp$end.p[i]-2):(cov9.cp$end.p[i])])
  }else{
    cov9.cp$depth.cer.p[i] <- mean(cer9.p$depth[(cov9.cp$start.p[i]):(cov9.cp$start.p[i]+2)])
    cov9.cp$depth.IXpar.p[i] <- mean(IXpar9.p$depth[(cov9.cp$start.p[i]):(cov9.cp$start.p[i]+2)])
    cov9.cp$depth.Xpar.p[i] <- mean(Xpar9.p$depth[(cov9.cp$start.p[i]):(cov9.cp$start.p[i]+2)])
    cov9.cp$depth.XVpar.p[i] <- mean(XVpar9.p$depth[(cov9.cp$start.p[i]):(cov9.cp$start.p[i]+2)])
  }
}

######test here#########

  cov9.cp$seqid <- "chr09"
test <- get_cov_gene(cov9.cp, "depth.cer.c2", "depth.cer.p2", cer, "chr09")

#########test here#########

#take a look at it
plot(cov9.cp$depth.cer.c, type="l")
lines(cov9.cp$depth.cer.p, col="pink", lty=2)
lines(cov9.cp$depth.IXpar.c, col="blue")
lines(cov9.cp$depth.IXpar.p, col="purple", lty=2)
lines(cov9.cp$depth.Xpar.c, col="green")
lines(cov9.cp$depth.Xpar.p, col="darkgreen", lty=2)
lines(cov9.cp$depth.XVpar.c, col="red")
lines(cov9.cp$depth.XVpar.p, col="darkred", lty=2)

#pull out only genes
row_odd <- seq_len(nrow(cov9.cp)) %% 2
cov9.cp.genes <- cov9.cp[row_odd==0,]

#plot just gene info
plot(cov9.cp.genes$depth.cer.c, type="l", ylim=c(0,50))
lines(cov9.cp.genes$depth.cer.p, col="pink", lty=2)
lines(cov9.cp.genes$depth.IXpar.c, col="blue")
lines(cov9.cp.genes$depth.IXpar.p, col="purple", lty=2)
lines(cov9.cp.genes$depth.Xpar.c, col="green")
lines(cov9.cp.genes$depth.Xpar.p, col="darkgreen", lty=2)
lines(cov9.cp.genes$depth.XVpar.c, col="red")
lines(cov9.cp.genes$depth.XVpar.p, col="darkred", lty=2)
#gets rid of some peaks but not much diff otherwise

#####only find region of interest - region A####
#region A - start
grep("YIL166C", cov9.cp$ID)
#2
abline(v=2, lty=2, col="forestgreen")
#end
grep("YIL156W", cov9.cp$ID)
cov9.cp[22:26,]
#26 real one
abline(v=26, lty=2, col="forestgreen")
#region with rec selected - start
grep("YIL151C", cov9.cp$ID)
#36
abline(v=36, lty=4, col="purple")
#end
grep("YIL169C", cov9.cp$ID)
#not there - just the end somewhere?

#make cut down data frame for just the region of interest (where rec selected)
cov9.cp.A <- cov9.cp[1:36,]
cov9.cp.A.genes <- cov9.cp.genes[1:18,]


##############Step 1: check for mismapping - A###########
range(cov9.cp.A.genes$depth.cer.c)
range(cov9.cp.A.genes$depth.cer.p)
range(cov9.cp.A.genes$depth.IXpar.c)
range(cov9.cp.A.genes$depth.IXpar.p)
range(cov9.cp.A.genes$depth.Xpar.c)
range(cov9.cp.A.genes$depth.Xpar.p)
range(cov9.cp.A.genes$depth.XVpar.c)
range(cov9.cp.A.genes$depth.XVpar.p)
#no mismapping

#plots
#all blocks
plot(cov9.cp.A$depth.cer.c, type="l", ylim=c(0,100))
lines(cov9.cp.A$depth.cer.p, col="red", lty=2)
lines(cov9.cp.A$depth.IXpar.c, col="black")
lines(cov9.cp.A$depth.IXpar.p, col="red", lty=2)
lines(cov9.cp.A$depth.Xpar.c, col="black")
lines(cov9.cp.A$depth.Xpar.p, col="red", lty=2)
lines(cov9.cp.A$depth.XVpar.c, col="black")
lines(cov9.cp.A$depth.XVpar.p, col="red", lty=2)

#only genes
plot(cov9.cp.A.genes$depth.cer.c, type="l", ylim=c(0,100))
lines(cov9.cp.A.genes$depth.cer.p, col="red", lty=2)
lines(cov9.cp.A.genes$depth.IXpar.c, col="black")
lines(cov9.cp.A.genes$depth.IXpar.p, col="red", lty=2)
lines(cov9.cp.A.genes$depth.Xpar.c, col="black")
lines(cov9.cp.A.genes$depth.Xpar.p, col="red", lty=2)
lines(cov9.cp.A.genes$depth.XVpar.c, col="black")
lines(cov9.cp.A.genes$depth.XVpar.p, col="red", lty=2)
#pattern not exactly the same with different cerevisiae


#########find DNA loading correction - chr1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 14, 16########
#limit to only genes on chr interested in
genes_rest.c <- subset(genes.c, seqid%in%c("chr01","chr02","chr03","chr04","chr05","chr06","chr07","chr08","chr11","chr13","chr14","chr16"))
genes_rest.p <- subset(genes.p, seqid%in%c("chr01","chr02","chr03","chr04","chr05","chr06","chr07","chr08","chr11","chr13","chr14","chr16"))


#find which genes in both and in same order
c.geneind_rest <-c()
for (i in 1:length(genes_rest.c$ID)){
  c.geneind_rest <- c(c.geneind_rest, which(genes_rest.c$ID%in%genes_rest.p$ID[i]))
}

#in order
is.sorted(c.geneind_rest)

#make n17 (par) equivalent
p.geneind_rest <- c()
for (i in 1:length(c.geneind_rest)){
  id <- as.character(genes_rest.c$ID[c.geneind_rest[i]])
  p.geneind_rest <- c(p.geneind_rest, which(genes_rest.p$ID==id))
}
#in order
is.sorted(p.geneind_rest)

#pull out only genes that are in both
genes.c.both_rest <- genes_rest.c[c.geneind_rest,]
genes.p.both_rest <- genes_rest.p[p.geneind_rest,]

#make data frame with genes only
cov_rest.cp <- merge(select(genes.c.both_rest, seqid, start, end, strand, ID), select(genes.p.both_rest, seqid, start, end, strand, ID), by=c("seqid","ID"))
names(cov_rest.cp) <- c("seqid","ID","start.c","end.c","strand.c","start.p","end.p","strand.p")

#######here#########

cer9.c <- subset(cer, Chr=="W303.chr09")
cer9.p <- subset(cer, Chr=="N_17.chr09")
IXpar9.c <- subset(IXpar, Chr=="W303.chr09")
IXpar9.p <- subset(IXpar, Chr=="N_17.chr09")
Xpar9.c <- subset(Xpar, Chr=="W303.chr09")
Xpar9.p <- subset(Xpar, Chr=="N_17.chr09")
XVpar9.c <- subset(XVpar, Chr=="W303.chr09")
XVpar9.p <- subset(XVpar, Chr=="N_17.chr09")


test <- get_cov_gene(cov_rest.cp, "depth.cer.c", "depth.cer.p", cer, "chr01")





#get coverage for genes
cov_rest.cp$depth.cer.c <- 0
cov_rest.cp$depth.cer.p <- 0
cov_rest.cp$depth.IXpar.c <- 0
cov_rest.cp$depth.IXpar.p <- 0
cov_rest.cp$depth.Xpar.c <- 0
cov_rest.cp$depth.Xpar.p <- 0
cov_rest.cp$depth.XVpar.c <- 0
cov_rest.cp$depth.XVpar.p <- 0
for (i in 1:length(cov_rest.cp$ID)){
  if(cov_rest.cp$strand.c[i]=="-"){
    cov_rest.cp$depth.cer.c[i] <- mean(cer.c$depth[(cov9.cp$end.c[i]-2):(cov9.cp$end.c[i])])
    cov_rest.cp$depth.IXpar.c[i] <- mean(IXpar9.c$depth[(cov9.cp$end.c[i]-2):(cov9.cp$end.c[i])])
    cov_rest.cp$depth.Xpar.c[i] <- mean(Xpar9.c$depth[(cov9.cp$end.c[i]-2):(cov9.cp$end.c[i])])
    cov_rest.cp$depth.XVpar.c[i] <- mean(XVpar9.c$depth[(cov9.cp$end.c[i]-2):(cov9.cp$end.c[i])])
  }else{
    cov9.cp$depth.cer.c[i] <- mean(cer9.c$depth[(cov9.cp$start.c[i]):(cov9.cp$start.c[i]+2)])
    cov9.cp$depth.IXpar.c[i] <- mean(IXpar9.c$depth[(cov9.cp$start.c[i]):(cov9.cp$start.c[i]+2)])
    cov9.cp$depth.Xpar.c[i] <- mean(Xpar9.c$depth[(cov9.cp$start.c[i]):(cov9.cp$start.c[i]+2)])
    cov9.cp$depth.XVpar.c[i] <- mean(XVpar9.c$depth[(cov9.cp$start.c[i]):(cov9.cp$start.c[i]+2)])
  }
  if(cov9.cp$strand.p[i]=="-"){
    cov9.cp$depth.cer.p[i] <- mean(cer9.p$depth[(cov9.cp$end.p[i]-2):(cov9.cp$end.p[i])])
    cov9.cp$depth.IXpar.p[i] <- mean(IXpar9.p$depth[(cov9.cp$end.p[i]-2):(cov9.cp$end.p[i])])
    cov9.cp$depth.Xpar.p[i] <- mean(Xpar9.p$depth[(cov9.cp$end.p[i]-2):(cov9.cp$end.p[i])])
    cov9.cp$depth.XVpar.p[i] <- mean(XVpar9.p$depth[(cov9.cp$end.p[i]-2):(cov9.cp$end.p[i])])
  }else{
    cov9.cp$depth.cer.p[i] <- mean(cer9.p$depth[(cov9.cp$start.p[i]):(cov9.cp$start.p[i]+2)])
    cov9.cp$depth.IXpar.p[i] <- mean(IXpar9.p$depth[(cov9.cp$start.p[i]):(cov9.cp$start.p[i]+2)])
    cov9.cp$depth.Xpar.p[i] <- mean(Xpar9.p$depth[(cov9.cp$start.p[i]):(cov9.cp$start.p[i]+2)])
    cov9.cp$depth.XVpar.p[i] <- mean(XVpar9.p$depth[(cov9.cp$start.p[i]):(cov9.cp$start.p[i]+2)])
  }
}



######here######

#for now just use global sum
cer_glob_cov <- sum(cer$depth)
IX_glob_cov <- sum(IXpar$depth)
X_glob_cov <- sum(Xpar$depth)
XV_glob_cov <- sum(XVpar$depth)

#pull out only chr of interest
#'for A - chr9
#'do by species
cer9.c <- subset(cer, Chr=="W303.chr09")
cer9.p <- subset(cer, Chr=="N_17.chr09")
IXpar9.c <- subset(IXpar, Chr=="W303.chr09")
IXpar9.p <- subset(IXpar, Chr=="N_17.chr09")
Xpar9.c <- subset(Xpar, Chr=="W303.chr09")
Xpar9.p <- subset(Xpar, Chr=="N_17.chr09")
XVpar9.c <- subset(XVpar, Chr=="W303.chr09")
XVpar9.p <- subset(XVpar, Chr=="N_17.chr09")



plot(cov9.cp$depth.cer.c, type="l")
lines(cov9.cp$depth.cer.p, col="pink", lty=2)
lines(cov9.cp$depth.IXpar.c, col="blue")
lines(cov9.cp$depth.IXpar.p, col="purple", lty=2)
lines(cov9.cp$depth.Xpar.c, col="green")
lines(cov9.cp$depth.Xpar.p, col="darkgreen", lty=2)
lines(cov9.cp$depth.XVpar.c, col="red")
lines(cov9.cp$depth.XVpar.p, col="darkred", lty=2)
#region A
grep("YIL166C", cov9.cp$ID)
#2
abline(v=2, lty=2, col="forestgreen")
grep("YIL156W", cov9.cp$ID)
cov9.cp[22:26,]
#26 real one
abline(v=26, lty=2, col="forestgreen")
#region with rec selected
grep("YIL151C", cov9.cp$ID)
#36
abline(v=36, lty=4, col="purple")

grep("YIL169C", cov9.cp$ID)
#not there - just the end somewhere?

#pull out only genes
row_odd <- seq_len(nrow(cov9.cp)) %% 2
cov9.cp.genes <- cov9.cp[row_odd==0,]

#plot just gene info
plot(cov9.cp.genes$depth.cer.c, type="l", ylim=c(0,50))
lines(cov9.cp.genes$depth.cer.p, col="pink", lty=2)
lines(cov9.cp.genes$depth.IXpar.c, col="blue")
lines(cov9.cp.genes$depth.IXpar.p, col="purple", lty=2)
lines(cov9.cp.genes$depth.Xpar.c, col="green")
lines(cov9.cp.genes$depth.Xpar.p, col="darkgreen", lty=2)
lines(cov9.cp.genes$depth.XVpar.c, col="red")
lines(cov9.cp.genes$depth.XVpar.p, col="darkred", lty=2)
#gets rid of some peaks but not much diff otherwise






#########find DNA loading correction + global expected val - A##########
#find factor to multiply each gene by
loadcorr_A.c <- IX_glob_cov/cer_glob_cov
loadcorr_A.p <- cer_glob_cov/IX_glob_cov

#multiply each gene
cov9.cp.A.genes$dnacor.cer.c <- loadcorr_A.c*cov9.cp.A.genes$depth.cer.c
cov9.cp.A.genes$dnacor.IXpar.p <- loadcorr_A.p*cov9.cp.A.genes$depth.IXpar.p

#find global expected value
A_cov_both_par <- c(cov9.cp.A.genes$dnacor.cer.c,cov9.cp.A.genes$dnacor.IXpar.p)
mean(A_cov_both_par)
#34.21748



#########find genes for region C - chr10########
#limit to only genes on chr interested in
genes10.c <- subset(genes.c, seqid=="chr10")
genes10.p <- subset(genes.p, seqid=="chr10")


#find which genes in both and in same order
c.geneind10 <-c()
for (i in 1:length(genes10.c$ID)){
  c.geneind10 <- c(c.geneind10, which(genes10.c$ID%in%genes10.p$ID[i]))
}

#in order
is.sorted(c.geneind10)

#make n17 (par) equivalent
p.geneind10 <- c()
for (i in 1:length(c.geneind10)){
  id <- as.character(genes10.c$ID[c.geneind10[i]])
  p.geneind10 <- c(p.geneind10, which(genes10.p$ID==id))
}
#in order
is.sorted(p.geneind10)

#pull out only genes that are in both
genes.c.both10 <- genes10.c[c.geneind10,]
genes.p.both10 <- genes10.p[p.geneind10,]

#make data frame with genes and between genes regions
cer10.c <- subset(cer, Chr=="W303.chr10")
cer10.p <- subset(cer, Chr=="N_17.chr10")
ID <- 1
start.c <- 1
end.c <- c()
start.p <- 1
end.p <- c()
for (i in 1:length(genes.c.both10$X)){
  ID <- c(ID, as.character(genes.c.both10$ID[i]), i+1)
  start.c <- c(start.c, genes.c.both10$start[i], (genes.c.both10$end[i]+1))
  end.c <- c(end.c, (genes.c.both10$start[i]-1), genes.c.both10$end[i])
  start.p <- c(start.p, genes.p.both10$start[i], (genes.p.both10$end[i]+1))
  end.p <- c(end.p, (genes.p.both10$start[i]-1), genes.p.both10$end[i])
}
end.c <- c(end.c, max(cer10.c$locus))
end.p <- c(end.p, max(cer10.p$locus))
cov10.cp <- data.frame(ID, start.c, end.c, start.p, end.p)

#########figure out coverage of each gene for each sample - region C######
#pull out only chr of interest
#'for C - chr10
#'do by species
cer10.c <- subset(cer, Chr=="W303.chr10")
cer10.p <- subset(cer, Chr=="N_17.chr10")
IXpar10.c <- subset(IXpar, Chr=="W303.chr10")
IXpar10.p <- subset(IXpar, Chr=="N_17.chr10")
Xpar10.c <- subset(Xpar, Chr=="W303.chr10")
Xpar10.p <- subset(Xpar, Chr=="N_17.chr10")
XVpar10.c <- subset(XVpar, Chr=="W303.chr10")
XVpar10.p <- subset(XVpar, Chr=="N_17.chr10")


cov10.cp$depth.cer.c <- 0
cov10.cp$depth.cer.p <- 0
cov10.cp$depth.IXpar.c <- 0
cov10.cp$depth.IXpar.p <- 0
cov10.cp$depth.Xpar.c <- 0
cov10.cp$depth.Xpar.p <- 0
cov10.cp$depth.XVpar.c <- 0
cov10.cp$depth.XVpar.p <- 0
for (i in 1:length(cov10.cp$ID)){
  cov10.cp$depth.cer.c[i] <- mean(cer10.c$depth[cov10.cp$start.c[i]:cov10.cp$end.c[i]])
  cov10.cp$depth.cer.p[i] <- mean(cer10.p$depth[cov10.cp$start.p[i]:cov10.cp$end.p[i]])
  cov10.cp$depth.IXpar.c[i] <- mean(IXpar10.c$depth[cov10.cp$start.c[i]:cov10.cp$end.c[i]])
  cov10.cp$depth.IXpar.p[i] <- mean(IXpar10.p$depth[cov10.cp$start.p[i]:cov10.cp$end.p[i]])
  cov10.cp$depth.Xpar.c[i] <- mean(Xpar10.c$depth[cov10.cp$start.c[i]:cov10.cp$end.c[i]])
  cov10.cp$depth.Xpar.p[i] <- mean(Xpar10.p$depth[cov10.cp$start.p[i]:cov10.cp$end.p[i]])
  cov10.cp$depth.XVpar.c[i] <- mean(XVpar10.c$depth[cov10.cp$start.c[i]:cov10.cp$end.c[i]])
  cov10.cp$depth.XVpar.p[i] <- mean(XVpar10.p$depth[cov10.cp$start.p[i]:cov10.cp$end.p[i]])
}


plot(cov10.cp$depth.cer.c, type="l")
lines(cov10.cp$depth.cer.p, col="pink", lty=2)
lines(cov10.cp$depth.IXpar.c, col="blue")
lines(cov10.cp$depth.IXpar.p, col="purple", lty=2)
lines(cov10.cp$depth.Xpar.c, col="green")
lines(cov10.cp$depth.Xpar.p, col="darkgreen", lty=2)
lines(cov10.cp$depth.XVpar.c, col="red")
lines(cov10.cp$depth.XVpar.p, col="darkred", lty=2)
#region C
grep("YJL218W", cov10.cp$ID)
#2
abline(v=2, lty=2, col="forestgreen")
grep("YJL165C", cov10.cp$ID)
#108
abline(v=108, lty=2, col="forestgreen")
#region with rec selected
grep("YJL164C", cov10.cp$ID)
#110
abline(v=110, lty=4, col="purple")

grep("YJL219W", cov10.cp$ID)
#not there - just the end somewhere?


#pull out only genes
row_odd <- seq_len(nrow(cov10.cp)) %% 2
cov10.cp.genes <- cov10.cp[row_odd==0,]


#plot just gene info
plot(cov10.cp.genes$depth.cer.c, type="l", ylim=c(0,50))
lines(cov10.cp.genes$depth.cer.p, col="pink", lty=2)
lines(cov10.cp.genes$depth.IXpar.c, col="blue")
lines(cov10.cp.genes$depth.IXpar.p, col="purple", lty=2)
lines(cov10.cp.genes$depth.Xpar.c, col="green")
lines(cov10.cp.genes$depth.Xpar.p, col="darkgreen", lty=2)
lines(cov10.cp.genes$depth.XVpar.c, col="red")
lines(cov10.cp.genes$depth.XVpar.p, col="darkred", lty=2)
#gets rid of some peaks but not much diff otherwise


#####only find region of interest - region C####
cov10.cp.C <- cov10.cp[1:110,]
cov10.cp.C.genes <- cov10.cp.genes[1:55,]

#all blocks
plot(cov10.cp.C$depth.cer.c, type="l", ylim=c(0,100))
lines(cov10.cp.C$depth.cer.p, col="pink", lty=2)
lines(cov10.cp.C$depth.IXpar.c, col="blue")
lines(cov10.cp.C$depth.IXpar.p, col="purple", lty=2)
lines(cov10.cp.C$depth.Xpar.c, col="green")
lines(cov10.cp.C$depth.Xpar.p, col="darkgreen", lty=2)
lines(cov10.cp.C$depth.XVpar.c, col="red")
lines(cov10.cp.C$depth.XVpar.p, col="darkred", lty=2)

#only genes
plot(cov10.cp.C.genes$depth.cer.c, type="l", ylim=c(0,100))
lines(cov10.cp.C.genes$depth.cer.p, col="pink", lty=2)
lines(cov10.cp.C.genes$depth.IXpar.c, col="blue")
lines(cov10.cp.C.genes$depth.IXpar.p, col="purple", lty=2)
lines(cov10.cp.C.genes$depth.Xpar.c, col="green")
lines(cov10.cp.C.genes$depth.Xpar.p, col="darkgreen", lty=2)
lines(cov10.cp.C.genes$depth.XVpar.c, col="red")
lines(cov10.cp.C.genes$depth.XVpar.p, col="darkred", lty=2)
#definitely less peaky with genes but still some variation
#pattern not exactly the same with both cerevisiae, but one is s288c?
#more similar pattern than A?
#also some mis-mapping this time

mod <- lm(cov10.cp.C.genes$depth.XVpar.c~cov10.cp.C.genes$depth.IXpar.c)
summary(mod)

#####calculate even coverage#####
#the problem seems to be in reliable coverage across the region
plot(cov10.cp.C.genes$depth.cer.c/cov10.cp.C.genes$depth.cer.c)
lines(cov10.cp.C.genes$depth.IXpar.c/cov10.cp.C.genes$depth.cer.c)
lines(cov10.cp.C.genes$depth.XVpar.c/cov10.cp.C.genes$depth.cer.c, col="green")
lines(cov10.cp.C.genes$depth.XVpar.c/cov10.cp.C.genes$depth.IXpar.c)
lines(cov10.cp.C.genes$depth.XVpar.c/30, col="red")


plot(cov10.cp.C.genes$depth.Xpar.p)
plot(cov10.cp.C.genes$depth.Xpar.c)
#one gene in this region obviously mismapping but that's about it


summary(lm(cov10.cp.C.genes$depth.IXpar.c~cov10.cp.C.genes$depth.XVpar.c))
plot(cov10.cp.C.genes$depth.IXpar.c~cov10.cp.C.genes$depth.XVpar.c)

#ctl = IXpar
IXparavg.c <- mean(cov10.cp.C.genes$depth.IXpar.c)
cov10.cp.C.genes$c_IXpar.c <- IXparavg.c/cov10.cp.C.genes$depth.IXpar.c
cov10.cp.C.genes$correctedXVpar.c <- cov10.cp.C.genes$c_IXpar.c*cov10.cp.C.genes$depth.XVpar.c
plot(cov10.cp.C.genes$depth.IXpar.c, type="l", ylim=c(0,100), lty=2)
lines(cov10.cp.C.genes$depth.XVpar.c, col="blue")
lines(cov10.cp.C.genes$correctedXVpar.c, col="red")

plot(cov10.cp.C.genes$depth.IXpar.c~cov10.cp.C.genes$depth.XVpar.c)

var(cov10.cp.C.genes$depth.XVpar.c)
var(cov10.cp.C.genes$correctedXVpar.c)


#########find genes for region F - chr15########
#limit to only genes on chr interested in
genes15.c <- subset(genes.c, seqid=="chr15")
genes15.p <- subset(genes.p, seqid=="chr15")


#find which genes in both and in same order
c.geneind15 <-c()
for (i in 1:length(genes15.c$ID)){
  c.geneind15 <- c(c.geneind15, which(genes15.c$ID%in%genes15.p$ID[i]))
}

#in order
is.sorted(c.geneind15)

#make n17 (par) equivalent
p.geneind15 <- c()
for (i in 1:length(c.geneind15)){
  id <- as.character(genes15.c$ID[c.geneind15[i]])
  p.geneind15 <- c(p.geneind15, which(genes15.p$ID==id))
}
#in order
is.sorted(p.geneind15)

#pull out only genes that are in both
genes.c.both15 <- genes15.c[c.geneind15,]
genes.p.both15 <- genes15.p[p.geneind15,]

#make data frame with genes and between genes regions
cer15.c <- subset(cer, Chr=="W303.chr15")
cer15.p <- subset(cer, Chr=="N_17.chr15")
ID <- 1
start.c <- 1
end.c <- c()
start.p <- 1
end.p <- c()
for (i in 1:length(genes.c.both15$X)){
  ID <- c(ID, as.character(genes.c.both15$ID[i]), i+1)
  start.c <- c(start.c, genes.c.both15$start[i], (genes.c.both15$end[i]+1))
  end.c <- c(end.c, (genes.c.both15$start[i]-1), genes.c.both15$end[i])
  start.p <- c(start.p, genes.p.both15$start[i], (genes.p.both15$end[i]+1))
  end.p <- c(end.p, (genes.p.both15$start[i]-1), genes.p.both15$end[i])
}
end.c <- c(end.c, max(cer15.c$locus))
end.p <- c(end.p, max(cer15.p$locus))
cov15.cp <- data.frame(ID, start.c, end.c, start.p, end.p)

#########figure out coverage of each gene for each sample - region F######
#pull out only chr of interest
#'for F - chr15
#'do by species
cer15.c <- subset(cer, Chr=="W303.chr15")
cer15.p <- subset(cer, Chr=="N_17.chr15")
IXpar15.c <- subset(IXpar, Chr=="W303.chr15")
IXpar15.p <- subset(IXpar, Chr=="N_17.chr15")
Xpar15.c <- subset(Xpar, Chr=="W303.chr15")
Xpar15.p <- subset(Xpar, Chr=="N_17.chr15")
XVpar15.c <- subset(XVpar, Chr=="W303.chr15")
XVpar15.p <- subset(XVpar, Chr=="N_17.chr15")


cov15.cp$depth.cer.c <- 0
cov15.cp$depth.cer.p <- 0
cov15.cp$depth.IXpar.c <- 0
cov15.cp$depth.IXpar.p <- 0
cov15.cp$depth.Xpar.c <- 0
cov15.cp$depth.Xpar.p <- 0
cov15.cp$depth.XVpar.c <- 0
cov15.cp$depth.XVpar.p <- 0
for (i in 1:length(cov15.cp$ID)){
  cov15.cp$depth.cer.c[i] <- mean(cer15.c$depth[cov15.cp$start.c[i]:cov15.cp$end.c[i]])
  cov15.cp$depth.cer.p[i] <- mean(cer15.p$depth[cov15.cp$start.p[i]:cov15.cp$end.p[i]])
  cov15.cp$depth.IXpar.c[i] <- mean(IXpar15.c$depth[cov15.cp$start.c[i]:cov15.cp$end.c[i]])
  cov15.cp$depth.IXpar.p[i] <- mean(IXpar15.p$depth[cov15.cp$start.p[i]:cov15.cp$end.p[i]])
  cov15.cp$depth.Xpar.c[i] <- mean(Xpar15.c$depth[cov15.cp$start.c[i]:cov15.cp$end.c[i]])
  cov15.cp$depth.Xpar.p[i] <- mean(Xpar15.p$depth[cov15.cp$start.p[i]:cov15.cp$end.p[i]])
  cov15.cp$depth.XVpar.c[i] <- mean(XVpar15.c$depth[cov15.cp$start.c[i]:cov15.cp$end.c[i]])
  cov15.cp$depth.XVpar.p[i] <- mean(XVpar15.p$depth[cov15.cp$start.p[i]:cov15.cp$end.p[i]])
}


plot(cov15.cp$depth.cer.c, type="l")
lines(cov15.cp$depth.cer.p, col="pink", lty=2)
lines(cov15.cp$depth.IXpar.c, col="blue")
lines(cov15.cp$depth.IXpar.p, col="purple", lty=2)
lines(cov15.cp$depth.Xpar.c, col="green")
lines(cov15.cp$depth.Xpar.p, col="darkgreen", lty=2)
lines(cov15.cp$depth.XVpar.c, col="red")
lines(cov15.cp$depth.XVpar.p, col="darkred", lty=2)
#region F
grep("YOL138C", cov15.cp$ID)
#48
abline(v=48, lty=2, col="forestgreen")
grep("YOL024W", cov15.cp$ID)
#280
abline(v=280, lty=2, col="forestgreen")
#region with rec selected
grep("YOL020W", cov15.cp$ID)
#288
abline(v=288, lty=4, col="purple")

grep("YOL141W", cov10.cp$ID)
#not there - just the end somewhere?


#pull out only genes
row_odd <- seq_len(nrow(cov15.cp)) %% 2
cov15.cp.genes <- cov15.cp[row_odd==0,]


#plot just gene info
plot(cov15.cp.genes$depth.cer.c, type="l", ylim=c(0,50))
lines(cov15.cp.genes$depth.cer.p, col="pink", lty=2)
lines(cov15.cp.genes$depth.IXpar.c, col="blue")
lines(cov15.cp.genes$depth.IXpar.p, col="purple", lty=2)
lines(cov15.cp.genes$depth.Xpar.c, col="green")
lines(cov15.cp.genes$depth.Xpar.p, col="darkgreen", lty=2)
lines(cov15.cp.genes$depth.XVpar.c, col="red")
lines(cov15.cp.genes$depth.XVpar.p, col="darkred", lty=2)
#gets rid of some peaks but not much diff otherwise


#####only find region of interest - region C####
cov10.cp.C <- cov10.cp[1:110,]
cov10.cp.C.genes <- cov10.cp.genes[1:55,]

#all blocks
plot(cov10.cp.C$depth.cer.c, type="l", ylim=c(0,100))
lines(cov10.cp.C$depth.cer.p, col="pink", lty=2)
lines(cov10.cp.C$depth.IXpar.c, col="blue")
lines(cov10.cp.C$depth.IXpar.p, col="purple", lty=2)
lines(cov10.cp.C$depth.Xpar.c, col="green")
lines(cov10.cp.C$depth.Xpar.p, col="darkgreen", lty=2)
lines(cov10.cp.C$depth.XVpar.c, col="red")
lines(cov10.cp.C$depth.XVpar.p, col="darkred", lty=2)

#only genes
plot(cov10.cp.C.genes$depth.cer.c, type="l", ylim=c(0,100))
lines(cov10.cp.C.genes$depth.cer.p, col="pink", lty=2)
lines(cov10.cp.C.genes$depth.IXpar.c, col="blue")
lines(cov10.cp.C.genes$depth.IXpar.p, col="purple", lty=2)
lines(cov10.cp.C.genes$depth.Xpar.c, col="green")
lines(cov10.cp.C.genes$depth.Xpar.p, col="darkgreen", lty=2)
lines(cov10.cp.C.genes$depth.XVpar.c, col="red")
lines(cov10.cp.C.genes$depth.XVpar.p, col="darkred", lty=2)
#definitely less peaky with genes but still some variation
#pattern not exactly the same with both cerevisiae, but one is s288c?
#more similar pattern than A?
#also some mis-mapping this time

#####calculate even coverage#####
#the problem seems to be in reliable coverage across the region
plot(cov10.cp.C.genes$depth.cer.c/cov10.cp.C.genes$depth.cer.c)
lines(cov10.cp.C.genes$depth.IXpar.c/cov10.cp.C.genes$depth.cer.c)
lines(cov10.cp.C.genes$depth.XVpar.c/cov10.cp.C.genes$depth.cer.c, col="green")
lines(cov10.cp.C.genes$depth.XVpar.c/cov10.cp.C.genes$depth.IXpar.c)
lines(cov10.cp.C.genes$depth.XVpar.c/30, col="red")


plot(cov10.cp.C.genes$depth.Xpar.p)
plot(cov10.cp.C.genes$depth.Xpar.c)
#one gene in this region obviously mismapping but that's about it



#########old###################

#####find coverage ratio#####
#find proportion paradoxus
cov.cp.A$fp.a1 <- cov.cp.A$depth.a1.p/(cov.cp.A$depth.a1.c+cov.cp.A$depth.a1.p)
cov.cp.A$fp.a2 <- cov.cp.A$depth.a2.p/(cov.cp.A$depth.a2.c+cov.cp.A$depth.a2.p)
plot(cov.cp.A$fp.a1, type="l")
lines(cov.cp.A$fp.a2)

cov.cp.A.genes$fp.a1 <- cov.cp.A.genes$depth.a1.p/(cov.cp.A.genes$depth.a1.c+cov.cp.A.genes$depth.a1.p)
cov.cp.A.genes$fp.a2 <- cov.cp.A.genes$depth.a2.p/(cov.cp.A.genes$depth.a2.c+cov.cp.A.genes$depth.a2.p)
plot(cov.cp.A.genes$fp.a1, type="l", ylab="proportion paradoxus")
lines(cov.cp.A.genes$fp.a2, col="blue")
#looking at only genes is smoother - continue with this

a1.diffs <- c(0)
for (i in 2:length(cov.cp.A.genes$ID)){
  a1.diffs <- c(a1.diffs, cov.cp.A.genes$fp.a1[i]-cov.cp.A.genes$fp.a1[i-1])
}

a2.diffs <- c(0)
for (i in 2:length(cov.cp.A.genes$ID)){
  a2.diffs <- c(a2.diffs, cov.cp.A.genes$fp.a2[i]-cov.cp.A.genes$fp.a2[i-1])
}

#add to plot the diffs
#for a1 do *-1 so positive number if decrease in p
lines(a1.diffs*-1, col="black")
points(a1.diffs*-1, col="black")
#a1 decr in par most around 15
#a1 also decr in par around 10
lines((a2.diffs), col="blue")
points(a2.diffs, col="blue")
#a2 incr in par most around 15
#a2 slower incr in par around 10
abline(h=0, lty=2, col="red")



######find min of fp1 and fp2#####
#plot fpa1 +fpa2 to find min - subtract 0.5 to plot
lines(cov.cp.A.genes$fp.a1+cov.cp.A.genes$fp.a2-0.5, col="red")
#make vector to find min
fpoverall <- cov.cp.A.genes$fp.a1+cov.cp.A.genes$fp.a2
grep(min(fpoverall), fpoverall)
#11 is the min - 0.9542517

cov.cp.A.genes$totaldepth.a1 <- cov.cp.A.genes$depth.a1.c + cov.cp.A.genes$depth.a1.p
cov.cp.A.genes$totaldepth.a2 <- cov.cp.A.genes$depth.a2.c + cov.cp.A.genes$depth.a2.p

mean(cov.cp.A.genes$totaldepth.a1)
mean(cov.cp.A.genes$totaldepth.a2)
#########oldA########
#load in A1 coverage
A1 <- read.table("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/A1.coverage",
                                header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

A1 <- rename(A1, c(V1="Chr", V2="locus", V3="depth")) # renames the header

#load in A2 coverage
A2 <- read.table("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/A2.coverage",
                 header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

A2 <- rename(A2, c(V1="Chr", V2="locus", V3="depth")) # renames the header

#pull out only chr of interest
#'for A - chr9
#'do by species
A1_chr9.c <- subset(A1, Chr=="W303.chr09")
A1_chr9.p <- subset(A1, Chr=="N_17.chr09")
A2_chr9.c <- subset(A2, Chr=="W303.chr09")
A2_chr9.p <- subset(A2, Chr=="N_17.chr09")

######
#load in gene files
genes.c <- read.csv("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/W303.genes.csv", head=TRUE)
genes.p <- read.csv("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/N17.genes.csv", head=TRUE)

#limit to only genes on chr interested in
genes_chr9.c <- subset(genes.c, seqid=="chr09")
genes_chr9.p <- subset(genes.p, seqid=="chr09")


#find which genes in both and in same order
c.geneind <-c()
for (i in 1:length(genes_chr9.c$ID)){
  c.geneind <- c(c.geneind, which(genes_chr9.c$ID%in%genes_chr9.p$ID[i]))
}

#in order
is.sorted(c.geneind)

#make n17 (par) equivalent
p.geneind <- c()
for (i in 1:length(c.geneind)){
  id <- as.character(genes_chr9.c$ID[c.geneind[i]])
  p.geneind <- c(p.geneind, which(genes_chr9.p$ID==id))
}
#in order
is.sorted(p.geneind)



#pull out only genes that are in both
genes.c.both <- genes_chr9.c[c.geneind,]
genes.p.both <- genes_chr9.p[p.geneind,]



#make data frame with genes and between genes regions
ID <- 1
start.c <- 1
end.c <- c()
start.p <- 1
end.p <- c()
for (i in 1:length(genes.c.both$X)){
      ID <- c(ID, as.character(genes.c.both$ID[i]), i+1)
      start.c <- c(start.c, genes.c.both$start[i], (genes.c.both$end[i]+1))
      end.c <- c(end.c, (genes.c.both$start[i]-1), genes.c.both$end[i])
      start.p <- c(start.p, genes.p.both$start[i], (genes.p.both$end[i]+1))
      end.p <- c(end.p, (genes.p.both$start[i]-1), genes.p.both$end[i])
  }
    end.c <- c(end.c, max(A1_chr9.c$locus))
    end.p <- c(end.p, max(A1_chr9.p$locus))
cov.cp <- data.frame(ID, start.c, end.c, start.p, end.p)

###figure out coverage of each gene
cov.cp$depth.a1.c <- 0
cov.cp$depth.a1.p <- 0
cov.cp$depth.a2.c <- 0
cov.cp$depth.a2.p <- 0
for (i in 1:length(cov.cp$ID)){
       cov.cp$depth.a1.c[i] <- mean(A1_chr9.c$depth[cov.cp$start.c[i]:cov.cp$end.c[i]])
       cov.cp$depth.a1.p[i] <- mean(A1_chr9.p$depth[cov.cp$start.p[i]:cov.cp$end.p[i]])
       cov.cp$depth.a2.c[i] <- mean(A2_chr9.c$depth[cov.cp$start.c[i]:cov.cp$end.c[i]])
       cov.cp$depth.a2.p[i] <- mean(A2_chr9.p$depth[cov.cp$start.p[i]:cov.cp$end.p[i]])
     }


plot(cov.cp$depth.a1.c, type="l")
lines(cov.cp$depth.a1.p, col="black", lty=2)
lines(cov.cp$depth.a2.c, col="blue")
lines(cov.cp$depth.a2.p, col="blue", lty=2)
#region A
grep("YIL166C", cov.cp$ID)
#2
abline(v=2, lty=2, col="forestgreen")
grep("YIL156W", cov.cp$ID)
cov.cp[22:26,]
#26 real one
abline(v=26, lty=2, col="forestgreen")
#region with rec selected
grep("YIL151C", cov.cp$ID)
#36
abline(v=36, lty=4, col="purple")

grep("YIL169C", cov.cp$ID)
#not there - just the end somewhere?
grep("YIL043C", cov.cp$ID)
#258
abline(v=258, lty=4, col="red")

#pull out only genes
row_odd <- seq_len(nrow(cov.cp)) %% 2
cov.cp.genes <- cov.cp[row_odd==0,]
cov.cp.genes

#plot just gene info
plot(cov.cp.genes$depth.a1.c, type="l")
lines(cov.cp.genes$depth.a1.p, col="black", lty=2)
lines(cov.cp.genes$depth.a2.c, col="blue")
lines(cov.cp.genes$depth.a2.p, col="blue", lty=2)
#gets rid of weird peak but not much diff otherwise

#####only find region of interest####
cov.cp.A <- cov.cp[1:36,]
cov.cp.A.genes <- cov.cp.genes[1:18,]

#all blocks
plot(cov.cp.A$depth.a1.c, type="l")
lines(cov.cp.A$depth.a1.p, col="black", lty=2)
lines(cov.cp.A$depth.a2.c, col="blue")
lines(cov.cp.A$depth.a2.p, col="blue", lty=2)

#only genes
plot(cov.cp.A.genes$depth.a1.c, type="l")
lines(cov.cp.A.genes$depth.a1.p, col="black", lty=2)
lines(cov.cp.A.genes$depth.a2.c, col="blue")
lines(cov.cp.A.genes$depth.a2.p, col="blue", lty=2)


#####find coverage ratio#####
#find proportion paradoxus
cov.cp.A$fp.a1 <- cov.cp.A$depth.a1.p/(cov.cp.A$depth.a1.c+cov.cp.A$depth.a1.p)
cov.cp.A$fp.a2 <- cov.cp.A$depth.a2.p/(cov.cp.A$depth.a2.c+cov.cp.A$depth.a2.p)
plot(cov.cp.A$fp.a1, type="l")
lines(cov.cp.A$fp.a2)

cov.cp.A.genes$fp.a1 <- cov.cp.A.genes$depth.a1.p/(cov.cp.A.genes$depth.a1.c+cov.cp.A.genes$depth.a1.p)
cov.cp.A.genes$fp.a2 <- cov.cp.A.genes$depth.a2.p/(cov.cp.A.genes$depth.a2.c+cov.cp.A.genes$depth.a2.p)
plot(cov.cp.A.genes$fp.a1, type="l", ylab="proportion paradoxus")
lines(cov.cp.A.genes$fp.a2, col="blue")
#looking at only genes is smoother - continue with this

a1.diffs <- c(0)
for (i in 2:length(cov.cp.A.genes$ID)){
  a1.diffs <- c(a1.diffs, cov.cp.A.genes$fp.a1[i]-cov.cp.A.genes$fp.a1[i-1])
}

a2.diffs <- c(0)
for (i in 2:length(cov.cp.A.genes$ID)){
  a2.diffs <- c(a2.diffs, cov.cp.A.genes$fp.a2[i]-cov.cp.A.genes$fp.a2[i-1])
}

#add to plot the diffs
#for a1 do *-1 so positive number if decrease in p
lines(a1.diffs*-1, col="black")
points(a1.diffs*-1, col="black")
#a1 decr in par most around 15
#a1 also decr in par around 10
lines((a2.diffs), col="blue")
points(a2.diffs, col="blue")
#a2 incr in par most around 15
#a2 slower incr in par around 10
abline(h=0, lty=2, col="red")



######find min of fp1 and fp2#####
#plot fpa1 +fpa2 to find min - subtract 0.5 to plot
lines(cov.cp.A.genes$fp.a1+cov.cp.A.genes$fp.a2-0.5, col="red")
#make vector to find min
fpoverall <- cov.cp.A.genes$fp.a1+cov.cp.A.genes$fp.a2
grep(min(fpoverall), fpoverall)
#11 is the min - 0.9542517

cov.cp.A.genes$totaldepth.a1 <- cov.cp.A.genes$depth.a1.c + cov.cp.A.genes$depth.a1.p
cov.cp.A.genes$totaldepth.a2 <- cov.cp.A.genes$depth.a2.c + cov.cp.A.genes$depth.a2.p

mean(cov.cp.A.genes$totaldepth.a1)
mean(cov.cp.A.genes$totaldepth.a2)


######do it all again for C#####
#setwd("Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/")

#load in C1 coverage
C1 <- read.table("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/C1.coverage",
                 header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

C1 <- rename(C1, c(V1="Chr", V2="locus", V3="depth")) # renames the header

#load in C2 coverage
C2 <- read.table("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/C2.coverage",
                 header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

C2 <- rename(C2, c(V1="Chr", V2="locus", V3="depth")) # renames the header

#pull out only chr of interest
#'for C - chr10
#'do by species
C1_chr10.c <- subset(C1, Chr=="W303.chr10")
C1_chr10.p <- subset(C1, Chr=="N_17.chr10")
C2_chr10.c <- subset(C2, Chr=="W303.chr10")
C2_chr10.p <- subset(C2, Chr=="N_17.chr10")

######
#load in gene files
genes.c <- read.csv("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/W303.genes.csv", head=TRUE)
genes.p <- read.csv("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/N17.genes.csv", head=TRUE)

#limit to only genes on chr interested in
genes_chr10.c <- subset(genes.c, seqid=="chr10")
genes_chr10.p <- subset(genes.p, seqid=="chr10")


#find which genes in both and in same order
c.geneind <-c()
for (i in 1:length(genes_chr10.c$ID)){
  c.geneind <- c(c.geneind, which(genes_chr10.c$ID%in%genes_chr10.p$ID[i]))
}

#in order
is.sorted(c.geneind)

#make n17 (par) equivalent
p.geneind <- c()
for (i in 1:length(c.geneind)){
  id <- as.character(genes_chr10.c$ID[c.geneind[i]])
  p.geneind <- c(p.geneind, which(genes_chr10.p$ID==id))
}
#in order
is.sorted(p.geneind)



#pull out only genes that are in both
genes.c.both.10 <- genes_chr10.c[c.geneind,]
genes.p.both.10 <- genes_chr10.p[p.geneind,]

#exclude YJL217W - not in paradoxus ref genome
genes.c.both.10 <- genes.c.both.10[-which(genes.c.both.10$ID=="YJL217W"),]
genes.p.both.10 <- genes.p.both.10[-which(genes.p.both.10$ID=="YJL217W"),]

#exclude YJL214W  - weird mapping to par
genes.c.both.10 <- genes.c.both.10[-which(genes.c.both.10$ID=="YJL214W"),]
genes.p.both.10 <- genes.p.both.10[-which(genes.p.both.10$ID=="YJL214W"),]



#make data frame with genes and between genes regions
ID <- 1
start.c <- 1
end.c <- c()
start.p <- 1
end.p <- c()
for (i in 1:length(genes.c.both.10$X)){
  ID <- c(ID, as.character(genes.c.both.10$ID[i]), i+1)
  start.c <- c(start.c, genes.c.both.10$start[i], (genes.c.both.10$end[i]+1))
  end.c <- c(end.c, (genes.c.both.10$start[i]-1), genes.c.both.10$end[i])
  start.p <- c(start.p, genes.p.both.10$start[i], (genes.p.both.10$end[i]+1))
  end.p <- c(end.p, (genes.p.both.10$start[i]-1), genes.p.both.10$end[i])
}
end.c <- c(end.c, max(C1_chr10.c$locus))
end.p <- c(end.p, max(C1_chr10.p$locus))
cov.cp <- data.frame(ID, start.c, end.c, start.p, end.p)

###figure out coverage of each gene
cov.cp$depth.c1.c <- 0
cov.cp$depth.c1.p <- 0
cov.cp$depth.c2.c <- 0
cov.cp$depth.c2.p <- 0
#for (i in 1:length(cov.cp$ID)){
#  cov.cp$depth.c1.c[i] <- mean(C1_chr10.c$depth[cov.cp$start.c[i]:cov.cp$end.c[i]])
#  cov.cp$depth.c1.p[i] <- mean(C1_chr10.p$depth[cov.cp$start.p[i]:cov.cp$end.p[i]])
#  cov.cp$depth.c2.c[i] <- mean(C2_chr10.c$depth[cov.cp$start.c[i]:cov.cp$end.c[i]])
#  cov.cp$depth.c2.p[i] <- mean(C2_chr10.p$depth[cov.cp$start.p[i]:cov.cp$end.p[i]])
#}

for (i in 1:length(cov.cp$ID)){
  cov.cp$depth.c1.c[i] <- median(C1_chr10.c$depth[cov.cp$start.c[i]:cov.cp$end.c[i]])
  cov.cp$depth.c1.p[i] <- median(C1_chr10.p$depth[cov.cp$start.p[i]:cov.cp$end.p[i]])
  cov.cp$depth.c2.c[i] <- median(C2_chr10.c$depth[cov.cp$start.c[i]:cov.cp$end.c[i]])
  cov.cp$depth.c2.p[i] <- median(C2_chr10.p$depth[cov.cp$start.p[i]:cov.cp$end.p[i]])
}


plot(cov.cp$depth.c1.c, type="l")
lines(cov.cp$depth.c1.p, col="black", lty=2)
lines(cov.cp$depth.c2.c, col="blue")
lines(cov.cp$depth.c2.p, col="blue", lty=2)
#region C
grep("YJL218W", cov.cp$ID)
#2
abline(v=2, lty=2, col="forestgreen")
grep("YJL165C", cov.cp$ID)
#106
abline(v=106, lty=2, col="forestgreen")
#region with rec selected
grep("YJL164C", cov.cp$ID)
#108
abline(v=108, lty=4, col="purple")

grep("YJL219W", cov.cp$ID)
#not there - just the end somewhere?


#pull out only genes
row_odd <- seq_len(nrow(cov.cp)) %% 2
cov.cp.genes <- cov.cp[row_odd==0,]
cov.cp.genes

#plot just gene info
plot(cov.cp.genes$depth.c1.c, type="l")
lines(cov.cp.genes$depth.c1.p, col="black", lty=2)
lines(cov.cp.genes$depth.c2.c, col="blue")
lines(cov.cp.genes$depth.c2.p, col="blue", lty=2)
#smoother

#####only find region of interest####
cov.cp.C <- cov.cp[1:108,]
cov.cp.C.genes <- cov.cp.genes[1:54,]

#all blocks
plot(cov.cp.C$depth.c1.c, type="l")
lines(cov.cp.C$depth.c1.p, col="black", lty=2)
lines(cov.cp.C$depth.c2.c, col="blue")
lines(cov.cp.C$depth.c2.p, col="blue", lty=2)

#only genes
plot(cov.cp.C.genes$depth.c1.c, type="l")
lines(cov.cp.C.genes$depth.c1.p, col="black", lty=2)
lines(cov.cp.C.genes$depth.c2.c, col="blue")
lines(cov.cp.C.genes$depth.c2.p, col="blue", lty=2)


#####find coverage ratio#####
#find proportion paradoxus
cov.cp.C$fp.c1 <- cov.cp.C$depth.c1.p/(cov.cp.C$depth.c1.c+cov.cp.C$depth.c1.p)
cov.cp.C$fp.c2 <- cov.cp.C$depth.c2.p/(cov.cp.C$depth.c2.c+cov.cp.C$depth.c2.p)
plot(cov.cp.C$fp.c1, type="l")
lines(cov.cp.C$fp.c2)

cov.cp.C.genes$fp.c1 <- cov.cp.C.genes$depth.c1.p/(cov.cp.C.genes$depth.c1.c+cov.cp.C.genes$depth.c1.p)
cov.cp.C.genes$fp.c2 <- cov.cp.C.genes$depth.c2.p/(cov.cp.C.genes$depth.c2.c+cov.cp.C.genes$depth.c2.p)
plot(cov.cp.C.genes$fp.c1, type="l", ylab="proportion paradoxus")
lines(cov.cp.C.genes$fp.c2, col="blue")

cov.cp.C.genes$index <- c(1:54)
lo <- loess(cov.cp.C.genes$fp.c1~cov.cp.C.genes$index)
plot(cov.cp.C.genes$index, cov.cp.C.genes$fp.c1, type="l")
#lines(predict(lo), col='black', lwd=2)
lo2 <- loess(cov.cp.C.genes$fp.c2~cov.cp.C.genes$index)
lines(cov.cp.C.genes$index, cov.cp.C.genes$fp.c2, type="l", col="blue")
#lines(predict(lo2), col='blue', lwd=2)
#lines(predict(lo)+predict(lo2)-0.5, col="red")
abline(h=0.5)
#make vector to find min
fpoverall <- predict(lo)+predict(lo2)
which(fpoverall==min(fpoverall))
min(fpoverall)
#54 is the min = 0.9809415
#before end it's 21 - 0.9919789


#looking at only genes is smoother - continue with this

c1.diffs <- c(0)
for (i in 2:length(cov.cp.C.genes$ID)){
  c1.diffs <- c(c1.diffs, cov.cp.C.genes$fp.c1[i]-cov.cp.C.genes$fp.c1[i-1])
}

c2.diffs <- c(0)
for (i in 2:length(cov.cp.C.genes$ID)){
  c2.diffs <- c(c2.diffs, cov.cp.C.genes$fp.c2[i]-cov.cp.C.genes$fp.c2[i-1])
}

#add to plot the diffs
#for a1 do *-1 so positive number if decrease in p
lines(c1.diffs*-1, col="black", lty=2)
points(c1.diffs*-1, col="black")
#a1 decr in par most around 15
#a1 also decr in par around 10
lines((c2.diffs), col="blue", lty=2)
points(c2.diffs, col="blue")
#a2 incr in par most around 15
#a2 slower incr in par around 10
abline(h=0, lty=2, col="red")



######find min of fp1 and fp2#####
#plot fpa1 +fpa2 to find min - subtract 0.5 to plot
lines(cov.cp.C.genes$fp.c1+cov.cp.C.genes$fp.c2-0.5, col="red")
#make vector to find min
fpoverall <- cov.cp.C.genes$fp.c1+cov.cp.C.genes$fp.c2
which(fpoverall==min(fpoverall))
#3 is the min = 0.8254203

cov.cp.C.genes$totaldepth.c1 <- cov.cp.C.genes$depth.c1.c + cov.cp.C.genes$depth.c1.p
cov.cp.C.genes$totaldepth.c2 <- cov.cp.C.genes$depth.c2.c + cov.cp.C.genes$depth.c2.p

mean(cov.cp.C.genes$totaldepth.c1)
#1218.79
mean(cov.cp.C.genes$totaldepth.c2)
#1147.195

cov.cp.C.genes[3,]
# ID start.c end.c start.p end.p depth.c1.c depth.c1.p depth.c2.c depth.c2.p     fp.c1      fp.c2
#6 YJL214W   30577 32286   15749 17456   216.6275    640.822    1088.01   92.12354 0.7473583 0.07806196
#totaldepth.c1 totaldepth.c2
#6      857.4495      1180.133
#no n's in par or cer
plot(C1_chr10.c[30577:32286,]$depth, ylim=c(0,2000))
points(C1_chr10.p[15749:17456,]$depth, col="red")
points(C2_chr10.c[30577:32286,]$depth)
points(C2_chr10.p[15749:17456,]$depth, col="red")

g31c<-C1_chr10.c[30577:32286,]
g31p <- C1_chr10.p[15749:17456,]
g32c<-C2_chr10.c[30577:32286,]
g32p <- C2_chr10.p[15749:17456,]

median(g31c$depth)*0.25
g31c[which(g31c$depth<55),]
median(g31p$depth)*0.25
g31p[which(g31p$depth<130.875),]$locus
#16484-16502, 16909-17071
median(g32c$depth)*0.25
which(g32c$depth<270.25)
median(g32p$depth)*0.25
g32p[which(g32p$depth<18.75),]$locus
#16479-16502, 16909-17094














lines((cov.cp.C.genes$totaldepth.c1/1215.644)-0.5, col="green", lty=3)
lines((cov.cp.C.genes$totaldepth.c2/1144.886)-0.5, col="purple", lty=3)

plot(cov.cp.C.genes$totaldepth.c1, col="green", ylim=c(0,2000), type="l")
plot(cov.cp.C.genes$totaldepth.c2, col="purple", ylim=c(0,2000), type="l")

tail(cov.cp.C.genes)
c <- C2_chr10.c[25663:114847,]
plot(c$depth, type="l", ylim=c(0,2000))
abline(h=0)
p <- C2_chr10.p[2529:102835,]
plot(p$depth, type="l", ylim=c(0,2000))
abline(h=0)


C1_chr10.c <- subset(C1, Chr=="W303.chr10")
C1_chr10.p <- subset(C1, Chr=="N_17.chr10")
C2_chr10.c <- subset(C2, Chr=="W303.chr10")
C2_chr10.p <- subset(C2, Chr=="N_17.chr10")


#####try for C W303 and N17 separate######
#load in C1W coverage
C1w <- read.table("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/C1W.coverage",
                 header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

C1w <- rename(C1w, c(V1="Chr", V2="locus", V3="depth")) # renames the header

#load in C1N coverage
C1n <- read.table("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/C1N.coverage",
                  header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

C1n <- rename(C1n, c(V1="Chr", V2="locus", V3="depth")) # renames the header

#load in C2W coverage - w303
C2w <- read.table("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/C2W.coverage",
                 header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

C2w <- rename(C2w, c(V1="Chr", V2="locus", V3="depth")) # renames the header

#load in C2N coverage - n17
C2n <- read.table("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/C2N.coverage",
                  header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

C2n <- rename(C2n, c(V1="Chr", V2="locus", V3="depth")) # renames the header

#pull out only chr of interest
#'for C - chr10
#'do by species
C1w_chr10.c <- subset(C1w, Chr=="W303.chr10")
C1n_chr10.p <- subset(C1n, Chr=="N_17.chr10")
C2w_chr10.c <- subset(C2w, Chr=="W303.chr10")
C2n_chr10.p <- subset(C2n, Chr=="N_17.chr10")

######
#find which genes in both and in same order
c.geneind10 <-c()
for (i in 1:length(genes_chr10.c$ID)){
  c.geneind10 <- c(c.geneind10, which(genes_chr10.c$ID%in%genes_chr10.p$ID[i]))
}

#in order
is.sorted(c.geneind10)

#make n17 (par) equivalent
p.geneind10 <- c()
for (i in 1:length(c.geneind10)){
  id <- as.character(genes_chr10.c$ID[c.geneind10[i]])
  p.geneind10 <- c(p.geneind10, which(genes_chr10.p$ID==id))
}
#in order
is.sorted(p.geneind10)



#pull out only genes that are in both
genes.c.both10 <- genes_chr10.c[c.geneind10,]
genes.p.both10 <- genes_chr10.p[p.geneind10,]

#exclude YJL217W - not in paradoxus ref genome
genes.c.both10 <- genes.c.both10[-which(genes.c.both10$ID=="YJL217W"),]
genes.p.both10 <- genes.p.both10[-which(genes.p.both10$ID=="YJL217W"),]





#make data frame with genes and between genes regions
ID <- 1
start.c <- 1
end.c <- c()
start.p <- 1
end.p <- c()
for (i in 1:length(genes.c.both10$X)){
  ID <- c(ID, as.character(genes.c.both10$ID[i]), i+1)
  start.c <- c(start.c, genes.c.both10$start[i], (genes.c.both10$end[i]+1))
  end.c <- c(end.c, (genes.c.both10$start[i]-1), genes.c.both10$end[i])
  start.p <- c(start.p, genes.p.both10$start[i], (genes.p.both10$end[i]+1))
  end.p <- c(end.p, (genes.p.both10$start[i]-1), genes.p.both10$end[i])
}
end.c <- c(end.c, max(C1w_chr10.c$locus))
end.p <- c(end.p, max(C1n_chr10.p$locus))
cov.cp10 <- data.frame(ID, start.c, end.c, start.p, end.p)

###figure out coverage of each gene
cov.cp10$depth.c1.c <- 0
cov.cp10$depth.c1.p <- 0
cov.cp10$depth.c2.c <- 0
cov.cp10$depth.c2.p <- 0
for (i in 1:length(cov.cp10$ID)){
  cov.cp10$depth.c1.c[i] <- mean(C1w_chr10.c$depth[cov.cp10$start.c[i]:cov.cp10$end.c[i]])
  cov.cp10$depth.c1.p[i] <- mean(C1n_chr10.p$depth[cov.cp10$start.p[i]:cov.cp10$end.p[i]])
  cov.cp10$depth.c2.c[i] <- mean(C2w_chr10.c$depth[cov.cp10$start.c[i]:cov.cp10$end.c[i]])
  cov.cp10$depth.c2.p[i] <- mean(C2n_chr10.p$depth[cov.cp10$start.p[i]:cov.cp10$end.p[i]])
}


plot(cov.cp10$depth.c1.c, type="l")
lines(cov.cp10$depth.c1.p, col="black", lty=2)
lines(cov.cp10$depth.c2.c, col="blue")
lines(cov.cp10$depth.c2.p, col="blue", lty=2)
#region C
grep("YJL218W", cov.cp10$ID)
#2
abline(v=2, lty=2, col="forestgreen")
grep("YJL165C", cov.cp10$ID)
#106
abline(v=106, lty=2, col="forestgreen")
#region with rec selected
grep("YJL164C", cov.cp10$ID)
#108
abline(v=108, lty=4, col="purple")

grep("YJL219W", cov.cp10$ID)
#not there - just the end somewhere?


#pull out only genes
row_odd <- seq_len(nrow(cov.cp10)) %% 2
cov.cp10.genes <- cov.cp10[row_odd==0,]
cov.cp10.genes

#plot just gene info
plot(cov.cp10.genes$depth.c1.c, type="l")
lines(cov.cp10.genes$depth.c1.p, col="black", lty=2)
lines(cov.cp10.genes$depth.c2.c, col="blue")
lines(cov.cp10.genes$depth.c2.p, col="blue", lty=2)
#gets rid of lot of noise but not much diff otherwise

#####only find region of interest####
cov.cp10.C <- cov.cp10[1:108,]
cov.cp10.C.genes <- cov.cp10.genes[1:54,]

#all blocks
plot(cov.cp10.C$depth.c1.c, type="l")
lines(cov.cp10.C$depth.c1.p, col="black", lty=2)
lines(cov.cp10.C$depth.c2.c, col="blue")
lines(cov.cp10.C$depth.c2.p, col="blue", lty=2)

#only genes
plot(cov.cp10.C.genes$depth.c1.c, type="l")
lines(cov.cp10.C.genes$depth.c1.p, col="black", lty=2)
lines(cov.cp10.C.genes$depth.c2.c, col="blue")
lines(cov.cp10.C.genes$depth.c2.p, col="blue", lty=2)


#####find coverage ratio#####
#find proportion paradoxus
cov.cp10.C$fp.c1 <- cov.cp10.C$depth.c1.p/(cov.cp10.C$depth.c1.c+cov.cp10.C$depth.c1.p)
cov.cp10.C$fp.c2 <- cov.cp10.C$depth.c2.p/(cov.cp10.C$depth.c2.c+cov.cp10.C$depth.c2.p)
plot(cov.cp10.C$fp.c1, type="l")
lines(cov.cp10.C$fp.c2, col="blue")

cov.cp10.C.genes$fp.c1 <- cov.cp10.C.genes$depth.c1.p/(cov.cp10.C.genes$depth.c1.c+cov.cp10.C.genes$depth.c1.p)
cov.cp10.C.genes$fp.c2 <- cov.cp10.C.genes$depth.c2.p/(cov.cp10.C.genes$depth.c2.c+cov.cp10.C.genes$depth.c2.p)
plot(cov.cp10.C.genes$fp.c1, type="l", ylab="proportion paradoxus",ylim=c(0,1))
lines(cov.cp10.C.genes$fp.c2, col="blue")
#looking at only genes is smoother - continue with this

c1.diffs <- c(0)
for (i in 2:length(cov.cp10.C.genes$ID)){
  c1.diffs <- c(c1.diffs, cov.cp10.C.genes$fp.c1[i]-cov.cp10.C.genes$fp.c1[i-1])
}

c2.diffs <- c(0)
for (i in 2:length(cov.cp10.C.genes$ID)){
  c2.diffs <- c(c2.diffs, cov.cp10.C.genes$fp.c2[i]-cov.cp10.C.genes$fp.c2[i-1])
}

#add to plot the diffs
#for c1 do *-1 so positive number if decrease in p
lines(c1.diffs*-1, col="black")
points(c1.diffs*-1, col="black")
#a1 decr in par most around 15
#a1 also decr in par around 10
lines((c2.diffs), col="blue")
points(c2.diffs, col="blue")
#a2 incr in par most around 15
#a2 slower incr in par around 10
abline(h=0, lty=2, col="red")



######find min of fp1 and fp2#####
#plot fpc1 +fpc2 to find min - subtract 0.5 to plot
lines(cov.cp10.C.genes$fp.c1+cov.cp10.C.genes$fp.c2-0.5, col="red")
#make vector to find min
fpoverall <- cov.cp10.C.genes$fp.c1+cov.cp10.C.genes$fp.c2
grep(min(fpoverall), fpoverall)
#3 is the min - 0.7098976

cov.cp10.C.genes$totaldepth.c1 <- cov.cp10.C.genes$depth.c1.c + cov.cp10.C.genes$depth.c1.p
cov.cp10.C.genes$totaldepth.c2 <- cov.cp10.C.genes$depth.c2.c + cov.cp10.C.genes$depth.c2.p

mean(cov.cp10.C.genes$totaldepth.c1)
#2241.07
mean(cov.cp10.C.genes$totaldepth.c2)
#2097.208
#reads mapping to both probably

lines((cov.cp10.C.genes$totaldepth.c1/2241.07)-0.5, col="green", lty=3)
lines((cov.cp10.C.genes$totaldepth.c2/2097.208)-0.5, col="purple", lty=3)


#####################


cov.wn <- list()
x<-1
for (j in unique(genes.c.both$seqid)){
  subc <- subset(genes.c.both, seqid==j)
  subp <- subset(genes.p.both, seqid==j)
  chr <- as.character(subc$seqid[1])
  ID <- 1
  start.c <- 1
  end.c <- c()
  start.p <- 1
  end.p <- c()
  for (i in 1:length(subc$X)){
    chr <- c(chr, as.character(subc$seqid[i]), as.character(subc$seqid[i]))
    ID <- c(ID, as.character(subc$ID[i]), i+1)
    start.c <- c(start.c, subc$start[i], (subc$end[i]+1))
    end.c <- c(end.c, (subc$start[i]-1), subc$end[i])
    start.p <- c(start.p, subp$start[i], (subp$end[i]+1))
    end.p <- c(end.p, (subp$start[i]-1), subp$end[i])
  }
  end.c <- c(end.c, max(sample007.w[grep(chr[1], sample007.w$Chr),]$locus))
  end.p <- c(end.p, max(sample007.n[grep(chr[1], sample007.n$Chr),]$locus))
  cov.wn[[x]] <- data.frame(chr, ID, start.w, end.w, start.n, end.n)
  x <- x+1
}

#########


###figure out coverage of each gene
for (i in 1:length(cov.wn)){

}
head(cov.wn[[1]])
head(sample007.w)

cov.wn[[1]]$chr

grep(cov.wn[[1]]$chr[1], sample007.w)


cov.wn$depth.w <- 0
cov.wn$depth.n <- 0
for (i in 1:length(cov.wn$ID)){
  cov.wn$depth.w[i] <- mean(sample007_chr2.w$depth[cov.wn$start.w[i]:cov.wn$end.w[i]])
  cov.wn$depth.n[i] <- mean(sample007_chr2.n$depth[cov.wn$start.n[i]:cov.wn$end.n[i]])
}

View(cov.wn)

head(sample007.w)

######################################
##figure out coverage between genes
genes_chr1.w$idepth <- 0
genes_chr1.w$iozdepth <- 0
genes_chr1.w$ipdepth <- 0
for (i in 1:(length(genes_chr1.w$seqid))-1){
  genes_chr1.w$idepth[i] <- sum(sample007_chr1.w$depth[which(genes_chr1.w[i,]$end<sample007_chr1.w$locus & sample007_chr1.w$locus<genes_chr1.w[(i+1),]$start)])
  genes_chr1.w$iozdepth[i] <- sum(oz.sample007_chr1.w$depth[which(genes_chr1.w[i,]$end<oz.sample007_chr1.w$locus & oz.sample007_chr1.w$locus<genes_chr1.w[(i+1),]$start)])
  gene <- as.character(genes_chr1.w$ID[i])
  gene1 <- as.character(genes_chr1.w$ID[(i+1)])
  genes_chr1.w$ipdepth[i] <- sum(sample007_chr1.n$depth[which(genes_chr1.n[genes_chr1.n$ID==gene,]$end<sample007_chr1.n$locus & sample007_chr1.n$locus<genes_chr1.n[genes_chr1.n$ID==gene1,]$start)])
}


131746
133653


length(sample007_chr1.w[which(131746<sample007_chr1.w$locus & sample007_chr1.w$locus<133653),]$depth)
oz.sample007_chr1.w[which(131746<oz.sample007_chr1.w$locus & oz.sample007_chr1.w$locus<133653),]

a <- merge(sample007_chr1.w[which(131746<sample007_chr1.w$locus & sample007_chr1.w$locus<133653),],
      oz.sample007_chr1.w[which(131746<oz.sample007_chr1.w$locus & oz.sample007_chr1.w$locus<133653),], by="locus")
sum(a$depth.x)/(133653-131746)
sum(a$depth.y)/(133653-131746)
###actually, fewer reads mapping to some areas in mine, but reads mapping to other areas as well so equalling out?

sample007_chr1.n[which((genes_chr1.n[genes_chr1.n$ID==gene,]$start+700)<sample007_chr1.n$locus & sample007_chr1.n$locus<(genes_chr1.n[genes_chr1.n$ID==gene,]$start+800)),]
i <- 14

mean(sample007_chr1.w$depth[c(2,5,100)])

genes_chr1.w[genes_chr1.w$ID=="YAR035W",]
genes_chr1.n[genes_chr1.n$ID=="YAR035W",]


oz.sample007_chr1.w[206908>oz.sample007_chr1.w$locus & oz.sample007_chr1.w$locus>204845,]
sample007_chr1.w[206908>sample007_chr1.w$locus & sample007_chr1.w$locus>204845,]
sample007_chr1.n[169518>sample007_chr1.n$locus & sample007_chr1.n$locus>167437,]


directoz.w <- read.table("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/comparecov/ozdirect.coverage",
                             header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

head(directoz.w)

directoz.w <- rename(directoz.w, c(V1="Chr", V2="locus", V3="depth")) # renames the header
#can only do 1 chr at a time
directoz_chr1.w <- subset(directoz.w, Chr=="W303.chr01")

##same if calc coverage from his
sum(directoz.w==oz.sample007.w)
length(directoz.w$Chr)


midoz.w <- read.table("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/comparecov/007_woz.coverage",
                         header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

head(midoz.w)

midoz.w <- rename(midoz.w, c(V1="Chr", V2="locus", V3="depth")) # renames the header
#can only do 1 chr at a time
midoz_chr1.w <- subset(midoz.w, Chr=="W303.chr01")

head(midoz_chr1.w)
head(oz.sample007_chr1.w)
head(sample007_chr1.w)

sum(midoz.w==oz.sample007.w)
#not the same

sum(midoz.w==sample007.w)
length(midoz.w$Chr)
#same as the one I did
#so something in step where make sam/bam


######do region E from the beginning####
library(reshape)
is.sorted = Negate(is.unsorted)
#load in E1 coverage
E1 <- read.table("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/E1/E1.coverage",
                 header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

E1 <- rename(E1, c(V1="Chr", V2="locus", V3="depth")) # renames the header

#load in E2 coverage
E2 <- read.table("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/E2/E2.coverage",
                 header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

E2 <- rename(E2, c(V1="Chr", V2="locus", V3="depth")) # renames the header

#pull out only chr of interest
#'for E - chr15
#'do by species
E1_chr15.c <- subset(E1, Chr=="W303.chr15")
E1_chr15.p <- subset(E1, Chr=="N_17.chr15")
E2_chr15.c <- subset(E2, Chr=="W303.chr15")
E2_chr15.p <- subset(E2, Chr=="N_17.chr15")

######
#load in gene files
genes.c <- read.csv("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/W303.genes.csv", head=TRUE)
genes.p <- read.csv("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/N17.genes.csv", head=TRUE)

#limit to only genes on chr interested in
genes_chr15.c <- subset(genes.c, seqid=="chr15")
genes_chr15.p <- subset(genes.p, seqid=="chr15")


#find which genes in both and in same order
c.geneind <-c()
for (i in 1:length(genes_chr15.c$ID)){
  c.geneind <- c(c.geneind, which(genes_chr15.c$ID%in%genes_chr15.p$ID[i]))
}

#in order
is.sorted(c.geneind)

#make n17 (par) equivalent
p.geneind <- c()
for (i in 1:length(c.geneind)){
  id <- as.character(genes_chr15.c$ID[c.geneind[i]])
  p.geneind <- c(p.geneind, which(genes_chr15.p$ID==id))
}
#in order
is.sorted(p.geneind)

#pull out only genes that are in both
genes.c.both <- genes_chr15.c[c.geneind,]
genes.p.both <- genes_chr15.p[p.geneind,]



#make data frame with genes and between genes regions
ID <- 1
start.c <- 1
end.c <- c()
start.p <- 1
end.p <- c()
for (i in 1:length(genes.c.both$X)){
  ID <- c(ID, as.character(genes.c.both$ID[i]), i+1)
  start.c <- c(start.c, genes.c.both$start[i], (genes.c.both$end[i]+1))
  end.c <- c(end.c, (genes.c.both$start[i]-1), genes.c.both$end[i])
  start.p <- c(start.p, genes.p.both$start[i], (genes.p.both$end[i]+1))
  end.p <- c(end.p, (genes.p.both$start[i]-1), genes.p.both$end[i])
}
end.c <- c(end.c, max(E1_chr15.c$locus))
end.p <- c(end.p, max(E1_chr15.p$locus))
cov.cp <- data.frame(ID, start.c, end.c, start.p, end.p)

###figure out coverage of each gene
cov.cp$depth.e1.c <- 0
cov.cp$depth.e1.p <- 0
cov.cp$depth.e2.c <- 0
cov.cp$depth.e2.p <- 0
for (i in 1:length(cov.cp$ID)){
  cov.cp$depth.e1.c[i] <- mean(E1_chr15.c$depth[cov.cp$start.c[i]:cov.cp$end.c[i]])
  cov.cp$depth.e1.p[i] <- mean(E1_chr15.p$depth[cov.cp$start.p[i]:cov.cp$end.p[i]])
  cov.cp$depth.e2.c[i] <- mean(E2_chr15.c$depth[cov.cp$start.c[i]:cov.cp$end.c[i]])
  cov.cp$depth.e2.p[i] <- mean(E2_chr15.p$depth[cov.cp$start.p[i]:cov.cp$end.p[i]])
}


plot(cov.cp$depth.e1.c, type="l")
lines(cov.cp$depth.e1.p, col="black", lty=2)
lines(cov.cp$depth.e2.c, col="blue")
lines(cov.cp$depth.e2.p, col="blue", lty=2)
#region E
grep("YOL160W", cov.cp$ID)
#2
abline(v=2, lty=2, col="forestgreen")
grep("YOL157C", cov.cp$ID)
#10
abline(v=10, lty=2, col="forestgreen")
#region with rec selected
grep("YOL154W", cov.cp$ID)
#18
abline(v=18, lty=4, col="purple")

grep("YOL164W-A", cov.cp$ID)
#not there
grep("YOL164W-A", genes_chr15.c$ID)
#4 here
grep("YOL154W", genes_chr15.c$ID)
#17 here
#just use end somewhere for YOL164W-A
abline(v=1, lty=4, col="red")

#pull out only genes
row_odd <- seq_len(nrow(cov.cp)) %% 2
cov.cp.genes <- cov.cp[row_odd==0,]
cov.cp.genes

#plot just gene info
plot(cov.cp.genes$depth.e1.c, type="l")
lines(cov.cp.genes$depth.e1.p, col="black", lty=2)
lines(cov.cp.genes$depth.e2.c, col="blue")
lines(cov.cp.genes$depth.e2.p, col="blue", lty=2)
#gets rid of some weird peaks but not much diff otherwise

#####only find region of interest####
cov.cp.E <- cov.cp[1:18,]
cov.cp.E.genes <- cov.cp.genes[1:9,]

#all blocks
plot(cov.cp.E$depth.e1.c, type="l", ylim=c(0,200))
lines(cov.cp.E$depth.e1.p, col="black", lty=2)
lines(cov.cp.E$depth.e2.c, col="blue")
lines(cov.cp.E$depth.e2.p, col="blue", lty=2)

#only genes
plot(cov.cp.E.genes$depth.e1.c, type="l", ylim=c(-100,200), col="blue", ylab="Coverage")
lines(cov.cp.E.genes$depth.e1.p, col="blue", lty=2)
lines(cov.cp.E.genes$depth.e2.c, col="black")
lines(cov.cp.E.genes$depth.e2.p, col="black", lty=2)

#lines(cov.cp.E.genes$depth.e1.c+cov.cp.E.genes$depth.e2.c-0.5, col="red")
#lines(cov.cp.E.genes$depth.e1.p+cov.cp.E.genes$depth.e2.p-0.5, col="red")

#lines((cov.cp.E.genes$fp.e1+cov.cp.E.genes$fp.e2)*100, col="red")

cove1.diffs <- c(0)
for (i in 2:length(cov.cp.E.genes$ID)){
  cove1.diffs <- c(cove1.diffs, cov.cp.E.genes$depth.e1.p[i]-cov.cp.E.genes$depth.e1.p[i-1])
}

cove2.diffs <- c(0)
for (i in 2:length(cov.cp.E.genes$ID)){
  cove2.diffs <- c(cove2.diffs, cov.cp.E.genes$depth.e2.p[i]-cov.cp.E.genes$depth.e2.p[i-1])
}
lines(cove1.diffs-50, col="blue")
points(cove1.diffs-50, col="blue")

lines((cove2.diffs*-1)-50, col="black")
points((cove2.diffs*-1)-50, col="black")

abline(h=-50, col="red")

cove2.diffs
cove1.diffs



#####find coverage ratio#####
#find proportion paradoxus
cov.cp.E$fp.e1 <- cov.cp.E$depth.e1.p/(cov.cp.E$depth.e1.c+cov.cp.E$depth.e1.p)
cov.cp.E$fp.e2 <- cov.cp.E$depth.e2.p/(cov.cp.E$depth.e2.c+cov.cp.E$depth.e2.p)
plot(cov.cp.E$fp.e1, type="l", ylim=c(0,1))
lines(cov.cp.E$fp.e2)

cov.cp.E.genes$fp.e1 <- cov.cp.E.genes$depth.e1.p/(cov.cp.E.genes$depth.e1.c+cov.cp.E.genes$depth.e1.p)
cov.cp.E.genes$fp.e2 <- cov.cp.E.genes$depth.e2.p/(cov.cp.E.genes$depth.e2.c+cov.cp.E.genes$depth.e2.p)
plot(cov.cp.E.genes$fp.e1, type="l", ylab="proportion paradoxus", ylim=c(-1,1), col="blue")
lines(cov.cp.E.genes$fp.e2, col="black")
#looking at only genes is smoother - continue with this

e1.diffs <- c(0)
for (i in 2:length(cov.cp.E.genes$ID)){
  e1.diffs <- c(e1.diffs, cov.cp.E.genes$fp.e1[i]-cov.cp.E.genes$fp.e1[i-1])
}

e2.diffs <- c(0)
for (i in 2:length(cov.cp.E.genes$ID)){
  e2.diffs <- c(e2.diffs, cov.cp.E.genes$fp.e2[i]-cov.cp.E.genes$fp.e2[i-1])
}

#add to plot the diffs
lines(e1.diffs, col="blue")
points(e1.diffs, col="blue")
#e1 all weird
#e1 incr in par most around 6 but 5 weird
#e1 also incr in par around 7
#for e2 do *-1 so positive number if decrease in p
lines((e2.diffs*-1), col="black")
points((e2.diffs*-1), col="black")
#e2 decr in par most around 5 but weird part
#e2 faster decr in par around 7???
abline(h=0, lty=2, col="red")


####smoothing####
cov.cp.E.genes$index <- c(1:9)
lo <- loess(cov.cp.E.genes$fp.e1~cov.cp.E.genes$index)
plot(cov.cp.E.genes$index, cov.cp.E.genes$fp.e1, type="l", ylim=c(0,1))
lines(predict(lo), col='black', lwd=2)
lo2 <- loess(cov.cp.E.genes$fp.e2~cov.cp.E.genes$index)
lines(cov.cp.E.genes$index, cov.cp.E.genes$fp.e2, type="l", col="blue")
lines(predict(lo2), col='blue', lwd=2)
lines(predict(lo)+predict(lo2)-0.5, col="red")
abline(h=0.5)
#make vector to find min
fpoverall <- predict(lo)+predict(lo2)
which(fpoverall==min(fpoverall))
min(fpoverall)
#didn't work


######find min of fp1 and fp2#####
#plot fpa1 +fpa2 to find min - subtract 0.5 to plot
lines(cov.cp.E.genes$fp.e1+cov.cp.E.genes$fp.e2-0.5, col="red")
#make vector to find min
fpoverall <- cov.cp.E.genes$fp.e1+cov.cp.E.genes$fp.e2
grep(min(fpoverall), fpoverall)
#9 is the min - 1
#doesn't really work

cov.cp.E.genes$totaldepth.e1 <- cov.cp.E.genes$depth.e1.c + cov.cp.E.genes$depth.e1.p
cov.cp.E.genes$totaldepth.e2 <- cov.cp.E.genes$depth.e2.c + cov.cp.E.genes$depth.e2.p

mean(cov.cp.E.genes$totaldepth.e1)
mean(cov.cp.E.genes$totaldepth.e2)



######do region F####
#load in F1 coverage
F1 <- read.table("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/F1/F1.coverage",
                 header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

F1 <- rename(F1, c(V1="Chr", V2="locus", V3="depth")) # renames the header

#load in F2 coverage
F2 <- read.table("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/F2/F2.coverage",
                 header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

F2 <- rename(F2, c(V1="Chr", V2="locus", V3="depth")) # renames the header

#pull out only chr of interest
#'for F - chr15
#'do by species
F1_chr15.c <- subset(F1, Chr=="W303.chr15")
F1_chr15.p <- subset(F1, Chr=="N_17.chr15")
F2_chr15.c <- subset(F2, Chr=="W303.chr15")
F2_chr15.p <- subset(F2, Chr=="N_17.chr15")

F1_test.p <- F1_chr15.p[-grep("N", Nchr15),]
F2_test.p <- F2_chr15.p[-grep("N", Nchr15),]
F1_test.c <- F1_chr15.c[-grep("N", Wchr15),]
F2_test.c <- F2_chr15.c[-grep("N", Wchr15),]

ID start.c end.c start.p end.p depth.f1.c   depth.f1.p   depth.f2.c depth.f2.p        fp.f1     fp.f2
42 YOL141W   58476 60563   46776 48863   1737.278   0.02394636   0.01532567   1583.472 1.378365e-05 0.9999903


######
#load in gene files
#genes.c <- read.csv("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/W303.genes.csv", head=TRUE)
#genes.p <- read.csv("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/N17.genes.csv", head=TRUE)

#limit to only genes on chr interested in
genes_chr15.c <- subset(genes.c, seqid=="chr15")
genes_chr15.p <- subset(genes.p, seqid=="chr15")


#find which genes in both and in same order
c.geneind <-c()
for (i in 1:length(genes_chr15.c$ID)){
  c.geneind <- c(c.geneind, which(genes_chr15.c$ID%in%genes_chr15.p$ID[i]))
}

#in order
is.sorted(c.geneind)

#make n17 (par) equivalent
p.geneind <- c()
for (i in 1:length(c.geneind)){
  id <- as.character(genes_chr15.c$ID[c.geneind[i]])
  p.geneind <- c(p.geneind, which(genes_chr15.p$ID==id))
}
#in order
is.sorted(p.geneind)

#pull out only genes that are in both
genes.c.both <- genes_chr15.c[c.geneind,]
genes.p.both <- genes_chr15.p[p.geneind,]

#exclude YOL090W - something wrong
genes.c.both <- genes.c.both[-which(genes.c.both$ID=="YOL090W"),]
genes.p.both <- genes.p.both[-which(genes.p.both$ID=="YOL090W"),]



#make data frame with genes and between genes regions
ID <- 1
start.c <- 1
end.c <- c()
start.p <- 1
end.p <- c()
for (i in 1:length(genes.c.both$X)){
  ID <- c(ID, as.character(genes.c.both$ID[i]), i+1)
  start.c <- c(start.c, genes.c.both$start[i], (genes.c.both$end[i]+1))
  end.c <- c(end.c, (genes.c.both$start[i]-1), genes.c.both$end[i])
  start.p <- c(start.p, genes.p.both$start[i], (genes.p.both$end[i]+1))
  end.p <- c(end.p, (genes.p.both$start[i]-1), genes.p.both$end[i])
}
end.c <- c(end.c, max(F1_chr15.c$locus))
end.p <- c(end.p, max(F1_chr15.p$locus))
cov.cp <- data.frame(ID, start.c, end.c, start.p, end.p)

###figure out coverage of each gene
cov.cp$depth.f1.c <- 0
cov.cp$depth.f1.p <- 0
cov.cp$depth.f2.c <- 0
cov.cp$depth.f2.p <- 0
for (i in 1:length(cov.cp$ID)){
  cov.cp$depth.f1.c[i] <- median(F1_chr15.c$depth[cov.cp$start.c[i]:cov.cp$end.c[i]])
  cov.cp$depth.f1.p[i] <- median(F1_chr15.p$depth[cov.cp$start.p[i]:cov.cp$end.p[i]])
  cov.cp$depth.f2.c[i] <- median(F2_chr15.c$depth[cov.cp$start.c[i]:cov.cp$end.c[i]])
  cov.cp$depth.f2.p[i] <- median(F2_chr15.p$depth[cov.cp$start.p[i]:cov.cp$end.p[i]])
}


plot(cov.cp$depth.f1.c, type="l")
lines(cov.cp$depth.f1.p, col="black", lty=2)
lines(cov.cp$depth.f2.c, col="blue")
lines(cov.cp$depth.f2.p, col="blue", lty=2)
#region F
grep("YOL138C", cov.cp$ID)
#48
abline(v=48, lty=2, col="forestgreen")
grep("YOL024W", cov.cp$ID)
#280
abline(v=280, lty=2, col="forestgreen")
#region with rec selected
grep("YOL020W", cov.cp$ID)
#288
abline(v=288, lty=4, col="purple")

grep("YOL141W", cov.cp$ID)
#42
abline(v=42, lty=4, col="red")

#pull out only genes
row_odd <- seq_len(nrow(cov.cp)) %% 2
cov.cp.genes <- cov.cp[row_odd==0,]
cov.cp.genes

#plot just gene info
plot(cov.cp.genes$depth.f1.c, type="l")
lines(cov.cp.genes$depth.f1.p, col="black", lty=2)
lines(cov.cp.genes$depth.f2.c, col="blue")
lines(cov.cp.genes$depth.f2.p, col="blue", lty=2)
#gets rid of some weird peaks but not much diff otherwise

#####only find region of interest####
cov.cp.F <- cov.cp[42:288,]
cov.cp.F.genes <- cov.cp.genes[21:144,]

#all blocks
plot(cov.cp.F$depth.f1.c, type="l")
lines(cov.cp.F$depth.f1.p, col="black", lty=2)
lines(cov.cp.F$depth.f2.c, col="blue")
lines(cov.cp.F$depth.f2.p, col="blue", lty=2)

#only genes
plot(cov.cp.F.genes$depth.f1.c, type="l", col="blue", ylab="Coverage")
lines(cov.cp.F.genes$depth.f1.p, col="blue", lty=2)
lines(cov.cp.F.genes$depth.f2.c, col="black")
lines(cov.cp.F.genes$depth.f2.p, col="black", lty=2)

#lines(cov.cp.E.genes$depth.e1.c+cov.cp.E.genes$depth.e2.c-0.5, col="red")
#lines(cov.cp.E.genes$depth.e1.p+cov.cp.E.genes$depth.e2.p-0.5, col="red")

#lines((cov.cp.E.genes$fp.e1+cov.cp.E.genes$fp.e2)*100, col="red")

#covf1.diffs <- c(0)
#for (i in 2:length(cov.cp.E.genes$ID)){
#  cove1.diffs <- c(cove1.diffs, cov.cp.E.genes$depth.e1.p[i]-cov.cp.E.genes$depth.e1.p[i-1])
#}

#cove2.diffs <- c(0)
#for (i in 2:length(cov.cp.E.genes$ID)){
#  cove2.diffs <- c(cove2.diffs, cov.cp.E.genes$depth.e2.p[i]-cov.cp.E.genes$depth.e2.p[i-1])
#}
#lines(cove1.diffs-50, col="blue")
#points(cove1.diffs-50, col="blue")

#lines((cove2.diffs*-1)-50, col="black")
#points((cove2.diffs*-1)-50, col="black")

#abline(h=-50, col="red")

#cove2.diffs
#cove1.diffs



#####find coverage ratio#####
#find proportion paradoxus
cov.cp.F$fp.f1 <- cov.cp.F$depth.f1.p/(cov.cp.F$depth.f1.c+cov.cp.F$depth.f1.p)
cov.cp.F$fp.f2 <- cov.cp.F$depth.f2.p/(cov.cp.F$depth.f2.c+cov.cp.F$depth.f2.p)
plot(cov.cp.F$fp.f1, type="l", ylim=c(0,1))
lines(cov.cp.F$fp.f2)

cov.cp.F.genes$fp.f1 <- cov.cp.F.genes$depth.f1.p/(cov.cp.F.genes$depth.f1.c+cov.cp.F.genes$depth.f1.p)
cov.cp.F.genes$fp.f2 <- cov.cp.F.genes$depth.f2.p/(cov.cp.F.genes$depth.f2.c+cov.cp.F.genes$depth.f2.p)
plot(cov.cp.F.genes$fp.f1, type="l", ylab="proportion paradoxus", col="blue")
lines(cov.cp.F.genes$fp.f2, col="black")
#looking at only genes is smoother - continue with this

f1.diffs <- c(0)
for (i in 2:length(cov.cp.F.genes$ID)){
  f1.diffs <- c(f1.diffs, cov.cp.F.genes$fp.f1[i]-cov.cp.F.genes$fp.f1[i-1])
}

f2.diffs <- c(0)
for (i in 2:length(cov.cp.F.genes$ID)){
  f2.diffs <- c(f2.diffs, cov.cp.F.genes$fp.f2[i]-cov.cp.F.genes$fp.f2[i-1])
}

#add to plot the diffs
lines(f1.diffs, col="blue")
points(f1.diffs, col="blue")
#e1 all weird
#e1 incr in par most around 6 but 5 weird
#e1 also incr in par around 7
#for e2 do *-1 so positive number if decrease in p
lines((f2.diffs*-1), col="black")
points((f2.diffs*-1), col="black")
#e2 decr in par most around 5 but weird part
#e2 faster decr in par around 7???
abline(h=0, lty=2, col="red")



####smoothing####
cov.cp.F.genes$index <- c(1:124)
lo <- loess(cov.cp.F.genes$fp.f1~cov.cp.F.genes$index)
plot(cov.cp.F.genes$index, cov.cp.F.genes$fp.f1, type="l", ylim=c(0,1))
lines(predict(lo), col='black', lwd=2)
lo2 <- loess(cov.cp.F.genes$fp.f2~cov.cp.F.genes$index)
lines(cov.cp.F.genes$index, cov.cp.F.genes$fp.f2, type="l", col="blue")
lines(predict(lo2), col='blue', lwd=2)
lines(predict(lo)+predict(lo2)-0.5, col="red")
abline(h=0.5)
#make vector to find min
fpoverall <- predict(lo)+predict(lo2)
which(fpoverall==min(fpoverall))
min(fpoverall)
#37 - 0.9496481

#probs something very wrong wiht 52, get rid of - YOL090W

######find min of fp1 and fp2#####
#plot fpa1 +fpa2 to find min - subtract 0.5 to plot
lines(cov.cp.F.genes$fp.f1+cov.cp.F.genes$fp.f2-0.5, col="red")
#make vector to find min
fpoverall <- cov.cp.F.genes$fp.f1+cov.cp.F.genes$fp.f2
grep(min(fpoverall), fpoverall)
#52
#should exclude 52
grep(min(fpoverall[-52]), fpoverall[-52])
#91
#probs exclude that too
nfpoverall <- fpoverall[-52]
nfpoverall <- nfpoverall[-91]
grep(min(nfpoverall), nfpoverall)
min(nfpoverall)
#24 - 0.9015046


cov.cp.F.genes$totaldepth.f1 <- cov.cp.F.genes$depth.f1.c + cov.cp.F.genes$depth.f1.p
cov.cp.F.genes$totaldepth.f2 <- cov.cp.F.genes$depth.f2.c + cov.cp.F.genes$depth.f2.p

mean(cov.cp.F.genes$totaldepth.f1)
mean(cov.cp.F.genes$totaldepth.f2)



######do region E from the beginning removing N's####
##removing N's doesn't really do anything in this
library(reshape)
is.sorted = Negate(is.unsorted)
library(readr)
#load in N17 chr15
N17chr15 <- read_file("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/plan/sequences/n17genome.chr15.fa")
#make a vector
Nchr15 <- strsplit(N17chr15,"")[[1]]
#get rid of line breaks
Nchr15 <- Nchr15[-grep("\n",Nchr15)]
#get rid of header thing
Nchr15 <- Nchr15[-(1:11)]
#find out which are N
grep("N", Nchr15)

#load in W303 chr15
W303chr15 <- read_file("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/plan/sequences/w303genome.chr15.fa")
#make a vector
Wchr15 <- strsplit(W303chr15,"")[[1]]
#get rid of line breaks
Wchr15 <- Wchr15[-grep("\n",Wchr15)]
#get rid of header thing
Wchr15 <- Wchr15[-(1:11)]
#find out which are N
grep("N", Wchr15)




#load in E1 coverage
E1 <- read.table("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/E1/E1.coverage",
                 header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

E1 <- rename(E1, c(V1="Chr", V2="locus", V3="depth")) # renames the header

#load in E2 coverage
E2 <- read.table("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/E2/E2.coverage",
                 header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

E2 <- rename(E2, c(V1="Chr", V2="locus", V3="depth")) # renames the header

#pull out only chr of interest
#'for E - chr15
#'do by species
E1_chr15.c <- subset(E1, Chr=="W303.chr15")
E1_chr15.p <- subset(E1, Chr=="N_17.chr15")
E2_chr15.c <- subset(E2, Chr=="W303.chr15")
E2_chr15.p <- subset(E2, Chr=="N_17.chr15")

E1_test.p <- E1_chr15.p[-grep("N", Nchr15),]
E2_test.p <- E2_chr15.p[-grep("N", Nchr15),]
E1_test.c <- E1_chr15.c[-grep("N", Wchr15),]
E2_test.c <- E2_chr15.c[-grep("N", Wchr15),]

######
#load in gene files
genes.c <- read.csv("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/W303.genes.csv", head=TRUE)
genes.p <- read.csv("/Users/jasmineono/Documents/GreigLabPDoc/BDM_Mapping/Chr_Replacement/bioinformatics/N17.genes.csv", head=TRUE)

#limit to only genes on chr interested in
genes_chr15.c <- subset(genes.c, seqid=="chr15")
genes_chr15.p <- subset(genes.p, seqid=="chr15")


#find which genes in both and in same order
c.geneind <-c()
for (i in 1:length(genes_chr15.c$ID)){
  c.geneind <- c(c.geneind, which(genes_chr15.c$ID%in%genes_chr15.p$ID[i]))
}

#in order
is.sorted(c.geneind)

#make n17 (par) equivalent
p.geneind <- c()
for (i in 1:length(c.geneind)){
  id <- as.character(genes_chr15.c$ID[c.geneind[i]])
  p.geneind <- c(p.geneind, which(genes_chr15.p$ID==id))
}
#in order
is.sorted(p.geneind)

#pull out only genes that are in both
genes.c.both <- genes_chr15.c[c.geneind,]
genes.p.both <- genes_chr15.p[p.geneind,]



#make data frame with genes and between genes regions
ID <- 1
start.c <- 1
end.c <- c()
start.p <- 1
end.p <- c()
for (i in 1:length(genes.c.both$X)){
  ID <- c(ID, as.character(genes.c.both$ID[i]), i+1)
  start.c <- c(start.c, genes.c.both$start[i], (genes.c.both$end[i]+1))
  end.c <- c(end.c, (genes.c.both$start[i]-1), genes.c.both$end[i])
  start.p <- c(start.p, genes.p.both$start[i], (genes.p.both$end[i]+1))
  end.p <- c(end.p, (genes.p.both$start[i]-1), genes.p.both$end[i])
}
end.c <- c(end.c, max(E1_chr15.c$locus))
end.p <- c(end.p, max(E1_chr15.p$locus))
cov.cp <- data.frame(ID, start.c, end.c, start.p, end.p)

###figure out coverage of each gene
cov.cp$depth.e1.c <- 0
cov.cp$depth.e1.p <- 0
cov.cp$depth.e2.c <- 0
cov.cp$depth.e2.p <- 0
for (i in 1:length(cov.cp$ID)){
  cov.cp$depth.e1.c[i] <- mean(E1_test.c$depth[cov.cp$end.c[i]>=E1_test.c$locus&E1_test.c$locus>=cov.cp$start.c[i]])
  cov.cp$depth.e1.p[i] <- mean(E1_test.p$depth[cov.cp$end.p[i]>=E1_test.p$locus&E1_test.p$locus>=cov.cp$start.p[i]])
  cov.cp$depth.e2.c[i] <- mean(E2_test.c$depth[cov.cp$end.c[i]>=E2_test.c$locus&E2_test.c$locus>=cov.cp$start.c[i]])
  cov.cp$depth.e2.p[i] <- mean(E2_test.p$depth[cov.cp$end.p[i]>=E2_test.p$locus&E2_test.p$locus>=cov.cp$start.p[i]])
}


plot(cov.cp$depth.e1.c, type="l")
lines(cov.cp$depth.e1.p, col="black", lty=2)
lines(cov.cp$depth.e2.c, col="blue")
lines(cov.cp$depth.e2.p, col="blue", lty=2)
#region E
grep("YOL160W", cov.cp$ID)
#2
abline(v=2, lty=2, col="forestgreen")
grep("YOL157C", cov.cp$ID)
#10
abline(v=10, lty=2, col="forestgreen")
#region with rec selected
grep("YOL154W", cov.cp$ID)
#18
abline(v=18, lty=4, col="purple")

grep("YOL164W-A", cov.cp$ID)
#not there
grep("YOL164W-A", genes_chr15.c$ID)
#4 here
grep("YOL154W", genes_chr15.c$ID)
#17 here
#just use end somewhere for YOL164W-A
abline(v=1, lty=4, col="red")

#pull out only genes
row_odd <- seq_len(nrow(cov.cp)) %% 2
cov.cp.genes <- cov.cp[row_odd==0,]
cov.cp.genes

#plot just gene info
plot(cov.cp.genes$depth.e1.c, type="l")
lines(cov.cp.genes$depth.e1.p, col="black", lty=2)
lines(cov.cp.genes$depth.e2.c, col="blue")
lines(cov.cp.genes$depth.e2.p, col="blue", lty=2)
#gets rid of some weird peaks but not much diff otherwise

#####only find region of interest####
cov.cp.E <- cov.cp[1:18,]
cov.cp.E.genes <- cov.cp.genes[1:9,]

#all blocks
plot(cov.cp.E$depth.e1.c, type="l", ylim=c(0,200))
lines(cov.cp.E$depth.e1.p, col="black", lty=2)
lines(cov.cp.E$depth.e2.c, col="blue")
lines(cov.cp.E$depth.e2.p, col="blue", lty=2)

#only genes
plot(cov.cp.E.genes$depth.e1.c, type="l", ylim=c(-100,200), col="blue", ylab="Coverage")
lines(cov.cp.E.genes$depth.e1.p, col="blue", lty=2)
lines(cov.cp.E.genes$depth.e2.c, col="black")
lines(cov.cp.E.genes$depth.e2.p, col="black", lty=2)

#lines(cov.cp.E.genes$depth.e1.c+cov.cp.E.genes$depth.e2.c-0.5, col="red")
#lines(cov.cp.E.genes$depth.e1.p+cov.cp.E.genes$depth.e2.p-0.5, col="red")

#lines((cov.cp.E.genes$fp.e1+cov.cp.E.genes$fp.e2)*100, col="red")

cove1.diffs <- c(0)
for (i in 2:length(cov.cp.E.genes$ID)){
  cove1.diffs <- c(cove1.diffs, cov.cp.E.genes$depth.e1.p[i]-cov.cp.E.genes$depth.e1.p[i-1])
}

cove2.diffs <- c(0)
for (i in 2:length(cov.cp.E.genes$ID)){
  cove2.diffs <- c(cove2.diffs, cov.cp.E.genes$depth.e2.p[i]-cov.cp.E.genes$depth.e2.p[i-1])
}
lines(cove1.diffs-50, col="blue")
points(cove1.diffs-50, col="blue")

lines((cove2.diffs*-1)-50, col="black")
points((cove2.diffs*-1)-50, col="black")

abline(h=-50, col="red")

cove2.diffs
cove1.diffs



#####find coverage ratio#####
#find proportion paradoxus
cov.cp.E$fp.e1 <- cov.cp.E$depth.e1.p/(cov.cp.E$depth.e1.c+cov.cp.E$depth.e1.p)
cov.cp.E$fp.e2 <- cov.cp.E$depth.e2.p/(cov.cp.E$depth.e2.c+cov.cp.E$depth.e2.p)
plot(cov.cp.E$fp.e1, type="l", ylim=c(0,1))
lines(cov.cp.E$fp.e2)

cov.cp.E.genes$fp.e1 <- cov.cp.E.genes$depth.e1.p/(cov.cp.E.genes$depth.e1.c+cov.cp.E.genes$depth.e1.p)
cov.cp.E.genes$fp.e2 <- cov.cp.E.genes$depth.e2.p/(cov.cp.E.genes$depth.e2.c+cov.cp.E.genes$depth.e2.p)
plot(cov.cp.E.genes$fp.e1, type="l", ylab="proportion paradoxus", ylim=c(-1,1), col="blue")
lines(cov.cp.E.genes$fp.e2, col="black")
#looking at only genes is smoother - continue with this

e1.diffs <- c(0)
for (i in 2:length(cov.cp.E.genes$ID)){
  e1.diffs <- c(e1.diffs, cov.cp.E.genes$fp.e1[i]-cov.cp.E.genes$fp.e1[i-1])
}

e2.diffs <- c(0)
for (i in 2:length(cov.cp.E.genes$ID)){
  e2.diffs <- c(e2.diffs, cov.cp.E.genes$fp.e2[i]-cov.cp.E.genes$fp.e2[i-1])
}

#add to plot the diffs
lines(e1.diffs, col="blue")
points(e1.diffs, col="blue")
#e1 all weird
#e1 incr in par most around 6 but 5 weird
#e1 also incr in par around 7
#for e2 do *-1 so positive number if decrease in p
lines((e2.diffs*-1), col="black")
points((e2.diffs*-1), col="black")
#e2 decr in par most around 5 but weird part
#e2 faster decr in par around 7???
abline(h=0, lty=2, col="red")



######find min of fp1 and fp2#####
#plot fpa1 +fpa2 to find min - subtract 0.5 to plot
lines(cov.cp.E.genes$fp.e1+cov.cp.E.genes$fp.e2-0.5, col="red")
#make vector to find min
fpoverall <- cov.cp.E.genes$fp.e1+cov.cp.E.genes$fp.e2
grep(min(fpoverall), fpoverall)
#9 is the min - 1
#doesn't really work

cov.cp.E.genes$totaldepth.e1 <- cov.cp.E.genes$depth.e1.c + cov.cp.E.genes$depth.e1.p
cov.cp.E.genes$totaldepth.e2 <- cov.cp.E.genes$depth.e2.c + cov.cp.E.genes$depth.e2.p

mean(cov.cp.E.genes$totaldepth.e1)
mean(cov.cp.E.genes$totaldepth.e2)


#################old stuff#############

ctl = Xpar
Xparavg <- mean(cov9.cp.A.genes$depth.Xpar.c)
cov9.cp.A.genes$corrXpar <- Xparavg/cov9.cp.A.genes$depth.Xpar.c
cov9.cp.A.genes$correctedXVpar <- cov9.cp.A.genes$corrXpar*cov9.cp.A.genes$depth.XVpar.c
plot(cov9.cp.A.genes$depth.Xpar.c, type="l", ylim=c(0,50), lty=2)
lines(cov9.cp.A.genes$depth.XVpar.c, col="blue")
lines(cov9.cp.A.genes$correctedXVpar, col="red")

plot(cov9.cp.A.genes$depth.Xpar.c~cov9.cp.A.genes$depth.XVpar.c)
plot(cov9.cp.A.genes$depth.Xpar.c~cov9.cp.A.genes$correctedXVpar)

#####calculate even coverage - region A#####
#don't need to worry about genes mapping to other genome for this region
#the problem seems to be in reliable coverage across the region
plot(cov9.cp.A.genes$depth.cer.c/cov9.cp.A.genes$depth.cer.c)
lines(cov9.cp.A.genes$depth.IXpar.p/cov9.cp.A.genes$depth.IXpar.p)
lines(cov9.cp.A.genes$depth.Xpar.c/cov9.cp.A.genes$depth.cer.c, col="green")
lines(cov9.cp.A.genes$depth.Xpar.c/30, col="darkgreen")
lines(cov9.cp.A.genes$depth.XVpar.c/cov9.cp.A.genes$depth.cer.c, col="red")
lines(cov9.cp.A.genes$depth.XVpar.c/30, col="darkred")
lines(cov9.cp.A.genes$depth.XVpar.c/cov9.cp.A.genes$depth.Xpar.c, col="yellow")

#what would it be like if there was recombination between cer and IXpar?
cer.c <- cov9.cp.A.genes$depth.cer.c
cer.p <- cov9.cp.A.genes$depth.cer.p
IX.c <- cov9.cp.A.genes$depth.IXpar.c
IX.p <- cov9.cp.A.genes$depth.IXpar.p


#let's say 5% recombination event between 3rd and 4th genes
IX.c[4:18] <- IX.c[4:18]+cer.c[4:18]*0.05
cer.c[4:18]<-cer.c[4:18]-cer.c[4:18]*0.05

cer.p[4:18] <- cer.p[4:18]+IX.p[4:18]*0.05
IX.p[4:18] <- IX.p[4:18]-IX.p[4:18]*0.05

plot(cer.c/(cer.c+cer.p),type="l",ylim=c(0,1))
lines(IX.c/(IX.c+IX.p), col="blue")


IX.c[5:18] <- IX.c[5:18]+cer.c[5:18]*0.1
cer.c[5:18]<-cer.c[5:18]-cer.c[5:18]*0.1

cer.p[5:18] <- cer.p[5:18]+IX.p[5:18]*0.1
IX.p[5:18] <- IX.p[5:18]-IX.p[5:18]*0.1

IX.c[8:18] <- IX.c[8:18]+cer.c[8:18]*0.15
cer.c[8:18]<-cer.c[8:18]-cer.c[8:18]*0.15

cer.p[8:18] <- cer.p[8:18]+IX.p[8:18]*0.15
IX.p[8:18] <- IX.p[8:18]-IX.p[8:18]*0.15

IX.c[10:18] <- IX.c[10:18]+cer.c[10:18]*0.3
cer.c[10:18]<-cer.c[10:18]-cer.c[10:18]*0.3

cer.p[10:18] <- cer.p[10:18]+IX.p[10:18]*0.3
IX.p[10:18] <- IX.p[10:18]-IX.p[10:18]*0.3

plot(cer.c/(cer.c+cer.p),type="l",ylim=c(0,1))
lines(IX.c/(IX.c+IX.p), col="blue")

#gets a bit messy
