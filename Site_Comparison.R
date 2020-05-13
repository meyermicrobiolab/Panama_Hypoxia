##site comparison with SSID 

#pull libraries
library(dada2)
library(ShortRead)
library(ggplot2)
library(phyloseq)
library(vegan)
library(knitr)
library(ALDEx2)
library(CoDaSeq)
library(zCompositions)
library(igraph)
library(car)
library(grDevices)
library(propr)
library(cowplot)
library(randomcoloR)
library(dplyr)
library(reshape2)
library(tibble)
library(exactRankTests)
library(nlme)
library(data.table)
library(Rmisc)

# READ IN OTU data that has been filtered for very low abundance sequences. 
otu <- read.table("Panama_ps5_silva_nochloronomito_otu_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("Panama_ps5_silva_nochloronomito_taxa_table.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("Panama_ps5_silva_metadata.txt",sep="\t",header=T,row.names=1)

##make a phyloseq object to be manipulated
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))
ps

##Subsetting species SSID
ps=subset_samples(ps,Species =="Siderastrea siderea")
ps=subset_samples(ps,Site!="Punta")
ps
otu = as(otu_table(ps), "matrix")
taxon = as(tax_table(ps), "matrix")
metadata = as(sample_data(ps), "matrix")
write.table(otu,"SSID_OTU.txt",sep="\t",col.names=NA)
write.table(taxon,"SSID_taxa.txt",sep="\t",col.names=NA)
write.table(metadata,"SSID_metadata.txt",sep="\t",col.names=NA)
otu <- read.table("SSID_OTU.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("SSID_taxa.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("SSID_metadata.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))
ps

##need to subset samples by site too (since only need Finca and Tierra Oscura)

# First, replace 0 values with an estimate (because normalization is taking log, can't have 0)
# Also transposing here, need samples as rows
d.czm <- cmultRepl(t(otu), method="CZM", label=0)
# Perform the center-log-ratio (CLR) transformation 
d.clr <- codaSeq.clr(d.czm)
# transpose matrix of CLR transformed data for ordination and dendrogram
E.clr <- t(d.clr)
# plot compositional PCA biplot (perform a singular value decomposition)
d.pcx <- prcomp(E.clr)
# calculate percent variance explained for the axis labels
pc1 <- round(d.pcx$sdev[1]^2/sum(d.pcx$sdev^2),2)
pc2 <- round(d.pcx$sdev[2]^2/sum(d.pcx$sdev^2),2)
xlab <- paste("PC1: ", pc1, sep="")
ylab <- paste("PC2: ", pc2, sep="")
biplot(d.pcx, cex=c(0.6,0.4), var.axes=F,scale=1, xlab=xlab, ylab=ylab, ylabs=genus)
summary(d.pcx)
str(d.pcx)
screeplot(d.pcx)
# replot PCA with ggplot2 (showing samples only)
df_out <- as.data.frame(d.pcx$x)
theme_set(theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()))
cols<-c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pdf("PCA_SSID_site.pdf",width=8.5)
p<-ggplot(df_out,aes(x=PC1,y=PC2,fill=samples$Site,shape=samples$Species))
p<-p+geom_point(size=3)+
  theme(axis.title=element_text(size=14))+
  theme(axis.text=element_text(size=12))+
  theme(legend.title=element_text(size=14))+
  theme(legend.text=element_text(size=12))+
  scale_fill_manual(values=cols)+
  scale_shape_manual(values=c(21,22,24))+
  guides(fill = guide_legend(override.aes=list(shape=21)))
p + labs(x=xlab, y=ylab, fill="Site",shape="Species") + coord_fixed()
dev.off()
####### Use phyloseq/vegan to perform ANOSIM/PERMANOVA
# set metadata as factors for anosim
treat<-as.character(samples$Treatment)
species<-as.character(samples$Species)
origin<-as.character(samples$Origin)
logger<-as.character(samples$Logger.MiniDOT)
site<-as.character(samples$Site)
# set Aitchison distance as the distance measure for the anosim
dist.clr <- dist(E.clr)
# anosim between groups using Aitchison distance, ARE THE GROUPS DIFFERENT (look at between data)
ano.treat <- anosim(dist.clr, treat, permutations=999)
pdf("Panama_SSID_ANOSIM_site.pdf",width=8.5)
plot(ano.treat)
dev.off()
#permanove punta
perm<-adonis(dist.clr~species*site,as(sample_data(ps),"data.frame"))
print(perm)
