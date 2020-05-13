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
## save version information 
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

###### Quality-filter reads and create Amplicon Sequence Variant tables
# put parsed, adaptors & primers removed, unjoined (R1 and R2 separate) fastq files
# into directory for DADA2 & make sure the full path is updated in the next line:

###adjust path name
path <- "~/Desktop/PANAMA_HYPOXIA/cutadapt"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_cut.fastq.gz and SAMPLENAME_R2_cut.fastq.gz
# Samplename is everything before the first underscore
fnFs <- sort(list.files(path, pattern="_R1_cut.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_cut.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Examine quality profiles of forward and reverse reads
plotQualityProfile(fnFs[1:6])
plotQualityProfile(fnRs[1:6])

# Perform filtering and trimming
# Assign the filenames for the filtered fastq.gz files.
# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Filter the forward and reverse reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

# Learn the Error Rates, it TAKES TIME!
# Forward reads
errF <- learnErrors(filtFs, multithread=TRUE)
# Reverse reads
errR <- learnErrors(filtRs, multithread=TRUE)

# visualize the estimated error rates
plotErrors(errF, nominalQ=TRUE)

# Dereplicate the filtered fastq files
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Infer the sequence variants in each sample
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Inspecting the dada-class object returned by dada:
dadaFs[[1]]

# Merge the denoised forward and reverse reads:
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers) ## The sequences being tabled vary in length.
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Remove chimeric sequences:
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# Track reads through the pipeline
# As a final check of our progress, we’ll look at the number of reads that made it through each step in the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, "dada_read_stats.txt",sep="\t",col.names=NA)

# SAVE THIS FILE SO YOU DON'T HAVE TO REPEAT ALL OF THE ABOVE STEPS, adjust name
saveRDS(seqtab.nochim, file="~/Desktop/PANAMA_HYPOXIA/seqtab.nochim.rds")
# RELOAD THE SAVED INFO FROM HERE (if you have closed the project):
# seqtab.nochim <- readRDS("~/Desktop/PANAMA_HYPOXIA/seqtab.nochim.rds")
seqtab.nochim <- readRDS("~/Desktop/PANAMA_HYPOXIA/seqtab.nochim.rds")

# Assign taxonomy
# Make sure the appropriate database is available in the DADA2 directory
taxa <- assignTaxonomy(seqtab.nochim, "~/Desktop/PANAMA_HYPOXIA/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
# FIX the NAs in the taxa table
taxon <- as.data.frame(taxa,stringsAsFactors=FALSE)
taxon$Phylum[is.na(taxon$Phylum)] <- taxon$Kingdom[is.na(taxon$Phylum)]
taxon$Class[is.na(taxon$Class)] <- taxon$Phylum[is.na(taxon$Class)]
taxon$Order[is.na(taxon$Order)] <- taxon$Class[is.na(taxon$Order)]
taxon$Family[is.na(taxon$Family)] <- taxon$Order[is.na(taxon$Family)]
taxon$Genus[is.na(taxon$Genus)] <- taxon$Family[is.na(taxon$Genus)]
write.table(taxon,"Panama_silva_taxa_table.txt",sep="\t",col.names=NA)
write.table(seqtab.nochim, "Panama_silva_otu_table.txt",sep="\t",col.names=NA)

# Create phyloseq object from otu and taxonomy tables from dada2, along with the sample metadata.
otu <- read.table("Panama_silva_otu_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("Panama_silva_taxa_table.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("Panama_metadata.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))
ps
#24083 taxa and 57 samples

# remove chloroplasts and mitochondria and Eukaryota
get_taxa_unique(ps, "Family") #682
get_taxa_unique(ps, "Order") #412
get_taxa_unique(ps, "Kingdom") #4
ps <- subset_taxa(ps, Family !="Mitochondria")
ps <- subset_taxa(ps, Order !="Chloroplast")
ps <- subset_taxa(ps, Kingdom !="Eukaryota")
ps <- subset_taxa(ps, Kingdom !="NA")
get_taxa_unique(ps, "Family") #678
get_taxa_unique(ps, "Order") #409
get_taxa_unique(ps, "Kingdom") #2
ps
#22815 taxa and 57 samples

# filtered taxa with phyloseq, now export cleaned otu and taxa tables from phyloseq
otu = as(otu_table(ps), "matrix")
taxon = as(tax_table(ps), "matrix")
metadata = as(sample_data(ps), "matrix")
write.table(otu,"Panama_silva_nochloronomito_otu_table.txt",sep="\t",col.names=NA)
write.table(taxon,"Panama_silva_nochloronomito_taxa_table.txt",sep="\t",col.names=NA)

# remove control samples for plotting, remaining samples = 62
ps = subset_samples(ps, Site != "control") ##! is how toy exclude a sample
ps

#22815 taxa and 56 samples

# plot number of observed ASVs in coral samples
plot_richness(ps,x="Treatment",color="Species",measures=c("Observed"))

# look at data and chose filtering method for very low abundance ASVs
ntaxa(ps) #22815
ps5<-filter_taxa(ps, function(x) mean(x) >5, TRUE)
ntaxa(ps5) #879
get_taxa_unique(ps, "Genus") # 1328
get_taxa_unique(ps5, "Genus") #304

# filtered ASVs with very low abundance with phyloseq, now export otu and taxa tables from phyloseq for codaseq
otu = as(otu_table(ps5), "matrix")
taxon = as(tax_table(ps5), "matrix")
metadata = as(sample_data(ps5), "matrix")
write.table(otu,"Panama_ps5_silva_nochloronomito_otu_table.txt",sep="\t",col.names=NA)
write.table(taxon,"Panama_ps5_silva_nochloronomito_taxa_table.txt",sep="\t",col.names=NA)
write.table(metadata,"Panama_ps5_silva_metadata.txt",sep="\t",col.names=NA)

######### Perform center-log-ratio transformation on ASVs and calculate Aitchison Distance and principal components
# READ IN OTU data that has been filtered for very low abundance sequences; do not clear data here. Keep phyloseq object ps5 for anosim/permanova
otu <- read.table("Panama_ps5_silva_nochloronomito_otu_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("Panama_ps5_silva_nochloronomito_taxa_table.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("Panama_ps5_silva_metadata.txt",sep="\t",header=T,row.names=1)
genus<-as.character(taxon$Genus)

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
#pdf("PCA_species_treatment.pdf",width=8.5)
p<-ggplot(df_out,aes(x=PC1,y=PC2,fill=samples$Species,shape=samples$Treatment))
p<-p+geom_point(size=3)+
  theme(axis.title=element_text(size=14))+
  theme(axis.text=element_text(size=12))+
  theme(legend.title=element_text(size=14))+
  theme(legend.text=element_text(size=12))+
  scale_fill_manual(values=cols)+
  scale_shape_manual(values=c(21,22,24))+
  guides(fill = guide_legend(override.aes=list(shape=21)))
p + labs(x=xlab, y=ylab, fill="Species",shape="Treatment") + coord_fixed()
#dev.off()

pdf("PCA_origin_treatment.pdf",width=8.5)
p<-ggplot(df_out,aes(x=PC1,y=PC2,fill=samples$Origin,shape=samples$Treatment))
p<-p+geom_point(size=3)+
  theme(axis.title = element_text(size=14))+
  theme(axis.text=element_text(size=12))+
  theme(legend.title = element_text(size=14))+
  theme(legend.text = element_text(size=12))+
  scale_shape_manual(values=c(21,22,24))+
  guides(fill = guide_legend(override.aes=list(shape=21)))
p + labs(x=xlab, y=ylab, fill="Origin", shape="Treatment") + coord_fixed()
dev.off()

####### Use phyloseq/vegan to perform ANOSIM/PERMANOVA
# set metadata as factors for anosim
treat<-as.character(samples$Treatment)
species<-as.character(samples$Species)
origin<-as.character(samples$Origin)
logger<-as.character(samples$Logger.MiniDOT)
# set Aitchison distance as the distance measure for the anosim
dist.clr <- dist(E.clr)

# anosim between groups using Aitchison distance
ano.species <- anosim(dist.clr, species, permutations=999)
pdf("Panama_ANOSIM_Species.pdf",width=8.5)
plot(ano.species)
dev.off()

ano.origin <- anosim(dist.clr, origin, permutations=999)
pdf("Panama_ANOSIM_origin.pdf",width=8.5)
plot(ano.origin)
dev.off()

ano.logger <- anosim(dist.clr, logger, permutations=999)
pdf("Panama_ANOSIM_logger.pdf",width=8.5)
plot(ano.logger)
dev.off()

# permanova between groups using Aitchison distance
perm<-adonis(dist.clr~treat*species,as(sample_data(ps5),"data.frame"))
print(perm)

perm<-adonis(dist.clr~species*treat,as(sample_data(ps5),"data.frame"))
print(perm)

perm<-adonis(dist.clr~origin*treat,as(sample_data(ps5),"data.frame"))
print(perm)

perm<-adonis(dist.clr~treat*logger,as(sample_data(ps5),"data.frame"))
print(perm)

perm<-adonis(dist.clr~species*logger,as(sample_data(ps5),"data.frame"))
print(perm)


########## BARCHARTS
ps5
ps_ra<-transform_sample_counts(ps5, function(OTU) OTU/sum(OTU))
ps_ra_alam = subset_samples(ps_ra, Species == "Agaricia lamarcki")
ps_ra_ssid = subset_samples(ps_ra, Species == "Siderastrea siderea")
#figure out how many colors you need
get_taxa_unique(ps_ra, "Order") #108
get_taxa_unique(ps_ra, "Class") #38
#you can make n any number of colors you want; with as much difference between the colors as possible (distinct colors)
n <- 38
palette <- distinctColorPalette(n)  #you can rerun this line to get a new selection of colors
#save list of colors for palette that you like, this kind of works:
cat(capture.output(print(palette), file="palette38.txt"))
palette38<-c("#AC99E4","#BA4AB4","#4BA68D","#A240ED","#79DE65","#E7A24F","#A3A47A","#6EE295","#6AEB44","#CEF146","#D9CA79","#E2B7E9","#E0E8DB","#D580B8","#9DE9C1","#6DA0E2","#95B2AA","#6CE8ED","#E25D60","#DAEAB6","#7DAA5E","#CEEB86","#7080E4","#B4E7E3","#E840DB","#AE98B3","#E8D2DD","#BACCEC","#E8D84A","#B0735D","#53EECF","#664D75","#E5559B","#E8BD9F","#E297A4","#6D56CF","#DF87EB","#69B8D9")

p1=plot_bar(ps_ra_alam, fill="Class")+
  geom_bar(aes(fill=Class), stat="identity",position="stack")+theme(strip.text=element_text(face="bold"))+
  theme(axis.text.x=element_text(angle = 90))+scale_fill_manual(values=palette38)+
  facet_grid(.~Treatment,scales="free",space="free")+
  ggtitle("Agaricia lamarcki")+
  theme(plot.title = element_text(face="italic"))+
  theme(legend.position = "bottom")+
  theme(axis.title.x = element_blank())
p1

p2=plot_bar(ps_ra_ssid, fill="Class")+
  geom_bar(aes(fill=Class), stat="identity",position="stack")+theme(strip.text=element_text(face="bold"))+
  theme(axis.text.x=element_text(angle = 90))+scale_fill_manual(values=palette38)+
  facet_grid(.~Treatment,scales="free",space="free")+
  ggtitle("Siderastrea siderea")+
  theme(plot.title = element_text(face="italic"))+
  theme(legend.position="bottom")+
  theme(axis.title.x = element_blank())
p2

# adjust width and height until it looks right for double columns
pdf("Panama_Barcharts_Class.pdf",width=20, height=10)
plot_grid(p1,p2,labels=c("A","B"), ncol=2, nrow=1)
dev.off()

                               
################### try bubbleplots

pdf(file="bubbleplot.pdf",width=8.5)
ggplot(ps_ra,aes(Sample,Treatment))+
  geom_point(aes(size=Abundance),alpha=0.5)+
  scale_size(range = c(0, 100), name="Relative Abundance")+
  scale_y_discrete(limits = rev(levels(ps_ra$Treatment)))+
  theme_bw()+
  theme(axis.text.x=element_text(size=10,angle=90,hjust=1,vjust=0))+
  theme(axis.text.y=element_text(size=10))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  theme(legend.text = element_text(size=10))
   #+scale_colour_manual(values=myPalette1)
dev.off()


################################# ANCOM TEST OF DIFFERENTIALLY ABUNDANT FAMILIES
#ANCOM Function - compare across multiple treatments groups using a compositional appproach
#https://sites.google.com/site/siddharthamandal1985/research

###Need to run this first in order to run ANCOM on data
ancom.W = function(otu_data,var_data,
                   adjusted,repeated,
                   main.var,adj.formula,
                   repeat.var,long,rand.formula,
                   multcorr,sig){
  
  n_otu=dim(otu_data)[2]-1
  
  otu_ids=colnames(otu_data)[-1]
  
  if(repeated==F){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID",all.y=T),row.names=NULL)
    #data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var)],by="Sample.ID",all.y=T),row.names=NULL)
  }else if(repeated==T){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID"),row.names=NULL)
    # data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var,repeat.var)],by="Sample.ID"),row.names=NULL)
  }
  
  base.formula = paste0("lr ~ ",main.var)
  if(repeated==T){
    repeat.formula = paste0(base.formula," | ", repeat.var)
  }
  if(adjusted==T){
    adjusted.formula = paste0(base.formula," + ", adj.formula)
  }
  
  if( adjusted == F & repeated == F ){
    fformula  <- formula(base.formula)
  } else if( adjusted == F & repeated == T & long == T ){
    fformula  <- formula(base.formula)   
  }else if( adjusted == F & repeated == T & long == F ){
    fformula  <- formula(repeat.formula)   
  }else if( adjusted == T & repeated == F  ){
    fformula  <- formula(adjusted.formula)   
  }else if( adjusted == T & repeated == T  ){
    fformula  <- formula(adjusted.formula)   
  }else{
    stop("Problem with data. Dataset should contain OTU abundances, groups, 
         and optionally an ID for repeated measures.")
  }
  
  
  
  if( repeated==FALSE & adjusted == FALSE){
    if( length(unique(data_comp[,which(colnames(data_comp)==main.var)]))==2 ){
      tfun <- exactRankTests::wilcox.exact
    } else{
      tfun <- stats::kruskal.test
    }
  }else if( repeated==FALSE & adjusted == TRUE){
    tfun <- stats::aov
  }else if( repeated== TRUE & adjusted == FALSE & long == FALSE){
    tfun <- stats::friedman.test
  }else if( repeated== TRUE & adjusted == FALSE & long == TRUE){
    tfun <- nlme::lme
  }else if( repeated== TRUE & adjusted == TRUE){
    tfun <- nlme::lme
  }
  
  logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
  for(ii in 1:(n_otu-1)){
    for(jj in (ii+1):n_otu){
      data.pair <- data_comp[,which(colnames(data_comp)%in%otu_ids[c(ii,jj)])]
      lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
      
      lr_dat <- data.frame( lr=lr, data_comp,row.names=NULL )
      
      if(adjusted==FALSE&repeated==FALSE){  ## Wilcox, Kruskal Wallis
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==FALSE&repeated==TRUE&long==FALSE){ ## Friedman's 
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==TRUE&repeated==FALSE){ ## ANOVA
        model=tfun(formula=fformula, data = lr_dat,na.action=na.omit)   
        picker=which(gsub(" ","",row.names(summary(model)[[1]]))==main.var)  
        logratio.mat[ii,jj] <- summary(model)[[1]][["Pr(>F)"]][picker]
      }else if(repeated==TRUE&long==TRUE){ ## GEE
        model=tfun(fixed=fformula,data = lr_dat,
                   random = formula(rand.formula),
                   correlation=corAR1(),
                   na.action=na.omit)   
        picker=which(gsub(" ","",row.names(anova(model)))==main.var)
        logratio.mat[ii,jj] <- anova(model)[["p-value"]][picker]
      }
      
    }
  } 
  
  ind <- lower.tri(logratio.mat)
  logratio.mat[ind] <- t(logratio.mat)[ind]
  
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
  
  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  #########################################
  ### Code to extract surrogate p-value
  surr.pval <- apply(mc.pval,1,function(x){
    s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
    # s0=max(x[which(as.numeric(as.character(x))<alpha)])
    return(s0)
  })
  #########################################
  ### Conservative
  if(multcorr==1){
    W <- apply(b,1,function(x){
      subp <- length(which(x<sig))
    })
    ### Moderate
  } else if(multcorr==2){
    W <- apply(mc.pval,1,function(x){
      subp <- length(which(x<sig))
    })
    ### No correction
  } else if(multcorr==3){
    W <- apply(logratio.mat,1,function(x){
      subp <- length(which(x<sig))
    })
  }
  
  return(W)
  }



ANCOM.main = function(OTUdat,Vardat,
                      adjusted,repeated,
                      main.var,adj.formula,
                      repeat.var,longitudinal,
                      random.formula,
                      multcorr,sig,
                      prev.cut){
  
  p.zeroes=apply(OTUdat[,-1],2,function(x){
    s=length(which(x==0))/length(x)
  })
  
  zeroes.dist=data.frame(colnames(OTUdat)[-1],p.zeroes,row.names=NULL)
  colnames(zeroes.dist)=c("Taxon","Proportion_zero")
  
  zero.plot = ggplot(zeroes.dist, aes(x=Proportion_zero)) + 
    geom_histogram(binwidth=0.1,colour="black",fill="white") + 
    xlab("Proportion of zeroes") + ylab("Number of taxa") +
    theme_bw()
  
  #print(zero.plot)
  
  OTUdat.thinned=OTUdat
  OTUdat.thinned=OTUdat.thinned[,c(1,1+which(p.zeroes<prev.cut))]
  
  otu.names=colnames(OTUdat.thinned)[-1]
  
  W.detected   <- ancom.W(OTUdat.thinned,Vardat,
                          adjusted,repeated,
                          main.var,adj.formula,
                          repeat.var,longitudinal,random.formula,
                          multcorr,sig)
  
  W_stat       <- W.detected
  
  
  ### Bubble plot
  
  W_frame = data.frame(otu.names,W_stat,row.names=NULL)
  W_frame = W_frame[order(-W_frame$W_stat),]
  
  W_frame$detected_0.9=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.8=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.7=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.6=rep(FALSE,dim(W_frame)[1])
  
  W_frame$detected_0.9[which(W_frame$W_stat>0.9*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.8[which(W_frame$W_stat>0.8*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.7[which(W_frame$W_stat>0.7*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.6[which(W_frame$W_stat>0.6*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  
  final_results=list(W_frame,zero.plot)
  names(final_results)=c("W.taxa","PLot.zeroes")
  return(final_results)
}


#####ANCOM results
otu <- read.table("Panama_silva_nochloronomito_otu_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("Panama_silva_nochloronomito_taxa_table.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("Panama_metadata.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))
ps
#22815 taxa and 57 samples
# remove control samples for plotting, remaining samples = 56
ps = subset_samples(ps, Species != "control")
ps
#22815 taxa and 56 samples

#Separate coral species
ps_ssid <- subset_samples(ps, Species == "Siderastrea siderea")
ps_alam <- subset_samples(ps, Species == "Agaricia lamarcki")

############Siderastrea siderea
# aggregate ASVs by family
dat <- tax_glom(ps_ssid, taxrank = "Family")
#melt the data, so it's like a dataframe
datm <- psmelt(dat)
#Cast the new datatable with columns that are of interest
datc <- data.table::dcast(datm, Sample + Logger.MiniDOT + Treatment + Origin + Species ~ Family, value.var = 'Abundance', fun.aggregate = sum)
dim(datc) #dimensions of the table
#683 families
otud <- datc[,c(1,7:683)] #select the first column, and then all of the taxa columns  
colnames(otud)[1] <- "Sample.ID" #rename the first column to Sample.ID - this is to match ANCOM syntax
metadat <- sample_data(ps_ssid) #get the sample data
metadat <- as.data.frame(as.matrix(metadat)) #make into into a matrix
# at this point, my sample id numbers are the row names, not a separate column, move row names to column with dplylr
metadat <- tibble::rownames_to_column(metadat, "Sample.ID") #make sure the sample names column is called Sample.ID
names(otud) <- make.names(names(otud)) #get the names from the table for 
otu_test <- otud #rename otud to otu_test, for syntax in ANCOM
metadat <- select(metadat, c("Sample.ID","Species","Treatment","Origin")) # use select to only use treatment columns of interest
map_test <- metadat #rename map_TEst
Vardat <- map_test #specify that this for Vardat - ANCOM syntax
#### ANCOM test - not adjusted, more than 2 levels = Kruskal Wallis
comparison_test_treat=ANCOM.main(OTUdat=otu_test, #calling the OTU table
                                 Vardat=map_test, #calling the metadata
                                 adjusted=FALSE, #true if covariates are to be included for adjustment
                                 repeated=FALSE, #repeated measure
                                 main.var="Treatment", #main variable or fator
                                 adj.formula= NULL, #other factors to include
                                 repeat.var=FALSE, #repeated measure
                                 long = FALSE, #longitudinal study
                                 multcorr=2,
                                 sig=0.05, #significance level
                                 prev.cut=0.90) #OTUs with proportion of zeroes greater than prev.cut are not included in the analysis
res <- comparison_test_treat$W.taxa #taxa that significantly vary across factor level of interest
write.table(res,"SSID_ANCOM_family_KruskallWallis_Treatment.txt",sep="\t",col.names=NA)
res2 <- res[which(res$detected_0.7==TRUE),] 

sig_sites <- glue::glue_collapse(droplevels(factor(res2$otu.names)), sep = ", ") #this is to get a list of the families that are different
print(sig_sites)
#Family_XII, Marinifilaceae, Pseudomonadaceae, Parvibaculaceae, Nitrincolaceae, Marinobacteraceae, Francisellaceae, Halomonadaceae, Shewanellaceae
#Calculate relative abundance
datc_relabund <-  sweep(datc[,7:683], 1, rowSums(datc[,7:683]), '/')
datc_relnames <- cbind(datc[,1:6],datc_relabund)

#only selet the significant families
sig_dis <- select(datc_relnames, Sample, Species, Treatment, Family_XII, Marinifilaceae, Pseudomonadaceae, Parvibaculaceae, Nitrincolaceae, Marinobacteraceae, Francisellaceae, Halomonadaceae, Shewanellaceae )
sig_long <- melt(sig_dis, id.vars=c("Sample","Species","Treatment"),variable.name="Family",value.name="Proportion")
sum_sig <- Rmisc::summarySE(sig_long, measurevar = "Proportion", groupvars = c("Treatment","Family"), na.rm=TRUE)

cols<-c("Natural population"="#56B4E9","Open"="#E69F00","Tent"="#CC79A7")
sum_sig$Treatment<-factor(sum_sig$Treatment, levels=c("Natural population","Open","Tent"))
#pdf("ANCOM_Families_Treatment_covarSpecies.pdf",width=8.5)
fams <- ggplot(sum_sig, aes(x=Family, y=Proportion))+
  geom_point(size=4,aes(color=Treatment))+
  scale_colour_manual(values=cols)+
  coord_flip()+
  theme_bw()+
  theme(axis.text.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=14))+
  theme(axis.title.x=element_text(size=14))+
  theme(axis.title.y=element_text(size=14))+
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  geom_errorbar(aes(ymin=Proportion-se, ymax=Proportion+se), width=.1)+
  ggtitle("Siderastrea siderea")+
  theme(plot.title = element_text(face="italic"))+
  theme(legend.title = element_blank())+
  theme(legend.text = element_text(size=12))
fams
#dev.off()


############ Agaracia lamarcki
# aggregate ASVs by family
dat <- tax_glom(ps_alam, taxrank = "Family")
#melt the data, so it's like a dataframe
datm <- psmelt(dat)
#Cast the new datatable with columns that are of interest
datc <- data.table::dcast(datm, Sample + Logger.MiniDOT + Treatment + Origin + Species ~ Family, value.var = 'Abundance', fun.aggregate = sum)
dim(datc) #dimensions of the table
#683 families
otud <- datc[,c(1,7:683)] #select the first column, and then all of the taxa columns  
colnames(otud)[1] <- "Sample.ID" #rename the first column to Sample.ID - this is to match ANCOM syntax
metadat <- sample_data(ps_alam) #get the sample data
metadat <- as.data.frame(as.matrix(metadat)) #make into into a matrix
# at this point, my sample id numbers are the row names, not a separate column, move row names to column with dplylr
metadat <- tibble::rownames_to_column(metadat, "Sample.ID") #make sure the sample names column is called Sample.ID
names(otud) <- make.names(names(otud)) #get the names from the table for 
otu_test <- otud #rename otud to otu_test, for syntax in ANCOM
metadat <- select(metadat, c("Sample.ID","Species","Treatment","Origin")) # use select to only use treatment columns of interest
map_test <- metadat #rename map_TEst
Vardat <- map_test #specify that this for Vardat - ANCOM syntax
#### ANCOM test - not adjusted, more than 2 levels = Kruskal Wallis
comparison_test_treat=ANCOM.main(OTUdat=otu_test, #calling the OTU table
                                 Vardat=map_test, #calling the metadata
                                 adjusted=FALSE, #true if covariates are to be included for adjustment
                                 repeated=FALSE, #repeated measure
                                 main.var="Treatment", #main variable or fator
                                 adj.formula= NULL, #other factors to include
                                 repeat.var=FALSE, #repeated measure
                                 long = FALSE, #longitudinal study
                                 multcorr=2,
                                 sig=0.05, #significance level
                                 prev.cut=0.90) #OTUs with proportion of zeroes greater than prev.cut are not included in the analysis
resA <- comparison_test_treat$W.taxa #taxa that significantly vary across factor level of interest
write.table(res,"ALAM_ANCOM_family_KruskallWallis_Treatment.txt",sep="\t",col.names=NA)
resA2 <- res[which(resA$detected_0.7==TRUE),] 

sig_sites <- glue::glue_collapse(droplevels(factor(resA2$otu.names)), sep = ", ") #this is to get a list of the families that are different
print(sig_sites)
#Family_XII, Marinifilaceae, Pseudomonadaceae, Parvibaculaceae, Nitrincolaceae
datc_relabund <-  sweep(datc[,7:683], 1, rowSums(datc[,7:683]), '/')
datc_relnames <- cbind(datc[,1:6],datc_relabund)

#only selet the significant families
sig_dis <- select(datc_relnames, Sample, Species, Treatment, Family_XII, Marinifilaceae, Pseudomonadaceae, Parvibaculaceae, Nitrincolaceae )
sig_long <- melt(sig_dis, id.vars=c("Sample","Species","Treatment"),variable.name="Family",value.name="Proportion")
library(Rmisc)
sum_sig <- Rmisc::summarySE(sig_long, measurevar = "Proportion", groupvars = c("Treatment","Family"), na.rm=TRUE)

#cols<-c("lesion"="#D55E00","near"="#E69F00","far"="#999999","undiseased"="#000000")
#sum_sig$Treatment<-factor(sum_sig$Treatment, levels=c("lesion","near","far","undiseased"))
pdf("ANCOM_Families_Treatment_covarSpecies.pdf",width=8.5)
fams <- ggplot(sum_sig, aes(x=Family, y=Proportion))+
  geom_point(size=4,aes(color=Treatment))+
  #scale_colour_manual(values=cols)+
  coord_flip()+
  theme_bw()+
  theme(axis.text.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=14))+
  theme(axis.title.x=element_text(size=14))+
  theme(axis.title.y=element_text(size=14))+
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  geom_errorbar(aes(ymin=Proportion-se, ymax=Proportion+se), width=.1)+
  ggtitle("Agaricia lamarcki")+
  theme(plot.title = element_text(face="italic"))+
  theme(legend.title = element_blank())+
  theme(legend.text = element_text(size=12))
fams
dev.off()




############ BOTH SPECIES TOGETHER
# aggregate ASVs by family
dat <- tax_glom(ps, taxrank = "Family")
#melt the data, so it's like a dataframe
datm <- psmelt(dat)
#Cast the new datatable with columns that are of interest
datc <- data.table::dcast(datm, Sample + Logger.MiniDOT + Treatment + Origin + Species ~ Family, value.var = 'Abundance', fun.aggregate = sum)
dim(datc) #dimensions of the table
#683 families
otud <- datc[,c(1,7:683)] #select the first column, and then all of the taxa columns  
colnames(otud)[1] <- "Sample.ID" #rename the first column to Sample.ID - this is to match ANCOM syntax
metadat <- sample_data(ps) #get the sample data
metadat <- as.data.frame(as.matrix(metadat)) #make into into a matrix
# at this point, my sample id numbers are the row names, not a separate column, move row names to column with dplylr
metadat <- tibble::rownames_to_column(metadat, "Sample.ID") #make sure the sample names column is called Sample.ID
names(otud) <- make.names(names(otud)) #get the names from the table for 
otu_test <- otud #rename otud to otu_test, for syntax in ANCOM
metadat <- select(metadat, c("Sample.ID","Species","Treatment","Origin")) # use select to only use treatment columns of interest
map_test <- metadat #rename map_TEst
Vardat <- map_test #specify that this for Vardat - ANCOM syntax
#### ANCOM test - not adjusted, more than 2 levels = Kruskal Wallis
comparison_test_treat=ANCOM.main(OTUdat=otu_test, #calling the OTU table
                                 Vardat=map_test, #calling the metadata
                                 adjusted=FALSE, #true if covariates are to be included for adjustment
                                 repeated=FALSE, #repeated measure
                                 main.var="Treatment", #main variable or fator
                                 adj.formula= NULL, #other factors to include
                                 repeat.var=FALSE, #repeated measure
                                 long = FALSE, #longitudinal study
                                 multcorr=2,
                                 sig=0.05, #significance level
                                 prev.cut=0.90) #OTUs with proportion of zeroes greater than prev.cut are not included in the analysis
resALL <- comparison_test_treat$W.taxa #taxa that significantly vary across factor level of interest
write.table(res,"ANCOM_family_KruskallWallis_Treatment.txt",sep="\t",col.names=NA)
resALL2 <- res[which(resALL$detected_0.7==TRUE),] 

#### ANCOM test - Adjusted by Coral species, ANOVA
comparison_test_treat=ANCOM.main(OTUdat=otu_test, #calling the OTU table
                                 Vardat=map_test, #calling the metadata
                                 adjusted=TRUE, #true if covariates are to be included for adjustment
                                 repeated=FALSE, #repeated measure
                                 main.var="Treatment", #main variable or fator
                                 adj.formula= "Species", #other factors to include
                                 repeat.var=FALSE, #repeated measure
                                 long = FALSE, #longitudinal study
                                 multcorr=2,
                                 sig=0.05, #significance level
                                 prev.cut=0.90) #OTUs with proportion of zeroes greater than prev.cut are not included in the analysis
resALL3 <- comparison_test_treat$W.taxa #taxa that significantly vary across factor level of interest
write.table(resALL3,"SSID_ANCOM_family_ANOVA_Treatment_covarSpecies.txt",sep="\t",col.names=NA)
resALL4 <- res[which(resALL3$detected_0.9==TRUE),] 

sig_sites <- glue::glue_collapse(droplevels(factor(resALL4$otu.names)), sep = ", ") #this is to get a list of the families that are different
print(sig_sites)
#Family_XII, Marinifilaceae, Pseudomonadaceae, Parvibaculaceae, Nitrincolaceae, Marinobacteraceae, Francisellaceae, Halomonadaceae, Shewanellaceae, Midichloriaceae, Schekmanbacteria, Desulfovibrionaceae, Bacteroidales, Anaerolineaceae
#Calculate relative abundance
datc_relabund <-  sweep(datc[,7:683], 1, rowSums(datc[,7:683]), '/')
datc_relnames <- cbind(datc[,1:6],datc_relabund)

#only selet the significant families
sig_dis <- select(datc_relnames, Sample, Species, Treatment, Family_XII, Marinifilaceae, Pseudomonadaceae, Parvibaculaceae, Nitrincolaceae, Marinobacteraceae, Francisellaceae, Halomonadaceae, Shewanellaceae, Midichloriaceae, Schekmanbacteria, Desulfovibrionaceae, Bacteroidales, Anaerolineaceae )
sig_long <- melt(sig_dis, id.vars=c("Sample","Species","Treatment"),variable.name="Family",value.name="Proportion")
sum_sig <- Rmisc::summarySE(sig_long, measurevar = "Proportion", groupvars = c("Species","Treatment","Family"), na.rm=TRUE)

#position=position_jitter(height=.00001)

pdf("ANCOM_Families_Treatment_covarSpecies.pdf",width=11)
cols<-c("Natural population"="#56B4E9","Open"="#E69F00","Tent"="#CC79A7")
sum_sig$Treatment<-factor(sum_sig$Treatment, levels=c("Natural population","Open","Tent"))
fams <- ggplot(sum_sig, aes(x=Family, y=Proportion))+
  geom_point(size=4,aes(color=Treatment))+
  scale_colour_manual(values=cols)+
  facet_grid(.~Species)+
  coord_flip()+
  theme_bw()+
  theme(plot.margin=unit(c(1,1,1,1),"cm"))+
  theme(axis.text.x=element_text(size=10))+
  theme(axis.text.y=element_text(size=14))+
  theme(axis.title.x=element_text(size=14))+
  theme(axis.title.y=element_blank())+
  theme(strip.text.x = element_text(size=14, face="italic"))+
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  geom_errorbar(aes(ymin=Proportion-se, ymax=Proportion+se), width=.1)+
  theme(legend.title = element_blank())+
  theme(legend.text = element_text(size=12))
fams
dev.off()

##########LOG SCALE
pdf("ANCOM_Families_Treatment_covarSpecies_log.pdf",width=11)
require(scales)
cols<-c("Natural population"="#56B4E9","Open"="#E69F00","Tent"="#CC79A7")
sum_sig$Treatment<-factor(sum_sig$Treatment, levels=c("Natural population","Open","Tent"))
fams <- ggplot(sum_sig, aes(x=Family, y=Proportion))+
  geom_point(size=4,aes(color=Treatment))+
  scale_colour_manual(values=cols)+
  facet_grid(.~Species)+
  coord_flip()+
  theme_bw()+
  theme(plot.margin=unit(c(1,1,1,1),"cm"))+
  theme(axis.text.x=element_text(size=10))+
  theme(axis.text.y=element_text(size=14))+
  theme(axis.title.x=element_text(size=14))+
  theme(axis.title.y=element_blank())+
  theme(strip.text.x = element_text(size=14, face="italic"))+
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  geom_errorbar(aes(ymin=Proportion+0.001-se, ymax=Proportion+0.001+se), width=.1)+
  theme(legend.title = element_blank())+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  theme(legend.text = element_text(size=12))
fams
dev.off()



###### Bar charts with one coral species at a time, finer resolution than Class
get_taxa_unique(ps_ra_alam, "Genus") #304
# get rid of ASVs with no counts
ps_ra_alam <- prune_taxa(taxa_sums(ps_ra_alam) != 0, ps_ra_alam)
get_taxa_unique(ps_ra_alam, "Genus") #285
get_taxa_unique(ps_ra_alam, "Family") #169

n <- 169
palette169 <- distinctColorPalette(n)

pdf("Alam_Families.pdf",width=11, height=20)
p1=plot_bar(ps_ra_alam, fill="Family")+
  geom_bar(aes(fill=Family), stat="identity",position="stack")+theme(strip.text=element_text(face="bold"))+
  theme(axis.text.x=element_text(angle = 90))+scale_fill_manual(values=palette169)+
  facet_grid(.~Treatment,scales="free",space="free")+
  ggtitle("Agaricia lamarcki")+
  theme(plot.title = element_text(face="italic"))+
  theme(legend.position = "bottom")+
  theme(axis.title.x = element_blank())
p1
dev.off()

get_taxa_unique(ps_ra_ssid, "Genus") #304
# get rid of ASVs with no counts
ps_ra_ssid <- prune_taxa(taxa_sums(ps_ra_ssid) != 0, ps_ra_ssid)
get_taxa_unique(ps_ra_ssid, "Genus") #297
get_taxa_unique(ps_ra_ssid, "Family") #174
n <- 174
palette174 <- distinctColorPalette(n)

pdf("Ssid_Families.pdf",width=11, height=20)
p2=plot_bar(ps_ra_ssid, fill="Family")+
  geom_bar(aes(fill=Family), stat="identity",position="stack")+theme(strip.text=element_text(face="bold"))+
  theme(axis.text.x=element_text(angle = 90))+scale_fill_manual(values=palette174)+
  facet_grid(.~Treatment,scales="free",space="free")+
  ggtitle("Diploria labyrinthiformis")+
  theme(plot.title = element_text(face="italic"))+
  theme(legend.position="bottom")+
  theme(axis.title.x = element_blank())
p2
dev.off()


