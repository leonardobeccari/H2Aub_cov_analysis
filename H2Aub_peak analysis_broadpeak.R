### load required libraries
require(rtracklayer)
library(chipenrich)
library(GenomicRanges)
library(IRanges)
library(AnnotationHub)
library(eulerr)
library(UpSetR)


### Select working directory and destination folder
setwd(choose.dir())
path_res=paste0(choose.dir(), "/")

### Import MACS2 broadpeak files 
extraCols_broadPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric")

infiles <- list.files(path = choose.dir(), pattern = ".bed", )
infiles=grep(pattern = "MACS_H2A", infiles, value = TRUE, invert = FALSE, fixed = TRUE)
infiles
sample.names=gsub(".bed", "", infiles)
sample.names=gsub("MACS_", "", sample.names)
sample.names

Alldata=GRangesList()

for(i in 1:length(infiles)) {                              # Head of for-loop
  Alldata[[length(Alldata)+1]]= import(infiles[i], format = "BED", extraCols = extraCols_broadPeak)
}
names(Alldata)=sample.names

# Keep only autosomes and merge intervals separated by less than 250bp 

for(i in 1:length(Alldata)){
  Alldata[[i]]=keepStandardChromosomes(Alldata[[i]], species = "Mus_musculus" , pruning.mode = "coarse")
  Alldata[[i]]=dropSeqlevels(Alldata[[i]], "chrM" , pruning.mode = "coarse")
  
  Alldata[[i]]=reduce(Alldata[[i]], min.gapwidth=250L)
}

# Check summary of all datasets
Nb.of.regions=vector()
for(i in 1:length(Alldata)){
 x=length(Alldata[[i]])
 Nb.of.regions=c(Nb.of.regions, x)
}
names(Nb.of.regions)=sample.names
print(Nb.of.regions)

### Extract gene promoter coordinates 
ah=AnnotationHub()
grs <- query(ah, "GRanges")
mm10=query(grs, "mm10")
refseq=query(mm10, "ucsc")
genes=subset(refseq, title=="Ensembl Genes")[[1]]
head(genes)

prom=promoters(genes,upstream= 5000, downstream = 2000)
prom=keepStandardChromosomes(prom, species = "Mus_musculus" , pruning.mode = "coarse")
prom=dropSeqlevels(prom, "chrM" , pruning.mode = "coarse")

### Filter intervals that overlap in at least two replicates

## WT saline
int_WTsal_r1_r4=intersect(Alldata$H2Aub_WT_sal_r1, Alldata$H2Aub_WT_sal_r4)
int_WTsal_r1_r7=intersect(Alldata$H2Aub_WT_sal_r1, Alldata$H2Aub_WT_sal_r7)
int_WTsal_r4_r7=intersect(Alldata$H2Aub_WT_sal_r4, Alldata$H2Aub_WT_sal_r7)

all_H2Aub_WT_sal=Reduce(union, list(int_WTsal_r1_r4, int_WTsal_r1_r7, int_WTsal_r4_r7 ))
H2Aub_WT_sal=subsetByOverlaps(prom, all_H2Aub_WT_sal, minoverlap = 250)

## WT cocaine
int_WTcoc_r1_r4=intersect(Alldata$H2Aub_WT_coc_r1, Alldata$H2Aub_WT_coc_r4)
int_WTcoc_r1_r7=intersect(Alldata$H2Aub_WT_coc_r1, Alldata$H2Aub_WT_coc_r7)
int_WTcoc_r4_r7=intersect(Alldata$H2Aub_WT_coc_r4, Alldata$H2Aub_WT_coc_r7)

all_H2Aub_WT_coc=Reduce(union, list(int_WTcoc_r1_r4, int_WTcoc_r1_r7, int_WTcoc_r4_r7 ))
H2Aub_WT_coc=subsetByOverlaps(prom, all_H2Aub_WT_coc, minoverlap = 250)

## Maged1 cKO saline
int_cKOsal_r1_r2=intersect(Alldata$H2Aub_KO_sal_rep1, Alldata$H2Aub_KO_sal_rep2)
int_cKOsal_r1_r3=intersect(Alldata$H2Aub_KO_sal_rep1, Alldata$H2Aub_KO_sal_rep3)
int_cKOsal_r2_r3=intersect(Alldata$H2Aub_KO_sal_rep2, Alldata$H2Aub_KO_sal_rep3)

all_H2Aub_KO_sal=Reduce(union, list(int_cKOsal_r1_r2, int_cKOsal_r1_r3, int_cKOsal_r2_r3 ))
H2Aub_KO_sal=subsetByOverlaps(prom, all_H2Aub_KO_sal, minoverlap = 250)

## Maged1 cKO cocaine
all_H2Aub_KO_coc=intersect(Alldata$H2Aub_KO_coc_rep1, Alldata$H2Aub_KO_coc_rep2)
H2Aub_KO_coc=subsetByOverlaps(prom, all_H2Aub_KO_coc, minoverlap = 250)

# Create a summary of number of elements
summary_elements=cbind(c(length(Alldata$H2Aub_WT_sal_r1), length(Alldata$H2Aub_WT_sal_r4), length(Alldata$H2Aub_WT_sal_r7), length(all_H2Aub_WT_sal), length(H2Aub_WT_sal)),
                       c(length(Alldata$H2Aub_WT_coc_r1), length(Alldata$H2Aub_WT_coc_r4), length(Alldata$H2Aub_WT_coc_r7), length(all_H2Aub_WT_coc), length(H2Aub_WT_coc)), 
                       c(length(Alldata$H2Aub_KO_sal_rep1), length(Alldata$H2Aub_KO_sal_rep2), length(Alldata$H2Aub_KO_sal_rep3),length(all_H2Aub_KO_sal), length(H2Aub_KO_sal)),
                       c(length(Alldata$H2Aub_KO_coc_rep1), length(Alldata$H2Aub_KO_coc_rep2), "na",length(all_H2Aub_KO_coc), length(H2Aub_KO_coc)))
colnames(summary_elements)=c("H2Aub_WT_sal","H2Aub_WT_coc", "H2Aub_KO_sal","H2Aub_KO_coc")
rownames(summary_elements)=c("rep1", "rep2","rep3", "All commons", "Commons in prom")
print(summary_elements)
write.table(summary_elements, paste0(path_res, "summary_H2A_elements.txt"), row.names = TRUE, col.names = TRUE, sep="\t"  )

### Plot distance of H2Aub+ elements from gene TSSs
pdf(paste0(path_res, "dist_tss_h2a_wt_sal.pdf"))
plot_dist_to_tss(peaks = as.data.frame(all_H2Aub_WT_sal), genome = 'mm10')
dev.off()
pdf(paste0(path_res, "dist_tss_h2a_wt_coc.pdf"))
plot_dist_to_tss(peaks = as.data.frame(all_H2Aub_WT_coc), genome = 'mm10')
dev.off()
pdf(paste0(path_res, "dist_tss_h2a_ko_sal.pdf"))
plot_dist_to_tss(peaks = as.data.frame(all_H2Aub_KO_sal), genome = 'mm10')
dev.off()
pdf(paste0(path_res, "dist_tss_h2a_ko_coc.pdf"))
plot_dist_to_tss(peaks = as.data.frame(all_H2Aub_KO_coc), genome = 'mm10')
dev.off()

ASSPEAK_H2AWTsal= chipenrich(as.data.frame(all_H2Aub_WT_sal), genome = 'mm10', 
             locusdef = "nearest_tss", qc_plots = FALSE, randomization = 'complete',
             out_name = NULL, n_cores = 1)

ASSPEAK_H2AWTcoc= chipenrich(as.data.frame(all_H2Aub_WT_coc), genome = 'mm10', 
                             locusdef = "nearest_tss", qc_plots = FALSE, randomization = 'complete',
                             out_name = NULL, n_cores = 1)

ASSPEAK_H2AKOsal= chipenrich(as.data.frame(all_H2Aub_KO_sal), genome = 'mm10', 
                             locusdef = "nearest_tss", qc_plots = FALSE, randomization = 'complete',
                             out_name = NULL, n_cores = 1)

ASSPEAK_H2AKOcoc= chipenrich(as.data.frame(all_H2Aub_KO_coc), genome = 'mm10', 
                             locusdef = "nearest_tss", qc_plots = FALSE, randomization = 'complete',
                             out_name = NULL, n_cores = 1)

TSSdist_ASSPEAK_H2AWTsal=as.vector(c(
        length(which(abs(ASSPEAK_H2AWTsal$peaks$dist_to_tss) < 5000)),
        length(which(abs(ASSPEAK_H2AWTsal$peaks$dist_to_tss) > 5000 & abs(ASSPEAK_H2AWTsal$peaks$dist_to_tss)< 10000)),
        length(which(abs(ASSPEAK_H2AWTsal$peaks$dist_to_tss) > 10000 & abs(ASSPEAK_H2AWTsal$peaks$dist_to_tss)< 50000)),
        length(which(abs(ASSPEAK_H2AWTsal$peaks$dist_to_tss) > 50000 & abs(ASSPEAK_H2AWTsal$peaks$dist_to_tss)< 100000)),
        length(which(abs(ASSPEAK_H2AWTsal$peaks$dist_to_tss) > 100000))
        ))

TSSdist_ASSPEAK_H2AWTcoc=as.vector(c(
  length(which(abs(ASSPEAK_H2AWTcoc$peaks$dist_to_tss) < 5000)),
  length(which(abs(ASSPEAK_H2AWTcoc$peaks$dist_to_tss) > 5000 & abs(ASSPEAK_H2AWTcoc$peaks$dist_to_tss)< 10000)),
  length(which(abs(ASSPEAK_H2AWTcoc$peaks$dist_to_tss) > 10000 & abs(ASSPEAK_H2AWTcoc$peaks$dist_to_tss)< 50000)),
  length(which(abs(ASSPEAK_H2AWTcoc$peaks$dist_to_tss) > 50000 & abs(ASSPEAK_H2AWTcoc$peaks$dist_to_tss)< 100000)),
  length(which(abs(ASSPEAK_H2AWTcoc$peaks$dist_to_tss) > 100000))
))

TSSdist_ASSPEAK_H2AKOsal=as.vector(c(
  length(which(abs(ASSPEAK_H2AKOsal$peaks$dist_to_tss) < 5000)),
  length(which(abs(ASSPEAK_H2AKOsal$peaks$dist_to_tss) > 5000 & abs(ASSPEAK_H2AKOsal$peaks$dist_to_tss)< 10000)),
  length(which(abs(ASSPEAK_H2AKOsal$peaks$dist_to_tss) > 10000 & abs(ASSPEAK_H2AKOsal$peaks$dist_to_tss)< 50000)),
  length(which(abs(ASSPEAK_H2AKOsal$peaks$dist_to_tss) > 50000 & abs(ASSPEAK_H2AKOsal$peaks$dist_to_tss)< 100000)),
  length(which(abs(ASSPEAK_H2AKOsal$peaks$dist_to_tss) > 100000))
))

TSSdist_ASSPEAK_H2AKOcoc=as.vector(c(
  length(which(abs(ASSPEAK_H2AKOcoc$peaks$dist_to_tss) < 5000)),
  length(which(abs(ASSPEAK_H2AKOcoc$peaks$dist_to_tss) > 5000 & abs(ASSPEAK_H2AKOcoc$peaks$dist_to_tss)< 10000)),
  length(which(abs(ASSPEAK_H2AKOcoc$peaks$dist_to_tss) > 10000 & abs(ASSPEAK_H2AKOcoc$peaks$dist_to_tss)< 50000)),
  length(which(abs(ASSPEAK_H2AKOcoc$peaks$dist_to_tss) > 50000 & abs(ASSPEAK_H2AKOcoc$peaks$dist_to_tss)< 100000)),
  length(which(abs(ASSPEAK_H2AKOcoc$peaks$dist_to_tss) > 100000))
))


TSSDISTbas=data.frame("WTsal"=TSSdist_ASSPEAK_H2AWTsal, "WTcoc"=TSSdist_ASSPEAK_H2AWTcoc,
                      "KOsal"=TSSdist_ASSPEAK_H2AKOsal, "KOcoc"=TSSdist_ASSPEAK_H2AWTcoc)
TSSDIST_rel=TSSDISTbas
TSSDIST_rel$WTsal=TSSDISTbas$WTsal/sum(TSSDISTbas$WTsal)
TSSDIST_rel$WTcoc=TSSDISTbas$WTcoc/sum(TSSDISTbas$WTcoc)
TSSDIST_rel$KOsal=TSSDISTbas$KOsal/sum(TSSDISTbas$KOsal)
TSSDIST_rel$KOcoc=TSSDISTbas$KOcoc/sum(TSSDISTbas$KOcoc)

dist.names=c("<5kb", "5-10kb", "10-50kb", "50-100kb", ">100kb")

rownames(TSSDIST_rel)=dist.names
rownames(TSSDISTbas)=dist.names

write.table(TSSDISTbas, "absDISTtoTSS.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep="\t")

### Generate and export lists of H2Aub+ genes  

gene.info=read.table(choose.files(), header = TRUE,  sep = "\t", stringsAsFactors = FALSE) ##choose biomart mouse mm10 file mart_mm10.txt

Transcr_WT_sal=as.vector(H2Aub_WT_sal$name)
genes_WT_sal=unique(as.vector(gene.info[which(gene.info$Transcript.stable.ID %in% Transcr_WT_sal), 2]))

Transcr_WT_coc=as.vector(H2Aub_WT_coc$name)
genes_WT_coc=unique(as.vector(gene.info[which(gene.info$Transcript.stable.ID %in% Transcr_WT_coc), 2]))

Transcr_KO_sal=as.vector(H2Aub_KO_sal$name)
genes_KO_sal=unique(as.vector(gene.info[which(gene.info$Transcript.stable.ID %in% Transcr_KO_sal), 2]))

Transcr_KO_coc=as.vector(H2Aub_KO_coc$name)
genes_KO_coc=unique(as.vector(gene.info[which(gene.info$Transcript.stable.ID %in% Transcr_KO_coc), 2]))

# Create gene summary
summary_genes=cbind(c(length(genes_WT_sal), length(genes_WT_coc)),c(length(genes_KO_sal), length(genes_KO_coc)))
colnames(summary_genes)=c("WT", "KO")
rownames(summary_genes)=c("saline", "cocaine")
print(summary_genes)
write.table(summary_genes, paste0(path_res, "summary_H2Aub_genes.txt"), row.names = TRUE, col.names = TRUE, sep="\t"  )

# Export gene lists
write.table(as.data.frame(genes_WT_sal), paste0(path_res, "genelist_H2Aub_WTsal_allgenes.txt"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(as.data.frame(genes_WT_coc), paste0(path_res, "genelist_H2Aub_WTcoc_allgenes.txt"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(as.data.frame(genes_KO_sal), paste0(path_res, "genelist_H2Aub_KOsal_allgenes.txt"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(as.data.frame(genes_KO_coc), paste0(path_res, "genelist_H2Aub_KOcoc_allgenes.txt"), sep="\t", row.names=F, col.names=F, quote=F)

### Euler and Upset plots

listInput <-list(H2Asal_WT= genes_WT_sal, H2Acoc_WT= genes_WT_coc,
                 H2Asal_KO= genes_KO_sal, H2Acoc_KO= genes_KO_coc)

WT_coc.vs.sal=euler(fromList(listInput[c(1,2)]), shape="circle")
pdf(paste0(path_res, "euler_WT_coc_vs_sal.pdf"))
plot(WT_coc.vs.sal, quantities = TRUE)
dev.off()

KO_coc.vs.sal=euler(fromList(listInput[c(3,4)]), shape="circle")
pdf("euler_KO_coc_vs_sal.pdf")
plot(KO_coc.vs.sal, quantities = TRUE)
dev.off()

SAL_wt.vs.KO=euler(fromList(listInput[c(1,3)]), shape="circle")
pdf(paste0(path_res, "euler_SAL_wt_vs_ko.pdf"))
plot(SAL_wt.vs.KO, quantities = TRUE)
dev.off()

COC_wt.vs.KO=euler(fromList(listInput[c(2,4)]), shape="circle")
pdf(paste0(path_res, "euler_COC_wt_vs_ko.pdf"))
plot(COC_wt.vs.KO, quantities = TRUE)
dev.off()

pdf(paste0(path_res, "upsetR_H2Aub.pdf"))
upset(fromList(listInput), order.by = "freq")
dev.off()


# Export gene lists of dataset intersections
WT_common_sal_coc=intersect(genes_WT_sal, genes_WT_coc)
WT_saline_spec=setdiff(genes_WT_sal, WT_common_sal_coc)
WT_cocaine_spec=setdiff(genes_WT_coc, WT_common_sal_coc)

KO_common_sal_coc=intersect(genes_KO_sal, genes_KO_coc)
KO_saline_spec=setdiff(genes_KO_sal, KO_common_sal_coc)
KO_cocaine_spec=setdiff(genes_KO_coc, KO_common_sal_coc)

SAL_common_wt_ko=intersect(genes_WT_sal, genes_KO_sal)
SAL_wt_spec=setdiff(genes_WT_sal, SAL_common_wt_ko)
SAL_ko_spec=setdiff(genes_KO_sal, SAL_common_wt_ko)

COC_common_wt_ko=intersect(genes_WT_coc, genes_KO_coc)
COC_wt_spec=setdiff(genes_WT_coc, COC_common_wt_ko)
COC_ko_spec=setdiff(genes_KO_coc, COC_common_wt_ko)

H2A_coc_spec_MD=setdiff(WT_cocaine_spec, KO_cocaine_spec)
H2A_sal_spec_MD=setdiff(WT_saline_spec, KO_saline_spec)




write.table(as.data.frame(WT_cocaine_spec), paste0(path_res, "genelist_H2A_cocaine_spec_WT.txt"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(as.data.frame(WT_saline_spec), paste0(path_res, "genelist_H2A_WT_saline_spec.txt"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(as.data.frame(WT_common_sal_coc), paste0(path_res, "genelist_H2A_WT_common_sal_coc.txt"), sep="\t", row.names=F, col.names=F, quote=F)
            
write.table(as.data.frame(KO_cocaine_spec), paste0(path_res, "genelist_H2A_KO_cocaine_spec.txt"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(as.data.frame(KO_saline_spec), paste0(path_res, "genelist_H2A_KO_saline_spec.txt"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(as.data.frame(KO_common_sal_coc), paste0(path_res, "genelist_H2A_KO_common_sal_coc.txt"), sep="\t", row.names=F, col.names=F, quote=F)

write.table(as.data.frame(H2A_coc_spec_MD), paste0(path_res, "genelist_H2A_cocaine_spec_MD.txt"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(as.data.frame(H2A_sal_spec_MD), paste0(path_res, "genelist_H2A_saline_spec_MD.txt"), sep="\t", row.names=F, col.names=F, quote=F)

write.table(as.data.frame(SAL_common_wt_ko), paste0(path_res, "genelist_common_WT_KO_saline.txt"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(as.data.frame(SAL_wt_spec), paste0(path_res, "genelist_H2A_WT_spec_saline.txt"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(as.data.frame(SAL_ko_spec), paste0(path_res, "genelist_H2A_KO_spec_saline.txt"), sep="\t", row.names=F, col.names=F, quote=F)

write.table(as.data.frame(COC_common_wt_ko), paste0(path_res, "genelist_common_WT_KO_cocaine.txt"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(as.data.frame(COC_wt_spec), paste0(path_res, "genelist_H2A_WT_spec_cocaine.txt"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(as.data.frame(COC_ko_spec), paste0(path_res, "genelist_H2A_KO_spec_cocaine.txt"), sep="\t", row.names=F, col.names=F, quote=F)




