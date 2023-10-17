### Load libraries
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(gprofiler2)
library(org.Mm.eg.db)
library(GOSemSim)
library(compEpiTools)
BiocManager::install("compEpiTools")
BiocManager::install("DOSE")

BiocManager::install("org.Mm.eg.db")

### define result paths
# Create desired destination folder for gene lists and GO:BP tales and graphs
setwd(choose.dir())
genelists=paste0(choose.dir(), "\\" )
GOresults=paste0(choose.dir(), "\\" )

### GP table functions
exp_GPtable=function(x) {
  df.name <- deparse(substitute(x))
  GP=as.data.frame(x$result)
  GP=GP[,c("query", "source", "term_id",
           "term_name", "p_value", "query_size", 
           "intersection_size", "term_size", 
           "effective_domain_size", "intersection_size")]
  GP$GeneRatio = paste0(GP$intersection_size,  "/", GP$query_size)
  GP$BgRatio = paste0(GP$term_size, "/", GP$effective_domain_size)
  names(GP) = c("Cluster", "Category", "ID", "Description", "p.adjust", 
                "query_size", "Count", "term_size", "effective_domain_size", 
                "geneID", "GeneRatio", "BgRatio")
  GP$geneID = gsub(",", "/", GP$geneID)
  write.table(GP, paste0(GOresults, "GPtable_", df.name, ".txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  print(paste0("data ", df.name, "has been processed and exported"))
}

GP.plot.table=function(x) {
  GP=as.data.frame(x$result)
  GP=GP[,c("query", "source", "term_id",
                                   "term_name", "p_value", "query_size", 
                                   "intersection_size", "term_size", 
                                   "effective_domain_size", "intersection_size")]
  GP$GeneRatio=paste0(GP$intersection_size,  "/", GP$query_size)
  GP$BgRatio=paste0(GP$term_size, "/", GP$effective_domain_size)
  names(GP)=c("Cluster", "Category", "ID", "Description", "p.adjust", 
                "query_size", "Count", "term_size", "effective_domain_size", 
                "geneID", "GeneRatio", "BgRatio")
  GP$geneID =gsub(",", "/", GP$geneID)
  row.names(GP)=paste0(GP$Cluster,GP$ID)
  GP_cluster=new("compareClusterResult", compareClusterResult = GP)
  return(GP_cluster)
}


#####################################################################################################
                                       ## RNAseq data##
#####################################################################################################

############################ Cocaine vs saline comparison ###########################################

### Defining geneset
# Upregulated genes are defined as those whose Log2 -CPM ratio saline vs cocaine is <0.5
# Downregulated genes are defined as those whose Log2 -CPM ratio saline vs cocaine is >0.5

All_RNAsea=read.table(choose.files(),header = TRUE, sep = "\t", stringsAsFactors = FALSE )

SvC_WT=All_RNAsea[which(All_RNAsea$contrast =="Ctrl.Saline-Ctrl.Cocaine" & abs(All_RNAsea$log2FC)>= 0.5)  ,]
UP_SvC_wt=SvC_WT[which(SvC_WT$up.down== "down-regulated") ,]
Down_SvC_wt=SvC_WT[which(SvC_WT$up.down== "up-regulated") ,]
UP_SvC_wt=UP_SvC_wt[order(UP_SvC_wt$log2FC, decreasing = FALSE),]
Down_SvC_wt=Down_SvC_wt[order(Down_SvC_wt$log2FC, decreasing = TRUE),]

SvC_KO=All_RNAsea[which(All_RNAsea$contrast =="KO.Saline-KO.Cocaine" & abs(All_RNAsea$log2FC)>= 0.5)  ,]
UP_SvC_KO=SvC_KO[which(SvC_KO$up.down== "down-regulated") ,]
Down_SvC_KO=SvC_KO[which(SvC_KO$up.down== "up-regulated") ,]
UP_SvC_KO=UP_SvC_KO[order(UP_SvC_KO$log2FC, decreasing = FALSE),]
Down_SvC_KO=Down_SvC_KO[order(Down_SvC_KO$log2FC, decreasing = TRUE),]

Down_MD=setdiff(Down_SvC_wt$target,Down_SvC_KO$target)# Maged1 dependent (MD) downregulated genes are deefned as those that are downregulates in Cocaine vs Saline conditions 
# in WT mice but that are not differentially expressed (FDR >0,5 or Log2 |FC|<0,5) in Cocaine vs Saline in Maged1 cKO mice
Up_MD=setdiff(UP_SvC_wt$target,UP_SvC_KO$target)# Maged1 dependent(MD)  upregulated genes are deefned as those that are upregulates in Cocaine vs Saline conditions 
# in WT mice but that are not differentially expressed (FDR >0,5 or Log2 |FC|<0,5) in Cocaine vs Saline in Maged1 cKO mice

### gProfiler analysis

UP_SvC_wt_names = gconvert(UP_SvC_wt$target)
Down_SvC_wt_names = gconvert(Down_SvC_wt$target)
up_names_SvC_KO = gconvert(UP_SvC_KO$target)
down_names_SvC_KO = gconvert(Down_SvC_KO$target)
up_MD_names = gconvert(Up_MD)
down_MD_names = gconvert(Down_MD) 

GP_Up_SvC_wt_RNA=gost(UP_SvC_wt_names$name, organism ="mmusculus" , ordered_query = TRUE)
GP_Down_SvC_wt_RNA=gost(Down_SvC_wt_names$name, organism ="mmusculus" , ordered_query = TRUE)
GP_Up_SVC_RNA_KO=gost(up_names_SvC_KO$name, organism ="mmusculus" , ordered_query = TRUE)
GP_Down_SvC_KO_RNA=gost(down_names_SvC_KO$name, organism ="mmusculus" , ordered_query = TRUE)
GP_up_MD=gost(up_MD_names$name, organism ="mmusculus" , ordered_query = TRUE)
GP_down_MD=gost(down_MD_names$name, organism ="mmusculus" , ordered_query = TRUE)

### Export gene lists and GP tables
write.table(UP_SvC_wt$target, paste0(genelists, "genelist_Cocaine_upreg_genes_WT.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t" )
write.table(Down_SvC_wt$target, paste0(genelists, "genelist_Cocaine_downreg_genes_WT.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t" )
write.table(UP_SvC_KO$target, paste0(genelists, "genelist_Cocaine_upreg_genes_KO.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t" )
write.table(Down_SvC_KO$target, paste0(genelists, "genelist_Cocaine_downreg_genes_KO.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t" )
write.table(Up_MD, paste0(genelists, "genelist_Cocaine_upreg_genes_MD.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t" )
write.table(Down_MD, paste0(genelists, "genelist_Cocaine_downreg_genes_MD.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t" )

exp_GPtable(x=GP_Up_SvC_wt_RNA)
exp_GPtable(x=GP_Down_SvC_wt_RNA)
exp_GPtable(x=GP_Up_SVC_RNA_KO)
exp_GPtable(x=GP_Down_SvC_KO_RNA)
exp_GPtable(x=GP_up_MD)
exp_GPtable(x=GP_down_MD)

### GP plots 

# Up/down regulated genes in cocaine vs saline condition in WT mice
RNA_multi_SvC_wt=gost(list("up-regulated" = UP_SvC_wt_names$name,
                         "down-regulated" = Down_SvC_wt_names$name),organism ="mmusculus" ,ordered_query = TRUE)
pdf(paste0(GOresults, "gostplot_WT_sal_vs_coc_RNA.pdf"),
    width = 16, height = 4 )
gostplot(RNA_multi_SvC_wt, interactive = FALSE,)
dev.off()

# Up/down regulated genes in cocaine vs saline condition in Maged1 cKO mice
RNA_multi_SvC_KO=gost(list("up-regulated" = up_names_SvC_KO$name,
                         "down-regulated" = down_names_SvC_KO$name),organism ="mmusculus" ,ordered_query = TRUE)
pdf(paste0(GOresults,"gostplot_KO_sal_vs_coc_RNA.pdf"),
    width = 16, height = 4 )
gostplot(RNA_multi_SvC_KO, interactive = FALSE)
dev.off()

# Up/down regulated genes in cocaine vs saline condition in WT mice but not in Maged1 cKO animals
RNA_multi_sal_vscoc_MD=gost(list("up-regulated" = up_MD_names$name,
                                   "down-regulated" = down_MD_names$name), organism ="mmusculus" ,ordered_query = TRUE)
pdf(paste0(GOresults, "gostplot_WT_sal_vs_coc_wt spec_RNA.pdf"),
    width = 16, height = 4 )
gostplot(RNA_multi_sal_vscoc_MD, interactive = FALSE,)
dev.off()

### Enrichplot dotplot

sel20_DOWN_coc_GO=c("GO:0001578","GO:0007018","GO:0006928","GO:0000226","GO:0120036",
                    "GO:0030030","GO:1903530","GO:0007010","GO:0051046","GO:0015844",
                    "GO:0023061","GO:0006810","GO:0051047","GO:0032940","GO:0006811",
                    "GO:0006996","GO:0051937","GO:0099111","GO:0006812","GO:0015872")
sel20_UP_coc_GO=c("GO:0031175","GO:0007010","GO:0031344","GO:0008366","GO:0006915",
                  "GO:0033043","GO:0032291","GO:0030030","GO:0012501","GO:0033554",
                  "GO:0061564", "GO:0008219","GO:0097191","GO:0010941","GO:0007015",
                  "GO:0030029", "GO:0099565","GO:0032970","GO:0099601","GO:0006950")


RNA_multi_sal_vscoc_MD$result=RNA_multi_sal_vscoc_MD$result[
  c(which(RNA_multi_sal_vscoc_MD$result$query == "down-regulated" &  RNA_multi_sal_vscoc_MD$result$term_id %in% sel20_DOWN_coc_GO),
    which(RNA_multi_sal_vscoc_MD$result$query == "up-regulated" &  RNA_multi_sal_vscoc_MD$result$term_id  %in% sel20_UP_coc_GO)),]


DP_sal_vscoc_MD=GP.plot.table(RNA_multi_sal_vscoc_MD)
DP_sal_vscoc_MD@compareClusterResult$Cluster=gsub("down-regulated", "DOWN",DP_sal_vscoc_MD@compareClusterResult$Cluster )
DP_sal_vscoc_MD@compareClusterResult$Cluster=gsub("up-regulated", "UP",DP_sal_vscoc_MD@compareClusterResult$Cluster )

DP_sal_vscoc_MD@compareClusterResult$ID=paste0(DP_sal_vscoc_MD@compareClusterResult$Cluster, "_", DP_sal_vscoc_MD@compareClusterResult$ID)
DP_sal_vscoc_MD@compareClusterResult$Description=paste0(DP_sal_vscoc_MD@compareClusterResult$Cluster, "_", DP_sal_vscoc_MD@compareClusterResult$Description)


svg(paste0(GOresults,"enrichplot_GO_sal_vs_coc_WT_RNA_MD_final.svg"))
enrichplot::dotplot(DP_sal_vscoc_MD, showCategory=40, by = "Cluster" ,size = "geneRatio")+ scale_colour_gradient(low="red", high="blue", trans="log10")
dev.off()



#####################################################################################################
                                         ##ChIP-M data##
#####################################################################################################

### The H2Aub+ gene lists are the results of the CHIPH2Aub.R script analysis

############################ Cocaine vs saline comparison ###########################################

### Importing geneset
ChIP_coc_spec_WT=as.vector(read.table(choose.files(),header = FALSE, sep = "\t", stringsAsFactors = FALSE )[,1])
ChIP_sal_spec_WT=as.vector(read.table(choose.files(),header = FALSE, sep = "\t", stringsAsFactors = FALSE )[,1])
ChIP_coc_spec_KO=as.vector(read.table(choose.files(),header = FALSE, sep = "\t", stringsAsFactors = FALSE )[,1])
ChIP_sal_spec_KO=as.vector(read.table(choose.files(),header = FALSE, sep = "\t", stringsAsFactors = FALSE )[,1])

ChIP_coc_spec_MD=setdiff(ChIP_coc_spec_WT, ChIP_coc_spec_KO )# Maged1 dependent (MD) Cocaine H2Aub+ specific genes are defined as those that are H2Aub+ in Cocaine but not Saline conditions 
# in WT mice but that are not differentially enriched in Maged1 cKO mice
ChIP_sal_spec_MD=setdiff(ChIP_sal_spec_WT, ChIP_sal_spec_KO )# Maged1 dependent (MD) Saline H2Aub+ specific genes are defined as those that are H2Aub+ in Saline but not Cocaine conditions 
# in WT mice but that are not differentially enriched in Maged1 cKO mice

### GeneProfiler analysis
ChIP_coc_spec_WT.names=gconvert(ChIP_coc_spec_WT)
ChIP_sal_spec_WT.names=gconvert(ChIP_sal_spec_WT)
ChIP_coc_spec_KO.names=gconvert(ChIP_coc_spec_KO)
ChIP_sal_spec_KO.names=gconvert(ChIP_sal_spec_KO)
ChIP_coc_spec_MD.names=gconvert(ChIP_coc_spec_MD)
ChIP_sal_spec_MD.names=gconvert(ChIP_sal_spec_MD)

GP_ChIP_coc_spec_WT=gost(ChIP_coc_spec_WT.names$name, organism ="mmusculus" , ordered_query = FALSE)
GP_ChIP_sal_spec_WT=gost(ChIP_sal_spec_WT.names$name, organism ="mmusculus" , ordered_query = FALSE)
GP_ChIP_coc_spec_KO=gost(ChIP_coc_spec_KO.names$name, organism ="mmusculus" , ordered_query = FALSE)
GP_ChIP_sal_spec_KO=gost(ChIP_sal_spec_KO.names$name, organism ="mmusculus" , ordered_query = FALSE)
GP_ChIP_coc_spec_MD=gost(ChIP_coc_spec_MD.names$name, organism ="mmusculus" , ordered_query = FALSE)
GP_ChIP_sal_spec_MD=gost(ChIP_sal_spec_MD.names$name, organism ="mmusculus" , ordered_query = FALSE)

### Export GP tables
exp_GPtable(x=GP_ChIP_coc_spec_WT)
exp_GPtable(x=GP_ChIP_sal_spec_WT)
exp_GPtable(x=GP_ChIP_coc_spec_KO)
exp_GPtable(x=GP_ChIP_sal_spec_KO)
exp_GPtable(x=GP_ChIP_coc_spec_MD)
exp_GPtable(x=GP_ChIP_sal_spec_MD)

### GP plots
ChIP_multi_sal_vs_coc_WT = gost(list("Cocaine specific" = ChIP_coc_spec_WT.names$name,
                           "Saline specific" = ChIP_sal_spec_WT.names$name), organism ="mmusculus" ,ordered_query = FALSE)
pdf(paste0(GOresults,"gostplot_ChIP_sal_vs_coc_WT.pdf"),  width = 16, height = 4 )
gostplot(ChIP_multi_sal_vs_coc_WT, interactive = FALSE)
dev.off()

ChIP_multi_sal_vs_coc_KO = gost(list("Cocaine specific" = ChIP_coc_spec_KO.names$name,
                                "Saline specific" = ChIP_sal_spec_KO.names$name), organism ="mmusculus" ,ordered_query = FALSE)
pdf(paste0(GOresults,"gostplot_ChIP_sal_vs_coc_KO.pdf"),  width = 16, height = 4 )
gostplot(ChIP_multi_sal_vs_coc_KO, interactive = FALSE)
dev.off()

ChIP_multi_sal_vs_coc_MD = gost(list("Cocaine specific" = ChIP_coc_spec_MD.names$name,
                               "Saline specific" = ChIP_sal_spec_MD.names$name), organism ="mmusculus" ,ordered_query = FALSE)
pdf(paste0(GOresults,"gostplot_ChIP_WT_vs_KO_sal.pdf"),  width = 16, height = 4 )
gostplot(ChIP_multi_sal_vs_coc_MD, interactive = FALSE)
dev.off()

### Enrichplot dotplots
sel20_H2ACoc_MD=c("GO:0006810","GO:0006928","GO:0006811","GO:0030030","GO:0120036",
                  "GO:0006812","GO:0031175","GO:0034220","GO:0023061","GO:0098655", 
                  "GO:0032940","GO:0051046","GO:0007010","GO:0034765","GO:0098916",
                  "GO:0050803","GO:0097485","GO:0007409","GO:0007411","GO:0007416")




  ChIP_multi_sal_vs_coc_MD$result=ChIP_multi_sal_vs_coc_MD$result[
  c(which(ChIP_multi_sal_vs_coc_MD$result$query == "Cocaine specific" &  ChIP_multi_sal_vs_coc_MD$result$term_id %in% sel20_H2ACoc_MD),
    which(ChIP_multi_sal_vs_coc_MD$result$query == "Saline specific"&  ChIP_multi_sal_vs_coc_MD$result$source == "GO:BP")),]

DP_sal_vscoc_MD_chip=GP.plot.table(ChIP_multi_sal_vs_coc_MD)

DP_sal_vscoc_MD_chip@compareClusterResult$ID=paste0(DP_sal_vscoc_MD_chip@compareClusterResult$Cluster, "_", DP_sal_vscoc_MD_chip@compareClusterResult$ID)
DP_sal_vscoc_MD_chip@compareClusterResult$Description=paste0(DP_sal_vscoc_MD_chip@compareClusterResult$Cluster, "_", DP_sal_vscoc_MD_chip@compareClusterResult$Description)

svg(paste0(GOresults,"enrichplot_GO_sal_vs_coc_WT_ChIP_final.svg"))
enrichplot::dotplot(DP_sal_vscoc_MD_chip, showCategory=40, by = "Cluster" ,size = "geneRatio", font.size=8 )+scale_colour_gradient( low="red", high="blue", trans="log10") 
dev.off()




#####################################################################################################
## Analysis enrichment GOterms H2Aub total vs MD #
#####################################################################################################

GOBP_ChIP_coc_spec_WT=as.vector(GP_ChIP_coc_spec_WT$result[which(GP_ChIP_coc_spec_WT$result$source=="GO:BP"),9])
GOBP_ChIP_coc_spec_WT_s=ID_GP_H2A_s=simplifyGOterms(GOBP_ChIP_coc_spec_WT, maxOverlap= 0.7, ontology="BP", go2allEGs=org.Mm.egGO2ALLEGS)

GOBP_ChIP_coc_spec_MD=as.vector(GP_ChIP_coc_spec_MD$result[which(GP_ChIP_coc_spec_MD$result$source=="GO:BP"),9])
GOBP_ChIP_coc_spec_MD_s=simplifyGOterms(GOBP_ChIP_coc_spec_MD, maxOverlap= 0.7, ontology="BP", go2allEGs=org.Mm.egGO2ALLEGS)

GOBP_RNA_coc_spec_MD=as.vector(GP_down_MD$result[which(GP_down_MD$result$source=="GO:BP"),9])
GOBP_RNA_coc_spec_MD_s=simplifyGOterms(GOBP_RNA_coc_spec_MD, maxOverlap= 0.7, ontology="BP", go2allEGs=org.Mm.egGO2ALLEGS)

GOBP_RNA_coc_spec_WT=as.vector(GP_Down_SvC_wt_RNA$result[which(GP_Down_SvC_wt_RNA$result$source=="GO:BP"),9])
GOBP_RNA_coc_spec_WT_s=simplifyGOterms(GOBP_RNA_coc_spec_WT, maxOverlap= 0.7, ontology="BP", go2allEGs=org.Mm.egGO2ALLEGS)



GOBP_common_MD_s=intersect(GOBP_ChIP_coc_spec_MD_s, GOBP_RNA_coc_spec_MD_s)
GOBP_common_MD=intersect(GOBP_ChIP_coc_spec_MD, GOBP_RNA_coc_spec_MD)

common_MD_goterm_names=GP_ChIP_coc_spec_MD$result[which(GP_ChIP_coc_spec_MD$result$term_id %in% GOBP_common_MD_s), c(9,11)]



GOBP_coc_spec_WT_MD=intersect(sel20_H2ACoc_MD, intersect(GOBP_ChIP_coc_spec_WT,GOBP_ChIP_coc_spec_MD))




GOBP_ChIP_sal_spec_WT=GP_ChIP_sal_spec_WT$result[which(GP_ChIP_sal_spec_WT$result$source=="GO:BP"),9]
GOBP_ChIP_sal_spec_MD=GP_ChIP_sal_spec_MD$result[which(GP_ChIP_sal_spec_MD$result$source=="GO:BP"),9]

GOBP_sal_spec_WT_MD_S=intersect(GOBP_ChIP_sal_spec_WT,GOBP_ChIP_sal_spec_MD)
GOBP_sal_spec_WT_MD=intersect(GOBP_ChIP_sal_spec_WT,GOBP_ChIP_sal_spec_MD)


SEL_GP_ChIP_coc_spec_WT=GP_ChIP_coc_spec_WT$result[which(GP_ChIP_coc_spec_WT$result$term_id %in% GOBP_coc_spec_WT_MD),]
SEL_GP_ChIP_coc_spec_MD=GP_ChIP_coc_spec_MD$result[which(GP_ChIP_coc_spec_MD$result$term_id %in% GOBP_coc_spec_WT_MD),]


SEL_GP_ChIP_coc_spec_WT$E=(SEL_GP_ChIP_coc_spec_WT$intersection_size/SEL_GP_ChIP_coc_spec_WT$query_size)/(SEL_GP_ChIP_coc_spec_WT$term_size/SEL_GP_ChIP_coc_spec_WT$effective_domain_size)
SEL_GP_ChIP_coc_spec_MD$E=(SEL_GP_ChIP_coc_spec_MD$intersection_size/SEL_GP_ChIP_coc_spec_MD$query_size)/(SEL_GP_ChIP_coc_spec_MD$term_size/SEL_GP_ChIP_coc_spec_MD$effective_domain_size)

E_coc_spec_vs_MD=cbind(SEL_GP_ChIP_coc_spec_WT$term_id, SEL_GP_ChIP_coc_spec_WT$term_name, SEL_GP_ChIP_coc_spec_WT$E, SEL_GP_ChIP_coc_spec_MD$E)


SEL_GP_ChIP_sal_spec_WT=GP_ChIP_sal_spec_WT$result[which(GP_ChIP_sal_spec_WT$result$term_id %in% GOBP_sal_spec_WT_MD),]
colnames(SEL_GP_ChIP_sal_spec_WT)=c("term_ID", "term name", "E coc spec", "E coc MD")

length(intersect(ChIP_coc_spec_WT, Down_SvC_wt$target))

######################################


Evalues_RNA_DOWN_MD =RNA_multi_sal_vscoc_MD$result[which(RNA_multi_sal_vscoc_MD$result$query == "down-regulated" 
                                                                                   &  RNA_multi_sal_vscoc_MD$result$term_id %in% sel20_DOWN_coc_GO),]

Evalues_RNA_DOWN_MD$Evalue=(Evalues_RNA_DOWN_MD$intersection_size/Evalues_RNA_DOWN_MD$query_size)/
                          (Evalues_RNA_DOWN_MD$term_size/Evalues_RNA_DOWN_MD$effective_domain_size)




Evalues_RNA_DOWN_WT =RNA_multi_SvC_wt$result[which(RNA_multi_SvC_wt$result$query == "down-regulated" 
                                                   &  RNA_multi_SvC_wt$result$term_id %in% sel20_DOWN_coc_GO),]

Evalues_RNA_DOWN_WT$Evalue=(Evalues_RNA_DOWN_WT$intersection_size/Evalues_RNA_DOWN_WT$query_size)/
  (Evalues_RNA_DOWN_WT$term_size/Evalues_RNA_DOWN_WT$effective_domain_size)


common_Evalues=intersect(Evalues_RNA_DOWN_MD$term_id, Evalues_RNA_DOWN_WT$term_id)


Evalues_RNA_DOWN_MD=Evalues_RNA_DOWN_MD[which(Evalues_RNA_DOWN_MD$term_id %in% common_Evalues),]
Evalues_RNA_DOWN_WT=Evalues_RNA_DOWN_WT[which(Evalues_RNA_DOWN_WT$term_id %in% common_Evalues),]




Evalues_RNA_DOWN=data.frame("GO_ID"=Evalues_RNA_DOWN_MD$term_id, "GO_name"=Evalues_RNA_DOWN_MD$term_name, 
                               "EvalueMD"=Evalues_RNA_DOWN_MD$Evalue,  "EvalueWT"=Evalues_RNA_DOWN_WT$Evalue)


################################################


Evalues_ChIP_DOWN_MD =ChIP_multi_sal_vs_coc_MD$result[which(ChIP_multi_sal_vs_coc_MD$result$query == "Cocaine specific"
                                                         &  ChIP_multi_sal_vs_coc_MD$result$term_id %in% sel20_DOWN_coc_GO),]

Evalues_ChIP_DOWN_MD$Evalue=(Evalues_ChIP_DOWN_MD$intersection_size/Evalues_ChIP_DOWN_MD$query_size)/
  (Evalues_ChIP_DOWN_MD$term_size/Evalues_ChIP_DOWN_MD$effective_domain_size)


Evalues_ChIP_DOWN_WT =ChIP_multi_sal_vs_coc_WT$result[which(ChIP_multi_sal_vs_coc_WT$result$query == "Cocaine specific"
                                                            &  ChIP_multi_sal_vs_coc_WT$result$term_id %in% sel20_DOWN_coc_GO),]

Evalues_ChIP_DOWN_WT$Evalue=(Evalues_ChIP_DOWN_WT$intersection_size/Evalues_ChIP_DOWN_WT$query_size)/
  (Evalues_ChIP_DOWN_WT$term_size/Evalues_ChIP_DOWN_WT$effective_domain_size)



common_Evalues=intersect(Evalues_ChIP_DOWN_MD$term_id, Evalues_ChIP_DOWN_WT$term_id)


Evalues_ChIP_DOWN_MD=Evalues_ChIP_DOWN_MD[which(Evalues_ChIP_DOWN_MD$term_id %in% common_Evalues),]
Evalues_ChIP_DOWN_WT=Evalues_ChIP_DOWN_WT[which(Evalues_ChIP_DOWN_WT$term_id %in% common_Evalues),]






Evalues_ChIP_DOWN=data.frame("GO_ID"=Evalues_ChIP_DOWN_MD$term_id, "GO_name"=Evalues_ChIP_DOWN_MD$term_name, 
                            "EvalueMD"=Evalues_ChIP_DOWN_MD$Evalue,  "EvalueWT"=Evalues_ChIP_DOWN_WT$Evalue)
write.table(Evalues_ChIP_DOWN, "Evalues_ChIP_DOWN.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t" )

write.table(Evalues_RNA_DOWN, "Evalues_RNA_DOWN.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t" )
















