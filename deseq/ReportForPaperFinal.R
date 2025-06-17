#' ---
#' title: Deseq comparisons
#' author: George Savva
#' date: June 2025
#' ---

#' Code to generate figure 4 and associated senstivity analyses, tables etc.

#' This is an R source file that can be run with 'source' or spun into a report with rmarkdown::render

#' Don't output messages or warnings in the report (comment this out if you want the warnings and messages)
knitr::opts_chunk$set(message = FALSE, warning = FALSE)

#' Load packages
library(data.table) # Data wrangling
library(here)       # Locating files
library(DESeq2)     # Differential abundance
library(ggplot2)    # Plotting
library(GGally)     # Pairs plot
library(patchwork)  # Multi-panel plots

#' I ran the analyses with and without pairing, and with and without the excluded replicate pair.
#'
#' This iterates over all of the combinations of these options
#' so generates the graphs and tables for different values of DROPONE and PAIRED.
#' 
#' The main analysis in the paper has 'PAIRED=TRUE' and 'DROPONE=FALSE'
#'
#' You can't split a for loop across RMD chunks, so for this file I've hard-coded the primary analysis parameters.

PAIRED=TRUE
DROPONE=TRUE
# for(PAIRED in c(FALSE,TRUE)){
#   for(DROPONE in c(FALSE,TRUE)){

#' Read and encode the data, create the master DDS object

  Plant_data_3_reps <- read.csv(here("RNA_seq data","Plants","3_rep_analysis","Annotated_counts_in_plants_2.csv"), row.names = 1)
  # Extract gene IDs and names.
  plantGeneNames <- Plant_data_3_reps |> as.data.table(keep.rownames = "Gene") |> _[, .(Gene,NamePlant=Name,locus_tagPlant=locus_tag)]
  gene_ids <- rownames(Plant_data_3_reps)
  counts <- as.matrix(Plant_data_3_reps[, c("PK1", "PK2", "PK3", "VS1", "VS2", "VS3")])
  # Replace NA values with 0
  counts[is.na(counts)] <- 0
  # Remove genes with all zero counts from the plant dataset
  plant_counts <- counts[rowSums(counts) > 0, ]
  # Create metadata
  plant_metadata <- data.frame(
    row.names = colnames(counts),
    condition = factor(c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment")),
    pair = factor(c(1, 2, 3, 1, 2, 3))
  )

  # Load the in vitro data
  invitrocounts =fread("RNA_seq data/In Vitro/Annotated_counts_in_vitro2.csv")
  IVGeneNames <- invitrocounts[, .(Gene=Geneid,NameIV=Name,locus_tagIV=Genomeid)]
  # Extract gene IDs and names
  counts <- as.matrix(invitrocounts[, c("PK1", "PK2", "PK3","VS1", "VS2", "VS3")])
  colnames(counts) <- paste0("IV_",colnames(counts))
  rownames(counts) <- invitrocounts$Geneid
  # Replace NA values with 0
  counts[is.na(counts)] <- 0
  # Remove genes with all zero counts
  iv_counts <- counts[rowSums(counts) > 0, ]
  # Create metadata
  iv_metadata <- data.frame(
    row.names = colnames(counts),
    condition = factor(c("Control", "Control","Control", "Treatment","Treatment", "Treatment")),
    pair = factor(c(4,5,6,4,5,6))
  )
  
  ##### CONCATENATE INTO A MASTER ANALYSIS DATASET
  metadata <- rbind(plant_metadata,iv_metadata)[, c("condition","pair")]
  metadata$environment = rep(c("Plant","IV"),each=6) |> factor()
  metadata$pairwithingroup = factor(rep(rep(1:3,2),2))
  
  #### THERE IS ONE MORE ROW IN THE IV DATA?  "nbis-gene-1746"  GETS DROPPED FROM THE PLANT DATASET
  counts <- merge.data.frame(plant_counts,iv_counts,by=0)
  rownames(counts) <- counts$Row.names
  counts$Row.names <- NULL

  # Create DESeq dataset object with a simple design which gets overwritten later
  dds3 <- DESeqDataSetFromMatrix(countData = counts, 
                                 colData   = metadata,
                                 design    = ~ condition + environment)
  
  
#' Check that the counts of genes are similar magnitude across all reps and conditions.  
#' This gives us some reassurance that the genes are lined up correctly and are 
#' measuring roughly the same thing!
#' 
#' This does not depend on the analysis options, its just descriptive.
#' 
allCounts <- ggpairs(counts[1:1000,],upper = list(continuous = wrap(ggally_cor, method = "spearman"))) + 
  scale_x_log10() + 
  scale_y_log10()
#' We only need the pairwise counts file once as it's descriptive.
ggsave(plot=allCounts, filename = "pairwiseCounts.png",width=10,height=10)


#' Now we set up the design that we actually want.
#' Following the guidance from the DESeq vignette on pairing and nesting
  if(PAIRED) {
    formula_main <- ~ environment + environment:condition + environment:pairwithingroup
    formula_short <- ~ condition + pairwithingroup # can't use 'pair' here because it will complain about missing levels.
  }
  if(PAIRED==FALSE){
    formula_main <- ~ environment + environment:condition
    formula_short <- ~ condition
  }

  design(dds3) <- formula_main

  if(DROPONE){
    dds3 <- dds3[,dds3$pair!=1]
    colData(dds3)$pairwithingroup <- droplevels(colData(dds3)$pairwithingroup) 
    colData(dds3)$pairwithingroup <- factor(c(1,2,1,2,1,2,3,1,2,3))
  }

#' These DDS objects are for running the DESeq in the plant or IV environments independently  
#' Checking that the fact we're running all the samples in the same model doesn't cause any weirdness, given the differences
#' in variance between the pairs.

  dds3_plantonly <- dds3[, dds3$environment=="Plant"]
  dds3_invitroonly <- dds3[, dds3$environment=="IV"]
  
  design(dds3_plantonly) <- formula_short
  design(dds3_invitroonly) <- formula_short
  
#' Need to make sure there aren't any empty levels
  colData(dds3_plantonly)$pairwithingroup <- droplevels(colData(dds3_plantonly)$pairwithingroup)
  colData(dds3_invitroonly)$pairwithingroup <- droplevels(colData(dds3_invitroonly)$pairwithingroup)
  
#' Now we can get generate the deseq and get the results
  
## Generate the deseq results
  if(DROPONE==TRUE & PAIRED==TRUE){
    # Comes from the deseq vignette on paired tests in unbalanced groups
    # If we don't do this we have unidentifiable matrix
    m1 <- model.matrix( formula_main, colData(dds3))
    m1 |> unname()
    all.zero <- apply(m1, 2, function(x) all(x==0))
    all.zero
    idx <- which(all.zero)
    m1 <- m1[,-idx]
    dds3 <- dds3 |> DESeq(full = m1) 
  } else { dds3 <- dds3 |> DESeq() } 

  dds3_plantonly <- dds3_plantonly |> DESeq()
  dds3_invitroonly <- dds3_invitroonly |> DESeq()
  
#' Now Extract the interaction effects and main effects (along with shrunk estimators, p-values etc)
  
  interactionEffects <- results(dds3, contrast=list("environmentIV.conditionTreatment","environmentPlant.conditionTreatment")) |> as.data.table(keep.rownames="Gene")
  plantEffects <- results(dds3, name = "environmentPlant.conditionTreatment") |> as.data.table(keep.rownames="Gene")
  invitroEffects <- results(dds3, name = "environmentIV.conditionTreatment") |> as.data.table(keep.rownames="Gene")
  
  dds3.results.plant.shrunk <- lfcShrink(dds3, coef="environmentPlant.conditionTreatment") |> as.data.table(keep.rownames="Gene")
  dds3.results.invitro.shrunk <- lfcShrink(dds3, coef="environmentIV.conditionTreatment")  |> as.data.table(keep.rownames="Gene")
  
  plantOnlyEffects <- results(dds3_plantonly,contrast = c("condition", "Treatment", "Control")) |> as.data.table(keep.rownames="Gene")
  invitroOnlyEffects <- results(dds3_invitroonly,contrast = c("condition", "Treatment", "Control")) |> as.data.table(keep.rownames="Gene")
  
#' Combine the different effects into one table

  interactionEffects=interactionEffects[, .(Gene,INT_L2FC=log2FoldChange, INT_pvalue=pvalue, INT_padj=padj )]
  plantEffects=plantEffects[, .(Gene,Plant_L2FC=log2FoldChange, Plant_pvalue=pvalue, Plant_padj=padj )]
  invitroEffects=invitroEffects[, .(Gene,IV_L2FC=log2FoldChange, IV_pvalue=pvalue, IV_padj=padj )]
  
  invitroEffectsShrunk = dds3.results.invitro.shrunk[ , .(Gene,IVShrunk_L2FC=log2FoldChange, IVShrunk_pvalue=pvalue, IVShrunk_padj=padj ,IVShrunk_baseMean=baseMean)]
  plantEffectsShrunk = dds3.results.plant.shrunk[ , .(Gene,PlantShrunk_L2FC=log2FoldChange, PlantShrunk_pvalue=pvalue, PlantShrunk_padj=padj,PlantShrunk_baseMean=baseMean )]
  
  plantOnlyEffects=plantOnlyEffects[, .(Gene,PlantOnly_baseMean=baseMean,PlantOnly_L2FC=log2FoldChange, PlantOnly_pvalue=pvalue, PlantOnly_padj=padj )]
  invitroOnlyEffects=invitroOnlyEffects[, .(Gene,InvitroOnly_baseMean=baseMean,InvitroOnly_L2FC=log2FoldChange, InvitroOnly_pvalue=pvalue, InvitroOnly_padj=padj )]
  
  combinedEffects <- interactionEffects[plantEffects,on="Gene"][invitroEffects,on="Gene"][plantGeneNames,on="Gene"] |>
    _[IVGeneNames,on="Gene"][invitroOnlyEffects,on="Gene"][plantOnlyEffects,on="Gene"] |>
    _[invitroEffectsShrunk,on="Gene"][plantEffectsShrunk,on="Gene"]

  
## Just checking that all of the names and tags still line up across groups
  combinedEffects[,all.equal(NamePlant,NameIV)]
  combinedEffects[,all.equal(locus_tagIV,locus_tagPlant)]
## All good!

## Get the numbers of genes that are considered up- or down- regulated by each method:

combinedEffectsLong <- combinedEffects |> melt(id.vars = c("Gene","NamePlant","locus_tagPlant","NameIV","locus_tagIV")) |> _[, c("Environment","score") := tstrsplit(variable,"_",fixed=TRUE)] |> 
  dcast(Gene + Environment ~ score, value.var="value") |> 
  _[ , regulated:="Neither"] |>
  _[ L2FC>=1 & padj<=0.05 , regulated:="Up"] |>
  _[ L2FC<= -1 & padj<=0.05 , regulated:="Down"]

fwrite(file = sprintf("numbersOfDAGenes_DROPONE_%s_PAIRED_%s.csv",DROPONE,PAIRED), x=combinedEffectsLong[, table(regulated,Environment)])


    
fig4c <- ggplot(combinedEffects) + aes(x=IV_L2FC, y=Plant_L2FC,col=INT_padj<0.05)+ geom_point(size=1) + 
  scale_color_manual(values=c("black","red"),na.value = NA,name="Interaction effect",
                     labels=c("FALSE"="Not significant (adj.p > 0.05)",
                              "TRUE"="Significant (adj.p < 0.05)")) + coord_fixed() + 
  geom_hline(yintercept=c(-1,1),lty="dotted") + 
  geom_vline(xintercept=c(-1,1),lty="dotted") + theme_bw() + 
  geom_abline(lty="dashed") + stat_smooth(method="lm",col="black")+
  labs(x="Log 2 fold change (in vitro)", y="Log 2 fold change (in planta)")

fig4cShrunk <- ggplot(combinedEffects) + aes(x=IVShrunk_L2FC, y=PlantShrunk_L2FC,col=INT_padj<0.05)+ geom_point(size=1) + 
  scale_color_manual(values=c("black","red"),na.value = NA) + coord_fixed() + 
  geom_hline(yintercept=c(-1,1),lty="dotted") + 
  geom_vline(xintercept=c(-1,1),lty="dotted") + theme_bw() + 
  geom_abline() + theme(legend.position="none") + 
  labs(x="Shrunk Log 2 fold change (in vitro)", y="Shrunk Log 2 fold change (in planta)")


(fig4a <- ggplot(combinedEffectsLong[Environment%in%c("IVShrunk","PlantShrunk")]) + 
  aes(x=baseMean,y=L2FC,col=regulated) + 
  geom_point() + 
  facet_wrap(~Environment,labeller=labeller(.cols=c("IVShrunk"="In vitro","PlantShrunk"="In planta")))+
  scale_x_log10() + 
  geom_hline(yintercept=c(-1,1),lty="dotted") + 
  theme_bw() + theme(strip.background = element_rect(fill = NA))+
  scale_color_manual(values=c("Down"="darkgreen","Neither"="grey",Up="orange"),name="Regulation",labels=c("Down"="Downregulated","Neither"="Not significant","Up"="Upregulated")) + 
  labs(x="Mean expression", y ="Log2 fold change"))


#' Make figure 4
fig4a / fig4cShrunk + plot_layout(heights=c(1,1)) + plot_annotation(tag_levels = "A")

ggsave(sprintf("Figure4.png_DROPONE_%s_PAIRED_%s.png",DROPONE,PAIRED),width=8,height=8) 


#' Spearman correlation between L2FC in the different environments
combinedEffects[, cor.test(IV_L2FC,Plant_L2FC,use = "complete.obs",method="spearman")]

dds3.results.plant.shrunk <- lfcShrink(dds3, coef="environmentPlant.conditionTreatment")
dds3.results.invitro.shrunk <- lfcShrink(dds3, coef="environmentIV.conditionTreatment")

#' Generate the estimate of LFC within each pair by extracting the normalised counts from the deseq model
normalisedCounts_all <- counts(dds3, normalized=TRUE) |> as.data.frame() |> 
  as.data.table(keep.rownames = "Gene") |> setNames(c("Gene",paste(sep = "_",colData(dds3)$environment,colData(dds3)$pair,colData(dds3)$condition))) |> 
  melt(id.vars="Gene") |> _[, c("Env","Rep","Cond") := tstrsplit(variable,"_",fixed=TRUE)] |> 
  dcast(Env+Rep+Gene ~ Cond, value.var="value") |> _[, LFC:=log2(Treatment/Control)] |> _[ , Pair:=paste0(Env,Rep)] |> 
  _[ !is.finite(LFC) , LFC := NA] |>
  dcast(Gene ~ Pair , value.var="LFC")

#' Make the comparison of the LFCs across pairs.
lfccomp <- GGally::ggpairs(normalisedCounts_all[,-1]) + ggtitle("Pairwise correlation between estimated log2-fold change arising within each pair of samples") + theme_bw()

ggsave(plot = lfccomp,sprintf("LFCcomparison_DROPONE_%s_PAIRED_%s.png",DROPONE,PAIRED),width=8,height=8) 
fwrite(file = sprintf("allDeseqOutput_DROPONE_%s_PAIRED_%s.csv",DROPONE,PAIRED),x = combinedEffects)

#   }
# }     # THESE ARE THE CLOSE BRACES FOR THE FOR LOOP I


#' Do some model diagnostics
plot(sizeFactors(dds3), colSums(counts(dds3)))
plotDispEsts(dds3)
vsd <- vst(dds3)
vsn::meanSdPlot(assay(vsd)) 




#' Now reload the primary analysis results (only if you need to)
# 
# combinedEffects<- fread("allDeseqOutput_DROPONE_TRUE_PAIRED_TRUE.csv")
# combinedEffectsLong <- combinedEffects |> melt(id.vars = c("Gene","NamePlant","locus_tagPlant","NameIV","locus_tagIV")) |> _[, c("Environment","score") := tstrsplit(variable,"_",fixed=TRUE)] |> 
#   dcast(Gene + Environment ~ score, value.var="value") |> 
#   _[ , regulated:="Neither"] |>
#   _[ L2FC>=1 & padj<=0.05 , regulated:="Up"] |>
#   _[ L2FC<= -1 & padj<=0.05 , regulated:="Down"]
# 

#' Make the Excel workbook to keep all of the tables.

combinedEffects[order(IV_pvalue),.(locus_tagIV,NameIV,IVShrunk_L2FC,InvitroOnly_L2FC,IV_L2FC,IV_pvalue,IVShrunk_pvalue,IV_padj,IVShrunk_padj)]

pformat = function(x) ifelse(x==0,"<1e-300",sprintf("%0.3e",x))
lformat = function(x) sprintf("%0.3f",x)
combinedEffectsInVitroReports <- combinedEffects[order(-IVShrunk_L2FC),.(
                                        "Locus Tag"=locus_tagIV,
                                        "Name"=NameIV,
                                        "Base mean"=InvitroOnly_baseMean,
                                        "Log2 Fold change (raw)"=lformat(IV_L2FC),
                                        "Log2 Fold change (shrunk)"=lformat(IVShrunk_L2FC),
                                        "P-value"=pformat(IVShrunk_pvalue),
                                        "P-adj."=pformat(IVShrunk_padj))]
head(combinedEffectsInVitroReports,40)
tail(combinedEffectsInVitroReports,40)[order(.N:1)]

combinedEffectsInPlantaReports <- combinedEffects[!is.na(PlantShrunk_padj)][order(-PlantShrunk_L2FC),.("Base mean"=PlantOnly_baseMean,
                                                                         "Locus Tag"=locus_tagPlant,
                                                                         "Name"=NamePlant,
                                                                         "Log2 Fold change (raw)"=lformat(Plant_L2FC),
                                                                         "Log2 Fold change (shrunk)"=lformat(PlantShrunk_L2FC),
                                                                         "P-value"=pformat(PlantShrunk_pvalue),
                                                                         "P-adj."=pformat(PlantShrunk_padj))]
library(openxlsx)

wb=createWorkbook()
addWorksheet(wb,"Summary")
addWorksheet(wb,"Intersection")
addWorksheet(wb,"In planta upregulated")
addWorksheet(wb,"In planta downregulated")
addWorksheet(wb,"In vitro upregulated")
addWorksheet(wb,"In vitro downregulated")
addWorksheet(wb,"All results")
writeData(wb,"Summary",combinedEffectsLong[, table(regulated,Environment)])
writeData(wb,"In planta upregulated",head(combinedEffectsInPlantaReports,40))
writeData(wb,"In planta downregulated",tail(combinedEffectsInPlantaReports,40)[order(.N:1)])
writeData(wb,"In vitro upregulated",head(combinedEffectsInVitroReports,40))
writeData(wb,"In vitro downregulated",tail(combinedEffectsInVitroReports,40)[order(.N:1)])
writeData(wb,"All results",combinedEffects)

pairedRegulated <- combinedEffectsLong[Environment %in% c("Plant","IV"), .(Gene,Environment,regulated)] |> 
  dcast(Gene~Environment) |>
  _[,table(IV,Plant)]

pairedRegulated

discordantGenes <- combinedEffectsLong[Environment %in% c("Plant","IV"), .(Gene,Environment,regulated)] |> 
  dcast(Gene~Environment) |> 
  _[Plant!="Neither" & IV!="Neither"] |>
  _[order(IV,Plant)]

#' These aren't all discordant, it's just the list of the genes that are significant in both models.  
#' Turns out most are discordant.
discordantGenes <- IVGeneNames[discordantGenes , on="Gene"]

discordantGenes

writeData(wb,"Intersection",pairedRegulated)
writeData(wb,"Intersection",discordantGenes,startRow = 10)


saveWorkbook(wb, "DeseqResults.xlsx",overwrite=TRUE)


combinedEffectsLong[Environment %in% c("Plant","IV"), .(Gene,Environment,regulated)] |> 
  dcast(Gene~Environment) |>
  _[,table(IV,Plant)] |> chisq.test()

## Kendall's tau for the concordance between environments.
combinedEffectsLong[Environment %in% c("Plant","IV"), .(Gene,Environment,regulated)] |> dcast(Gene~Environment) |>
  _[,table(IV,Plant)] |> 
  DescTools::KendallTauB(conf.level=0.95)


