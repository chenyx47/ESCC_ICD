#Protein Omics Data - ESCC
##Data Input and Running----
library(GSVA)
library(openxlsx)
library(dplyr)
library(ConsensusClusterPlus)
library(pheatmap)
setwd("RunningHub")
source("Script/tools.R")

Protein_Omics_ClinicalData<-
  read.xlsx("SourceData/Proteomic_Data.xlsx",
            sheet=2)
Survival_Data<-
  Protein_Omics_ClinicalData[c(1,11:14)]%>%
  mutate(TumorID=paste0(.[,1],"T"))

Survival_Data_addStage<-
  Protein_Omics_ClinicalData[c(1,9,11:14)]%>%
  mutate(TumorID=paste0(.[,1],"T"))%>%
  mutate(pTNM.stage=gsub("A|B|C","",pTNM.stage))
saveRDS(
  Survival_Data_addStage,
  "tmpData/Survival_Data_addStage.RDS"
)

Protein_Omics_ProteinData<-
  read.xlsx("SourceData/Proteomic_Data.xlsx",
            sheet=3)
Protein_Omics_ProteinData_Tumor<-
  Protein_Omics_ProteinData[
    c(2,grep("T",colnames(Protein_Omics_ProteinData)))
    ]
Protein_Omics_ProteinData_Tumor<-
  aggregate(.~Gene.symbol,Protein_Omics_ProteinData_Tumor,mean)

Protein_Omics_ProteinData_Tumor_Matrix<-
  FeaturetoRow(Protein_Omics_ProteinData_Tumor)
Protein_Omics_ProteinData_Tumor_Matrix<-
  apply(Protein_Omics_ProteinData_Tumor_Matrix,1,scale)%>%
  t()
colnames(Protein_Omics_ProteinData_Tumor_Matrix)<-
  colnames(Protein_Omics_ProteinData_Tumor%>%FeaturetoRow())
##Consensus Cluster
# GeneList<-
#   readRDS("tmpData/CellDeath_GeneList_List_Integrate.RDS")
# Protein_Omics_ProteinData_Tumor_GSVA<-
#   gsva(Protein_Omics_ProteinData_Tumor_Matrix%>%
#          as.matrix(),
#        GeneList,method='ssgsea')
# 
# results = ConsensusClusterPlus(Protein_Omics_ProteinData_Tumor_GSVA[c(1,4,5),]%>%
#                                  as.matrix(),
#                                maxK=6,
#                                reps=1000,
#                                pItem=0.8,
#                                pFeature=1,
#                                title="Result/Figure4/ESCC_CD",
#                                clusterAlg="km",
#                                distance="euclidean",
#                                seed=2021,
#                                plot="png")
# clusterNum=2    
# cluster=results[[clusterNum]][["consensusClass"]]%>%as.data.frame()
# CC_cluster<-
#   merge(cluster,
#         Survival_Data,
#         by.x=0,
#         by.y="TumorID")
# names(CC_cluster)[2]<-"Cluster"
# 
# saveRDS(list(ExprMatrix=Protein_Omics_ProteinData_Tumor,
#              Cell_Death_GSVA=Protein_Omics_ProteinData_Tumor_GSVA,
#              CC_cluster=CC_cluster),
#         "tmpData/Proteomics_Data_ESCC_CDdata.RDS")

#####################################################################################################################
##Figure1
#####################################################################################################################
load("SourceData/TCGA_ESCC_Info.Rdata") 
pro.immune.score <- gsva(as.matrix(escc.pro.exp),
                         immune.sig.list) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var = "sampleID")

pro.cdd.score2 <- gsva(as.matrix(escc.pro.exp),
                       cdd.genelist
) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var = "sampleID")


pro.cdd.immune.res2 <- data.frame()
for(cdd.one in names(pro.cdd.score2)[2:6]){
  for(path.one in names(pro.immune.score)[2:12]){
    pro.cdd.immu.cor <- cor.test(pro.cdd.score2[, cdd.one], 
                                 pro.immune.score[, path.one], 
                                 method = "spearman")
    
    pro.cdd.immune.res0 <- data.frame(cdd.type = cdd.one,
                                      feature = path.one,
                                      pvalue = pro.cdd.immu.cor$p.value,
                                      Rho = pro.cdd.immu.cor$estimate[[1]])
    pro.cdd.immune.res2 <- rbind(pro.cdd.immune.res2, pro.cdd.immune.res0)
  }
}

pro.cdd.immune.res2.sig <- filter(pro.cdd.immune.res2,
                                  grepl("SYSTEM", feature) == FALSE,
                                  pvalue < 0.05, Rho >0.2) %>%
  mutate(p.size = pvalue2size(pvalue),
         cdd.type = strsplit2(cdd.type, "_")[,1],
         cdd.type = tolower(cdd.type),
         feature = tolower(feature))


aggregate(.~cdd.type, pro.cdd.immune.res2.sig[, c(1,3,4)], mean)

pro.cdd.immune.res2.plot <- 
  ggplot(pro.cdd.immune.res2.sig, aes(cdd.type, feature, fill = Rho))+
  # geom_tile(data = pro.cdd.immune.res2.sig,
  #           mapping = aes(x = cdd.type, y = factor(feature, unique(feature))),
  #           fill = "white", alpha = 1)+
  theme_bw()+
  geom_point(data = pro.cdd.immune.res2.sig, 
             aes(x = cdd.type, y = feature, color = Rho, size = p.size), inherit.aes = F)+
  scale_color_gradient2(low = 'royalblue2', high = 'orangered', mid = 'white')+
  theme(axis.text.x = element_text(size = 16, face = "plain", angle=30,hjust=1, vjust=1, colour="black"),
        axis.text.y = element_text(size = 16, face = "plain", hjust=1, colour="black"))+
  xlab("")+ylab("")+labs(title="")+
  theme(legend.title=element_text(face="plain", colour="black", size=12))+
  theme(legend.text=element_text(face="plain", colour="black", size=12))+
  theme(legend.title = element_text(color="black", size = 16, face="plain"))

pro.cdd.immune.res2.plot
ggsave(plot = pro.cdd.immune.res2.plot,
       "Result/Figure1/Figure1C.pdf",
       width = 9, height = 3.5)

####################################################################################################################
##Figure3###########################################################################################################
####################################################################################################################

Protein_Omics_ProteinData_Tumor_CD<-
  readRDS("tmpData/Proteomics_Data_ESCC_CDdata.RDS")
Protein_Omics_ProteinData_Tumor_GSVA<-
  Protein_Omics_ProteinData_Tumor_CD$Cell_Death_GSVA
CC_cluster<-
  Protein_Omics_ProteinData_Tumor_CD$CC_cluster
  
CC_cluster$`Overall.survival.(days)`<-
  CC_cluster$`Overall.survival.(days)`/30

CC_cluster$`Disease.free.survival.(days)`<-
  CC_cluster$`Disease.free.survival.(days)`/30

Survfit<-survfit(Surv(`Overall.survival.(days)`,
                      `Death.at.follow-up`=="Yes")~Cluster,CC_cluster)
coxph(Surv(`Overall.survival.(days)`,
           `Death.at.follow-up`=="Yes")~Cluster,CC_cluster)%>%summary()
Plot<-
  ggsurvplot(Survfit,CC_cluster,# surv.median.line ="hv",
             pval = "HR = 0.66\nLog-Rank Pvalue = 0.08",
             legend=c(0.7,0.9),
             legend.title="",
             risk.table = TRUE,
             legend.labs=c("ICD-High Cluster","ICD-Low Cluster"),
             pval.size = 5,
             font.x = c(15, "plain", "black"), # feature of curve xlab title
             font.y = c(15, "plain", "black"), # feature of curve ylab title
             font.tickslab = c(15, "plain", "black"),
             palette = c("#A73030FF",  "#003C67FF"),
             font.title=c(15, "plain", "black"))
pdf("Result/Figure3/Figure3F.pdf",
    width = 6,height = 6)
print(Plot)
dev.off()

Survfit<-survfit(Surv(`Disease.free.survival.(days)`,
                      `Disease.free.status`=="Yes")~Cluster,CC_cluster)
coxph(Surv(`Disease.free.survival.(days)`,
           `Disease.free.status`=="Yes")~Cluster,CC_cluster)%>%summary()
Plot<-ggsurvplot(Survfit,CC_cluster,
           pval = "HR = 0.60\nLog-Rank Pvalue = 0.03",
           legend=c(0.7,0.9),
           legend.title="",
           risk.table = TRUE,
           legend.labs=c("ICD-High Cluster","ICD-Low Cluster"),
           pval.size = 5,
           font.x = c(15, "plain", "black"), # feature of curve xlab title
           font.y = c(15, "plain", "black"), # feature of curve ylab title
           font.tickslab = c(15, "plain", "black"),
           palette = c("#A73030FF",  "#003C67FF"),
           font.title=c(15, "plain", "black"))
pdf("Result/Figure3/Figure3E.pdf",
    width = 6,height = 6)
print(Plot)
dev.off()

##GSVA pathway----
dir.create("Result/FigureS1")
CC_cluster_GSVA<-
  merge(CC_cluster[c(1,2)],
        t(Protein_Omics_ProteinData_Tumor_GSVA),
        by.x="Row.names",by.y=0)
names(CC_cluster_GSVA)[2]<-"Cluster"
names(CC_cluster_GSVA)[c(3:7)]<-
  c("PYROPTOSIS","APOPTOSIS","AUTOPHAGY","FERROPTOSIS","NECROPTOSIS")

CC_cluster_GSVA_melt<-
  reshape2::melt(CC_cluster_GSVA[-1],
                 id="Cluster")
CC_cluster_GSVA_melt$Cluster<-
  ifelse(CC_cluster_GSVA_melt$Cluster==1,
         "ICD-High","ICD-Low")

Plot<-
  ggplot(CC_cluster_GSVA_melt,
         aes(x=Cluster,y=value,fill=Cluster))+geom_boxplot()+
  facet_grid(.~variable)+
  geom_signif(comparisons = list(c("ICD-High","ICD-Low")))+
  scale_fill_manual(values = c("#A73030FF",  "#003C67FF"))+
  theme_classic()+
  ylab("GSVA Score")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18),
        axis.text.y = element_text(size=15),
        strip.background = element_rect(fill = 'white', colour = NULL,
                                        size = rel(5), linetype = 0),
        strip.text = element_text(size=10))+
  guides(fill=FALSE);Plot
pdf("Result/FigureS1/FigureS1C1.pdf",
    width = 10,height = 4)
print(Plot)
dev.off()


CC_cluster_GSVA<-
  CC_cluster_GSVA[order(CC_cluster_GSVA$Cluster),]

ann_col<-list(Cluster=c("ICD-High"="#A73030FF",
                        "ICD-Low"="#003C67FF"))
colData1<-CC_cluster_GSVA[c(1,2)]%>%FeaturetoRow()
colnames(CC_cluster_GSVA)[c(3:7)]<-
  gsub("_.*","",colnames(CC_cluster_GSVA)[c(3:7)])%>%
  toupper()

p1<-pheatmap::pheatmap(CC_cluster_GSVA[-2]%>%FeaturetoRow()%>%t(),
                       color=colorRampPalette(c('#3C5488FF','#3C5488FF','white',
                                                '#DC0000FF','#DC0000FF'), bias=1)(50), border_color=NA,
                       annotation_col = colData1,
                       annotation_colors = ann_col,
                       show_rownames = T,
                       show_colnames = F,
                       scale="row",
                       cluster_rows = F,
                       cluster_cols = F
)
pdf("Result/FigureS1/FigureS1C2.pdf",
    width = 10,height = 3)
print(p1)
dev.off()

##PCA test----
library(factoextra)
library(FactoMineR)
pca2 <- PCA(CC_cluster_GSVA[-2]%>%FeaturetoRow(),
            graph = TRUE)
Plot<-
  fviz_pca_ind(pca2,
               geom.ind = "point", # show points only (nbut not "text")
               col.ind = as.character(CC_cluster_GSVA$Cluster)%>%recode("1"="ICD-High","2"="ICD-Low"), # color by groups
               palette = c("#A73030FF",  "#003C67FF"),
               addEllipses = TRUE, # Concentration ellipses
               legend.title = "Groups"
  );Plot

pdf("Result/FigureS1/FigureS1B.pdf",
    width = 7,height = 6)
print(Plot)
dev.off()


######################################################################################################################
##Figure4#############################################################################################################
######################################################################################################################
##xCell----
calculateHeatmap_sig<-function(Heatmapinput,GroupingData,
                               HighGroup="High",
                               LowGroup="Low"){
  Groups<-unique(GroupingData[,1])
  Heatmapinput<-t(Heatmapinput)%>%as.data.frame()
  Heatmapinput$SampleID<-rownames(Heatmapinput)
  Heatmapinput<-melt(Heatmapinput)
  Heatmapinput$Group<-ifelse(Heatmapinput$SampleID%in%rownames(GroupingData)[GroupingData[,1]==Groups[1]],
                             Groups[1],Groups[2])
  Result_frame<-data.frame(
    Items=unique(Heatmapinput$variable),
    High_mean=rep(NA,length(unique(Heatmapinput$variable))),
    Low_mean=rep(NA,length(unique(Heatmapinput$variable))),
    P_value=rep(NA,length(unique(Heatmapinput$variable)))
  )
  for (i in 1:length(unique(Heatmapinput$variable))) {
    Heatmapinput_chosen<-Heatmapinput[Heatmapinput$variable==unique(Heatmapinput$variable)[i],]
    Wilcox<-wilcox.test(value~Group,Heatmapinput_chosen)
    Result_frame[i,4]<-Wilcox$p.value
    Result_frame[i,2]<-median(Heatmapinput_chosen$value[Heatmapinput_chosen$Group==HighGroup])
    Result_frame[i,3]<-median(Heatmapinput_chosen$value[Heatmapinput_chosen$Group==LowGroup])
  }
  return(Result_frame)
}

ImmuneCellMatrix<-
  read.delim("tmpData/ProteinOmics_xCellResult.txt")
colnames(ImmuneCellMatrix)<-
  gsub("X","",colnames(ImmuneCellMatrix))
ImmuneCellMatrix<-
  merge(CC_cluster[c(1,2)],t(ImmuneCellMatrix),
        by.y=0,by.x="Row.names")
ImmuneCellMatrix<-
  ImmuneCellMatrix[
    order(ImmuneCellMatrix$Cluster),
  ]

SelectCellType<-
  colnames(ImmuneCellMatrix)[
    grep("DC|T|Macro|Mas|Plasma|B|NK|Neutro|Fibro",colnames(ImmuneCellMatrix))
    ]

SigMatrix<-
  calculateHeatmap_sig(ImmuneCellMatrix[-2]%>%FeaturetoRow()%>%t(),
                       ImmuneCellMatrix[c(1,2)]%>%FeaturetoRow(),
                       HighGroup = "1",LowGroup = "2")
###Figure4C----
Boxplot_ForMainNPC<-
  ImmuneCellMatrix%>%FeaturetoRow()%>%
  .[c("Cluster","DC","aDC","Monocytes","Macrophages")]%>%
  reshape2::melt(.,id="Cluster")
Boxplot_ForMainNPC$Cluster<-
  case_when(Boxplot_ForMainNPC$Cluster==1~"ICD-High",
            Boxplot_ForMainNPC$Cluster==2~"ICD-Low")

Plot<-ggplot(Boxplot_ForMainNPC,aes(x=Cluster,y=value,
                              fill=Cluster))+
  geom_violin()+
  geom_boxplot(fill="white",width=0.2)+
  facet_grid(~variable)+
  geom_signif(comparisons = list(c("ICD-High","ICD-Low")))+
  scale_fill_manual(values = c("#A73030FF",  "#003C67FF"))+
  theme_classic()+
  ylab("GSVA Score")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18),
        axis.text.y = element_text(size=15),
        strip.background = element_rect(fill = 'white', colour = NULL,
                                        size = rel(5), linetype = 0),
        strip.text = element_text(size=15))+
  guides(fill=FALSE);Plot
pdf("Result/Figure4/Figure4C.pdf",
    width = 8,height = 4)
print(Plot)
dev.off()
###FigureS2A----
dir.create("Result/FigureS2/")
ImmuneCellMatrix$Cluster<-
  case_when(ImmuneCellMatrix$Cluster==1~"ICD-High",
            ImmuneCellMatrix$Cluster==2~"ICD-Low")

ann_col<-list(Cluster=c("ICD-High"="#A73030FF",
                        "ICD-Low"="#003C67FF"))
colData1<-ImmuneCellMatrix[c(1,2)]%>%FeaturetoRow()
colnames(ImmuneCellMatrix)[c(3:66)]<-
  gsub("_.*","",colnames(ImmuneCellMatrix)[c(3:66)])%>%
  toupper()

p1<-pheatmap::pheatmap(ImmuneCellMatrix[-2]%>%FeaturetoRow()%>%t(),
                       color=colorRampPalette(c('#3C5488FF','#3C5488FF','white',
                                                '#DC0000FF','#DC0000FF'), bias=1)(50), border_color=NA,
                       annotation_col = colData1,
                       annotation_colors = ann_col,
                       show_rownames = T,
                       show_colnames = F,
                       scale="row",
                       cluster_rows = F,
                       cluster_cols = F
)
pdf("Result/FigureS2/FigureS2A.pdf",
    width = 10,height = 8)
print(p1)
dev.off()

##Immune Signal----
ImmuneLandscapeMatrix<-
  read.delim("tmpData/ProteinOmics_IMsignalResult.txt")
colnames(ImmuneLandscapeMatrix)<-
  gsub("X","",colnames(ImmuneLandscapeMatrix))
ImmuneLandscapeMatrix<-
  merge(CC_cluster[c(1,2)],t(ImmuneLandscapeMatrix),
        by.y=0,by.x="Row.names")
Boxplot_ForMainSig<-
  ImmuneLandscapeMatrix%>%FeaturetoRow()%>%
  .[c("Cluster","LIexpression_score","CSF1_response","TGFB_score_21050467","Module3_IFN_score","CHANG_CORE_SERUM_RESPONSE_UP")]%>%
  reshape2::melt(.,id="Cluster")
Boxplot_ForMainSig$variable<-
  recode(Boxplot_ForMainSig$variable,
         "LIexpression_score"="Leukocyte Infiltration",
         "CSF1_response"="Macrophages",
         "TGFB_score_21050467"="TGF-beta",
         "Module3_IFN_score"="IFN-gamma",
         "CHANG_CORE_SERUM_RESPONSE_UP"="Wound healing")

Boxplot_ForMainSig$Cluster<-
  case_when(Boxplot_ForMainSig$Cluster==1~"ICD-High",
            Boxplot_ForMainSig$Cluster==2~"ICD-Low")

Plot<-ggplot(Boxplot_ForMainSig,aes(x=Cluster,y=value,
                                    fill=Cluster))+
  geom_violin()+
  geom_boxplot(fill="white",width=0.2)+
  facet_grid(~variable)+
  geom_signif(comparisons = list(c("ICD-High","ICD-Low")))+
  scale_fill_manual(values = c("#A73030FF",  "#003C67FF"))+
  theme_classic()+
  ylab("GSVA Score")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18),
        axis.text.y = element_text(size=15),
        strip.background = element_rect(fill = 'white', colour = NULL,
                                        size = rel(5), linetype = 0),
        strip.text = element_text(size=15))+
  guides(fill=FALSE);Plot
pdf("Result/Figure4/Figure4E.pdf",
    width = 8,height = 4)
print(Plot)
dev.off()

