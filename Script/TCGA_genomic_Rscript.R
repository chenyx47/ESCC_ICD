#TCGA Data Anlaysis----
library(data.table)
library(GSVA)
library(limma)
library(tidyverse)
library(maftools)
library(ggsignif)
library(ggplot2)
library(survival)
library(reshape2)
library(tidyverse)
library(maftools)
library(ggsignif)

#####################################################################################################################
##Figure2
#####################################################################################################################
###Mutation----

load("SourceData/TCGA_ESCC_Info.Rdata")
dir.create("Result/Figure2")

col = c("#8a8acb", "#2bb3f3", "#56b881", "#50667f", "#41afa5", 
        "#f9886c", "#e55e5e", "#ed6498")
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins',
               'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')

pdf("Result/Figure2/Figure2A.pdf",
    width = 8, height = 5.5)
oncoplot(escc.maf, colors = col)
dev.off()

###CNV----

CNV_matrix<-read.table("SourceData/TCGA-ESCA.gistic.tsv",sep="\t",header = T)
gene_probe<-read.table("SourceData/gencode.v22.annotation.gene.probeMap",sep="\t",header = T)
colnames(CNV_matrix)[1]<-"id"
gene_probe<-gene_probe[,c(1:2)]

CNV_matrix<-merge(gene_probe,CNV_matrix,by="id")
CNV_matrix<-CNV_matrix[which(!(duplicated(CNV_matrix$gene))),]

rownames(CNV_matrix)<-CNV_matrix$gene
CNV_matrix<-CNV_matrix[,c(-1,-2)]
colnames(CNV_matrix)<-gsub("\\.","-",colnames(CNV_matrix))

pheano_table<-read.table("6.CNV/pheano.txt",sep="\t",header = T)
CNV_matrix<-CNV_matrix[,which(colnames(CNV_matrix) %in% pheano_table$submitter_id.samples[which(pheano_table$primary_diagnosis.diagnoses=="Squamous cell carcinoma, NOS")])]

CNV_percentage_table<-data.frame(Gene=rownames(CNV_matrix), Loss=rowSums(CNV_matrix<0)/ncol(CNV_matrix),Gain=rowSums(CNV_matrix>0)/ncol(CNV_matrix))


gene_table<-readRDS("SourceData/CellDeath_GeneList_List_Integrate.RDS")

Pyroptosis_gene_list<-gene_table$Pyroptosis_Gene
Ferroptosis_gene_list<-gene_table$Ferroptosis_Gene_KEGG
Necroptosis_gene_list<-gene_table$Necroptosis_GeneList_GO_NCD
total_gene_list<-c(Pyroptosis_gene_list,Ferroptosis_gene_list,Necroptosis_gene_list)


gene_table<-list(Pyroptosis=Pyroptosis_gene_list,
                 Ferroptosis=Ferroptosis_gene_list,
                 Necroptosis=Necroptosis_gene_list,Total=total_gene_list)

type_list=c("Pyroptosis",
            "Ferroptosis",
            "Necroptosis","Total")

pdf("Result/Figure2/Figure2BCD.pdf",width = 8,height = 3)
for (i in 1:length(type_list)){
  type<-type_list[i]
  gene_list<-gene_table[[type]]
  print(gene_list)
  CNV_percentage_table1<-CNV_percentage_table[which(CNV_percentage_table$Gene %in% gene_list),]
  
  CNV_percentage_table1<-CNV_percentage_table1[order(CNV_percentage_table1$Loss,decreasing = T),]
  Loss_table<-head(CNV_percentage_table1,10)
  
  CNV_percentage_table1<-CNV_percentage_table1[order(CNV_percentage_table1$Gain,decreasing = T),]
  Gain_table<-head(CNV_percentage_table1,10)
  
  result_table<-rbind(Loss_table,Gain_table)
  
  result_table1<-data.frame(Gene=result_table$Gene,CNV=result_table$Loss,Type="Loss")
  result_table2<-data.frame(Gene=result_table$Gene,CNV=result_table$Gain,Type="Gain")
  
  result_table<-rbind(result_table1,result_table2)
  result_table$Flag<-paste(result_table$Gene,result_table$Type,sep="-")
  result_table<-result_table[which(!(duplicated(result_table$Flag))),]
  result_table$Gene<-factor(result_table$Gene,levels=unique(c(Loss_table$Gene,Gain_table$Gene)))
  
  
  p<-ggplot(data=result_table, aes(x=Gene, y=CNV, fill=Type)) + 
    geom_col(position = 'stack', width = 0.6)+
    theme_classic()+
    scale_fill_manual(values = c("#A73030FF","#003C67FF"))+
    labs(x = 'Gene name', y ="CNV frequency",title=paste("CNV frequency of",type)) +
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
    theme(axis.text = element_text(size = 11), axis.title = element_text(size = 11), legend.text = element_text(size = 11, face='italic'),
          plot.title = element_text(size = 12,vjust = 0.5,hjust=0.5),
          axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
          axis.title.x = element_blank())
  print(p)
}

dev.off()

################################################################################################
##Figure6
################################################################################################
###compare CNV of different groups ----
load("SourceData/TCGA_ESCC_Info.Rdata")

escc.gistic.11 <- filter(escc.gistic, `11q13.3`>=0) %>%
  mutate(amp.del = ifelse(`11q13.3` >0, "Amp", "Normal"))

escc.gistic.11.plot <- 
  ggplot(escc.gistic.11,
         aes(x = amp.del, y = SelectedGene, fill = amp.del))+
  geom_violin()+
  geom_boxplot(width=0.4)+
  theme_bw()+
  guides(fill="none")+
  # geom_dotplot(binaxis='y', stackdir='center', binwidth = 2.6,
  #              stackratio=1.5, dotsize=0.05)+
  scale_fill_manual(values = c("#fe5000", "#84bd00"))+
  ylab("Value of Signature")+xlab("")+labs(color = "")+
  theme(axis.text.x = element_text(size = 15, face = "plain", angle = 0, 
                                   hjust = 0.5, vjust = 1, colour="black"), #,
        axis.text.y = element_text(size = 15, face = "plain", hjust=1, colour="black"),
        axis.title.x = element_text(size = 10, face = "plain", colour="black"),
        axis.title.y = element_text(size = 16, face = "plain", colour="black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "bottom",
        strip.text = element_text(color="black", face = "plain", size = 16),
        strip.background = element_blank())+
  geom_signif(comparisons = list(c("Amp", "Normal")), test = "wilcox.test", 
              textsize = 5, family = "ArialMT")+
  coord_cartesian(ylim = c(-1, 1))

escc.gistic.11.plot
ggsave(escc.gistic.11.plot, 
       filename = "Result/Figure6/Figure6C.pdf",
       width = 3, height = 3.3)


escc.gistic.9 <- mutate(escc.gistic,
                        amp.del = ifelse(`9p21.3` < -0.5, "Del", 
                                         ifelse(`9p21.3` > 0.5, "Non.Del",
                                                ifelse(`9p21.3` == 0, "Non.Del", "Out")))) %>%
  filter(amp.del != "Out")

escc.gistic.9.plot <- 
  ggplot(escc.gistic.9,
         aes(x = amp.del, y = SelectedGene, fill = amp.del))+
  geom_violin()+
  geom_boxplot(width=0.4)+
  theme_bw()+
  guides(fill="none")+
  # geom_dotplot(binaxis='y', stackdir='center', binwidth = 2.6,
  #              stackratio=1.5, dotsize=0.05)+
  scale_fill_manual(values = c("#00aeff", "#84bd00"))+
  ylab("Value of Signature")+xlab("")+labs(color = "")+
  theme(axis.text.x = element_text(size = 15, face = "plain", angle = 0, 
                                   hjust = 0.5, vjust = 1, colour="black"), #,
        axis.text.y = element_text(size = 15, face = "plain", hjust=1, colour="black"),
        axis.title.x = element_text(size = 10, face = "plain", colour="black"),
        axis.title.y = element_text(size = 16, face = "plain", colour="black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "bottom",
        strip.text = element_text(color="black", face = "plain", size = 16),
        strip.background = element_blank())+
  geom_signif(comparisons = list(c("Non.Del", "Del")),
              # c("Del", "Normal"),
              # c("Amp", "Normal")), 
              test = "wilcox.test", 
              textsize = 5, family = "ArialMT",step_increase = 0.1)+
  coord_cartesian(ylim = c(-1, 1))

escc.gistic.9.plot
ggsave(escc.gistic.9.plot, 
       filename = "Result/Figure6/Figure6D.pdf",
       width = 3, height = 3.3)

