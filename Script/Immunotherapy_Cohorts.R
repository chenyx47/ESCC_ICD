#ImmuneCohort----
library(survival)
library(survminer)
IO_All_Sets_Input<-
  readRDS("SourceData/IO_Cohort.RDS")

Select_Cohort<-
  c("PRJEB23709_combo_df","nm_combo_df","van_combo_df","cpm_atz_combo_df")
GeneList<-
  readRDS("tmpData/XGboostGeneList.RDS")

dir.create("Result/Figure7",recursive = TRUE)
for(i in Select_Cohort){
  Select_IO_Dataset<-
    IO_All_Sets_Input[[i]]
  Select_IO_Dataset_SelectGene<-
    Select_IO_Dataset[
      c("ID","OS","status","response","BOR")
    ]
  Select_IO_Dataset_SelectGene_GSVA<-
    try(gsva(Select_IO_Dataset[c(1,(which(names(Select_IO_Dataset)=="BOR")+1):ncol(Select_IO_Dataset))]%>%
               FeaturetoRow()%>%t(),
             list(SelectGene=GeneList$V1))%>%
          t()%>%as.data.frame()%>%
          merge(.,Select_IO_Dataset_SelectGene,
                by.x=0,by.y="ID"))
  
  BestCut<-surv_cutpoint(Select_IO_Dataset_SelectGene_GSVA,
                         time = "OS",event = "status",variables = "SelectGene",
                         minprop = 0.2)%>%
    .$cutpoint%>%.[1,1]
  Select_IO_Dataset_SelectGene_GSVA$SelectGene_Status<-
    ifelse(Select_IO_Dataset_SelectGene_GSVA$SelectGene>BestCut,
           "ICD-High","ICD-Low")%>%factor(.,levels = c("ICD-Low","ICD-High"))
  CoxResult<-
    coxph(Surv(OS,status)~SelectGene_Status,
          Select_IO_Dataset_SelectGene_GSVA)%>%
    summary()
  Survfit<-
    survfit(Surv(OS,status)~SelectGene_Status,
            Select_IO_Dataset_SelectGene_GSVA)
  Plot<-
    ggsurvplot(Survfit,Select_IO_Dataset_SelectGene_GSVA,
               legend=c(0.7,0.9),
               pval = paste("HR = ",
                            CoxResult$coefficients[1,2]%>%round(.,digits = 2),
                            "\nLog-Rank Pvalue = ",
                            CoxResult$coefficients[1,5]%>%round(.,digits = 2)),
               legend.title="",
               risk.table = TRUE,
               legend.labs=c("ICD-Low Cluster","ICD-High Cluster"),
               pval.size = 5,
               font.x = c(15, "plain", "black"), # feature of curve xlab title
               font.y = c(15, "plain", "black"), # feature of curve ylab title
               font.tickslab = c(15, "plain", "black"),
               palette = c("#003C67FF",  "#A73030FF"),
               font.title=c(15, "plain", "black"))
  pdf(paste0("Result/Figure7/",i,"ICD_Signal_OS.pdf"),
      width = 6,height = 6)
  print(Plot)
  dev.off()
}
