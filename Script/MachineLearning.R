##Machine Learning----
dir.create("Result/Figure5")
library(xgboost)

CellDeath_GeneList_List_Integrate<-
  readRDS("tmpData/CellDeath_GeneList_List_Integrate.RDS")

# exp<-tumor_filter$exp_tumor_merge_remove
# pheno<-ESCC_CC_CD_GEO_GSVA_3Path
# gene<-CellDeath_GeneList_List_Integrate[c(1,4,5)]
# gene<-unlist(gene)
# gene<-unique(gene)
# 
# exp_filter<-exp[rownames(exp)%in%gene,]
# exp_filter<-t(exp_filter)
# label_train=c(pheno$.-1)
# traindata <- list(data=exp_filter,label=label_train)
# dtrain <- xgb.DMatrix(data = traindata$data, label = traindata$label)
# numberOfClasses <- length(unique(traindata$label))
# names <-  colnames(dtrain)
# bst_model <- xgb.train(data = dtrain,
#                        nrounds = 1000,
#                        print_every_n = 10,
#                        verbose = TRUE,)
# importance_matrix= xgb.importance( model = bst_model)
# save(importance_matrix,file="./tmpData/XGboost_feature_importance.Rdata")
load()
pdf("Result/Figure5/Figure5A.pdf")
xgb.ggplot.importance(importance_matrix)
dev.off()
