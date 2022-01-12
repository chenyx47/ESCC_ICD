selectAtype<-function(frame1,frame1select1,frame1select2,selectelement,frame2){
  #frame1[,frame1select1]<-gsub("-",".",frame1[,frame1select1])
  a<-frame1[,frame1select1][frame1[frame1select2]==selectelement]
  b<-intersect(as.character(a),names(frame2)[-1])
  c<-frame2[b]
  d<-cbind(frame2[1],c)
  return(d)
}
#单位标准化

#remove batch effect


#首列变行
FeaturetoRow<-function(a){
  b<-a
  row.names(b)<-a[,1]
  b<-b[-1]
  return(b)
}

#取交集基因函数
extractgene<-function(a,c){
  a<-FeaturetoRow(a)
  a1<-a[c,]
  X<-data.frame(X=row.names(a1))
  a1<-cbind(X,a1)
  return(a1)
}

#标准化

#提取方差top1000基因
extractVary<-function(a,top){
 a<-FeaturetoRow(a)[-1] 
 a.scale<-t(scale(t(a)))
 d<-apply(a.scale,1,sd)
 order<-order(d,decreasing = T)[1:top]
 total.matrixtop<-total.matrix.scale[order,]
 return(total.matrixtop)
}

#GSEA symbol transforming
gsea_visualized<-function(gseadata){
  gseadir<-read.delim("E:/data/gene set analysis/GSEA_Name_Map.txt")
  gseadirnew<-gseadir[gseadir$UP%in%gseadata[,1],]
  names(gseadata)[1]<-"UP"
  newgseadata<-merge(gseadirnew,gseadata,by="UP")
  return(newgseadata)
}
