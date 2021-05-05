require(Matrix)

###load normalized expression data###
load("starting_data/normdat1_tpm.rda")  ###genes in rows, cells in columns

###load cluster identities###
load("starting_data/clusters1.rda")    ###order of cells should be the same as in normdat1_tpm

###set display order for clusters###
cluster_order=c("Cluster1","Cluster2","Cluster3","Cluster4")
cluster_order_id=match(clusters1,cluster_order)

###generate annotation matrix###
annomat=data.frame(sample_name=names(clusters1),cluster_id=cluster_order_id,cluster_label=paste0(cluster_order_id,"_",clusters1),stringsAsFactors = F)

###save generated files - note that the expression matrix is re-saved, in case any changes are made to it###
save(annomat,normdat1_tpm,file="generated_files/anno_expr_mat.rda")

###calculate fractions of cells expressing each gene in each cluster###
#load("generated_files/anno_exp_mat.rda")
allcl=unique(annomat$cluster_label[order(annomat$cluster_id)])
fracmat=matrix(0,nrow=nrow(normdat1_tpm),ncol=length(allcl))
rownames(fracmat)=rownames(normdat1_tpm)
colnames(fracmat)=allcl
meanmat=fracmat
for (ii in allcl) {
  fracmat[,ii]=apply(normdat1_tpm[,annomat$sample_name[annomat$cluster_label==ii]]>0,1,mean)
  meanmat[,ii]=apply(normdat1_tpm[,annomat$sample_name[annomat$cluster_label==ii]],1,mean)
  print(c(ii," done"))
}
save(fracmat,file="generated_files/fracmat.rda")

