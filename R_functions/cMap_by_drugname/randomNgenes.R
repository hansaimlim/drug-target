random_N_genes_foreach_drug <-
function(N,outputdirectory){
x=1
colnum=dim(drugRL)[2]
while(x<=colnum){
randomindex=sample(1:11863,N)
drugname=colnames(drugRL)[x]
drugname=gsub("/","",drugname)
drugname=gsub("[()]","_",drugname)
drugname=gsub("[+]","p",drugname)
drugname=gsub("-_","m_",drugname)
col=drugRL[,x]
col.sorted=sort(col)
col.rand=col.sorted[randomindex]
row=c(col.rand)
write.table(row,file=paste(outputdirectory,drugname,".txt",sep=""),
quote=FALSE, sep="\t")
x<-x+1
}
}
