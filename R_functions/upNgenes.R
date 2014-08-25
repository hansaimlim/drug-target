up_N_genes_foreach_drug <-
function(N,outputdirectory){
x=1
colnum=dim(drugRL)[2]
while(x<=colnum){
drugname=colnames(drugRL)[x]
drugname=gsub("/","",drugname)
drugname=gsub("[()]","_",drugname)
drugname=gsub("[+]","p",drugname)
drugname=gsub("-_","m_",drugname)
col=drugRL[,x]
col.sorted=sort(col)
col.top=col.sorted[1:N]
row=c(col.top)
write.table(row,file=paste(outputdirectory,drugname,".txt",sep=""),
quote=FALSE, sep="\t")
x<-x+1
}
}
