down_N_genes_foreach_drug <-
function(N,outputdirectory){
x=1
D=11863-N+1
colnum=dim(drugRL)[2]
while(x<=colnum){
drugname=colnames(drugRL)[x]
drugname=gsub("/","",drugname)
drugname=gsub("[()]","_",drugname)
drugname=gsub("[+]","p",drugname)
drugname=gsub("-_","m_",drugname)
col=drugRL[,x]
col.sorted=sort(col)
col.bottom=col.sorted[D:11863]
row=c(col.bottom)
write.table(row,file=paste(outputdirectory,drugname,".txt",sep=""),
quote=FALSE, sep="\t")
x<-x+1
}
}
