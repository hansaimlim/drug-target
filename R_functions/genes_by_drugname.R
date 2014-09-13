get_all_genes_by_drugname <-
function(drug){
x=which(colnames(drugRL) == drug)
col=drugRL[,x]
col.sorted=sort(col)
row=c(colnames(drugRL)[x], col.sorted)
write.table(row,
file=paste("cMap_",drug,".txt",sep=""),
append=TRUE, quote=FALSE, sep="\t")
}
