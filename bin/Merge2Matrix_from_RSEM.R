#!/usr/bin/Rscript
args = commandArgs(T)
if(length(args) != 4 ){
	write("\nUsage:",stdout())
	write("\tRscript  Merge2Matrix_from_StringTie.R  <IN_abundance.list (ID Path)>  <OUT_ReadCount.matrix.txt>  <OUT_TPM.matrix.txt>  <OUT_FPKM.matrix.txt>",stdout())
	write("\nAuthor & Date:",stdout())
	write("\t\tqinshishang@genomicsc.cn, 2018-03-15",stdout())
        write("\n",stdout())
	q(status=1)
}

List = read.table(args[1], header = F, stringsAsFactors = F)
In = read.table( List[1,2], header = T, stringsAsFactors = F, sep = "\t")

rownames(In) <- In[, 1]
In = In[, -1]

all_rowname = rownames(In)
all_colname = c(List[1,1])
ReadCount = as.matrix( In[all_rowname,"expected_count"] )
TPM = as.matrix( In[all_rowname,"TPM"] )
FPKM = as.matrix( In[all_rowname,"FPKM"] )

for (i in 2:nrow(List)) {
	all_colname = c( all_colname, List[i,1] )
	In = read.table( List[i,2], header = T, stringsAsFactors = F, sep = "\t")
	rownames(In) <- In[, 1]
	In = In[, -1]
	ReadCount_tmp = as.matrix( In[all_rowname,"expected_count"] )
	ReadCount = cbind(ReadCount, ReadCount_tmp)
        TPM_tmp = as.matrix( In[all_rowname,"TPM"] )
        TPM = cbind (TPM, TPM_tmp)
	FPKM_tmp = as.matrix( In[all_rowname,"FPKM"] )
	FPKM = cbind ( FPKM, FPKM_tmp)
}

rownames(ReadCount) = all_rowname
colnames(ReadCount) = all_colname
rownames(TPM) = all_rowname
colnames(TPM) = all_colname
rownames(FPKM) = all_rowname
colnames(FPKM) = all_colname

write.table(ReadCount,  args[2], append = F, quote = F, sep = "\t", col.names=NA)
write.table(TPM,  args[3], append = F, quote = F, sep = "\t", col.names=NA)
write.table(FPKM, args[4], append = F, quote = F, sep = "\t", col.names=NA)

