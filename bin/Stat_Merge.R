#!/usr/bin/Rscript
args = commandArgs(T)
if( length(args) < 2 ){
	write("\nUsage:\n\tRscript  Stat_Merge.R  <Stat1.txt>  <Stat2.txt>  <Stat3.txt>  ...  >OUT_merge_stat.txt\n", stdout())
	q(status=1)
}

out = read.table( args[1], header=T, row.names=1, stringsAsFactor=F, check.names=F );
for (i in 2:length(args)) {
	this = read.table( args[i], header=T, row.names=1, stringsAsFactor=F, check.names=F );
	this = this[rownames(out), ,drop=F];
	out = cbind(out,this);
}

write.table(out, stdout(), row.names = T, col.names = NA, sep = "\t", quote = FALSE, append = F)
