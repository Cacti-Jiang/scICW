#!/usr/bin/Rscript
args = commandArgs(T)
if( length(args) != 2 &&  length(args) != 3 ){
	write("\nUsage:\n\tRscript exp.stat.R  <TPM.noERCC.mat>  [ReadsCount.ERCC.mat(Optional)]  <OUT.stat>\n", stdout())
	write("\t(If 'ReadsCount.ERCC.mat' is given, the script will also calculate Readscount, ERCC ReadsCount and ERCC ratio.)\n", stdout())
	q(status=1)
}

if ( length(args) == 2 ){
	TPM <- read.table(args[1], sep = "\t", header = T, row.names=1, check.names = F, stringsAsFactors = F)
	TPM_stat <- data.frame(	TPM0=(apply(TPM, 2, function(x) sum( x>0 ))),
							TPM1=(apply(TPM, 2, function(x) sum( x>1 )))  )
	write.table(TPM_stat, file=args[2], row.names = T, col.names = NA, sep = "\t", quote = FALSE, append = F)
}else{
	TPM <- read.table(args[1], sep = "\t", header = T, row.names=1, check.names = F, stringsAsFactors = F)
	RC  <- read.table(args[2], sep = "\t", header = T, row.names=1, check.names = F, stringsAsFactors = F)
	TPM_stat <- data.frame(	TPM0=(apply(TPM, 2, function(x) sum( x>0 ))),
							TPM1=(apply(TPM, 2, function(x) sum( x>1 )))  )
	
	sumr <- function(x, gtype = "ERCC"){
		if(gtype=="ERCC"){
			result <- sum(x[grep("ERCC", rownames(RC))])
		}else{
			result <- sum(x[grep("ERCC", rownames(RC), invert=T)])
		}
		# RSEM may generate non-int expected_count
		result = as.integer(result)
		return(result)
	}									
	RC_ercc <- data.frame( ReadsCountERCC = (apply(RC, 2, function(x) sumr(x, gtype="ERCC"))))
	RC_noercc <- data.frame( ReadsCount = (apply(RC, 2, function(x) sumr(x, gtype="NoERCC"))))
	RC_merge <- cbind(RC_noercc, RC_ercc)
	RC_merge <- RC_merge[rownames(TPM_stat), ]

	Merge_stat <- cbind(TPM_stat, RC_merge)
	Merge_stat$ERCCRatio <-  sprintf( "%.2f", round(Merge_stat$ReadsCountERCC*100/(Merge_stat$ReadsCountERCC + Merge_stat$ReadsCount), 4))
	write.table(Merge_stat, file=args[3], row.names = T, col.names = NA, sep = "\t", quote = FALSE, append = F)
}
