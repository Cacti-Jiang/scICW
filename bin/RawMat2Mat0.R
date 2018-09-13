#!/usr/bin/Rscript
##Usage: Rscript RawMat2Mat0.R <raw.mat> <percentage> <out.prefix>
args = commandArgs(T)

raw <- read.table(args[1], header = T, check.names=F, stringsAsFactors = T)
pct <- as.numeric(args[2])

remain <- apply(raw, 1, function(x) sum(!(is.na(x)))> pct*ncol(raw))
raw <- raw[remain, ]
raw[is.na(raw)]<-0
raw <- raw[apply(raw, 1, function(x) sum(x==0)<=0.8*ncol(raw)), ]
#write.table(raw, paste0(args[3], ".tri.mat"), row.names = T, col.names = NA, sep = " ", quote = F)

raw[raw==2]<-1
raw <- raw[apply(raw, 1, function(x) sum(x==1)<=0.8*ncol(raw)), ]
write.table(raw, paste0(args[3], ".bi.mat"), row.names = T, col.names = NA, sep = " ", quote = F)

