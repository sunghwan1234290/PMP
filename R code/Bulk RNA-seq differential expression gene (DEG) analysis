library(edgeR)
library(gdata)

#setup variables
workDir = "./"
count.file = "./PMP_fresh_sample_RawCount.txt"
log2cpm.file = "./PMP_fresh_sample_log2CPM.txt"
result.file = "SigTable.Nomal vs PMP.txt"
result.file.bonf = "SigTable.bonf.Nomal vs PMP.txt"

p_cutoff = 0.01
fc_cutoff = 1.5

#Loading raw data
setwd(workDir)
raw.data <- read.delim(count.file)
log2cpm.data <- read.delim(log2cpm.file)

#target sample selection 
probeIds = raw.data[,1]
exprMtrx = raw.data[,-c(1)]

#divide samples
control.grp = exprMtrx[,c(1:7)] #nomal
test.grp= exprMtrx[,c(8:19)] #PMP

d <- cbind(control.grp, test.grp)
grouptag = c(rep(1, ncol(control.grp)),rep(2, ncol(test.grp)))

#get DGE list using edgeR
rownames(d) <- raw.data[, 1]
d <- DGEList(counts = d, group = grouptag, genes=raw.data[,1])
dim(d)

#Filtering very low expressed genes
#Since the smallest group size is three, we keep genes that achieve 
#at least one count per million (cpm) in at least three samples.
keep <- rowSums(cpm(d)>1) >= 3
d <- d[keep,]
dim(d)

#normalization
d <- calcNormFactors(d)

#producing MDS plot
plotMDS(d, main=paste("Nomal vs PMP"))

#Estimating the dispersion
d <- estimateCommonDisp(d, verbose=TRUE)
d <- estimateTagwiseDisp(d)
plotBCV(d)

# differentially expressed genes
et <- exactTest(d)
toptagList <- topTags(et, n=nrow(raw.data))$table
pValues = toptagList[,4]
fcValues = toptagList[,2]
fdrValues = toptagList[,5]

#Bonferroni correction
bonf.pValues = pValues * nrow(exprMtrx)
finalStatTable = cbind(toptagList, bonf.pValues)
sigList = finalStatTable[pValues<p_cutoff & (fcValues > log2(fc_cutoff) | fcValues < -log2(fc_cutoff)),]
sigList.bonf = finalStatTable[bonf.pValues<p_cutoff & (fcValues > log2(fc_cutoff) | fcValues < -log2(fc_cutoff)),]

#file writing
mergedData <- merge(sigList, log2cpm.data, by.x = 1, by.y = 1)
write.table(mergedData, result.file, sep="\t", row.names = F, quote = F)
mergedData <- merge(sigList.bonf, log2cpm.data, by.x = 1, by.y = 1)
write.table(mergedData, result.file.bonf, sep="\t", row.names = F, quote = F)
