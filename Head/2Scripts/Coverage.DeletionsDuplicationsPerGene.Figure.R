rm(list=ls(all=TRUE))

Genes = read.table("../../Body/3Results/Coverage.DeletionsDuplicationsPerGene.Table.txt", header = TRUE)

pdf("../../Body/4Figures/Coverage.DeletionsDuplicationsPerGene.Figure.R01.pdf")
Genes = Genes[!is.na(Genes$DuplCov) & !is.na(Genes$DelCov),]

XMin = -10 # log2(0.001)
XMax = max(log2(Genes$DelCov))
YMin = -10
YMax = max(log2(Genes$DuplCov))

par(mfrow=c(2,2))
plot(log2(Genes$DelCov),log2(Genes$DuplCov), pch = 16, xlim = c(XMin,XMax), ylim = c(YMin,YMax), col = rgb(0.1,0.1,0.1,0.1))
cor.test(log2(Genes$DelCov),log2(Genes$DuplCov), method = 'spearman')
# par(new=TRUE)
plot(log2(Genes[Genes$GeneType == 'protein_coding',]$DelCov),log2(Genes[Genes$GeneType == 'protein_coding',]$DuplCov),  pch = 16, xlim = c(XMin,XMax), ylim = c(YMin,YMax), col = rgb(0.5,0.1,0.1,0.1))
cor.test(log2(Genes[Genes$GeneType == 'protein_coding',]$DelCov),log2(Genes[Genes$GeneType == 'protein_coding',]$DuplCov), method = 'spearman')

dev.off()