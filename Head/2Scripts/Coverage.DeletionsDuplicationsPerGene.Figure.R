rm(list=ls(all=TRUE))

Genes = read.table("../../Body/3Results/Coverage.DeletionsDuplicationsPerGene.Table.txt", header = TRUE)

pdf("../../Body/4Figures/Coverage.DeletionsDuplicationsPerGene.Figure.R01.pdf")
Genes = Genes[!is.na(Genes$DuplCov) & !is.na(Genes$DelCov),]
Genes = Genes[Genes$DuplCov > 0 & Genes$DelCov > 0,]  # at least one is higher than zero

XMin = -10 # log2(0.001)
XMax = max(log2(Genes$DelCov))
YMin = -10
YMax = max(log2(Genes$DuplCov))

par(mfrow=c(2,3))
plot(log2(Genes$DelCov),log2(Genes$DuplCov), pch = 16, xlim = c(XMin,XMax), ylim = c(YMin,YMax), col = rgb(0.1,0.1,0.1,0.1), xlab = 'deletions', ylab = 'duplications', main = 'all genes')
#a <- lm(log2(Genes$DuplCov) ~ log2(Genes$DelCov))
#abline(a, col = 'grey')
cor.test(log2(Genes$DelCov),log2(Genes$DuplCov), method = 'spearman')
summary(log2(Genes$DelCov))

plot(log2(Genes[Genes$GeneType == 'protein_coding',]$DelCov),log2(Genes[Genes$GeneType == 'protein_coding',]$DuplCov),  pch = 16, xlim = c(XMin,XMax), ylim = c(YMin,YMax), col = rgb(0.5,0.1,0.1,0.1), xlab = 'deletions', ylab = 'duplications', main = 'protein-coding genes')
#a <- lm(log2(Genes[Genes$GeneType == 'protein_coding',]$DelCov) ~ log2(Genes[Genes$GeneType == 'protein_coding',]$DuplCov))
#abline(a, col = 'red')
cor.test(log2(Genes[Genes$GeneType == 'protein_coding',]$DelCov),log2(Genes[Genes$GeneType == 'protein_coding',]$DuplCov), method = 'spearman')

GWHS = Genes[!is.na(Genes$GenomeWideHaploinsufficiencyScore) & Genes$GenomeWideHaploinsufficiencyScore > 0,]
plot(log2(GWHS[GWHS$GenomeWideHaploinsufficiencyScore > median(GWHS$GenomeWideHaploinsufficiencyScore),]$DelCov),log2(GWHS[GWHS$GenomeWideHaploinsufficiencyScore > median(GWHS$GenomeWideHaploinsufficiencyScore),]$DuplCov),  pch = 16, xlim = c(XMin,XMax), ylim = c(YMin,YMax), col = rgb(0.1,0.3,0.1,0.1), xlab = 'deletions', ylab = 'duplications', main = 'haploinsufficient genes')
#a <- lm(log2(GWHS[GWHS$GenomeWideHaploinsufficiencyScore > median(GWHS$GenomeWideHaploinsufficiencyScore),]$DelCov) ~ log2(GWHS[GWHS$GenomeWideHaploinsufficiencyScore > median(GWHS$GenomeWideHaploinsufficiencyScore),]$DuplCov))
#abline(a, col = 'green')
cor.test(log2(GWHS[GWHS$GenomeWideHaploinsufficiencyScore > median(GWHS$GenomeWideHaploinsufficiencyScore),]$DelCov),log2(GWHS[GWHS$GenomeWideHaploinsufficiencyScore > median(GWHS$GenomeWideHaploinsufficiencyScore),]$DuplCov), method = 'spearman')

dev.off()

