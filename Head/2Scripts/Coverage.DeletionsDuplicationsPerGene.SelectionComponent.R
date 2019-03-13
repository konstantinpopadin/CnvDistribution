rm(list=ls(all=TRUE))

Genes = read.table("../../Body/3Results/Coverage.DeletionsDuplicationsPerGene.Table.txt", header = TRUE)

pdf("../../Body/4Figures/Coverage.DeletionsDuplicationsPerGene.SelectionComponent.R01.pdf")

# if there are no DupCov or DelCov - it means they are zero covered...
####### 1: deletions are more deleterious than duplications

VecOfParameters = c('KnKsMouse','GenomeWideHaploinsufficiencyScore','FunctionalIndispensabilityScore','ProbOfBeingLofIntolerant','ResidualVariationIntoleranceScore')
for (param in 1:length(VecOfParameters))
{  # param = 1
if (param == 1) {TEMP = data.frame(Genes$KnKsMouse,Genes$DelCov,Genes$DuplCov)}
if (param == 2) {TEMP = data.frame(Genes$GenomeWideHaploinsufficiencyScore,Genes$DelCov,Genes$DuplCov)}
if (param == 3) {TEMP = data.frame(Genes$FunctionalIndispensabilityScore,Genes$DelCov,Genes$DuplCov)}  
if (param == 4) {TEMP = data.frame(Genes$ProbOfBeingLofIntolerant,Genes$DelCov,Genes$DuplCov)}    
if (param == 5) {TEMP = data.frame(Genes$ResidualVariationIntoleranceScore,Genes$DelCov,Genes$DuplCov)}    

TEMP[,2:3][is.na(TEMP[,2:3] )] = 0 

DelVecRho=c(); DelVecP=c(); DuplVecRho=c(); DuplVecP=c()
for (i in 1:500)
  {
  temp = TEMP[sample(nrow(TEMP),5000),]

  Rho = as.numeric(cor.test(temp[,1],temp[,2], method = 'spearman')[4])
  P = as.numeric(cor.test(temp[,1],temp[,2], method = 'spearman')[3])
  DelVecRho = c(DelVecRho,Rho); DelVecP = c(DelVecP,P)

  Rho = as.numeric(cor.test(temp[,1],temp[,3], method = 'spearman')[4])
  P = as.numeric(cor.test(temp[,1],temp[,3], method = 'spearman')[3])
  DuplVecRho = c(DuplVecRho,Rho); DuplVecP = c(DuplVecP,P)
  }

if (param == 1) {KnKsMouseDelVecRho = DelVecRho; KnKsMouseDuplVecRho = DuplVecRho; }
if (param == 2) {GenomeWideHaploinsufficiencyScoreDelVecRho = DelVecRho; GenomeWideHaploinsufficiencyScoreDuplVecRho = DuplVecRho; }
if (param == 3) {FunctionalIndispensabilityScoreDelVecRho = DelVecRho; FunctionalIndispensabilityScoreDuplVecRho = DuplVecRho; }
if (param == 4) {ProbOfBeingLofIntolerantDelVecRho = DelVecRho; ProbOfBeingLofIntolerantDuplVecRho = DuplVecRho; }
if (param == 5) {ResidualVariationIntoleranceScoreDelVecRho = DelVecRho; ResidualVariationIntoleranceScoreDuplVecRho = DuplVecRho; }
}

boxplot(KnKsMouseDelVecRho,KnKsMouseDuplVecRho,GenomeWideHaploinsufficiencyScoreDelVecRho,GenomeWideHaploinsufficiencyScoreDuplVecRho,FunctionalIndispensabilityScoreDelVecRho,FunctionalIndispensabilityScoreDuplVecRho,ProbOfBeingLofIntolerantDelVecRho,ProbOfBeingLofIntolerantDuplVecRho,ResidualVariationIntoleranceScoreDelVecRho,ResidualVariationIntoleranceScoreDuplVecRho, col = c('green','blue'), outline = FALSE, notch = TRUE)
abline(h=0, col = 'red')

####### 2: selection and frequency - doesn't work for now !!! ???

dev.off()
