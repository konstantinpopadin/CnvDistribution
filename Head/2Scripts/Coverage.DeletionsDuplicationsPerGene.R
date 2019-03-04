rm(list=ls(all=TRUE))

COV = read.table("../../Body/2Derived/CoverageCnvForWindows.txt", header = TRUE)
table(COV$Ploidy)

Genes = read.table("../../Body/1Raw/gencode.v25.annotation.gtf.Genes.Shet.pLI.FIS.RVIS.GHIS.KnKs.GC.BrainSpecificRanking.Branch", header = TRUE)
Genes = Genes[!Genes$GeneChr %in% c('X','Y','M'),]
# Genes = Genes[Genes$GeneType == 'protein_coding',]
Genes$DelCov = 0
Genes$DuplCov = 0

### count average coverage by deletions (pl = 1) and duplications (pl = 3) for each gene: 
for (i in 1:nrow(Genes))
{ # i = 10
  chr = Genes$GeneChr[i]
  start = Genes$GeneStart[i]
  end = Genes$GeneEnd[i]
  CovTemp = COV[COV$Chr == chr & COV$StartofWindow <= end & COV$EndofWindow >= start,]  
  
  # count DelCov and DupCov: 
  DEL = CovTemp[CovTemp$Ploidy == 1,]
  DEL$Length = DEL$EndofWindow - DEL$StartofWindow
  DEL$COV = DEL$Length*DEL$coverage
  DelCov = sum(DEL$COV) / sum(DEL$Length)
  
  DUPL = CovTemp[CovTemp$Ploidy == 3,]
  DUPL$Length = DUPL$EndofWindow - DUPL$StartofWindow
  DUPL$COV = DUPL$Length*DUPL$coverage
  DuplCov = sum(DUPL$COV) / sum(DUPL$Length)
  
  Genes$DelCov[i]  = DelCov
  Genes$DuplCov[i] = DuplCov
}

write.table(Genes, "../../Body/3Results/Coverage.DeletionsDuplicationsPerGene.Table.txt", row.names = FALSE, quote = FALSE)

