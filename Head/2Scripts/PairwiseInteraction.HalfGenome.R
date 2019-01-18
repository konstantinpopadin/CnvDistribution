rm(list=ls(all=TRUE))

unzip("../../Body/1Raw/004.AllCnvsUnrelatedEuropeans.Overlap.GeneProperties.DeletionCarriers.zip", exdir = "../../Body/1Raw/")
Carriers = read.table("../../Body/1Raw/004.AllCnvsUnrelatedEuropeans.Overlap.GeneProperties.DeletionCarriers.txt", header = TRUE)
if (file.exists("../../Body/1Raw/004.AllCnvsUnrelatedEuropeans.Overlap.GeneProperties.DeletionCarriers.txt")) file.remove("../../Body/1Raw/004.AllCnvsUnrelatedEuropeans.Overlap.GeneProperties.DeletionCarriers.txt")

Cover = read.table("../../Body/1Raw/004.AllCnvsUnrelatedEuropeans.Overlap.GeneProperties.DeletionsChrStartEnd.CnvCoverage.txt", header = TRUE)
Cover$CnvId = paste(Cover$CHR,Cover$START,Cover$END, sep = '_')

Carriers = merge(Carriers,Cover, by = 'CnvId') # size is the same - good

### pdf

pdf("../../Body/4Figures/PairwiseInteraction.HalfGenome.01.pdf", height = 30, width = 30)

##### delete too common deletions
summary(Carriers$TotalCov)
Carriers = Carriers[Carriers$TotalCov < 1290,] # medain

##### delete too loaded individuals
Contam = data.frame(table(Carriers$SampleNameShort))
names(Contam)=c('SampleName','Contam')
Contam=Contam[order(Contam$Contam, decreasing = TRUE),]
HighlyLoadedSamples = Contam[Contam$Contam > 8,]$SampleName; length(HighlyLoadedSamples)
Carriers = Carriers[!Carriers$SampleNameShort %in% HighlyLoadedSamples,]

#### 
for (constraints in 1:3)
{
  if (constraints == 2) {Carriers = Carriers[Carriers$NumberOfProtCodGenes > 0,]}
  if (constraints == 3) {Carriers = Carriers[Carriers$NumberOfProtCodGenesWithHighSHet > 0,]}
AllCarriers  = length(unique(Carriers$SampleNameShort)) # 115798
CarriersFirst  = unique(Carriers[Carriers$Chromosome <= 10,]$SampleNameShort); CarriersFirstFr = length(CarriersFirst)/AllCarriers
CarriersSecond  = unique(Carriers[Carriers$Chromosome > 10,]$SampleNameShort); CarriersSecondFr = length(CarriersSecond)/AllCarriers
ObservedOverlap = length(intersect(CarriersFirst,CarriersSecond))/AllCarriers; ObservedOverlap

#### permutations (unlink chromosome and PatientId)
PermVec = ObservedOverlap # I initialize with real value, just to be sure, that range (when I plot hist) will be broad enough
for (perm in 1:1000)
{
Carriers$ChromosomePermut = sample(Carriers$Chromosome)
CarriersFirst  = unique(Carriers[Carriers$ChromosomePermut <= 10,]$SampleNameShort); CarriersFirstFr = length(CarriersFirst)/AllCarriers
CarriersSecond  = unique(Carriers[Carriers$ChromosomePermut > 10,]$SampleNameShort); CarriersSecondFr = length(CarriersSecond)/AllCarriers
PermVec = c(PermVec,length(intersect(CarriersFirst,CarriersSecond))/AllCarriers)
}
par(mfrow=c(1,1))
hist(PermVec, breaks = 100)
abline(v = ObservedOverlap,col = 'red', lwd = 3)

### batch by batch:
BatchVec = unique(Carriers$Batch); length(BatchVec)
par(mfrow=c(5,5))
  for (batch in 1:length(BatchVec))
    { # batch = 1
      Temp = Carriers[Carriers$Batch == BatchVec[batch],]
      AllCarriers  = length(unique(Temp$SampleNameShort)) # 115798
      CarriersFirst  = unique(Temp[Temp$Chromosome <= 10,]$SampleNameShort); CarriersFirstFr = length(CarriersFirst)/AllCarriers
      CarriersSecond  = unique(Temp[Temp$Chromosome > 10,]$SampleNameShort); CarriersSecondFr = length(CarriersSecond)/AllCarriers
      ObservedOverlap = length(intersect(CarriersFirst,CarriersSecond))/AllCarriers; ObservedOverlap
      
      # permutations (unlink chromosome and PatientId)
      PermVec = ObservedOverlap # I initialize with real value, just to be sure, that range (when I plot hist) will be broad enough
      for (perm in 1:100)
      {
        Temp$ChromosomePermut = sample(Temp$Chromosome)
        CarriersFirst  = unique(Temp[Temp$ChromosomePermut <= 10,]$SampleNameShort); CarriersFirstFr = length(CarriersFirst)/AllCarriers
        CarriersSecond  = unique(Temp[Temp$ChromosomePermut > 10,]$SampleNameShort); CarriersSecondFr = length(CarriersSecond)/AllCarriers
        PermVec = c(PermVec,length(intersect(CarriersFirst,CarriersSecond))/AllCarriers)
      }
      
      hist(PermVec, breaks = 100)
      abline(v = ObservedOverlap,col = 'red', lwd = 3)
    }
}
dev.off()



