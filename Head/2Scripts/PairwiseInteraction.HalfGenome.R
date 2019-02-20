rm(list=ls(all=TRUE))

unzip("../../Body/1Raw/004.AllCnvsUnrelatedEuropeans.Overlap.GeneProperties.DeletionCarriers.zip", exdir = "../../Body/1Raw/")
Carriers = read.table("../../Body/1Raw/004.AllCnvsUnrelatedEuropeans.Overlap.GeneProperties.DeletionCarriers.txt", header = TRUE)
if (file.exists("../../Body/1Raw/004.AllCnvsUnrelatedEuropeans.Overlap.GeneProperties.DeletionCarriers.txt")) file.remove("../../Body/1Raw/004.AllCnvsUnrelatedEuropeans.Overlap.GeneProperties.DeletionCarriers.txt")

Cover = read.table("../../Body/1Raw/004.AllCnvsUnrelatedEuropeans.Overlap.GeneProperties.DeletionsChrStartEnd.CnvCoverage.txt", header = TRUE)
Cover$CnvId = paste(Cover$CHR,Cover$START,Cover$END, sep = '_')

Carriers = merge(Carriers,Cover, by = 'CnvId') # size is the same - good

### pdf
pdf("../../Body/4Figures/PairwiseInteraction.HalfGenome.01.pdf", height = 40, width = 20)

##### delete too common deletions
summary(Carriers$TotalCov)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.99    78.98  1290.96 13806.20 41529.40 44914.92 

# Carriers = Carriers[Carriers$TotalCov < 1290,] # medain

### loading of individuals by deletions
Contam = data.frame(table(Carriers$SampleNameShort))
names(Contam)=c('SampleName','Contam')
Contam=Contam[order(Contam$Contam, decreasing = TRUE),]
summary(Contam$Contam)  # median and 3rd Qu == 1! max = 16. So only rare individuals have > 1 CNV. 
# HighlyLoadedSamples = Contam[Contam$Contam > 8,]$SampleName; length(HighlyLoadedSamples)
# Carriers = Carriers[!Carriers$SampleNameShort %in% HighlyLoadedSamples,]

par(mfrow=c(5,3))
## go by decile to more and more rare:
for (i in (0:9))
{ # i = 1
  for (l in (1:3))
  {
    if (l == 1) {Carriers1 = Carriers[Carriers$TotalCov < quantile(Carriers$TotalCov,1-i/10),]}
    if (l == 2) {Carriers1 = Carriers[Carriers$TotalCov < quantile(Carriers$TotalCov,1-i/10),]; Carriers1 = Carriers1[Carriers1$Length_bp < median(Carriers1$Length_bp),]}
    if (l == 3) {Carriers1 = Carriers[Carriers$TotalCov < quantile(Carriers$TotalCov,1-i/10),]; Carriers1 = Carriers1[Carriers1$Length_bp >= median(Carriers1$Length_bp),]}
    
    AllCarriers  = length(unique(Carriers1$SampleNameShort));  AllCarriers # 115798
    CarriersFirst  = unique(Carriers1[Carriers1$Chromosome <= 10,]$SampleNameShort); CarriersFirstFr = length(CarriersFirst)/AllCarriers; CarriersFirstFr
    CarriersSecond  = unique(Carriers1[Carriers1$Chromosome > 10,]$SampleNameShort); CarriersSecondFr = length(CarriersSecond)/AllCarriers; CarriersSecondFr
    ObservedOverlap = length(intersect(CarriersFirst,CarriersSecond))/AllCarriers; ObservedOverlap 
    MedianLength =  median(Carriers1$Length_bp) 
    # observed overlap (0.111) is less than expected: (0.728*0.382=0.27) because our casesa are not independent = we don't have zero carriers... 
    # the correct and simple way to build a NULL expectation => to permut data.

    #### permutations (unlink chromosome and PatientId)
    PermVec = ObservedOverlap # I initialize with real value, just to be sure, that range (when I plot hist) will be broad enough
    for (perm in 1:500)
    {
    Carriers1$ChromosomePermut = sample(Carriers1$Chromosome)
    CarriersFirst  = unique(Carriers1[Carriers1$ChromosomePermut <= 10,]$SampleNameShort); CarriersFirstFr = length(CarriersFirst)/AllCarriers
    CarriersSecond  = unique(Carriers1[Carriers1$ChromosomePermut > 10,]$SampleNameShort); CarriersSecondFr = length(CarriersSecond)/AllCarriers
    PermVec = c(PermVec,length(intersect(CarriersFirst,CarriersSecond))/AllCarriers)
    }
    hist(PermVec, breaks = 100, main = paste(as.numeric(quantile(Carriers$TotalCov,1-i/10)), MedianLength, sep = ' '), col = 'grey')
    abline(v = ObservedOverlap,col = 'red', lwd = 5)
  }  
}

dev.off()







### batch by batch:
#BatchVec = unique(Carriers$Batch); length(BatchVec)
#par(mfrow=c(5,5))
#  for (batch in 1:length(BatchVec))
#    { # batch = 1
#      Temp = Carriers[Carriers$Batch == BatchVec[batch],]
#      AllCarriers  = length(unique(Temp$SampleNameShort)) # 115798
#      CarriersFirst  = unique(Temp[Temp$Chromosome <= 10,]$SampleNameShort); CarriersFirstFr = length(CarriersFirst)/AllCarriers
#      CarriersSecond  = unique(Temp[Temp$Chromosome > 10,]$SampleNameShort); CarriersSecondFr = length(CarriersSecond)/AllCarriers
#      ObservedOverlap = length(intersect(CarriersFirst,CarriersSecond))/AllCarriers; ObservedOverlap
#      
#      # permutations (unlink chromosome and PatientId)
#      PermVec = ObservedOverlap # I initialize with real value, just to be sure, that range (when I plot hist) will be broad enough
#      for (perm in 1:100)
#      {
#        Temp$ChromosomePermut = sample(Temp$Chromosome)
#        CarriersFirst  = unique(Temp[Temp$ChromosomePermut <= 10,]$SampleNameShort); CarriersFirstFr = length(CarriersFirst)/AllCarriers
#        CarriersSecond  = unique(Temp[Temp$ChromosomePermut > 10,]$SampleNameShort); CarriersSecondFr = length(CarriersSecond)/AllCarriers
#        PermVec = c(PermVec,length(intersect(CarriersFirst,CarriersSecond))/AllCarriers)
#      }
#      
#     hist(PermVec, breaks = 100)
#      abline(v = ObservedOverlap,col = 'red', lwd = 3)
#    }
#}




