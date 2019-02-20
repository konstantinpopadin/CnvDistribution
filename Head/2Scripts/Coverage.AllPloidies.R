rm(list=ls(all=TRUE))

untar("../../Body/1Raw/002.AllCnvsUnrelatedEuropeans.tar.xz", exdir = "../../Body/1Raw/")
Carriers = read.table("../../Body/1Raw/002.AllCnvsUnrelatedEuropeans.txt", header = TRUE)
if (file.exists("../../Body/1Raw/002.AllCnvsUnrelatedEuropeans.txt")) file.remove("../../Body/1Raw/002.AllCnvsUnrelatedEuropeans.txt")

### number of individuals:
length(unique(Carriers$SampleNameShort)) # 329592 


### choose CNV class 
table(Carriers$Copy_Number)
# 0       1       3       4 
# 3214 2151574  157775    4077 

VecOfPloidies = c('all','0','1','3','4')
for (pl in 1:length(VecOfPloidies))
{ # pl = 2
  if (VecOfPloidies[pl] == 'all') {CNV = Carriers[abs(Carriers$Quality_Score) > 0.8,]}
  if (VecOfPloidies[pl] == 0)     {CNV = Carriers[abs(Carriers$Quality_Score) > 0.8 & Carriers$Copy_Number == 0,]}
  if (VecOfPloidies[pl] == 1)     {CNV = Carriers[abs(Carriers$Quality_Score) > 0.8 & Carriers$Copy_Number == 1,]}
  if (VecOfPloidies[pl] == 3)     {CNV = Carriers[abs(Carriers$Quality_Score) > 0.8 & Carriers$Copy_Number == 3,]}
  if (VecOfPloidies[pl] == 4)     {CNV = Carriers[abs(Carriers$Quality_Score) > 0.8 & Carriers$Copy_Number == 4,]} # there are no high quality CNVs with ploidy == 4
  
  ### calculate coverage for each window
  for (j in 1:22) # KP!!! 
  { # j = 1 
    Cnv = CNV[CNV$Chromosome == j,]
    if (nrow(Cnv) > 0)
    {
      VecOfUniqueBreakPoints = sort(unique(c(Cnv$Start_Position_bp,Cnv$End_Position_bp)))
      VecOfStart = VecOfUniqueBreakPoints[(-length(VecOfUniqueBreakPoints))] # delete last
      VecEnd = VecOfUniqueBreakPoints[-1]  # delete first
      DataWind = data.frame(VecOfStart,VecEnd)
      DataWind$Coverage = 0
      DataWind$CHR = j
      DataWind$Ploidy = VecOfPloidies[pl]
      for (i in 1 : nrow(DataWind))
      {
        Start = DataWind$VecOfStart[i]
        End = DataWind$VecEnd[i]
        coverage = nrow(Cnv[Cnv$Start_Position_bp <= Start & Cnv$End_Position_bp >= End,]) 
        DataWind$Coverage[i] = coverage
      }
      if (pl == 1 & j == 1) {DataDrivenWindows = DataWind}                           
      if (pl != 1 | j != 1) {DataDrivenWindows = rbind(DataDrivenWindows,DataWind)}  
    }
  }
}  

names(DataDrivenWindows) <- c('StartofWindow', 'EndofWindow', 'coverage', 'Chr', 'Ploidy')
table(DataDrivenWindows$Ploidy)
table(DataDrivenWindows$Chr)
write.table(DataDrivenWindows,"../../Body/2Derived/CoverageCnvForWindows.txt", row.names = FALSE, quote = FALSE)
