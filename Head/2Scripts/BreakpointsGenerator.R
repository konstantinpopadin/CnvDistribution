rm(list=ls(all=TRUE))

untar("../../Body/1Raw/002.AllCnvsUnrelatedEuropeans.tar.xz", exdir = "../../Body/1Raw/")
Carriers = read.table("../../Body/1Raw/002.AllCnvsUnrelatedEuropeans.txt", header = TRUE)
if (file.exists("../../Body/1Raw/002.AllCnvsUnrelatedEuropeans.txt")) file.remove("../../Body/1Raw/002.AllCnvsUnrelatedEuropeans.txt")

Carriers = Carriers[abs(Carriers$Quality_Score) > 0.8,]
Carriers = Carriers[grep("Chromosome|Start_Position_bp|End_Position_bp",colnames(Carriers))]
Carriers = unique(Carriers) # 44 621 rows
Carriers = Carriers[order(Carriers$Chromosome,Carriers$Start_Position_bp),]

write.table(Carriers,"../../Body/2Derived/UniqueBreakpointsOfAllHighQualityCNVs.txt", row.names = FALSE, quote = FALSE)
