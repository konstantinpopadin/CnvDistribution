rm(list=ls(all=TRUE))

untar("../../Body/1Raw/europeanID_pca.tar.xz", exdir = "../../Body/1Raw/")
PCA = read.table("../../Body/1Raw/europeanID_pca.txt", header = FALSE)
if (file.exists("../../Body/1Raw/europeanID_pca.txt")) file.remove("../../Body/1Raw/europeanID_pca.txt")

PCA = PCA[,c(1:3)] #PCA = PCA[,c(1,40,41)]
names(PCA)=c('SampleName','PC1','PC2')
nrow(PCA) # 409702
AllVec = PCA$SampleName

### 
Related = read.table("../../Body/1Raw/id.related", header = FALSE)
nrow(Related) # 39305
RelatedVec = Related$V1

###
length(intersect(AllVec,RelatedVec)) # just 30152 != 39305

pdf("../../Body/4Figures/TouchingOfPCA.01.R.pdf")
plot(PCA$PC1,PCA$PC2, col = rgb(0.1,0.1,0.1,0.1), main = 'All')
plot(PCA[PCA$SampleName %in% RelatedVec,]$PC1,PCA[PCA$SampleName %in% RelatedVec,]$PC2, col = rgb(0.1,0.1,0.1,0.1), main = 'related only')
plot(PCA[!PCA$SampleName %in% RelatedVec,]$PC1,PCA[!PCA$SampleName %in% RelatedVec,]$PC2, col = rgb(0.1,0.1,0.1,0.1), main = 'ALL minus related')
par(mfrow=c(2,1))
hist(PCA[!PCA$SampleName %in% RelatedVec,]$PC1, breaks = 1000)
hist(PCA[!PCA$SampleName %in% RelatedVec,]$PC2, breaks = 1000)
dev.off()










