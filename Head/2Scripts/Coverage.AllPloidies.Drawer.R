rm(list=ls(all=TRUE))

Cov = read.table("../../Body/2Derived/CoverageCnvForWindows.txt", header = TRUE)
Cov$coverage = log10(Cov$coverage)

pdf("../../Body/4Figures/Coverage.AllPloidies.Drawer.R01.pdf", height = 600, width = 400) # Dima - why this -? , onefile = FALSE)
a = 0
b = max(Cov$EndofWindow)
d = 0
f = max(Cov$coverage)  

## 329592 unrelated individuals => 1% is 329592/100 = 3295.92; log10(329592/100) = 3.517977;  0.1% = log10(329592/1000) = 2.517977; 0.01% = log10(329592/10000) = 1.517977

par(mfrow=c(22,1))
for (chr in 1:22)
{
  #### chr=1
  Windows = Cov[Cov$Chr == chr,]
  WindowsAllCHR = Windows[Windows$Ploidy == 'all',]
  Windows0CHR = Windows[Windows$Ploidy == 0,]
  Windows1CHR = Windows[Windows$Ploidy == 1,]
  Windows3CHR = Windows[Windows$Ploidy == 3,]
  #a = min(WindowsAllCHR$StartofWindow)
  #b = max(WindowsAllCHR$EndofWindow)
  # d = min(WindowsAllCHR$coverage)
  #d = 0
  #f = max(WindowsAllCHR$coverage)  
  if (chr <  22) {plot(NA, xlim=c(a,b), ylim=c(d,f), xlab='Location of Windows', ylab='Coverage', cex.axis = 10, xaxt='n')}
  if (chr == 22) 
  {
    plot(NA, xlim=c(a,b), ylim=c(d,f), xlab='Location of Windows', ylab='Coverage', cex.axis = 10)
    text(200000000,3, 'deletions', col = 'green', cex = 50)
    text(200000000,4, 'duplications', col = 'blue', cex = 50)
  }
  
  title(paste(chr), line = -30, cex.main = 50, adj = 0.01)
  # for (i in 1:nrow(WindowsAllCHR))   {segments(WindowsAllCHR$StartofWindow[i],WindowsAllCHR$coverage[i],WindowsAllCHR$EndofWindow[i],WindowsAllCHR$coverage[i],col = 'dark grey', lwd = 30)}
  for (i in 1:nrow(Windows0CHR))     {segments(Windows0CHR$StartofWindow[i],Windows0CHR$coverage[i],Windows0CHR$EndofWindow[i],Windows0CHR$coverage[i],col = 'red', lwd = 10)}
  for (i in 1:nrow(Windows1CHR))     {segments(Windows1CHR$StartofWindow[i],Windows1CHR$coverage[i],Windows1CHR$EndofWindow[i],Windows1CHR$coverage[i],col = 'green', lwd = 10)}
  for (i in 1:nrow(Windows3CHR))     {segments(Windows3CHR$StartofWindow[i],Windows3CHR$coverage[i],Windows3CHR$EndofWindow[i],Windows3CHR$coverage[i],col = 'blue', lwd = 10)}
  abline(h=3.517977, col = 'light grey')
  abline(h=2.517977, col = 'light grey')
  abline(h=1.517977, col = 'light grey')
}

### PLOT WindowsAll, A, B, C in separate pages for each chromosome

dev.off()
