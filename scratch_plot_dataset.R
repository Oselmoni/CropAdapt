LA = WR  = data.frame()

lf  = list.files('DATA/META/')

for (fil in lf) {
    
  meta = read.csv(paste0('DATA/META/',fil))
  meta$SP =  strsplit(fil, '_')[[1]][1]
  
  if (strsplit(fil, '_')[[1]][2]=='la.csv') {
   
    LA = rbind(LA, meta)
    
  }  else {
    
    meta = read.csv(paste0('DATA/META/',fil))
    meta$SP =  strsplit(fil, '_')[[1]][1]
    
    WR = rbind(WR, meta)
    
  }
  
}

library(rnaturalearth)
land=countries110

par(mfrow=c(2,1))
par(mar=c(1,1,1,1))
plot(NA, xlim=c(-180,180), ylim=c(-50,+70), main='Landraces', axes=F)
plot(land, add=T, col='lightgrey', border=NA)
points(LA$LON, LA$LAT, col=as.factor(LA$SP), pch=16, cex=0.25)
legend('bottomleft', unique(LA$SP), col=1:length(unique(LA$SP)), pch=16, cex=0.5)

plot(NA, xlim=c(-180,180), ylim=c(-50,+70), main='Wild Relatives', axes=F)
plot(land, add=T, col='lightgrey', border=NA)
points(WR$LON, WR$LAT, col=as.factor(WR$SP), pch=16, cex=0.25)
legend('bottomleft', unique(WR$SP), col=1:length(unique(WR$SP)), pch=16, cex=0.5)



