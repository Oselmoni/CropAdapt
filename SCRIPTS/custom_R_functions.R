## custom function to rapidly get p-value from many kendall taus
kendall_pvalue <- function(T, n) {
  # Z-score under H0: tau = 0
  z <- T * sqrt(9 * n * (n - 1) / (2 * (2 * n + 5)))
  
  p <- 2 * (1 - pnorm(abs(z)))
  
  return(p)
}



## custom plot
cplot = function(x,y, col=1, pch=16, cex=1, main='', xlab='', ylab='', signif=2, adj=0) {
  
  
  plot(x,y, col=col, pch=pch, cex=cex,axes=F, xlab=NA, ylab=NA)
  
  # add axes
  axis(1, at=seq(min(x), max(x), length.out=3), labels=signif(seq(min(x), max(x), length.out=3),signif))
  axis(2, at=seq(min(y), max(y), length.out=3), labels=signif(seq(min(y), max(y), length.out=3),signif))
  
  # add labels
  title(xlab=xlab, ylab=ylab, line=2)
  
  # add title
  title(main=main, adj=adj)
  
  
}




# custom plot for GEA display 
plotGEAgeo = function(coord, col, env, eVar, main='', colBRK, df=100) {
  
  ### Load environmental raster
  if (  file.exists(paste0('DATA/ENV/ForPlotting/',envVars[eVar,"variableRaw"])) == F ) { # if layer doesn't exist... create for the first time, adjust CRS and resolution
    R = rast(paste0(paste0(envVars[eVar,'folder'],'/',envVars[eVar,"variableRaw"]))) # load raw raster
    R = project(R, y='epsg:4326') # reproject
    canvas = R;res(canvas) = c(0.5,0.5) # reduce spatial resolution
    R = resample(R, canvas)
    writeRaster(R, filename=paste0('DATA/ENV/ForPlotting/',envVars[eVar,"variableRaw"]))
  } else { # if layer is already created --> laad directly 
    R = rast(paste0('DATA/ENV/ForPlotting/',envVars[eVar,"variableRaw"]))
  }
  
  
  
  ### Set extent of plotting map area
  dX = diff(range(coord[,1]))
  dY = diff(range(coord[,2]))
  
  ### Set boundaries of plotted area, so that overall plot is a square
  if (dX<dY) {
    minY = min(coord[,2])-dY*0.05
    maxY = max(coord[,2])+dY*0.05
    
    delta = ((dY*1.1)-dX)/2
    minX = min(coord[,1])-delta
    maxX = max(coord[,1])+delta
    
  } else {
    minX = min(coord[,1])-dX*0.05
    maxX = max(coord[,1])+dX*0.05
    
    delta = ((dX*1.1)-dY)/2
    minY = min(coord[,2])-delta
    maxY = max(coord[,2])+delta
  }
  
  
  # crop area of interest
  R = crop(R, ext(c(minX, maxX, minY, maxY)))
  
  ## load topography for bg
  TOPO = rast('DATA/ENV/GMRT/GMRTv4_4_0_20251215topo.tif')
  TOPO = crop(TOPO, R)
  TOPO[TOPO<(-30)] = NA
  TOPO = terrain(TOPO, v='TRI')
  
  # resample
  R =resample(R, TOPO)
  
  
  ## rasterize land 
  LAND=rasterize(ne_countries(scale='large'), R)
  
  ## remove water pixel from raster
  R[is.na(LAND)] = NA
  TOPO[is.na(LAND)] = NA  
  
  
  ###
  ### plot background map
  ###
  par(mar=c(0,0,0,0))
  plot(NA, xlim=c(ceiling(minX),floor(maxX)), ylim=c(ceiling(minY),floor(maxY)), xaxs='i', yaxs='i', axes=F)
  plot(R, col=colorRampPalette(c('#7BD0F5','#FD7790'))(10), breaks=c(0,seq(quantile(env, 0.05, na.rm=T), quantile(env, 0.95, na.rm=T), length.out=9),Inf), add=T, legend=F)
  plot(TOPO, col=adjustcolor(colorRampPalette(c('grey90','grey20'))(20), 0.2), add=T)
  #plot(ne_coastline(scale = 'large'), add=T, col='grey30', border=NA)
  
  box()
  ###
  ### add points
  ###
  
  DT=(maxX-minX)/df
  
  
  rows = seq(minX, maxX, by=DT)
  cols = seq(minY, maxY, by=DT)
  
  availPOS = data.frame('LON'=rep(rows, each=length(cols)), 'LAT'=rep(cols, times=length(rows)))
  
  ### group together coordinates within the same distance
  geoCL = cutree(hclust(dist(coord[,c('LON','LAT')]), method = 'single'), h=DT)
  
  
  ## for every geographic cluster of points (from smallest to largest)
  for (geo in names(sort(table(geoCL)))) {
    
    
    coords_geo = coord[geoCL==geo,,drop=F]
    
    
    
    # find center of cluster
    centerLON = mean(coords_geo$LON)
    centerLAT = mean(coords_geo$LAT)
    
    
    # calculate distance from center
    availPOS$DC = sqrt((availPOS$LON-centerLON)^2+(availPOS$LAT-centerLAT)^2)
    
    
    
    # for every point, find a position from matrix
    for (i in 1:nrow(coords_geo)) {
      
      # find closest available point
      sel.pos = which.min(availPOS$DC)
      
      # add point coordinate
      coords_geo$LONplot[i] = availPOS$LON[sel.pos]
      coords_geo$LATplot[i] = availPOS$LAT[sel.pos]
      
      # removce chosen point
      availPOS = availPOS[-sel.pos,]
      
    }
    
    ### draw lines
    for (i in 1:nrow(coords_geo)) {
      lines(c(centerLON,coords_geo$LONplot[i]),
            c(centerLAT,coords_geo$LATplot[i]),
      )
    }
    
    ### draw points
    points(coords_geo$LONplot, coords_geo$LATplot, bg=COLMAF[which(geoCL==geo)], pch=21, lwd=.5, cex=1.5)
    
    
  }
  
  
  box()
  title(main=main, line=-2)
}










### custom manhattan plot
manhattanPlot = function(chr, p, pos, sig, main='', chrL, col='green3') {
  
  
  # set chromosome colors
  chrCol = rep(c('grey80','grey50'), length.out=length(unique(chr)))
  
  # calculate cumulative position
  chrN = as.numeric(as.factor(chr))
  cpos = pos[chrN==1]
  for (ch in 2:length(unique(chrN))) {
    
    cpos = c(cpos, max(cpos)+pos[chrN==ch])
    
  }
  
  ### Plot manhattan plot
  plot(cpos, p, pch=16, col=chrCol[chrN], axes=F, xlab='', ylab='', cex=0.5)
  
  # add circle to significant genes
  points(cpos[sig], p[sig], col=col, pch=21, lwd=2)
  
  # add chromosome labs
  axis(1, at=unlist(by(cpos, chr, function(x) {min(x)+(diff(range(x))/2)})), labels = chrL, lwd = 0, las=2, line=-1, chrCol, cex.axis=.5)
  axis(2, at=seq(par('usr')[3], par('usr')[4], length.out=4), labels = round(seq(par('usr')[3], par('usr')[4], length.out=4) ), cex.axis=0.75, lwd=0.5)
  title(ylab='-log(emp p-value)', line=2)
  
}



### Custom boxplot GEA
plotGEAbp = function(gt, env, envLab, COLCS) {
  
  
  # plot canvas
  plot(NA, xlim=c(-0.5,2.5), ylim=range(env,na.rm=T), axes=F, yaxs='i')
  
  # add colorscale in background
  #env_brks = seq(min(env), max(env), length.out=10)
  #for (i in 1:10) {  rect(par('usr')[1],  env_brks[i], par('usr')[2], env_brks[i+1] , border=NA, col=COLCS[i]) }
  
  #rect(par('usr')[1],  par('usr')[3], par('usr')[2], par('usr')[4] , border=NA, col=adjustcolor('white',0.2))
  
  # add gt distribution
  points(gt+runif(length(gt), -0.3,0.3), env, pch=16, col=adjustcolor(1,0.1), cex=.75)
  
  # add distribution values for every gt
  lines(c(0,0), quantile(env[gt==0], na.rm=T)[c(2,4)], lwd=2);points(0, median(env[gt==0], na.rm=T), cex=1.5, pch=21, bg=COLBOX[1])
  lines(c(1,1), quantile(env[gt==1], na.rm=T)[c(2,4)], lwd=2);points(1, median(env[gt==1], na.rm=T), cex=1.5, pch=21, bg=COLBOX[2])
  lines(c(2,2), quantile(env[gt==2], na.rm=T)[c(2,4)], lwd=2);points(2, median(env[gt==2], na.rm=T), cex=1.5, pch=21, bg=COLBOX[3])
  
  # Add axes
  axis(1, at=c(0,1,2), cex.axis=1, lwd=0.5)
  axis(2, at=seq(min(env, na.rm=T),max(env, na.rm=T), length.out=3), labels=signif(seq(min(env, na.rm=T),max(env, na.rm=T), length.out=3),3), cex.axis=0.85, lwd=0.5)
  box(lwd=0.5)
  title(ylab=envLab, xlab='Genotype', line=2)
}


### Stouffer p-value normalization
stP = function(ps) {
  z <- qnorm(1 - ps)
  Zg <- sum(z) / sqrt(length(z))
  p_gene <- 1 - pnorm(Zg) 
  return(p_gene) }




# custom plot for PCA display 
plotPCAgeo = function(coord, col, main='', df=100, mainT=T, cexP=1) {
  
  ### add color to df
  coord$col = col
  
  ### for display purposed, exclude samples located more than 5000 km away frm any other
  DIST = as.matrix(dist(coord[,c('LON','LAT')]))
  diag(DIST) = NA
  minDIST = apply(DIST, 1, min, na.rm=T)
  coord = coord[minDIST<50,]
  
  ### Set extent of plotting map area
  dX = diff(range(coord[,1]))
  dY = diff(range(coord[,2]))
  
  ### Set boundaries of plotted area, so that overall plot is a square
  if (dX<dY) {
    minY = min(coord[,2])-dY*0.05
    maxY = max(coord[,2])+dY*0.05
    
    delta = ((dY*1.1)-dX)/2
    minX = min(coord[,1])-delta
    maxX = max(coord[,1])+delta
    
  } else {
    minX = min(coord[,1])-dX*0.05
    maxX = max(coord[,1])+dX*0.05
    
    delta = ((dX*1.1)-dY)/2
    minY = min(coord[,2])-delta
    maxY = max(coord[,2])+delta
  }
  
  
  
  ## load topography for bg
  TOPO = rast('DATA/ENV/GMRT/GMRTv4_4_0_20251215topo.tif')
  TOPO = crop(TOPO, ext(c(minX, maxX, minY, maxY)))
  TOPO[TOPO<(-30)] = NA
  TOPO = terrain(TOPO, v='TRI')
  
  
  ## rasterize land 
  LAND=rasterize(ne_countries(scale='large'), TOPO)
  
  ## remove water pixel from raster
  TOPO[is.na(LAND)] = NA  
  
  ## set scalbar sizes
  max_range = max(c(dX,dY))
  sbsize =   1000 # by default, sb 1000 km
  if (max_range<20) { sbsize = 500} # if small areas, reduce sbsize
  if (max_range<10) { sbsize = 250} # if small areas, reduce sbsize
  if (max_range<5) { sbsize = 100}
  
  ###
  ### plot background map
  ###
  par(mar=c(1,1,3,1))
  plot(NA, xlim=c(ceiling(minX),floor(maxX)), ylim=c(ceiling(minY),floor(maxY)), xaxs='i', yaxs='i', axes=F)
  plot(ne_countries(scale = 'large'), col='grey80', border='NA', add=T)
  plot(TOPO, col=adjustcolor(colorRampPalette(c('grey80','grey20'))(20), 0.2), add=T, legend=F)
  sbar(sbsize/100, labels = paste0(sbsize, ' km'), lwd=1, xy=c(minX+diff(c(minX,maxX))*0.1, minY+diff(c(minY,maxY))*0.1))
  
  ###
  ### add points
  ###
  
  DT=(maxX-minX)/df
  
  
  rows = seq(minX, maxX, by=DT)
  cols = seq(minY, maxY, by=DT)
  
  availPOS = data.frame('LON'=rep(rows, each=length(cols)), 'LAT'=rep(cols, times=length(rows)))
  
  ### group together coordinates within the same distance
  geoCL = cutree(hclust(dist(coord[,c('LON','LAT')]), method = 'single'), h=DT)
  
  
  ## for every geographic cluster of points (from smallest to largest)
  for (geo in names(sort(table(geoCL)))) {
    
    
    coords_geo = coord[geoCL==geo,,drop=F]
    
    
    
    # find center of cluster
    centerLON = mean(coords_geo$LON)
    centerLAT = mean(coords_geo$LAT)
    
    
    # calculate distance from center
    availPOS$DC = sqrt((availPOS$LON-centerLON)^2+(availPOS$LAT-centerLAT)^2)
    
    
    
    # for every point, find a position from matrix
    for (i in 1:nrow(coords_geo)) {
      
      # find closest available point
      sel.pos = which.min(availPOS$DC)
      
      # add point coordinate
      coords_geo$LONplot[i] = availPOS$LON[sel.pos]
      coords_geo$LATplot[i] = availPOS$LAT[sel.pos]
      
      # removce chosen point
      availPOS = availPOS[-sel.pos,]
      
    }
    
    ### draw lines
    for (i in 1:nrow(coords_geo)) {
      lines(c(centerLON,coords_geo$LONplot[i]),
            c(centerLAT,coords_geo$LATplot[i]),
      )
    }
    
    ### draw points
    points(coords_geo$LONplot, coords_geo$LATplot, bg=coords_geo$col, pch=21, lwd=.5, cex=cexP)
    
    
  }
  box()
  if (mainT==T) {title(main='     E) Geographic distribution', line=0.5, adj=0)}
  
  ### if small area, add reference map
  if (max(c(dX,dY))<200) {
    
    plot(NA, xlim=c(-180,+180), ylim=c(-90,90), xlab='', ylab='', axes=F)
    rect(-180, -90,180,90, col=adjustcolor('white', 0.8))
    plot(land, col=adjustcolor('grey80', 0.8), border=NA, add=T, axes=F)
    rect(minX, max(c(-90,minY)), maxX, min(c(maxY,90)), border='red', col=NA)
    
  }
  
  
}


### add scalebar

scalebar <- function(x, y, length_km) {
  
  #Compute how many degrees of longitude = length_km at this latitude
  # using the haversine formula, solved by trial (or direct approximation)
  R <- 6371  # Earth's radius in km
  
  # distance for 1 degree of longitude at given latitude
  deg2rad <- function(deg) deg * pi / 180
  lat_rad <- deg2rad(y)
  
  # haversine distance between (x, lat) and (x + 1, lat)
  hav_dist <- function(lon1, lon2, lat) {
    lon1 <- deg2rad(lon1); lon2 <- deg2rad(lon2); lat <- deg2rad(lat)
    dlon <- lon2 - lon1
    a <- cos(lat)^2 * sin(dlon / 2)^2
    2 * R * asin(sqrt(a))
  }
  
  km_per_deg_lon <- hav_dist(0, 1, y)
  length_deg <- length_km / km_per_deg_lon
  
  # draw the bar
  if (is.na(length_km)==F) {
  lines(c(x-length_deg/2, x+length_deg/2), c(y,y))
  text(x, y, pos=1, paste0(length_km, ' km'), cex=0.75)
  }
}
