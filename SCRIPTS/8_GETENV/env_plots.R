### Script to evaluate correlation between environmental variables and plot environmental variation over space

# Load libraries
library(corrplot)
library(rnaturalearth)

### Load dataset info
dslist = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt', header=T, sep='\t')
rownames(dslist) = dslist$dataset

### Load list of environmental variables
envVars = read.csv('DATA/ENV/envlist.csv')
rownames(envVars) = envVars$VariableID






# create container of environmental correlations across datasets
ENV.COR.DS = array(NA, dim=c(nrow(envVars),nrow(envVars),nrow(dslist)))
dimnames(ENV.COR.DS) = list(envVars$VariableID,envVars$VariableID,dslist$dataset)


# for every datasets: get correlation matrix
for (ds in dslist$dataset) {
  
  # load imputed env data  
  load(paste0('DATA/GEA_INPUT/meta_env_imputed/',ds,'.rda'))
  
  # calculate correlation
  COR.i = cor(iENV[,envVars$VariableID])
  diag(COR.i) = NA
  
  # add to container  
  ENV.COR.DS[,,ds] = COR.i
  
}


## Summarize correlation across all datasets
MEAN.COR = apply(ENV.COR.DS, 1:2, mean)
rownames(MEAN.COR) = envVars[rownames(MEAN.COR),'Description']
colnames(MEAN.COR) = envVars[colnames(MEAN.COR),'Description']


## Plot
png(filename = 'FIGURES/envVars/correlation.png', width = 10, height = 10, res = 500, units = 'in')
corrplot(MEAN.COR, na.label = ' ', axes=F)
dev.off()




#### Plot map for every environmental variable

## Create a container of LON, LAT and ENV variables for every descriptor across all datasets
ALL.ENV = data.frame()


# for every datasets: get lon, lat and env matrix
for (ds in dslist$dataset) {
  
  # load imputed env data  
  load(paste0('DATA/GEA_INPUT/meta_env_imputed/',ds,'.rda'))
  
  # add to container
  ALL.ENV = rbind(ALL.ENV, iENV)
}  

## Plot all environmental variables
land = ne_countries(scale='small') # get land shapefile

## create function to make plot
plot_env_geo = function(evars) {
  
  for (env in evars) {
    
    # set colors 
    COLCS = colorRampPalette(c('blue3','white','red3'))(10)
    
    # set color breaks
    evar = ALL.ENV[,env]
    COLBKS = c(min(evar)-0.000001, seq(quantile(evar, 0.1), quantile(evar, 0.9), length.out=9), max(evar)+0.000001)
    
    par(mar=c(.1,.1,1,.1))
    
    range(ALL.ENV$LON)
    
    plot(NA, xlim=c(-140, 144), ylim=c(-35,65),  col=0, axes=F) # plot canvas
    plot(land, add=T, col='lightgrey', border=NA) # add land  
    points(ALL.ENV$LON, ALL.ENV$LAT, col=COLCS[cut(evar,COLBKS)], pch=16, cex=0.5)  
    box()
    title(main=envVars[env,'Description'])
    legend('topleft', box.lwd = 0, bg=NA, title = envVars[env,'Unit'], legend = '', cex=0.8)
    legend('topleft', legend = c(signif(COLBKS[2],2), rep('', 8), signif(COLBKS[10],2)), pt.bg=COLCS, pt.lwd = 0, pt.cex = 3, pch=22, cex=0.8, box.lwd = 0, bg=NA, title = '')
    
  }
}


## Create outputs

png(filename = 'FIGURES/envVars/envplots1.png', width = 8, height = 11, res = 500, units = 'in',)
par(mfrow=c(7,2))
plot_env_geo(evars = envVars$VariableID[1:14])
dev.off()

png(filename = 'FIGURES/envVars/envplots2.png', width = 8, height = 11, res = 500, units = 'in',)
par(mfrow=c(7,2))
plot_env_geo(evars = envVars$VariableID[15:28])
dev.off()

png(filename = 'FIGURES/envVars/envplots3.png', width = 8, height = 11, res = 500, units = 'in',)
par(mfrow=c(7,2))
plot_env_geo(evars = envVars$VariableID[29:42])
dev.off()

png(filename = 'FIGURES/envVars/envplots4.png', width = 8, height = 11, res = 500, units = 'in',)
par(mfrow=c(7,2))
plot_env_geo(evars = envVars$VariableID[43:43])
dev.off()
