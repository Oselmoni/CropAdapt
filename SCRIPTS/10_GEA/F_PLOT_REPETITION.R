library(SNPRelate)
library(rnaturalearth)
library(terra)
library(sf)
library(visreg)
source('SCRIPTS/custom_R_functions.R')


### Load list of datasets
dslist = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt',header=T, sep='\t')$dataset
dslist = c(dslist[2:14], dslist[1])

### Get list of selected datasets
selectedDS = read.table('DATA/GEA_INPUT/GEA_selected_ds.txt', header=T, row.names=1, sep='\t')

### Load orthogroups list
load('DATA/SNP_ANNOTATION/orthogroups.rda')
rownames(orthogroups) = orthogroups$Orthogroup
### get land for plotting
land = ne_countries(scale='large')

### Load list of environmental variables
envVars = read.csv('DATA/ENV/envlist.csv')
rownames(envVars)= envVars$VariableID
envVars$ShortName =  gsub("\\\\n", "\n", envVars$ShortName) # maintain the line break

## Set plot colors
ds_cols = c('tomato_wr'='#f4cccc', 
            'barley_la'='#e6b8af', 'barley_wr'='#e6b8af',
            'commonbean_la'='#d0e0e3','commonbean_wr'='#d0e0e3',
            'corn_la'='#fff2cc','corn_wr'='#fff2cc',
            'rice_la'='#efefef','rice_wr'='#efefef',
            'sorghum_la'='#fce5cd',
            'soybean_la'='#d9ead3','soybean_wr'='#d9ead3',
            'teparybean_la'='#d9d2e9', 'teparybean_wr'='#d9d2e9'
)


### Load list of top candidate genes from LFMM picmin
load('DATA/GEA_OUTPUT/LFMM/TOPGENES.rda')
TOPGENES.LFMM = TOPGENES
TOPGENES.LFMM$GEA = 'LFMM'
### Load list of top candidate genes from KTAU picmin
load('DATA/GEA_OUTPUT/KTAU/TOPGENES.rda')
TOPGENES.KTAU = TOPGENES
TOPGENES.KTAU$GEA = 'KTAU'

### Create unique table of top genes
TOPGENES = rbind(TOPGENES.LFMM, TOPGENES.KTAU)

### Retain only best repetition per orthogroup
TOPGENES = do.call(rbind, by(TOPGENES, TOPGENES$locus, function(x) {x[which.min(x$pooled_q),]}))

### Sort by variable, then by significance
TOPGENES = TOPGENES[order(TOPGENES$ENVD, TOPGENES$pooled_q),]


### Extract GEA association for every orthogroup
GEA_OG=list()

for (og in TOPGENES$locus) {

  print(og)
  
  ## Isolate gene of interest
  tg = TOPGENES[TOPGENES$locus==og,]
  
  
  
  ####
  ###### Show Output Plots
  ####
  
  sortedPs = as.numeric(tg[,paste0('P(',dslist,')')])
  names(sortedPs) = dslist
  sortedPs = sort(sortedPs)
  top_ds = names(sortedPs[1:tg$n_est])

  
  for (ds in dslist) {
      
      
      if (is.na(tg[,paste0('P(',ds,')')])) { # if ortholog not genotyped: empty entry
          GEA_OG[[og]][[ds]][['sign']] = NA
        
      } else { # otherwise... retrive SNPs and env values
        
        
        ## retrieve WZA input & output
        load(paste0('DATA/GEA_OUTPUT/',GEAT,'/',ds,'/SNP_ANNOTATION_OL.rda'))
        load(paste0('DATA/GEA_OUTPUT/',GEAT,'/',ds,'/',tg$ENV,'.rda'))
        
        WZAin = SNP_ANNOTATION_OL
        WZAin$pval = WZA[[1]][WZAin$id]
        WZAout = WZA[[2]]
        rownames(WZAout) = WZAout$gene
        
        
        ## find SNPs on gene
        WZAin_gene = WZAin[WZAin$rnd_og==tg$locus,]
        
        
        ## find genotype variation
        ### Load GT matrix
        SNPS = snpgdsOpen(paste0('DATA/GEA_INPUT/GDS/gds_',ds,'.gds'), readonly = T, allow.duplicate = T)
        
        ### Load imputed GT matrix
        gti_files=list.files('DATA/GEA_INPUT/GTI/', pattern = paste0('GTI_',ds,'_K'),full.names = T)
        load(gti_files[regexpr(ds, gti_files)!=-1])
        
        
        ### Load environmental data
        load(paste0('DATA/GEA_INPUT/meta_env_imputed/',ds,'.rda'))
        
        ### Get GT and ENV data
        gt = GTI[,WZAin_gene$id,drop=F]
        env = iENV[,tg$ENV]
        
       
        ### Store in object
        GEA_OG[[og]][[ds]][['gt']] = gt
        GEA_OG[[og]][[ds]][['env']] = env
        GEA_OG[[og]][[ds]][['sign']] = ds%in%top_ds
        GEA_OG[[og]][[ds]][['envVar']] = tg$ENVD
        GEA_OG[[og]][[ds]][['gea']] = tg$GEA
        
        
        
      }
      
      
    } 
  }      


#### Plot repeated GEA results


### Create massive layout to include  all plots for every orthogroup vs. environmental variable
pdf(file = 'FIGURES/GEAs/REPETITIONS_PLOT.pdf', width = 6, height = 8 )
{
layout(matrix(c(rep(225, times=14),rep(227, times=4),
                rep(225, times=14),rep(227, times=4),
                rep(225, times=14),rep(227, times=4),
                rep(225, times=14),rep(227, times=4),
                1:14,rep(226, times=4),
                15:28,rep(226, times=4),
                29:42,rep(226, times=4),
                43:56,rep(226, times=4),
                57:70,rep(226, times=4),
                71:84,rep(226, times=4),
                85:98,rep(226, times=4),
                99:112,rep(226, times=4),
                113:126,rep(226, times=4),
                127:140,rep(226, times=4),
                141:154,rep(226, times=4),
                155:168,rep(226, times=4),
                169:182,rep(226, times=4),
                183:196,rep(226, times=4),
                197:210,rep(226, times=4),
                211:224,rep(226, times=4)
                ), nrow=20, byrow=T))
       
  
### Plot all plots of orthogroups SNP frequency v. environmental variable     
par(mar=c(0.1,0.1,0.1,0.1))
for (og in TOPGENES$locus) {
      for (ds in dslist) {
        
        if (is.na(GEA_OG[[og]][[ds]]$sign)) {
          
          plot(1,1,axes=F, xlab='', ylab='', col=NA) # if ortholog not genotyped: empty plot
          rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], col=ds_cols[ds], border=NA)
          rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], col=adjustcolor('white',0.5), border = NA)
          
        } else { # otherwise... retrive SNPs and env values
          
          ### Plot GEA
          gt = GEA_OG[[og]][[ds]]$gt
          env = GEA_OG[[og]][[ds]]$env
          
          # plot canvas
          plot(NA, ylim=c(0,1), xlim=range(env), axes=F)
          rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], col=ds_cols[ds], border = NA)

          
          # add regression line of every SNP
          for (g in 1:ncol(gt)) {
            
            ### find out if GEA is positive or negative
            GEA_sign = sign(cor(gt[,g],env, use='pairwise.complete.obs'))
            
            ### encode alleles as binary (p/a)
            af = (gt[,g]>=1)+0
            
            if (GEA_sign==-1) { af = (gt[,g]<=1)+0 }
            
            ### get regression line
            data = data.frame(af,env)
            mod = glm(af~env, data=data, family='binomial')  
            vr = visreg(mod, plot = F, scale='response')  

       
            if (GEA_OG[[og]][[ds]]$gea=='LFMM') {lines(vr$fit$env, vr$fit$visregFit, lwd=0.5)}
            if (GEA_OG[[og]][[ds]]$gea=='KTAU') {lines(vr$fit$env, vr$fit$visregFit, lwd=0.5, lty=2)}
            

          }
          
          ### If not significant contribution to repetition, add some white transparency
          if (GEA_OG[[og]][[ds]]$sign!=T) {rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], col=adjustcolor('white',0.5), border = NA)} 
          
          
         
          
          
        }
        
      
      } 
}      



### Plot count of repeated signal per crop type
COUNTS.OG = rep(0,14);names(COUNTS.OG)=dslist
for (i in 1:nrow(TOPGENES)) {
  tg_pvals = as.numeric(TOPGENES[i,paste0('P(',dslist,')')]) # get pvalues
  names(tg_pvals) = dslist
  
  # retain only top pvalues
  top_tg = sort(tg_pvals)[1:TOPGENES[i,'n_est']]
  
  COUNTS.OG[names(top_tg)] = COUNTS.OG[names(top_tg)] + 1
  
}
barplot(COUNTS.OG, col=ds_cols[dslist], width=0.5, space=1, xlim=c(1,length(COUNTS.OG))-0.25, names='', axes=F)
axis(2)

### Plot significance barplots for orthogroups
SIGN.OG = rev(-log(TOPGENES$pooled_q, base = 10))
barplot(SIGN.OG, width=0.5, space=1, ylim=c(1,length(SIGN.OG))-0.25, names='', horiz=TRUE, axes=F, col=0)
axis(1)
}
dev.off()






