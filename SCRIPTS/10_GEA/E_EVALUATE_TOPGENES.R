### Load topgenes and merge into single table

load('DATA/GEA_OUTPUT/KTAU/TOPGENES.rda')
TG_KTAU = TOPGENES
TG_KTAU$GEAT = 'KTAU'

load('DATA/GEA_OUTPUT/LFMM/TOPGENES.rda')
TG_LFMM = TOPGENES
TG_LFMM$GEAT = 'LFMM'

TOPGENES = rbind(TG_KTAU, TG_LFMM)

TOPGENES_U = do.call(rbind, by(TOPGENES, TOPGENES$locus, function(x) {return(x[which.min(x$pooled_q),])}))

### Check how many topgenes on unique genes, how many per env variable, how many per geat
length(unique(TOPGENES$locus))

table(TOPGENES$ENV)

median(TOPGENES_U$n_est)

sort(apply(TOPGENES_U[9:22], 2, function(x) { sum(is.na(x)==F)}))
