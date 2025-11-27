######## Prepare genomic regions

### Load gtf file
GTF = read.table('DATA/RAW/corn_la/genomic.gtf', sep='\t')


### keep only genes
GENES = GTF[GTF$V3=='transcript',c(1,4,5)]




### Create chromsome id
chrom = 1:10
names(chrom) = sort(names(table(GENES$V1)))[3:12]

### convert chromosome id
GENES$V1 = chrom[GENES$V1]


## remove genes outside main chromosomes 
GENES=na.omit(GENES)


### write down genes table
write.table(GENES, file='DATA/PROCESSED_VCF/corn_la/GENES.tsv', sep='\t', row.names=F, quote=F, col.names=F)










########## Prepare VCF file

# load map and call files
MAP = revcfRMAP = read.table('DATA/RAW/corn_la/map_152k.txt')
CALLS = read.table('DATA/RAW/corn_la/snp_numeric_mvp.txt')

# create empty VCF output
VCF = data.frame('#CHROM'=MAP$CHROM, 'POS'=MAP$POS, 'ID'=MAP$ID, 'REF'='<A0>', 'ALT'='<A1>', 'QUAL'='.', 'FILTER'='.', 'INFO'='.', 'FORMAT'='GT')


## convert genotypes in vcf format
CALLS_VCF = CALLS

CALLS_VCF[CALLS==0] = '0/0'
CALLS_VCF[CALLS==1] = '0/1'
CALLS_VCF[CALLS==2] = '1/1'

## join calls to vcf
VCF = cbind(VCF, CALLS_VCF)

# copy header as first row, to keep #CHROM in header 
VCF = rbind(c('#CHROM', colnames(VCF)[-1]), VCF)

# write VCF
write.table(VCF, 'DATA/RAW/corn_la/map_152k.vcf', row.names=F, col.names=F, quote=F, sep='\t')

## add header
header = "##fileformat=VCFv4.0"
vcf <- readLines("DATA/RAW/corn_la/map_152k.vcf")
writeLines(c(header, vcf), "DATA/RAW/corn_la/map_152k.vcf")

