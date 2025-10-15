### load metadata
meta = read.csv('DATA/META_RAW/miller_wr/Table1_main_text.csv')




### save output
write.csv(meta, 'DATA/META/millet_wr.csv', row.names = F)

