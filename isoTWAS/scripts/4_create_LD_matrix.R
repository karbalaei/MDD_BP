library('bigsnpr')
library(here)

address_1000G = "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/fusion_twas_LDREF_hg38/1000G.EUR."


for ( i in 1:22){
bed_file = paste0(address_1000G , i , ".bed")

data <- snp_readBed2 (bed_file , backingfile=tempfile())

genos = snp_attach(data)

G <- genos$genotypes

corr <- snp_cor(G)

rownames(corr) =  genos$map$marker.ID

colnames(corr) =  genos$map$marker.ID

Message(paste0("Saving Chr " , i ))

save(corr , file = here("isotwas" , "LD_matrix" , paste0("chromosome_" , i ,  ".rds")))

}