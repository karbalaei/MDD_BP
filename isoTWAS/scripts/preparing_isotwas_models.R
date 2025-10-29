suppressMessages(library("optparse"))
suppressMessages(library('plink2R'))
suppressMessages(library('methods'))
suppressMessages(library("readr"))
suppressMessages(library('dplyr'))
suppressMessages(library('stringr'))
suppressMessages(library('here'))
suppressMessages(library("jaffelab"))
suppressMessages(library('SummarizedExperiment'))

option_list = list(
  make_option("--region", action="store", default="Amygdala,sACC", type='character',
              help="Brain region to be analysed [required]"),
  make_option("--gene_index", action="store", default=NA, type='integer',
              help="Index of gene, from gene_list_isotwas file [required]"),
  make_option("--gene_id", action="store", default=NA, type='character',
              help="Ensemble gene id , from gene_list_isotwas file [required]") 
)

opt = parse_args(OptionParser(option_list=option_list))

set.seed(123)



message("Loading genotyping data")

plink_data <- read_plink(here(paste0("twas/" , opt$region, "_gene/" , "bim_files" ), paste0(opt$region , "_gene_" , opt$gene_index , "/filtered_snps_" , opt$region , "_gene_" , opt$gene_index)))

message("extracting genotyping matirx and removing NA")

snp_data <- plink_data$bim 

colnames(snp_data) = c("Chromosome","SNP" , "0" , "Position" , "ALT", "REF" )

load(file = here(paste0("isotwas/" , "Results/" , opt$region ) , paste0("Gene_" , opt$gene_index  , ".rds")))


gene_df = data.frame()

for ( i in 1 : length(model$Model) ){
df = model$Model[[i]]$Model %>% mutate(Feature = model$Model[[i]]$Transcript , R2 = model$Model[[i]]$R2 ,
                                       Gene = opt$gene_id ,Build = "hg38") 

  gene_df = rbind(df , gene_df)
  
}


#df_name <-  paste0(opt$gene_id , "_isotwas")
gene_df %>% left_join(snp_data[ , -3] , by = "SNP") -> df_ready

order_names = c("Feature" , "SNP" , "Chromosome" , "Position" , "Build" , "ALT" , "REF" , "Weight","R2" , "Gene")

df_ready = df_ready[ , order_names]
env <- parent.frame() # Get the current environment

assign(paste0(opt$gene_id , "_isotwas"), df_ready, envir = env) # Assign to a new name

  
save( list = paste0(opt$gene_id , "_isotwas") , file = here(paste0("isotwas/" , "Results/" , opt$region ) , paste0(opt$gene_id , "_isotwas" , ".rds")))



