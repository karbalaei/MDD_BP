suppressMessages(library("optparse"))
suppressMessages(library('plink2R'))
suppressMessages(library('glmnet'))
suppressMessages(library('methods'))
suppressMessages(library("readr"))
suppressMessages(library('dplyr'))
suppressMessages(library('stringr'))
suppressMessages(library('here'))
suppressMessages(library("jaffelab"))
suppressMessages(library('sessioninfo'))
suppressMessages(library('tibble'))
suppressMessages(library('isotwas'))

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


message("Loading expression matirx")

load(file = here("isotwas", paste0(opt$region, "_transcript/expression_matrix.rda")))

message("Loading gene transcript table")

gene_tx_list = read.table(file = here("isotwas", paste0(opt$region, "_transcript/gene_transcript_" , opt$region , ".txt")), header=T)

message("extracting transcripts")

tx_list = gene_tx_list %>% dplyr::filter(common_gene_id == opt$gene_id ) %>% pull(1)

isoform_mat = expression_clean %>% as.data.frame %>%  dplyr::select(starts_with (tx_list)) %>% data.matrix()

if (ncol(isoform_mat) < 2) {

message(" Isoforms are less than two. Process stopped!")

} else {

message("Loading genotyping data")

plink_data <- read_plink(here(paste0("twas/" , opt$region, "_gene/" , "bim_files" ), paste0(opt$region , "_gene_" , opt$gene_index , "/filtered_snps_" , opt$region , "_gene_" , opt$gene_index)))

message("extracting genotyping matirx and snp data and removing NA")

snp_mat <- plink_data$bed

snp_mat =  snp_mat[ , colSums(is.na(snp_mat)) == 0]

snp_data <- plink_data$bim 

colnames(snp_data) = c("Chromosome","SNP" , "0" , "Position" , "ALT", "REF" )

message("preparing bootstrap data")

boot_list = lapply(1:10,function(k){jitter(isoform_mat)})
isoform_mat_rep = rlist::list.rbind(boot_list)
rownames(isoform_mat_rep) = rep(rownames(isoform_mat),10)
colnames(isoform_mat_rep) =colnames(isoform_mat)
dim(isoform_mat_rep)

gene_df = data.frame()

message("running compute_isotwas function")

isotwas_model = compute_isotwas(X = snp_mat, 
                                Y = isoform_mat, 
                                Y.rep = isoform_mat_rep,
                                R = 10, 
                                id = rownames(isoform_mat_rep), 
                                omega_est = 'replicates', 
                                omega_nlambda = 5, 
                                method = c('mrce_lasso', 
                                           'multi_enet', 
                                           'univariate',
                                           'joinet',
                                           'spls'),
                                predict_nlambda = 10, 
                                family = 'gaussian', 
                                scale = FALSE, 
                                alpha = 0.5, 
                                nfolds = 5, 
                                verbose = TRUE, 
                                tx_names = tx_list, 
                                seed = 1789, 
                                run_all = FALSE, 
                                return_all = TRUE)
								
model = isotwas_model$isotwas_mod

saveRDS(model , file = here(paste0("isotwas/" , "Results/" , opt$region ) , paste0(opt$gene_id ,    "_raw_model.RDS")))

#save(model , file = here(paste0("isotwas/" , "Results/" , opt$region ) , paste0("Gene_" , opt$gene_index  , ".rds")))


#isotwas_model$isotwas_mod$R2 %>% mutate(Gene_name = opt$gene_id) %>% relocate(Gene_name) %>%
#write.table(file = here(paste0("isotwas/" , "Results/" , opt$region ) , paste0("Gene_" , opt$gene_index  , ".txt")) , sep = "\t" , row.names = F , quote = F)

for ( i in 1 : length(model$Model) ){
df = model$Model[[i]]$Model %>% mutate(Feature = model$Model[[i]]$Transcript , R2 = model$Model[[i]]$R2 , R2.P = model$Model[[i]]$P ,
                                       Gene = opt$gene_id ,Build = "hg38") 

  gene_df = rbind(df , gene_df)
  
}


#df_name <-  paste0(opt$gene_id , "_isotwas")
gene_df %>% left_join(snp_data[ , -3] , by = "SNP") -> df_ready

order_names = c("Feature" , "SNP" , "Chromosome" , "Position" , "Build" , "ALT" , "REF" , "Weight","R2" , "Gene")

df_ready = df_ready[ , order_names] %>% as.data.frame()

env <- parent.frame() # Get the current environment

assign(paste0(opt$gene_id , "_isoTWAS"), df_ready, envir = env) # Assign to a new name

  
save( list = paste0(opt$gene_id , "_isoTWAS") , file = here(paste0("isotwas/" , "Results/" , opt$region ) , paste0(opt$gene_id , "_isoTWAS" , ".RDS")))

}
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


