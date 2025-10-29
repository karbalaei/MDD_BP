suppressMessages(library("optparse"))
suppressMessages(library('glmnet'))
suppressMessages(library('methods'))
suppressMessages(library("readr"))
suppressMessages(library('dplyr'))
suppressMessages(library('stringr'))
suppressMessages(library('here'))
suppressMessages(library("jaffelab"))
suppressMessages(library('getopt'))
suppressMessages(library('tibble'))
suppressMessages(library('isotwas'))
suppressMessages(library('data.table'))


spec <- matrix(
  c(
    'region',
    'r',
    1,
    'character',
    'Either Amygdala or SACC',
    'gwas',
    'd',
    1,
    'character',
    'Either BP or MDD',    'help' ,
    'h',
    0,
    'logical',
    'Display help'
  ),
  byrow = TRUE,
  ncol = 5
)
opt <- getopt(spec)

set.seed(123)

set.seed(123)

out_df = data.frame(Gene = c(),
                    Feature = c(),
                    Z = c(),
                    P = c(),
                    permute.P = c(),
                    topSNP = c(),
                    topSNP.P = c())

source(here("isotwas","burdenTest.R"))


gwas <- fread(here("twas/hg38" , "bip2024_multianc_no23andMe_isotwas_hg38.txt"))

RDS_list = read.table(here("isotwas" , paste0(opt$region , "_RDS_list"))) %>% pull(1)

for ( i in 1:22){
  
  load(here("isotwas/LD_matrix" , paste0("chromosome_", i , ".rds")))
  
  gene_list = read.table(here("isotwas" , paste0(opt$region, "_transcript/" ,"gene_list_isotwas_" , opt$region , ".txt")) , header = T) %>%
    dplyr::filter(common_gene_id %in%  RDS_list) %>%
    dplyr::filter(seqnames == i ) %>% pull(2)
  
  message (paste("Working on Chromosome " , i  ))
  
  message (paste("There are " ,  length(gene_list)   , " genes on the list of this chromosome"))
  
  
  for (j in  1:length(gene_list)){
    
   load(here("isotwas/Results/" , opt$region , paste0(gene_list[j] , "_isoTWAS.RDS")))
    
    model = get(paste0(gene_list[j],"_isoTWAS")) %>% as.data.frame()
    
    message(paste("Working on Gene number " , j , "and its name is" ,   gene_list[j]  ))
    
    for (tx in unique(model$Feature)){
      
      sumstats.cur = subset(gwas,SNP %in% subset(model, Feature == tx)$SNP)
      tx_df = burdenTest(mod = subset(model, Feature == tx),
                                  ld = corr,
                                  gene = gene_list[j],
                                  sumStats = sumstats.cur,
                                  chr = 'Chromosome',
                                  pos = 'Position',
                                  a1 = 'A1',
                                  a2 = 'A2',
                                  a1_mod = 'ALT',
                                  a2_mod = 'REF',
                                  snpName = 'SNP',
                                  Z = 'Z',
                                  beta = NULL,
                                  se = NULL,
                                  featureName = 'Feature',
                                  R2cutoff = .01,
                                  alpha = 1e-3,
                                  nperms = 1000,
                                  usePos = F)
      out_df = rbind(out_df,tx_df)
      
    }
  }

}

save(out_df , file = here("isotwas" , "Results" , paste0(opt$region,"_BurdenTest_results_BP.RDS")))

write.table(out_df, file = here("isotwas" , "Results" , paste0(opt$region,"_BurdenTest_results_BP.txt")), sep = "\t", row.names = FALSE, quote = FALSE)




