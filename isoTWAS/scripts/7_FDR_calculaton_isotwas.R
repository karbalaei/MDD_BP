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
suppressMessages(library('tidyverse'))
suppressMessages(library('SummarizedExperiment'))

option_list = list(
  make_option("--region", action="store", default="Amygdala,sACC", type='character',
              help="Brain region to be analysed [required]"),
  make_option("--diag", action="store", default="BP,MDD", type='character',
              help="Diagnons[required]")
  )

opt = parse_args(OptionParser(option_list=option_list))


set.seed(123)


#opt$region = "Amygdala"
#opt$diag = "MDD"
#results_list <- list.files(path=here("isotwas/Results") ,  pattern = "*.txt")

out_df <-  read.table(here("isotwas/Results" , paste0(opt$region , "_BurdenTest_results_" ,opt$diag , ".txt")), header = T , sep = "\t")%>%
           filter( !is.na(as.numeric(Z)))

#out_df <-  
gene =   out_df %>% 
   #left_join((out_df %>%
    #group_by(Gene) %>%
  #dplyr::summarize( n = n_distinct(Feature))) , by = "Gene") %>% 
    mutate(across((all_of(c(3:5,7))), ~ as.numeric(.))) %>%
  #mutate(P = if_else(P == 1, 0.999, if_else(P == 0 , 1e-200 , P)))%>%
   #mutate(P_Z = 2*pnorm(-abs(Z)))%>%
    group_by(Gene) %>%
  # mutate(Screen.P = ACAT::ACAT((P_Z)))
    mutate(Screen.P = isotwas::p_screen((P)))

alpha1=.05
G = nrow(gene)
gene$Screen.P.Adjusted = p.adjust(gene$Screen.P,method = 'fdr')
R = length(unique(gene$Gene[gene$Screen.P.Adjusted < alpha1]))
alpha2 = (R*alpha1)/G
head(gene)


isoform_new = as.data.frame(matrix(nrow = 0,
                                   ncol = ncol(out_df)+2))

colnames(isoform_new) = c(colnames(out_df),'Screen.P','Confirmation.P')
gene = gene[order(gene$Screen.P),]

isoform_new = gene %>%
  group_by(Gene) %>%
  reframe(Feature = Feature,
            Confirmation.P = isotwas::p_confirm(P,alpha = alpha2))

isoform_new = merge(isoform_new,gene,by=c('Gene','Feature'))

isoform_new$Confirmation.P = ifelse(isoform_new$Screen.P.Adjusted < 0.05,
                                    isoform_new$Confirmation.P,
                                    1)
isoform_new = isoform_new[,c('Gene','Feature','Z','P','permute.P',
                             'topSNP','topSNP.P',
                             'Screen.P','Screen.P.Adjusted','Confirmation.P')]


load(here("data/rse_gene.Rdata"))

bm <-  rowRanges(rse_gene) %>% as.data.frame() %>% dplyr::select(1:3) %>% rownames_to_column("Gene") %>%
  setNames(c("Gene" , "Chromosome" , "Start" , "End")) %>% dplyr::filter(!Chromosome %in% c("chrX", "chrY" , "chrM"))%>% 
  mutate(Chromosome =as.numeric(str_remove_all(Chromosome , "chr")))

isoform_new = merge(bm,isoform_new,by='Gene')

isoform_new = isoform_new[order(isoform_new$Chromosome,
                                isoform_new$Start,
                                decreasing = F),]
isoform_sig = subset(isoform_new,
                     Screen.P.Adjusted < alpha1 &
                       Confirmation.P < alpha2 &
                       permute.P < 0.05)


head(isoform_new)

head(isoform_sig)


message( paste("number of statistacally meaningfull isofroms is" , nrow(isoform_sig)))

write.csv(isoform_new, here("isotwas/Results" , paste0(opt$region , "_all_results_" ,opt$diag , ".csv")), row.names = F , quote=F)

write.csv(isoform_sig, here("isotwas/Results" , paste0(opt$region , "_FDR_filtered_results_" ,opt$diag , ".csv")), row.names = F , quote=F)



keep.isoform = c()
if (nrow(isoform_sig) > 1){
  for (i in 1:(nrow(isoform_sig)-1)){
    if (isoform_sig$End[i] > isoform_sig$Start[i+1] - 1e6){
      keep.isoform = unique(c(keep.isoform,
                              c(isoform_sig$Feature[c(i,i+1)])))
    }
  }
}

isoform_sig = subset(isoform_sig,Feature %in% keep.isoform)

message( paste("number of genes not needed finemapping  is" , nrow(isoform_sig)))

write.csv(isoform_sig, here("isotwas/Results" , paste0(opt$region , "_FDR_filtered_results_" ,opt$diag , "_no_need_finemapping.csv")), row.names = F , quote=F)




