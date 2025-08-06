 ### Compare BPseq and DPseq models results #####
rm(list= ls())

library(here)
library(readr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(dplyr)
library(tibble)


here()

load(here("results" ,  "BP_MDD_All_results.RDS"))


scale_fill_Reza <- function(...){
  ggplot2:::manual_scale(
    'color', 
    values = setNames(c('#e41a1c', '#377eb8', '#984ea3', '#4daf4a'), c("Both Region", "Amygdala","sACC", "None" )), 
    ...
  )
}


Volcano_graphs <-  function(i , type , cutoff , xlim1 , xlim2 , ylim1 ) {
  
	 Amygdala = Outlist[[i]][[2]] %>% rownames_to_column("Id") 
         
         sACC = Outlist[[i]][[1]] %>% rownames_to_column("Id")

        df = merge(Amygdala, sACC , by = "Id")
  
  
 df_V_Amygdala <-   df%>%
    dplyr::mutate("DEp" = if_else( q_PrimaryDxBipolar.x < 0.05 , "Amygdala" , "None" ) , 
   Expression = case_when(PrimaryDxBipolar.x > 0  & q_PrimaryDxBipolar.x <0.05~ "Up-regulated",
                          PrimaryDxBipolar.x < 0  & q_PrimaryDxBipolar.x <0.05~ "Down-regulated",
                          TRUE ~ "Unchanged"))
  
 top_genes_Amygdala <- bind_rows(
   df_V_Amygdala %>% 
    dplyr::filter (Expression == 'Up-regulated') %>% 
     dplyr::arrange(q_PrimaryDxBipolar.x, desc(abs(PrimaryDxBipolar.x))) %>% 
     head(topgenes , n = 5),
   df_V_Amygdala %>% 
     dplyr::filter(Expression == 'Down-regulated') %>% 
     dplyr::arrange(q_PrimaryDxBipolar.x, desc(abs(PrimaryDxBipolar.x))) %>% 
     head(topgenes , n = 5)
 )
 
 
 plot1 = ggplot(df_V_Amygdala, aes(PrimaryDxBipolar.x, -log(q_PrimaryDxBipolar.x,10))) +
   geom_point(aes(color = Expression), size = 1) +
   xlab(expression("log"[2]*"FC")) + 
   ylab(expression("-log"[10]*"p-value")) +
   scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
   guides(colour = guide_legend(override.aes = list(size=1.5))) +
   geom_label_repel(data = top_genes_Amygdala,
                    mapping = aes(label = common_gene_id.x),
                    size = 2) + theme_classic()+
   theme(legend.position= "none") +
   xlim(xlim1 , xlim2) + ggtitle( paste0("Amygdala " , type)) +
   ylim(0 , ylim1)
   
 ggsave(plot = plot1 , filename = here("graphs", "volcano_plot" ,  paste0("Amygdala region " , type ,  " .jpeg")) , width = 8 , height = 8 )
 
 # pdf(here("graphs", "volcano_plot" ,  paste0("Amygdala region " , type ,  " .pdf")) , width = 8 , height = 8 )
 # 
 # Volcano_Amygdala
 # 
 # dev.off()
 
 
 
 df_V_sACC <-   df%>%
   dplyr::mutate("DEp" = if_else( q_PrimaryDxBipolar.y < 0.05 , "sACC" , "None" ) , 
                 Expression = case_when(PrimaryDxBipolar.y > 0  & q_PrimaryDxBipolar.y <0.05~ "Up-regulated",
                                        PrimaryDxBipolar.y < 0  & q_PrimaryDxBipolar.y <0.05~ "Down-regulated",
                                        TRUE ~ "Unchanged"))
 
 top_genes_sACC <- bind_rows(
   df_V_sACC %>% 
     dplyr::filter (Expression == 'Up-regulated') %>% 
     dplyr::arrange(q_PrimaryDxBipolar.y, desc(abs(PrimaryDxBipolar.y))) %>% 
     head(topgenes , n = 5),
   df_V_sACC %>% 
     dplyr::filter(Expression == 'Down-regulated') %>% 
     dplyr::arrange(q_PrimaryDxBipolar.y, desc(abs(PrimaryDxBipolar.y))) %>% 
     head(topgenes , n = 5)
 )
 
 
  plot2 = ggplot(df_V_sACC, aes(PrimaryDxBipolar.y, -log(q_PrimaryDxBipolar.y,10))) +
   geom_point(aes(color = Expression), size = 1) +
   xlab(expression("log"[2]*"FC")) + 
   ylab(expression("-log"[10]*"p-value")) +
   scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
   guides(colour = guide_legend(override.aes = list(size=1.5))) +
   geom_label_repel(data = top_genes_sACC,
                    mapping = aes(label = common_gene_id.y),
                    size = 2) + theme_classic()+
   theme(legend.position= "none") +
   xlim(xlim1 , xlim2) + ggtitle( paste0("sACC " , type)) +
   ylim(0 , ylim1)
 
 ggsave(plot = plot2 , filename =  here("graphs",  "volcano_plot" ,  paste0("sACC region " , type ,  " ..jpeg")) , width = 8 , height = 8 )
 
 # pdf(here("graphs", "volcano_plot" ,  paste0("sACC region " , type ,  " .pdf")) , width = 8 , height = 8 )
 # 
 # Volcano_sACC
 # 
 # dev.off()
 
}

Volcano_graphs(1 , "Gene" , 0.05 , -3 , 3 , 15)
Volcano_graphs(2 , "Exon" , 0.05 , -3 , 3  , 15)
Volcano_graphs(3 , "Junction" , 0.05 , -5 , 5  , 20)
Volcano_graphs(4 , "Transcript" , 0.05 , -3 , 3  , 30)


## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()
