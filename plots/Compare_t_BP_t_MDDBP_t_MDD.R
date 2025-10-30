#################################################################
### SECTION 1: LOAD LIBRARIES
#################################################################
message("Loading libraries...")
library(here)
library(tidyverse)
library(SummarizedExperiment)
library(jaffelab)
library(gridExtra)
library(ragg) # For high-quality plot saving
library(grid) 
#################################################################
### SECTION 2: LOAD DATA
#################################################################
message("Loading data objects...")

# Load Bipolar Disorder differential expression results
load(here("results", "BP results.RDS"), verbose = TRUE) # Assumes this creates 
load(here("results", "BP_MDD_All_results_test.RDS"), verbose = TRUE) # Assumes this creates 

Outlist_test <- map(Outlist_test, ~ {
  
    return(rev(.x))
  })

data_type <- c("gene", "exon", "jxn","tx")

files <- here("results", paste0("qSVA_MDD_",data_type,"_DEresults.rda"))
names(files) <- data_type

allOut <- lapply(files, function(x) get(load(x, verbose = TRUE)))

## sep results only
allOut <- transpose(allOut)
allOut <- allOut$sep

map_depth(allOut, 3, nrow)

Outlist_MDD <- lapply(allOut, function(main_element) {
  lapply(main_element, function(sub_object) {
    # For each sub-object, return only the first element (the first data frame)
    return(sub_object[[1]])
  })
})

map_depth(Outlist_MDD, 2, nrow)

message("Defining the main analysis function...")

Compare_models <- function(region, type, i, j) {

    message("... Performing general T-stat comparison analysis")
    BP_results <- Outlist_BP[[i]][[j]] %>% rownames_to_column("Id") %>% dplyr::filter(!(is.na(t)))
    MDD_results <- Outlist_MDD[[i]][[j]] %>% rownames_to_column("Id") %>% dplyr::filter(!(is.na(t)))
    MDDBP_results <- Outlist_test[[i]][[j]] %>% rownames_to_column("Id") %>% dplyr::filter(!(is.na(t)))



    colnames(BP_results) <- paste0(colnames(BP_results), "_BP")
    colnames(MDDBP_results) <- paste0(colnames(MDDBP_results), "_MDDBP")
    colnames(MDD_results) <- paste0(colnames(MDD_results), "_MDD")

    All_results <- BP_results %>%
            inner_join(MDDBP_results, by = c("Id_BP" = "Id_MDDBP")) %>%
      dplyr::mutate("Group" = if_else((adj.P.Val_BP < 0.05 & q_PrimaryDxBipolar_MDDBP < 0.05), "Both DEf",
        if_else((adj.P.Val_BP < 0.05 & q_PrimaryDxBipolar_MDDBP >= 0.05), "BP DEf",
          if_else((adj.P.Val_BP >= 0.05 & q_PrimaryDxBipolar_MDDBP < 0.05), "MDDBP DEf", "None"))))

    All_results2 <- MDD_results %>%
            inner_join(MDDBP_results, by = c("Id_MDD" = "Id_MDDBP")) %>%
      dplyr::mutate("Group" = if_else((adj.P.Val_MDD < 0.05 & q_PrimaryDxBipolar_MDDBP < 0.05), "Both DEf",
        if_else((adj.P.Val_MDD < 0.05 & q_PrimaryDxBipolar_MDDBP >= 0.05), "MDD DEf",
          if_else((adj.P.Val_MDD >= 0.05 & q_PrimaryDxBipolar_MDDBP < 0.05), "MDDBP DEf", "None")
        )
      ))

        All_results3 <- MDD_results %>%
            inner_join(BP_results, by = c("Id_MDD" = "Id_BP")) %>%
      dplyr::mutate("Group" = if_else((adj.P.Val_MDD < 0.05 & adj.P.Val_BP < 0.05), "Both DEf",
        if_else((adj.P.Val_MDD < 0.05 & adj.P.Val_BP >= 0.05), "MDD DEf",
          if_else((adj.P.Val_MDD >= 0.05 & adj.P.Val_BP < 0.05), "BP DEf", "None")
        )
      ))



scale_fill_p1 <- function(...){
  ggplot2:::manual_scale(
    'color', 
    values = setNames(c('#e41a1c', '#377eb8', '#984ea3', '#4daf4a'), c("Both DEf", "BP DEf","MDDBP DEf", "None" )), 
    ...
  )
}

scale_fill_p2 <- function(...){
  ggplot2:::manual_scale(
    'color', 
    values = setNames(c('#e41a1c', '#377eb8', '#984ea3', '#4daf4a'), c("Both DEf", "MDD DEf","MDDBP DEf", "None" )), 
    ...
  )
}

scale_fill_p3 <- function(...){
  ggplot2:::manual_scale(
    'color', 
    values = setNames(c('#e41a1c', '#377eb8', '#984ea3', '#4daf4a'), c("Both DEf", "MDD DEf","BP DEf", "None" )), 
    ...
  )
}



    p1 <- All_results %>% ggplot(aes(t_BP, t_MDDBP, color = Group)) + geom_point(alpha = 0.7) +
	theme_classic()+
    scale_fill_p1() +
     labs(title = "Compare BP vs. MDD/BP" ,
       x ="BP vs. control", y = "MDD vs. BP") +
  xlim(-10 , 10) + ylim(-10 , 10)+
  theme(        axis.text.x = element_text(face="bold",  
                                   size=14),
        axis.text.y = element_text(face="bold",  
                                   size=14) ,
        #axis.title.x = element_blank(),
        #axis.title.y = element_blank(),
        plot.title = element_text(size=20, face="bold.italic")
  )


    p2 <- All_results2 %>% ggplot(aes(t_MDD, t_MDDBP, color = Group)) + geom_point(alpha = 0.7) +
	theme_classic()+
    scale_fill_p2() +
     labs(title = "Compare MDD vs. MDD/BP" ,
       x ="MDD vs. control", y = "MDD vs. BP") +
  xlim(-10 , 10) + ylim(-10 , 10)+
  theme(        axis.text.x = element_text(face="bold",  
                                   size=14),
        axis.text.y = element_text(face="bold",  
                                   size=14) ,
        #axis.title.x = element_blank(),
        #axis.title.y = element_blank(),
        plot.title = element_text(size=20, face="bold.italic")
  )


    p3 <- All_results3 %>% ggplot(aes(t_MDD, t_BP, color = Group)) + geom_point(alpha = 0.7) +
	theme_classic()+
    scale_fill_p3() +
     labs(title = "Compare MDD vs. BP" ,
       x ="MDD vs. control", y = "BP vs. control") +
  xlim(-10 , 10) + ylim(-10 , 10)+
  theme(        axis.text.x = element_text(face="bold",  
                                   size=14),
        axis.text.y = element_text(face="bold",  
                                   size=14) ,
        #axis.title.x = element_blank(),
        #axis.title.y = element_blank(),
        plot.title = element_text(size=20, face="bold.italic")
  )

    plot_lists <- list(p1, p2 , p3)

  jpeg(here("graphs" , "twas_Enrichment results" ,   paste0("Compare t stat_" , type , "_" , region , ".jpeg")) , width = 18000, height = 6000 , res = 600)

    grid.arrange(grobs = plot_lists, ncol = 3)

 dev.off()

}

#################################################################
### SECTION 5: EXECUTION
#################################################################

# Define the plan for all analyses using character strings for object names
plot_params <- tibble::tribble(
  ~region,        ~type,        ~i, ~j,
  "Amygdala",  "Gene",       1,  1,
  "sACC",  "Gene",       1,  2,
  "Amygdalal",  "Exon",       2,  1,
  "sACC",  "Exon",       2,  2,
  "Amygdalal",  "Junction",   3,  1,
  "sACC",  "Junction",   3,  2,
  "Amygdalal",  "Transcript", 4,  1,
  "sACC",  "Transcript", 4,  2
)

pwalk(.l = plot_params , .f = Compare_models , .progress = T )

#################################################################
### SECTION 6: REPRODUCIBILITY
#################################################################
message("--- All analyses complete. ---")
Sys.time()
proc.time()
options(width = 120)
library(sessioninfo)
sessioninfo::session_info()