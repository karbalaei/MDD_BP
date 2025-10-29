### Library #####
library(here)
library("ggVennDiagram")
library(purrr)
library("ggvenn")
library(UpSetR)
library(VennDetail)
library(magrittr)
library(plyr)
library(tidyverse)
library(ggplot2)
library(reshape2)


load(file = here::here("results" ,"BP_MDD_All_results.RDS"))


#### Functions for extracting gene sets ####

# To select differentially expressed assigned genes 
get_signif_genes <- function(outFeature, colname = "common_feature_id", cutoff = 0.05, return_unique = FALSE){
  signif <- outFeature[[colname]][outFeature$q_PrimaryDxBipolar < cutoff]
  if(return_unique) signif <- unique(signif)
  signif <- signif[!is.na(signif)]
  return(signif)
}

# To select differentially expressed features 

get_signif_feature <- function(outFeature, cutoff = 0.05, return_unique = FALSE){
  signif <- rownames(outFeature)[outFeature$q_PrimaryDxBipolar < cutoff]
  if(return_unique) signif <- unique(signif)
  signif <- signif[!is.na(signif)]
  return(signif)
}


my_flatten <- function (x, use.names = TRUE, classes = "ANY") {
  #' Source taken from rlist::list.flatten
  len <- sum(rapply(x, function(x) 1L, classes = classes))
  y <- vector("list", len)
  i <- 0L
  items <- rapply(x, function(x) {
    i <<- i + 1L
    y[[i]] <<- x
    TRUE
  }, classes = classes)
  if (use.names && !is.null(nm <- names(items))) 
    names(y) <- nm
  y
}

#### Get signif genes####

signif_genes <- map_depth(Outlist, 2, ~get_signif_genes(.x, colname = "common_gene_id", return_unique = TRUE))
map_depth(signif_genes, 2, length)


signif_features <- map_depth(Outlist, 2, ~get_signif_feature(.x,  return_unique = TRUE))
map_depth(signif_features, 2, length)


signif_genes_flat <- my_flatten(map_depth(signif_genes, 2, unlist))
map_int(signif_genes_flat, length)


signif_feature_flat <- my_flatten(map_depth(signif_features, 2, unlist))

map_int(signif_feature_flat, length)


Amygdala_Genes_All <-   unique(list_c(signif_genes_flat[c(1 , 3, 5 , 7)]))
sACC_Genes_All <-   unique(list_c(signif_genes_flat[c(2 , 4, 6 , 8)]))

List_All <-  list( Amygdala = Amygdala_Genes_All ,  sACC = sACC_Genes_All)

### Venn diagram Genes ####


# Step 1: Define all plot parameters in a data frame or create a "manifest"

output_dir <- here::here("graphs", "Venn_diagram")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

plot_params <- tibble::tribble(
  ~plot_type, ~file_suffix,             ~data_indices,         ~extra_args,
  "upset",    "Upset graph All",        1:8,                   list(nsets = 8),
  "upset",    "Upset graph Amygdala",   c(1, 3, 5, 7),         list(nsets = 4),
  "upset",    "Upset graph sACC",       c(2, 4, 6, 8),         list(nsets = 4),
  "ggvenn",   "Venn_diagram_Amygdala Assigned genes", c(1, 3, 5, 7), list(set_name_size = 3.5, text_size = 3),
  "ggvenn",   "Venn_diagram_sACC Assigned genes",     c(2, 4, 6, 8), list(set_name_size = 3.5, text_size = 3),
  "ggvenn",   "Venn_diagram_gene",      c(1, 2),               list(set_name_size = 4, text_size = 4),
  "ggvenn",   "Venn_diagram_exon",      c(3, 4),               list(set_name_size = 4, text_size = 4),
  "ggvenn",   "Venn_diagram_junction",  c(5, 6),               list(set_name_size = 4, text_size = 4),
  "ggvenn",   "Venn_diagram_transcript",c(7, 8),               list(set_name_size = 4, text_size = 4)
)

#Step 2: Create a robust, server-safe (!) plotting function

create_plot_non_interactive <- function(plot_type, file_suffix, data_indices, extra_args) {
  plot_data <- signif_genes_flat[data_indices]
  
  if (plot_type == "ggvenn") {
    ggvenn_args <- c(list(plot_data), extra_args, list(stroke_size = 0.5, fill_color = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")))
    p <- rlang::exec(ggvenn, !!!ggvenn_args)
    ggsave(filename = file.path(output_dir, paste0(file_suffix, ".pdf")), plot = p, width = 8, height = 8)
    ggsave(filename = file.path(output_dir, paste0(file_suffix, ".jpg")), plot = p, width = 10, height = 8.33, dpi = 600)
    
  } else if (plot_type == "upset") {
    upset_args <- c(list(fromList(plot_data)), 
		extra_args, list(matrix.color = "#1b9e77", 
				main.bar.color = "#7570b3", 
				sets.bar.color = "#e7298a" , 
				mainbar.y.max = 3000
				))

    # PDF Version 
    pdf(file.path(output_dir, paste0(file_suffix, ".pdf")), width = 8, height = 8)
    print(do.call('upset', upset_args))
    dev.off()
    
    # JPEG Version 
    jpeg(file.path(output_dir, paste0(file_suffix, ".jpg")), width = 6000, height = 5000, res = 600)
    print(do.call('upset', upset_args))
    dev.off()
  }
  
  message("Generated plot: ", file_suffix)
}


#Step 3: Iterate and create all plots: iterate through each row of the parameter data frame and execute the plotting function
pwalk(plot_params, create_plot_non_interactive)

# Create the Venn diagram for all assigned features :

message("Generating final plot for List_All...")
ggvenn_all_plot <- ggvenn(
  List_All, stroke_size = 0.5, fill_color = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"),
  set_name_size = 3.5, text_size = 3
)
ggsave(
  filename = file.path(output_dir, "Venn_diagram_Amygdala sACC compare All Assigned genes.pdf"), 
  plot = ggvenn_all_plot, width = 8, height = 8
)
ggsave(
  filename = file.path(output_dir, "Venn_diagram_Amygdala sACC compare All Assigned genes.jpg"), 
  plot = ggvenn_all_plot, width = 10, height = 8.33, dpi = 600
)
message("All plots generated successfully.")

### Venn diagram features ####

# Step 1: Define all plot parameters in a data frame 
output_dir <- here::here("graphs", "Venn_diagram")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

plot_params_features <- tibble::tribble(
  ~file_suffix,             ~data_indices,         ~extra_args,
  "Venn_diagram_gene",      c(1, 2),               list(set_name_size = 4, text_size = 4),
  "Venn_diagram_exon",      c(3, 4),               list(set_name_size = 4, text_size = 4),
  "Venn_diagram_junction",  c(5, 6),               list(set_name_size = 4, text_size = 4),
  "Venn_diagram_transcript",c(7, 8),               list(set_name_size = 4, text_size = 4)
)

#Step 2: Create a robust, server-safe (!) plotting function
create_plot_non_interactive_feature <- function(file_suffix, data_indices, extra_args) {
  plot_data <- signif_feature_flat[data_indices]
  
    ggvenn_args <- c(list(plot_data), extra_args, list(stroke_size = 0.5, fill_color = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")))
    p <- rlang::exec(ggvenn, !!!ggvenn_args)
    ggsave(filename = file.path(output_dir, paste0(file_suffix, "_feature.pdf")), plot = p, width = 8, height = 8)
    ggsave(filename = file.path(output_dir, paste0(file_suffix, "_feature.jpg")), plot = p, width = 10, height = 8.33, dpi = 600)
    
  message("Generated plot: ", file_suffix)
}


# Step 3: Iterate and create all plots 
pwalk(plot_params_features, create_plot_non_interactive_feature)

### Bar plots #####

custom_order <- c("Gene" ,  "Exon", "Junction",  "Transcript" )

prepare_barplot_data <- function(data_list) {
  # Define the pairs to compare
  feature_pairs <- list(
    Gene = data_list[c(1, 2)],
    Exon = data_list[c(3, 4)],
    Junction = data_list[c(5, 6)],
    Transcript = data_list[c(7, 8)]
  )
  
  # Use map_dfr to iterate, calculate details, and row-bind into a single data frame
  map_dfr(feature_pairs, ~as.data.frame(detail(venndetail(.x))), .id = "Type") %>%
    # Give clear names to the two resulting columns
    `colnames<-`(c("Type", "Fre")) %>%
    # Group by the feature type
    group_by(Type) %>%
    # Within each group, add the Region name based on the known order from venndetail
    mutate(Region = rep(c("Common", "sACC", "Amygdala"), 4)) %>%
    # It's good practice to ungroup after the operation is done
    ungroup() %>%
    # Ensure the columns are in a logical order
    select(Type, Region, Fre)
}



create_and_save_barplot <- function(plot_data, file_suffix , Title) {
  # Create the ggplot object
  p <- plot_data %>%
    mutate(Type = fct_relevel(Type, custom_order)) %>%
    ggplot(aes(x = Type, y = Fre, fill = Region)) +
    geom_bar(stat = 'identity', position = position_dodge()) +
    geom_text(
      stat = 'identity',
      position = position_dodge(width = 0.8),
      aes(label = Fre),
      vjust = -1,
      size = 3
    ) +
    scale_fill_manual(values = c("#66c2a5", "grey", "#fc8d62")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
    labs(x = "", y = "No. significant features") + ggtitle(label = Title ) +
    ylim(0 , 3000)
  
  # Save PDF
  pdf(here::here("graphs", paste0(file_suffix, ".pdf")), width = 9, height = 6)
  print(p)
  dev.off()
  
  # Save JPEG
  jpeg(here::here("graphs", paste0(file_suffix, ".jpg")), width = 6000, height = 4000, res = 600)
  print(p)
  dev.off()
}


# Refactored Bar Plots Section 

# 1. Prepare data for both gene and feature plots
gene_barplot_data <- prepare_barplot_data(signif_genes_flat)
feature_barplot_data <- prepare_barplot_data(signif_feature_flat)

# 2. Create and save both plots
create_and_save_barplot(gene_barplot_data, "Barplot_DEf_gene" , "Assigned genes" )
create_and_save_barplot(feature_barplot_data, "Barplot_DEf_feature" , "Features")


#### Gene type ######


Amygdala_gene_type <-  Outlist$gene$amyg  %>% 
  dplyr::filter(q_PrimaryDxBipolar < 0.05) %>% 
  dplyr::select(c(common_gene_id ,gene_type )) %>% 
  mutate("Region" = "Amygdala")

sACC_gene_type <-  Outlist$gene$sacc %>% 
  dplyr::filter(q_PrimaryDxBipolar < 0.05) %>% 
  dplyr::select(c(common_gene_id ,gene_type )) %>% 
  mutate("Region" = "sACC")

Gene_type <-  rbind(Amygdala_gene_type , sACC_gene_type) 


unique(Gene_type$gene_type)

Common_genes <-  intersect(sACC_gene_type$common_gene_id , Amygdala_gene_type$common_gene_id)

Gene_type$Region <-replace (Gene_type$Region, Gene_type$common_gene_id %in% Common_genes, "Common") 

Gene_type <-   Gene_type %>% distinct()

Gene_type$Region <- factor(Gene_type$Region , levels = c("Amygdala"  , "Common" ,  "sACC" ))
Gene_type$Region %>% unique()

Other_RNA <-  c("snRNA","misc_RNA"  ,"sense_overlapping"  , "snoRNA" , "TEC" , "rRNA" , "Mt_rRNA" , "Mt_tRNA" , "miRNA")

pseudogene_RNA <-  c("transcribed_processed_pseudogene" , "processed_pseudogene" ,"transcribed_unitary_pseudogene")

lincRNA <- c("sense_intronic" , "lincRNA")

replace_RNA_vector <-  rep(c("Other_RNA"  , "pseudogene_RNA" , "lincRNA" ), c(length(Other_RNA),length(pseudogene_RNA) , length(lincRNA ) ))


Gene_type$gene_type <-  mapvalues(Gene_type$gene_type ,  c(Other_RNA , pseudogene_RNA , lincRNA ) , replace_RNA_vector )


unique(Gene_type$gene_type)

Gene_type$gene_type <-  factor(Gene_type$gene_type , levels = c( "protein_coding" , "antisense" , "pseudogene_RNA" , "lincRNA" , "Other_RNA" ))


### Final graph #####

gene_type_graph <- function() {

  dat = dcast (Gene_type, gene_type ~ Region, fun.aggregate = length)
  dat.melt = melt(dat, id.vars = "gene_type", measure.vars = c("Amygdala"  , "Common" ,  "sACC"))
  
  g <-  ggplot(dat.melt, aes(x = gene_type,y = value, fill = variable)) + 
    geom_bar(stat = "identity", position = position_dodge(width = .8), width = 0.7) +
    #ylim(0, 14) +
    geom_text(aes(label = value), position = position_dodge(width = .8), vjust = -0.5) +
    scale_fill_manual(values=c( "#66c2a5" , "grey" , "#fc8d62")) +
    theme_classic()  + theme(axis.text.x = element_text(size = 15 , angle = 45, vjust = 0.5, hjust=0.5),
                             axis.title.y = element_text(size = 15 , vjust = 3, hjust=0.5)   ,       axis.text.y = element_text(size = 15 )) +
    xlab("") + ylab("\nNo. significant features\n") +scale_x_discrete(breaks=c("protein_coding","pseudogene_RNA","lincRNA" , "Other_RNA" , "antisense"),
                                                                      labels=c("Coding", "Pseudogene", "lncRNA" , "Other RNA" , "Antisense")) +
    ylim( 0, 150)

  
  
  pdf(here::here("graphs" ,   "Barplot_gene_type.pdf") , width = 9 , height = 6)
  print(g)
  dev.off()
  
  # Save JPEG
  jpeg(here::here("graphs" ,   "Barplot_gene_type.jpg") , width = 6000 , height = 4000 , res = 600)
  print(g)
  dev.off()
  
}

gene_type_graph()

## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()

