#' Make gene level plot (Robust Version)
#'
#' @import ggplot2
#' @export

make_gene_plot <- function(gene_name,
                           clusterID=NULL,
                           cluster_list=NULL,
                           introns=NULL,
                           introns_to_plot=NULL,
                           counts, # no longer used
                           exons_table=NULL,
                           len=500,
                           min_exon_length=0.5,
                           main_title=NA,
                           snp_pos=NA,
                           snp=NA,
                           summary_func=colSums,
                           legend_title="Mean counts", region , 
                           debug=F){
  
  
  library(ggplot2)
  library(dplyr)
  library(foreach)
  library(reshape2)
  library(stringr)
  library(intervals)
  library(ggrepel)
  
  message(paste0("gene name:" , gene_name))

  exonMax <- 1000
  
  stopifnot( !is.null(exons_table) )
  stopifnot( !is.null(introns_to_plot))
  
  exons <- exons_table[ exons_table$gene_name == gene_name, ]
  if(nrow(exons) == 0) stop(paste("No exons found for gene_name:", gene_name))
  
  gene_start <- min(exons$start)
  gene_end <- max(exons$end)
  stopifnot( length( unique(exons$chr) ) == 1 )
  
  introns_to_plot$chr = gsub( "^chr", "", introns_to_plot$chr )
  myChr <- gsub( "^chr", "", unique(exons$chr) )
  
  myclusters <- introns_to_plot[ introns_to_plot$start >= gene_start & introns_to_plot$end <= gene_end & introns_to_plot$chr == myChr , ]
  if(nrow(myclusters) == 0){
    stop(paste("No junctions found within the genomic range of gene '", gene_name, "'."))
  }
  cluster_ids <- myclusters$clu
  
  exons <- exons[ !duplicated(paste( exons$start, exons$end)) ,]
  exons <- exons[ exons$end - exons$start <= exonMax,]
  exons <- exons[ order(exons$end),]
  
  all_junctions <- myclusters[, c("start","end", "clu")]
  length_transform <- function(g) log(g+1)
  
  m <- reshape2::melt(all_junctions, id.vars = "clu")
  s <- unique(m$value)
  
  if (!is.na(snp_pos)){
    SNP_pos <- as.numeric(str_split_fixed(snp_pos, ":", 2)[,2])
    s <- c(s, SNP_pos)
  }
  s <- sort(s)
  
  d <- s[2:length(s)] - s[1:length(s)-1]
  trans_d <- length_transform(d)
  coords <- c(0, cumsum(trans_d))
  names(coords) <- s
  total_length <- sum(trans_d)
  
  my_xlim <- c(-.05 * total_length, 1.05 * total_length)
  if(!is.na(snp_pos)){
    snp_coord <- coords[as.character(SNP_pos)]
  }
  
  invert_mapping <- function(pos, s_map, coords_map, xlim_map){
    if (pos %in% s_map) {
      coords_map[as.character(pos)]
    } else if (pos < min(s_map)) {
      xlim_map[1]
    } else if (pos > max(s_map)) {
      xlim_map[2]
    } else {
      w <- which(pos > s_map[1:(length(s_map)-1)] & pos < s_map[2:length(s_map)])
      if (length(w) == 0) return(NA_real_) 
      stopifnot(length(w) == 1)
      coords_map[w] + (coords_map[w+1] - coords_map[w]) * (pos - s_map[w]) / (s_map[w+1] - s_map[w])
    }
  }
  
  # --- CREATE EXON_DF WITH VAPPLY FOR TYPE SAFETY ---
  x_coords <- vapply(exons$start, FUN = invert_mapping, FUN.VALUE = numeric(1), s_map = s, coords_map = coords, xlim_map = my_xlim)
  xend_coords <- vapply(exons$end, FUN = invert_mapping, FUN.VALUE = numeric(1), s_map = s, coords_map = coords, xlim_map = my_xlim)
  
  exon_df <- data.frame(
    x = x_coords,
    xend = xend_coords,
    y = 0,
    yend = 0,
    label = exons$gene_name
  )
  
  if(any(is.na(exon_df$x)) || any(is.na(exon_df$xend))){
    # warning("Some exons could not be mapped to coordinates and were removed.") # Commented out to suppress warning
    exon_df <- exon_df[!is.na(exon_df$x) & !is.na(exon_df$xend), ]
  }
  
  small_indices <- which((exon_df$xend - exon_df$x) < min_exon_length)
  if (length(small_indices) > 0) {
    exon_df$xend[small_indices] <- exon_df$x[small_indices] + min_exon_length
  }
  my_xlim <- c(min(exon_df$x, na.rm=T), max(exon_df$xend, na.rm=T))

  YLIMN <- 10
  YLIMP <- -10
  curv <- 0.4
  
  allEdgesList <- foreach (i = 1:nrow(all_junctions)) %do% {
    start_coord <- coords[as.character(all_junctions$start[i])]
    end_coord <- coords[as.character(all_junctions$end[i])]
    if(is.na(start_coord) || is.na(end_coord)) return(NULL)
    l <- end_coord - start_coord
    direction <- if (i %% 2 == 1) 1 else -1 
    data.frame(start = start_coord, end = end_coord, startv = all_junctions$start[i],
               endv = all_junctions$end[i], clu = all_junctions$clu[i], Group = i,
               curve_direction = direction)
  }
  junctions <- do.call(rbind, allEdgesList)
  
  if(!is.null(introns)){
    if (is.data.frame(introns) && all(c("deltapsi", "clusterID", "start", "end") %in% names(introns))) {
      junctions$deltaPSI <- introns$deltapsi[match(paste(junctions$clu, junctions$startv, junctions$endv), paste(introns$clusterID, introns$start, introns$end))]
    } else {
      warning("'introns' argument is not a valid data frame. Skipping dPSI coloring.")
      introns <- NULL
    }
  }
  
  label_df <- junctions %>% group_by(clu) %>%
    summarise(start = min(start), end = max(end), middle = start + ((end - start) / 2), .groups = 'drop')
  
  if(!is.null(cluster_list)){
    label_df <- left_join(label_df, cluster_list %>% select(clusterID, FDR), by = c("clu" = "clusterID"))
    label_df$FDR[is.na(label_df$FDR)] <- "."
  } else {
    label_df$FDR <- "."
  }
  label_df$label <- gsub("\n.$", "", gsub("_", "\n", gsub("_[+-]", "", paste0(label_df$clu, "\n", label_df$FDR))))
  label_df <- label_df[order(label_df$middle),]
  label_df$labelY <- YLIMP - 0.2*YLIMP 
  if (nrow(label_df) > 1) { label_df$labelY[seq(1, nrow(label_df), 2)] <- YLIMN - 0.2*YLIMN }
  
   p <- ggplot()
  
  junctions_to_plot <- left_join(junctions, label_df %>% select(clu, FDR), by = "clu")
  sig_junctions <- junctions_to_plot %>% filter(FDR != ".")
  nonsig_junctions <- junctions_to_plot %>% filter(FDR == ".")
 
  # Non-significant junctions
  nonsig_up <- nonsig_junctions %>% filter(curve_direction == 1)
  nonsig_down <- nonsig_junctions %>% filter(curve_direction == -1)
  
  if (nrow(nonsig_up) > 0) {
    p <- p + geom_curve(data = nonsig_up, aes(x = start, y = 0, xend = end, yend = 0, group = Group),
                        curvature = curv, color = "grey80", size = 0.5)
  }
  if (nrow(nonsig_down) > 0) {
    p <- p + geom_curve(data = nonsig_down, aes(x = start, y = 0, xend = end, yend = 0, group = Group),
                        curvature = -curv, color = "grey80", size = 0.5)
  }
  
  # Significant junctions
  sig_up <- sig_junctions %>% filter(curve_direction == 1)
  sig_down <- sig_junctions %>% filter(curve_direction == -1)
  
  if (!is.null(introns) && "deltaPSI" %in% names(sig_junctions)) {
    # WITH dPSI coloring
    if (nrow(sig_up) > 0) {
      p <- p + geom_curve(data = sig_up, aes(x = start, y = 0, xend = end, yend = 0, group = Group, color = deltaPSI > 0),
                          curvature = curv, size = 0.8)
    }
    if (nrow(sig_down) > 0) {
      p <- p + geom_curve(data = sig_down, aes(x = start, y = 0, xend = end, yend = 0, group = Group, color = deltaPSI > 0),
                          curvature = -curv, size = 0.8)
    }
    # Important: The scale applies to the whole plot, so add it only once.
    # Note: R sorts logicals as FALSE, TRUE. So the first label/color is for FALSE (deltaPSI <= 0, "Down").
    p <- p + scale_color_manual("dPSI", values = c("FALSE" = "darkturquoise", "TRUE" = "firebrick2"), labels = c("Down", "Up"))
  } else {
    # WITHOUT dPSI coloring
    if (nrow(sig_up) > 0) {
      p <- p + geom_curve(data = sig_up, aes(x = start, y = 0, xend = end, yend = 0, group = Group),
                          curvature = curv, color = "#d66464", size = 0.8)
    }
    if (nrow(sig_down) > 0) {
      p <- p + geom_curve(data = sig_down, aes(x = start, y = 0, xend = end, yend = 0, group = Group),
                          curvature = -curv, color = "#d66464", size = 0.8)
    }
  }

  # Add central gene structure
  p <- p +
    geom_hline(yintercept = 0, size = 3, colour = "white") +
    geom_segment(data = exon_df, aes(x = x, y = y, xend = xend, yend = yend), size = 6, colour = 'black') +
    geom_hline(yintercept = 0, alpha = .9, size = 1)
  
  # Add cluster labels
  if (!is.null(cluster_list)) {
    p <- p +
      geom_segment(data = label_df, aes(x = start, xend = middle, y = 0, yend = labelY), colour = "grey", linetype = "dashed") +
      geom_segment(data = label_df, aes(x = end, xend = middle, y = 0, yend = labelY), colour = "grey", linetype = "dashed") +
      geom_text_repel(data = label_df, aes(x = middle, y = labelY, label = label), point.padding = NA, 
                      direction = "y", segment.alpha = 0, size = 4)
  }
  
  # Add final theme
  p <- p +
    ylim(YLIMP, YLIMN) +
    labs(title = paste(gene_name, collapse = "+"), x = NULL, y = NULL) +
    theme_bw(base_size = 15) +
    theme(plot.title = element_text(face = "bold.italic", size = 20, hjust = 0.5), panel.grid = element_blank(),
          panel.border = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
  
  ggsave( plot = p, filename = here("graphs" , "leafcutter" , region , paste0(clusterID , "_" , gene_name, ".pdf")) , width = 9, height = 6 )
}
