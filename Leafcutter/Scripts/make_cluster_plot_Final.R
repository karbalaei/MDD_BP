#here("graphs" , "BPseq_model" , "leafcutter" , region)

get_intron_meta=function(introns){
  intron_meta=do.call(rbind,strsplit(introns,":"))
  colnames(intron_meta)=c("chr","start","end","clu")
  intron_meta=as.data.frame(intron_meta,stringsAsFactors = F)
  intron_meta$start=as.numeric(intron_meta$start)
  intron_meta$end=as.numeric(intron_meta$end)
  intron_meta$middle=.5*(intron_meta$start+intron_meta$end)
  intron_meta
}

make_cluster_plot <- function(
    cluster_to_plot,
    main_title = NA,
    exons_table = NULL,
    meta = NULL,
    cluster_ids = NULL,
    counts = NULL,
    introns = NULL,
    snp_pos = NA,
    file_path = NULL){ # New argument to specify output file path
  
  
  library(ggplot2)
  library(dplyr)
  library(foreach)
  library(reshape2)
  library(stringr)
  library(intervals)
  library(ggrepel)
  library(gridExtra)
  library(grid)
 
  if( is.null(cluster_to_plot)){
    print("no cluster selected!")
    return(NULL)
  }
  
  # Ensure required inputs are not null
  stopifnot(!is.null(meta))
  stopifnot(!is.null(counts))
  stopifnot(!is.null(cluster_ids))
  stopifnot(!is.null(introns))
  
  meta$group=as.factor(meta$group)
  group_names=levels(meta$group)
  
  stopifnot(cluster_to_plot %in% cluster_ids)
  
  # --- FIX for "argument is not a matrix" error ---
  # Add drop = FALSE to ensure that even if one row is selected,
  # the result is still a matrix and can be transposed by t().
  y <- t(counts[ cluster_ids==cluster_to_plot, , drop = FALSE])
  # --- END OF FIX ---
  
  x <- meta$group
  length_transform <- function(g){ log(g+1) }
  introns_to_plot <- introns[ introns$clusterID == cluster_to_plot, ]
  summary_func <- colSums
  legend_title <- "Mean counts"
  alphabet <- c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z")
  
  junction_colour <- "red"
  cryptic_colour <- "pink"
  
  # convert colnames(y) into intron meta data
  intron_meta=get_intron_meta(colnames(y))
  
  intron_meta$verdict <- introns_to_plot$verdict[match(paste(intron_meta$start,intron_meta$end), paste(introns_to_plot$start, introns_to_plot$end) ) ]
  intron_meta$dPSI <- introns_to_plot$deltapsi[match(paste(intron_meta$start,intron_meta$end), paste(introns_to_plot$start, introns_to_plot$end) ) ]
  
  if (nrow(intron_meta) > length(alphabet)){
    stop("Not enough letters in `alphabet` for ranking. Please extend it.")
  }
  ranks <- alphabet[1:nrow(intron_meta)]
  
  absPSI <- intron_meta$dPSI[ order( abs(intron_meta$dPSI), decreasing=TRUE) ]
  intron_meta$rank <- ranks[match(intron_meta$dPSI, absPSI)]
  
  # make sure intron_meta has "chr" in front of chromosome name so it plays nice with the exon table
  if( all( ! grepl("chr", intron_meta$chr))){
    intron_meta$chr <- paste0("chr", as.character(intron_meta$chr))
  }
  
  if(!is.null(exons_table)) {
    if( all(! grepl("chr", exons_table$chr))){
      exons_table$chr <- paste0("chr", as.character(exons_table$chr))
    }
  }
  
  new_theme_empty <- theme_bw(base_size = 15 )
  new_theme_empty$panel.background = element_rect(fill="white", colour = "white")
  new_theme_empty$line <- element_blank()
  new_theme_empty$rect <- element_blank()
  new_theme_empty$strip.text <- element_blank()
  new_theme_empty$axis.text <- element_blank()
  
  groups=sort(unique(x))
  
  max_log=.5*ceiling(2*log10( 1+max( unlist( foreach (tis=groups) %do% { intron_meta$counts=summary_func(y[ tis==x,,drop=F]) } ) ) ))
  breaks=if (max_log <= 2.5) seq(0,max_log,by=0.5) else seq(0,ceiling(max_log),by=1)
  limits=c(0.0,max_log)
  
  intron_meta$id=as.factor(1:nrow(intron_meta)) # number each junction
  temp=intron_meta[,c("id","start","end")]
  m=reshape2::melt(temp, id.vars = "id") # melt to get list of all coordinates
  
  s=unique(m$value) # get unique start and end values
  if (!is.na(snp_pos)) s=c(s,snp_pos)
  s=sort(s)
  d=s[2:length(s)]-s[1:length(s)-1]
  trans_d <- length_transform(d) 
  coords <- c(0,cumsum(trans_d))
  names(coords)=s
  
  snp_coord=coords[as.character(snp_pos)]
  
  total_length=sum(trans_d)
  my_xlim=c(-.05*total_length,1.05*total_length)
  first_plot=T
  
 
  min_height=0
  max_height=0
  curv <- 0.1
  min_exon_length <- 0.5
  maxratio=0
  minratio=1.0
  yFactor = 0.65
  yConstant = -0.25
  labelTextSize=3.5
  curveMax = 10
  curveExponent = 2
  yOffset = 0
  centreLineWidth = 3
  
  mainPalette <- c("annotated" = junction_colour, "cryptic" = cryptic_colour)
  summary_func=function(a) apply( sweep(a,1,rowSums(a, na.rm=T),"/"),2, function(g) mean(g, na.rm=T) )
  
  plots <- foreach (tis=groups) %do% {
    intron_meta$counts=summary_func(y[ tis==x,,drop=F])
    maxratio=max(c(max(intron_meta$counts/sum(intron_meta$counts)), maxratio))
    minratio=min(c(min(intron_meta$counts/sum(intron_meta$counts)), minratio))
  }
  last_group=groups[length(groups)]
  
  plots <- list()
  for( fancyVar in 1:length(groups) ){
    
    intron_meta$counts=summary_func(y[ groups[fancyVar]==x,,drop=F])
    intron_meta$prop=intron_meta$counts
    
    group_sample_size=sum(groups[fancyVar]==x)
    
    allEdges=do.call(rbind,foreach (i=1:nrow(intron_meta)) %do% {
      if (i%%2==1) return(NULL)
      start=coords[ as.character(intron_meta$start[i]) ]
      end=coords[ as.character(intron_meta$end[i]) ]
      if(is.na(start) || is.na(end)) return(NULL)
      l=end-start
      edge = data.frame(start, end)
      edge$startv <- intron_meta$start[i]
      edge$endv <- intron_meta$end[i]
      edge$log10counts=intron_meta$counts[i]+1
      edge$label=paste0(format(intron_meta$prop[i],digits=2, scientific=FALSE),"^", intron_meta$rank[i])
      edge$clu<- intron_meta$clu[i]
      edge$Group <- i
      edge$xtext <-start+l/2
      edge$ytext <- -( ( l^(yFactor) / 2 ) + yConstant)
      edge$verdict <- ifelse( intron_meta$verdict[i] == "annotated", yes = "annotated", no ="cryptic")
      edge
    })
    
    allEdgesP=do.call(rbind,foreach (i=1:nrow(intron_meta)) %do% {
      if (i%%2==0) return(NULL)
      start=coords[ as.character(intron_meta$start[i]) ]
      end=coords[ as.character(intron_meta$end[i]) ]
      if(is.na(start) || is.na(end)) return(NULL)
      l=end-start
      edge = data.frame(start, end)
      edge$startv <- intron_meta$start[i]
      edge$endv <- intron_meta$end[i]
      edge$log10counts=intron_meta$counts[i]+1
      edge$label=paste0(format(intron_meta$prop[i],digits=2, scientific=FALSE),"^",intron_meta$rank[i])
      edge$clu<- intron_meta$clu[i]
      edge$Group <- i
      edge$xtext <-start+l/2
      edge$ytext <- l^(yFactor)/2+yConstant
      edge$SIZE <- intron_meta$prop[i]+1
      edge$verdict <- ifelse( intron_meta$verdict[i] == "annotated", yes = "annotated", no ="cryptic")
      edge
    })
    
    # Handle cases where all junctions are odd or even
    if(is.null(allEdgesP) && is.null(allEdges)) next
    YLIMP <- if (!is.null(allEdgesP) && nrow(allEdgesP) > 0) max(allEdgesP$ytext) * 1.25 else 1
    YLIMN <- if (!is.null(allEdges) && nrow(allEdges) > 0) min(allEdges$ytext) * 1.25 else -1
    
    g <- ggplot()
    if(!is.null(allEdgesP) && nrow(allEdgesP) > 0){
      g <- g + geom_curve(data=allEdgesP, aes(x = start, xend = xtext, y = yOffset, yend = ytext, group = Group, colour = verdict, size = curveMax * (log10counts - 1)^curveExponent ),
                          angle=90, curvature=-curv,lineend="round") +
        geom_curve(data=allEdgesP, aes(x = xtext, xend = end, y = ytext, yend = yOffset, group = Group, colour = verdict, size = curveMax * (log10counts - 1)^curveExponent ),
                   angle=90, curvature=-curv,lineend="round")
    }
    if(!is.null(allEdges) && nrow(allEdges) > 0){
      g <- g + geom_curve(data=allEdges, aes(x = start, xend = xtext, y = -yOffset, yend = ytext, group = Group, colour = verdict, size = curveMax * (log10counts - 1)^curveExponent ),
                          angle=90,curvature=curv,lineend="round") +
        geom_curve(data=allEdges, aes(x = xtext, xend = end, y = ytext, yend = -yOffset, group = Group, colour = verdict, size = curveMax * (log10counts - 1)^curveExponent ),
                   angle=90,curvature=curv,lineend="round")
    }
    
    g <- g + new_theme_empty +
      ylab(paste0(groups[fancyVar]," (n=",group_sample_size,")")) +
      xlab("") +
      xlim(my_xlim) +
      ggtitle(paste0(groups[fancyVar]," (n=",group_sample_size,")" ) ) +
      geom_hline(yintercept=0, size = centreLineWidth, colour = "white") +
      geom_hline(yintercept=0,alpha=.9, size=1) +
      ylim(YLIMN,YLIMP) +
      scale_size_continuous(limits=c(0,10),guide='none')
    
    if(!is.null(allEdgesP) && nrow(allEdgesP) > 0){
      g <- g + geom_label(data=allEdgesP,aes(x=xtext,y=0.95*ytext,label=label), size = labelTextSize, label.size = NA, parse=TRUE, fill = "white",colour = "black", label.r = unit(0.3,"lines"), label.padding = unit(0.3,"lines") )
    }
    if(!is.null(allEdges) && nrow(allEdges) > 0){
      g <- g + geom_label(data=allEdges,aes(x=xtext,y=0.95*ytext,label=label), size= labelTextSize, label.size = NA, parse=TRUE, fill = "white", colour = "black", label.r = unit(0.3,"lines"), label.padding = unit(0.3,"lines") )
    }
    
    if (!is.na(snp_coord)) {
      df=data.frame(x=snp_coord,xend=snp_coord,y=0,yend=max(c(allEdgesP$ytext,0))*1.1)
      g <- g + geom_segment(data=df,aes(x=x,y=y,xend=xend,yend=yend))
    }
    plots[[fancyVar]] <- g
  }
  
  if (length(plots) == 0) {
    message("No plots could be generated. Check input data.")
    return(NULL)
  }
  
   if (!is.null(exons_table)) {
    
    exons_chr <- exons_table[ exons_table$chr==intron_meta$chr[1], ] 
    if (nrow(exons_chr) > 0) {
      exons_here <- exons_chr[ ( min(s) <= exons_chr$start & exons_chr$start <= max(s) ) | ( min(s) <= exons_chr$end & exons_chr$end <= max(s) ), ]
      
      if( nrow(exons_here) > 0 ){
        exons_here <- unique(exons_here[ ( exons_here$end %in% intron_meta$start | exons_here$start %in% intron_meta$end ) & ( exons_here$end - exons_here$start <= 500 | exons_here$end == min(intron_meta$start) | exons_here$start == max(intron_meta$end) ), ])
      }
      if ( nrow( exons_here) > 0) {
        exons_here$gene_name=factor(exons_here$gene_name)
        gene_name_df <- data.frame( label=rev(levels(exons_here$gene_name)) )
        gene_name_df$N <- table(exons_here$gene_name)[ as.character(gene_name_df$label) ]
        gene_name_df <- gene_name_df[ order(gene_name_df$N, decreasing = TRUE), ]
        
        cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        gene_colors <- cbbPalette[ 1:length(gene_name_df$label) ]
        names(gene_colors) <- gene_name_df$label
        mainPalette <- c(gene_colors, mainPalette)
        
        invert_mapping=function(pos){
          if (pos %in% s) coords[as.character(pos)] else
            if (pos < min(s)) my_xlim[1] else
              if (pos > max(s)) my_xlim[2] else {
                w=which( pos < s[2:length(s)] & pos > s[1:(length(s)-1)] )
                if(length(w) == 0) return(NA)
                stopifnot(length(w)==1)
                coords[w] + (coords[w+1]-coords[w])*(pos - s[w])/(s[w+1]-s[w])
              }
        }
        
        exon_df <- data.frame( x=sapply(exons_here$start,invert_mapping),
                               xend=sapply(exons_here$end,invert_mapping),
                               y=0, yend=0, label = exons_here$gene_name)
        
        # check for NAs from invert_mapping and remove
        exon_df <- exon_df[ !is.na(exon_df$x) & !is.na(exon_df$xend), ]
        
        if(nrow(exon_df) > 0){
          exon_df$xend[which((exon_df$xend - exon_df$x) < min_exon_length)] <- exon_df$x[which((exon_df$xend - exon_df$x) < min_exon_length)] + min_exon_length
          exon_df <- exon_df[ !duplicated( paste(exon_df$x, exon_df$xend) ),]
          
          # Add exons to plots
          for (i in 1:length(plots) ){
            plots[[i]] <- plots[[i]] +
              geom_segment( data=exon_df, aes(x=x,y=y,xend=xend,yend=yend, colour = label), alpha=1, size=6)
          }
        }
      }
    }
  }
  
  if( all( !is.na(main_title) ) ){
    plots[[1]] <- plots[[1]] + ggtitle(main_title[1])
  }
  
 
  # Add color palette and manage legends
  # The legend should only appear on the last plot in the stack.
  for(i in 1:length(plots)) {
    plots[[i]] <- plots[[i]] + scale_colour_manual("", values = mainPalette)
    if (i < length(plots)) {
      # Hide legend on all but the last plot
      plots[[i]] <- plots[[i]] + guides(colour = "none")
    } else {
      # Style the legend on the last plot
      plots[[i]] <- plots[[i]] + theme(legend.position="bottom", legend.justification = 'right')
    }
  }
  
  # Arrange the plots into a single graphical object (grob)
  final_plot_object <- if (length(plots) > 1) {
    gridExtra::arrangeGrob(grobs = plots, ncol = 1)
  } else {
    plots[[1]]
  }
  
  # If a file_path is provided, save the plot to the file.
  # Otherwise, display the plot in the current graphics device (e.g., RStudio Plots pane).
  if (!is.null(file_path)) {
    # Dynamically set height based on number of plots
    plot_height <- 2 + 4 * length(plots)
    ggsave(filename = file_path, plot = final_plot_object, width = 8, height = plot_height, units = "in", dpi = 300)
    message(paste("Plot successfully saved to:", file_path))
  } else {
    # This will draw the plot to the screen if running interactively
    if (length(plots) > 1){
      grid::grid.draw(final_plot_object)
    } else {
      print(final_plot_object)
    }
  }
  
  # Invisibly return the plot object so it can be further manipulated if needed
  invisible(final_plot_object)
}
