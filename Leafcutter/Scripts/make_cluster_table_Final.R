tableTheme <- ttheme_minimal(
  core=list(bg_params = list(fill = c("whitesmoke","white"), col=NA)
  ),
  colhead=list(fg_params=list(col="black", fontface="bold"),
               bg_params = list(fill="white")),
  rowhead=list(fg_params=list(col="black"),
               bg_params = list(fill=c("white", "whitesmoke"))))

filter_intron_table <- function(introns, clu, toSave=FALSE , region){
  d <- dplyr::filter(introns, clusterID == clu) %>%
    dplyr::select( -clusterID, -gene, -ensemblID, -transcripts) %>%
    arrange( desc(abs(deltapsi)))
  if( !toSave ){
    d <- rename(d, "ΔPSI" = deltapsi )
  }else{
    d <- rename(d, "dPSI" = deltapsi ) # fudge as grid arrange doesn't like greek letters
  }
  row.names(d) <- letters[1:nrow(d)] # letters is just a:z
  
  mytable <- tableGrob(d, theme = tableTheme )
  mycols <- ncol(mytable)
  
  mytable$widths <- unit( c( 1/(3*mycols), rep(1/mycols, mycols-1) ), "npc")
  
  mytable <- gtable_add_grob(mytable,
                             grobs = segmentsGrob( # line across the bottom
                               x0 = unit(0,"npc"),
                               y0 = unit(0,"npc"),
                               x1 = unit(1,"npc"),
                               y1 = unit(0,"npc"),
                               gp = gpar(lwd = 2.0)),
                             t = 2, b = nrow(mytable), l = 1, r = mycols)
  
  mytable <- gtable_add_grob(mytable,
                             grobs = segmentsGrob( # line across the bottom
                               x0 = unit(0,"npc"),
                               y0 = unit(0,"npc"),
                               x1 = unit(1,"npc"),
                               y1 = unit(0,"npc"),
                               gp = gpar(lwd = 2.0)),
                             t = 1, b = 1, l = 1, r = mycols)
  #return(mytable)
  ggsave( plot = mytable, filename = here("graphs" , "leafcutter" , region , paste0(clu , "_table.pdf")) , width = 9, height = 6 )
}

