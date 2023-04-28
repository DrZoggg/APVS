
# checkGroups
plotMsigNetwork_galaxy <- function (ig, markGroups = NULL, genesetStat = NULL, nodeSF = 1, 
                                edgeSF = 1, lytFunc = "graphopt", lytParams = list()) 
{
  stopifnot(nodeSF > 0)
  stopifnot(edgeSF > 0)
  stopifnot(is.null(genesetStat) || !is.null(names(genesetStat)))
  # stopifnot(is.null(markGroups) || checkGroups(markGroups, 
  #                                              V(ig)$name))
  if (length(markGroups) > 12) {
    warning("Only the first 12 components will be plot")
    markGroups = markGroups[seq_len(12)]
  }
  ig = igraph::induced_subgraph(ig, V(ig)[igraph::degree(ig) > 
                                            0])
  if (!all(c("Category") %in% igraph::list.vertex.attributes(ig))) {
    V(ig)$Category = rep("custom", length(V(ig)))
  }
  colmap_nodes = RColorBrewer::brewer.pal(10, "Set3")
  names(colmap_nodes) = c("archived", "c1", "c2", 
                          "c3", "c4", "c5", "c6", "c7", 
                          "h", "custom")
  lytParams = c(list(graph = igraph::as.directed(ig), layout = lytFunc), 
                lytParams)
  p1 = do.call(ggraph::ggraph, lytParams) + ggraph::geom_edge_link(edge_width = 0.2 * 
                                                                     edgeSF, alpha = 1/log10(length(igraph::E(ig))), colour = "#66666666") + 
    ggraph::geom_node_point(aes(size = Size), colour = "#FFFFFF")
  if (is.null(genesetStat)) {
    p1 = p1 + ggraph::geom_node_point(aes(fill = Category, 
                                          size = Size), alpha = 0.75, shape = 21, stroke = 0.2 * 
                                        nodeSF) + ggplot2::scale_fill_manual(values = colmap_nodes) + 
      guides(fill = guide_legend(ncol = 4, override.aes = list(size = 5)))
  }
  else {
    p1$data$genesetStat = genesetStat[p1$data$name]
    p1 = p1 + ggraph::geom_node_point(aes(fill = genesetStat, 
                                          size = Size), alpha = 0.75, shape = 21, stroke = 0.2 * 
                                        nodeSF)
    if (!all(genesetStat >= 0)) {
      lims = c(-1, 1) * stats::quantile(abs(genesetStat), 
                                        0.99)
      palname = "cork"
      dir = 1
    }
    else {
      lims = stats::quantile(abs(genesetStat), c(0.01, 
                                                 0.99))
      palname = "tokyo"
      dir = -1
    }
    p1 = p1 + scico::scale_fill_scico(palette = palname, 
                                      na.value = "#AAAAAA", limits = lims, oob = scales::squish, 
                                      direction = dir)
  }
  p1 = p1 + guides(size = guide_legend(ncol = 2)) + ggplot2::scale_size_continuous(range = c(0.1, 
                                                                                             6) * nodeSF) + ggplot2::theme_void() + ggplot2::theme(legend.position = "bottom", 
                                                                                                                                                   plot.title = element_text(hjust = 0.5, size = rel(1.5)))
  if (!is.null(markGroups)) {
    # checkGroups(markGroups, p1$data$name)
    names(markGroups) = sapply(names(markGroups), function(x) {
      paste0(x, " (n = ", length(markGroups[[x]]), 
             ")")
    })
    hulldf = plyr::ldply(markGroups, function(x) {
      df = p1$data[p1$data$name %in% x, ]
      df = df[grDevices::chull(df$x, df$y), ]
    }, .id = "NodeGroup")
    p1 = p1 + ggforce::geom_shape(aes(x, y, colour = NodeGroup), 
                                  fill = NA, radius = ggplot2::unit(0.01, "npc"), 
                                  expand = ggplot2::unit(0.02, "npc"), data = hulldf) + 
      # ggplot2::scale_colour_brewer(palette = "Paired") + 
      scale_color_manual(values = alpha(colorRampPalette(	
        paletteer_d("ggthemes::Red_Blue_Brown") , 
                                                          alpha=TRUE)(n= 12 ),1)) + #边框跟填充色一致
      guides(colour = guide_legend(ncol = 4))
  }
  return(p1)
}




