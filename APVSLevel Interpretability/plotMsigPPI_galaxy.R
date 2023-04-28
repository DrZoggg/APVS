plotMsigPPI_galaxy <- function (ppidf, msigGsc, groups, geneStat = NULL, statName = "Gene-level statistic", 
                            threshConfidence = 0, threshFrequency = 0.25, threshStatistic = 0, 
                            threshUseAbsolute = TRUE, topN = 5, nodeSF = 1, edgeSF = 1, 
                            lytFunc = "graphopt", lytParams = list()) 
{
  # checkGroups(groups, names(msigGsc))
  stopifnot(threshConfidence >= 0)
  stopifnot(threshFrequency >= 0)
  stopifnot(nodeSF > 0)
  stopifnot(edgeSF > 0)
  stopifnot(is.null(geneStat) || !is.null(names(geneStat)))
  gr = computeMsigGroupPPI(ppidf, msigGsc, groups, geneStat, 
                           threshConfidence, threshFrequency, threshStatistic, threshUseAbsolute, 
                           topN)
  lytParams = c(list(graph = gr, layout = lytFunc), lytParams)
  p1 = do.call(ggraph::ggraph, lytParams) + ggraph::geom_edge_link(aes(colour = weight), 
                                                                   show.legend = FALSE, colour = "#666666", lwd = 0.5) + 
    ggraph::geom_node_point(aes(size = Degree, text = label), 
                            colour = "#FFFFFF") + ggraph::geom_node_text(aes(label = Label), 
                                                                         repel = TRUE) + ggraph::geom_node_point(aes(fill = geneStat, 
                                                                                                                     size = Degree), alpha = 0.75, shape = 21, stroke = 0.2 * 
                                                                                                                   nodeSF)
  if (is.null(geneStat)) {
    p1 = p1 + ggplot2::scale_fill_manual(values = c(A = 1), 
                                         na.value = "orchid")
  }
  else {
    if (threshUseAbsolute & !all(geneStat >= 0)) {
      lims = c(-1, 1) * stats::quantile(abs(geneStat), 
                                        0.99)
      palname = "vik"
    }
    else {
      lims = stats::quantile(abs(geneStat), c(0.01, 0.99))
      palname = "tokyo"
    }
    p1 = p1 + scico::scale_fill_scico(palette = palname, 
                                      na.value = "orchid", limits = lims, oob = scales::squish)
  }
  p1 = p1 + 
    ggplot2::facet_wrap(~Group, ncol = 4,scales = "free") + 
    guides(size = guide_legend(ncol = 2)) + ggplot2::scale_size_continuous(range = c(0.1, 
                                                                                     6) * nodeSF) + bhuvad_theme() + ggplot2::theme(legend.position = "bottom", 
                                                                                                                                    plot.title = element_text(hjust = 0.5, size = rel(1.5)), 
                                                                                                                                    panel.background = element_rect(fill = "#FFFFFF", 
                                                                                                                                                                    colour = "#000000"), strip.background = element_rect(colour = "#000000"), 
                                                                                                                                    axis.text = element_blank(), axis.title = element_blank())
  return(p1)
}