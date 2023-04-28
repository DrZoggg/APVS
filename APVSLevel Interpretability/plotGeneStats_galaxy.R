

plotGeneStats_galaxy <- function (geneStat, msigGsc, groups, statName = "Gene-level statistic", 
          topN = 5) 
{
  # checkGroups(groups, names(msigGsc))
  genefreq = plyr::ldply(groups, function(x) {
    gc = table(unlist(lapply(msigGsc[x], GSEABase::geneIds)))
    gc = data.frame(Gene = names(gc), Count = as.numeric(gc))
    return(gc)
  }, .id = "Group")
  genefreq$GeneStat = geneStat[genefreq$Gene]
  stopifnot(any(!is.na(genefreq$GeneStat)))
  genefreq = plyr::ddply(genefreq, "Group", function(x) {
    x = x[!is.na(x$GeneStat) & !is.infinite(x$GeneStat), 
    ]
    st = rank(abs(x$GeneStat)) * rank(x$Count)
    x$rank = rank(-st)
    return(x)
  })
  p1 = ggplot(genefreq, 
              aes(Count, GeneStat)) + 
    ggplot2::geom_jitter(data = genefreq[genefreq$rank > topN, ], size = 1.5,
                          colour = "gray70") + 
    ggplot2::geom_jitter(data = genefreq[genefreq$rank <= 
                                           topN, ], 
                         size = 2.8,
                         colour = '#2C8EAE'
                         ) + 
    ggrepel::geom_text_repel(aes(label = Gene),  data = genefreq[genefreq$rank <= topN, ]) + ggplot2::xlab("Gene-set count") + 
    ggplot2::ylab(statName) + 
    ggplot2::facet_wrap(~Group, ncol = 4,
                                                  scales = "free_x") + bhuvad_theme()+
    geom_hline(yintercept = 0, colour = '#F06563', lty = 2,lwd = .6)
  return(p1)
}


