plotMsigWordcloud_galaxy <- function (msigGsc, groups, weight = NULL, measure = c("tfidf", 
                                                      "tf"), version = msigdb::getMsigdbVersions(), org = c("auto", 
                                                                                                            "hs", "mm"), rmwords = getMsigBlacklist(), type = c("Name", 
                                                                                                                                                                "Short")) 
{
  # checkGroups(groups, names(msigGsc))
  measure = match.arg(measure)
  org = match.arg(org)
  type = match.arg(type)
  names(groups) = sapply(names(groups), function(x) {
    paste0(x, " (n = ", length(groups[[x]]), ")")
  })
  msigGsc_list = lapply(groups, function(x) msigGsc[x])
  worddf = plyr::ldply(msigGsc_list, function(x) {
    df = computeMsigWordFreq(x, weight, measure, version, 
                             org, rmwords)[[type]]
    df$freq = df$freq/max(df$freq)
    df = df[seq_len(min(30, nrow(df))), ]
    df$angle = sample(c(0, 90), nrow(df), replace = TRUE, 
                      prob = c(0.65, 0.35))
    return(df)
  }, .id = "NodeGroup")
  p1 = ggplot(worddf, aes(label = word, size = freq, color = freq, 
                          angle = angle)) + ggwordcloud::geom_text_wordcloud(rm_outside = TRUE, 
                                                                             shape = "circle", eccentricity = 0.65) + 
    ggplot2::facet_wrap(~NodeGroup, ncol = 4,
                                                                                                                                          scales = "free") + scico::scale_colour_scico(palette = "acton", 
                                                                                                                                                                                       direction = -1) + ggplot2::scale_size_area(max_size = 6/log10(1 + 
                                                                                                                                                                                                                                                       length(msigGsc_list))) + bhuvad_theme()
  return(p1)
}





