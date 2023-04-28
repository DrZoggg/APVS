





ggline_galaxy <- function(data = data, color,sizeline,sizepoint,id     ){
   p <-  ggplot(data, aes(x = group2, y = mean) ) +
  geom_line(color = color,size = sizeline ) +
  labs(y='APVSLevel (mean)',x='Status',title = id)+
  geom_point(shape=21, 
             color=color, fill=color, 
             size=sizepoint) +
  theme_bw(base_rect_size = 2)+
  theme(
    axis.text.x = element_blank(),
    # axis.text.x = element_text(size = 10,colour = 'black',face='bold',angle = 20,vjust = 0.6),
    axis.text.y = element_text(size = 10,colour = 'black',face='bold'),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12,colour = 'black',face='bold'),
    # axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    
    panel.grid.major   = element_line(color = "grey92", size = 2, linetype = "solid"),
    panel.background = element_rect(fill=  scales::alpha('white',.8)  ),
    panel.grid.minor    = element_line(color = "grey92", size = 1, linetype = "solid"),
    
    legend.position = c(0.8,0.7),
    legend.key.size   = unit(.8, "cm"),
    legend.background = element_blank(),
    legend.text = element_text(size = 12,colour = 'grey40',face='bold'),
    legend.title = element_text(size = 13,colour = 'grey40',face='bold'),
    #legend.direction = "horizontal",
    # plot.title = element_text(hjust = 0.5,size = 12,colour = 'black',face='bold')
    plot.title = element_blank()
  )+
  # scale_y_continuous(limits = c(min(ggdata_list[[id]]$mean)-0.1,
  #                               max(ggdata_list[[id]]$mean)+0.1)
  # )+
  scale_y_continuous(expand = c(0.1,0.08)
  )+
  scale_x_continuous(expand = c(0.08,0.08)
     )

} 



















