clean_merge <- function(df_kde, df_ics){
  col_df <- data.frame(y1 = paste0("Column",1:100), Competitor = seq(0,1,length.out = 100))
  
  df_kde$Facilitator <- seq(0,1,length.out = 100)
  df_kde <- df_kde %>% pivot_longer(cols = colnames(df_kde)[1:100],names_to = "y1")
  df_kde <- merge(df_kde, col_df)
  df_kde <- df_kde %>% select(-y1)
  
  df_ics <- df_ics %>% rename(Competitor = V1, Facilitator = V2, tip = V3)
  
  df_ics$Competitor <- round(df_ics$Competitor,8)
  df_ics$Facilitator <- round(df_ics$Facilitator,8)
  
  df_kde$Competitor <- round(df_kde$Competitor,8)
  df_kde$Facilitator <- round(df_kde$Facilitator,8)
  
  df_kde <- merge(df_kde, df_ics, all = TRUE)
  
  return(df_kde)
}

plot_dyn <- function(df_kde, scale_max){
  
  grid_size <- diff(sort(unique(df_kde$Competitor)))[1]
  
  tip_col <- "#F7DE37"
  notip_col <- "#461D39"
  
  plt1 <- df_kde %>% filter(tip != -1) %>% mutate(tip2 = ifelse(tip == 1, "tip", "no tip")) %>% 
    ggplot() + aes(x = Competitor, y = Facilitator, color = tip2, fill = tip2) + geom_tile() + theme_classic() +
    scale_color_manual("", values = c("tip" = tip_col, "no tip" = notip_col)) + 
    scale_fill_manual("", values = c("tip" = tip_col, "no tip" = notip_col)) + 
    theme(legend.position = 'top',
          legend.key.width = unit(1.2, "cm")) + 
    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + 
    xlab("Scaled competitor density") + 
    ylab("Scaled facilitator density")
  
  plt2 <- df_kde %>% mutate(val_scale = value*grid_size*grid_size) %>% filter(tip != -1) %>% 
    ggplot() + aes(x = Competitor, y = Facilitator, color = val_scale, fill = val_scale) + geom_tile() + theme_classic() + 
    scale_color_viridis_c("", option = 'turbo')+ 
    scale_fill_viridis_c("", option = 'turbo')+ 
    theme(legend.position = 'top',
          legend.key.width = unit(1.2, "cm"))+ 
    coord_cartesian(xlim = c(0,1), ylim = c(0,1))+ 
    xlab("Scaled competitor density") + 
    ylab("Scaled facilitator density")
  
  plt3 <- df_kde %>% filter(tip != -1) %>% 
    mutate(val_scale = value*grid_size*grid_size*tip) %>% 
    ggplot() + aes(x = Competitor, y = Facilitator, color = val_scale, fill = val_scale) + geom_tile() + theme_classic() + 
    scale_color_viridis_c("",option = 'turbo', limits = c(0.0, scale_max))+ 
    scale_fill_viridis_c("",option = 'turbo', limits = c(0.0, scale_max))+ 
    theme(legend.position = 'top',
          legend.key.width = unit(1.2, "cm"))+ 
    coord_cartesian(xlim = c(0,1), ylim = c(0,1))+ 
    xlab("Scaled competitor density") + 
    ylab("Scaled facilitator density")
  
  plt_obj <- plot_grid(plt1, plt2, plt3, nrow = 1, align = "hv")
  
  return(plt_obj)
}

plot_dyn_cor <- function(df_kde, scale_max, sigma, cor, digits){
  
  grid_size <- diff(sort(unique(df_kde$Competitor)))[1]
  
  tip_col <- "#F7DE37"
  notip_col <- "#461D39"
  
  ptip = round(calc_ptip(df_kde)*100, digits)
  
  plt1 <- df_kde %>% filter(tip != -1) %>% mutate(tip2 = ifelse(tip == 1, "tip", "no tip")) %>% 
    ggplot() + aes(x = Competitor, y = Facilitator, color = tip2, fill = tip2) + geom_tile() + theme_classic() +
    scale_color_manual("", values = c("tip" = tip_col, "no tip" = notip_col)) + 
    scale_fill_manual("", values = c("tip" = tip_col, "no tip" = notip_col)) + 
    theme(legend.position = 'top',
          legend.key.width = unit(1.2, "cm")) + 
    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + 
    xlab("Scaled competitor density") + 
    ylab("Scaled facilitator density") +
    #ggtitle(expression(paste0("Cor = ", cor, " , ", sigma, " = ", sigma, ", ptip = ", ptip))) +
    labs(title = bquote("Correlation" == .(cor) * "," ~ sigma == .(sigma) ~ ", ptip" == .(ptip) *"%"))
  
  plt2 <- df_kde %>% mutate(val_scale = value*grid_size*grid_size) %>% filter(tip != -1) %>% 
    ggplot() + aes(x = Competitor, y = Facilitator, color = val_scale, fill = val_scale) + geom_tile() + theme_classic() + 
    scale_color_viridis_c("", option = 'turbo')+ 
    scale_fill_viridis_c("", option = 'turbo')+ 
    theme(legend.position = 'top',
          legend.key.width = unit(1.2, "cm"))+ 
    coord_cartesian(xlim = c(0,1), ylim = c(0,1))+ 
    xlab("Scaled competitor density") + 
    ylab("Scaled facilitator density")
  
  plt3 <- df_kde %>% filter(tip != -1) %>% 
    mutate(val_scale = value*grid_size*grid_size*tip) %>% 
    ggplot() + aes(x = Competitor, y = Facilitator, color = val_scale, fill = val_scale) + geom_tile() + theme_classic() + 
    scale_color_viridis_c("",option = 'turbo', limits = c(0.0, scale_max))+ 
    scale_fill_viridis_c("",option = 'turbo', limits = c(0.0, scale_max))+ 
    theme(legend.position = 'top',
          legend.key.width = unit(1.2, "cm"))+ 
    coord_cartesian(xlim = c(0,1), ylim = c(0,1))+ 
    xlab("Scaled competitor density") + 
    ylab("Scaled facilitator density")
  
  plt_obj <- plot_grid(plt1, plt2, plt3, nrow = 1, align = "hv")
  
  return(plt_obj)
}


plot_dyn_p2 <- function(df_kde){
  
  grid_size <- diff(sort(unique(df_kde$Competitor)))[1]
  
  plt2 <- df_kde %>% mutate(val_scale = value*grid_size*grid_size) %>% filter(tip != -1, val_scale > 1e-7) %>% 
    ggplot() + aes(x = Competitor, y = Facilitator, color = val_scale , fill =val_scale) + geom_tile() + theme_classic() + 
    scale_color_viridis_c("Likelihood\nof\nbeing in\nstate space", option = 'turbo')+ 
    scale_fill_viridis_c("Likelihood\nof\nbeing in\nstate space", option = 'turbo')+ 
    theme(legend.position = 'top',
          legend.key.width = unit(1.5, "cm"))+ 
    coord_cartesian(xlim = c(0,1), ylim = c(0,1))
  
  return(plt2)
}


calc_ptip <- function(df_kde){
  
  grid_size <- diff(sort(unique(df_kde$Competitor)))[1]
  
  ptip <- df_kde %>% ungroup() %>% filter(tip != -1) %>% 
    mutate(val_scale = value*grid_size*grid_size*tip) %>% 
    summarize(ptip = sum(val_scale)) %>% pull(ptip)
  
  return(ptip)
}
