library(tidyverse)
library(gridExtra)
library(cowplot)
library(scico)

`%ni%` <- Negate(`%in%`)

case2_ics <- read.table("data/ICgrid/ICgrid_case2_alpha2theta40rF2p5to8p5in50.txt", sep = ",")

source("functions.R")

ptip_df <- read_csv("data/kde/case2_ptip.csv")

case2_0.05_wide <- read_csv("data/kde/case2_0.05_wide.csv")
case2_0.1_wide <- read_csv("data/kde/case2_0.1_wide.csv")
case2_0.15_wide <- read_csv("data/kde/case2_0.15_wide.csv")

case2_0.05 <- case2_0.05_wide %>% pivot_longer(cols = colnames(case2_0.05_wide)[5:25], names_to = "cor", values_to = "value")
case2_0.1 <- case2_0.1_wide %>% pivot_longer(cols = colnames(case2_0.05_wide)[5:25], names_to = "cor", values_to = "value")
case2_0.15 <- case2_0.15_wide %>% pivot_longer(cols = colnames(case2_0.05_wide)[5:25], names_to = "cor", values_to = "value")

case2_noise <- rbind(case2_0.05, case2_0.1, case2_0.15)

grid_size <- diff(sort(unique(case2_0.05$Competitor)))[1]

case2_0.05 %>% 
  ggplot() + aes(x = Competitor, y = Facilitator, color = value, fill = value) + geom_tile() + 
  theme_classic() + 
  facet_wrap(~cor)

plt1 <- case2_0.05 %>% group_by(Competitor, Facilitator, tip, cor) %>%
  summarize(mean_value = mean(value)) %>% 
  #filter(tip != -1) %>% 
  filter(cor %in% c(-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0)) %>% 
  ggplot() + aes(x = Competitor,y = Facilitator, z = mean_value, group = cor, color = cor) + geom_contour() + 
  theme_classic() +
  scale_color_scico(expression(rho), palette = 'roma')+ 
  coord_cartesian(xlim = c(0,0.5),  ylim = c(0,1.0)) + 
  xlab("Scaled competitor density") + ylab("Scaled facilitator density") +
  theme(legend.position = 'none')

plt2 <- case2_0.1 %>% group_by(Competitor, Facilitator, tip, cor) %>%
  summarize(mean_value = mean(value)) %>% 
  #filter(tip != -1) %>% 
  filter(cor %in% c(-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0)) %>% 
  ggplot() + aes(x = Competitor,y = Facilitator, z = mean_value, group = cor, color = cor) + geom_contour() + 
  theme_classic() +
  scale_color_scico(expression(rho), palette = 'roma')+ 
  coord_cartesian(xlim = c(0,0.5),  ylim = c(0,1.0)) + 
  xlab("Scaled competitor density") + ylab("Scaled facilitator density")+
  theme(legend.position = 'none')

plt3 <- case2_0.15 %>% group_by(Competitor, Facilitator, tip, cor) %>%
  summarize(mean_value = mean(value)) %>% 
  #filter(tip != -1) %>% 
  filter(cor %in% c(-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0)) %>% 
  ggplot() + aes(x = Competitor,y = Facilitator, z = mean_value, group = cor, color = cor) + geom_contour() + 
  theme_classic() +
  scale_color_scico(expression(rho), palette = 'roma')+ 
  coord_cartesian(xlim = c(0,0.5),  ylim = c(0,1.0)) + 
  xlab("Scaled competitor density") + ylab("Scaled facilitator density")+
  theme(legend.position = 'none')

plt4 <- case2_0.2 %>% group_by(Competitor, Facilitator, tip, cor) %>%
  summarize(mean_value = mean(value)) %>%
  #filter(tip != -1) %>% 
  filter(cor %in% c(-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0)) %>% 
  ggplot() + aes(x = Competitor,y = Facilitator, z = mean_value, group = cor, color = cor) + geom_contour() + 
  theme_classic() +
  scale_color_scico(expression(rho), palette = 'roma')+ 
  coord_cartesian(xlim = c(0,1.0),  ylim = c(0,1.0)) + 
  xlab("Scaled competitor density") + ylab("Scaled facilitator density")+
  theme(legend.position = 'right')

pdf("figures/kde/contour_rtip_case2.pdf", height = 4, width = 12)
plot_grid(plt1,plt2,plt3,plt4,nrow = 1, rel_widths = c(1,1,1,1.2))
dev.off()


ptip_baseline <- ptip_df %>% filter(cor == 0) %>% rename(baseline = ptip) %>% select(-cor)

ptip_df <- merge(ptip_df, ptip_baseline, all.x = TRUE)

ylims <- ptip_df %>% mutate(diff = abs(ptip - baseline)) %>% group_by(sigma,baseline) %>% 
  summarize(max_diff = max(diff)) %>% mutate(ymin = baseline - max_diff, ymax = baseline + max_diff)

plt1 <- ptip_df %>% filter(sigma == 0.05) %>% ggplot() + aes(x = cor, y = ptip) + geom_point(color = "#2b83ba") + theme_classic() + 
  xlab(expression("Correlation ("*rho*")")) + ylab("Probability of tipping (%)") + 
  geom_hline(yintercept = ptip_df[ptip_df$cor == 0 & ptip_df$sigma == 0.05,]$ptip, color = 'grey55', linetype = 'dashed') + 
  ggtitle(expression(sigma == 0.05)) +
  coord_cartesian(ylim = c(ylims[ylims$sigma == 0.05,]$ymin, ylims[ylims$sigma == 0.05,]$ymax))

plt2 <- ptip_df %>% filter(sigma == 0.1) %>% ggplot() + aes(x = cor, y = ptip) + geom_point(color = "#abdda4") + theme_classic() + 
  xlab(expression("Correlation ("*rho*")")) + ylab("Probability of tipping (%)") + 
  geom_hline(yintercept = ptip_df[ptip_df$cor == 0 & ptip_df$sigma == 0.1,]$ptip, color = 'grey55', linetype = 'dashed') + 
  ggtitle(expression(sigma == 0.1)) +
  coord_cartesian(ylim = c(ylims[ylims$sigma == 0.1,]$ymin, ylims[ylims$sigma == 0.1,]$ymax))

plt3 <- ptip_df %>% filter(sigma == 0.15) %>% ggplot() + aes(x = cor, y = ptip) + geom_point(color = "#fdae61") + theme_classic() + 
  xlab(expression("Correlation ("*rho*")")) + ylab("Probability of tipping (%)") + 
  geom_hline(yintercept = ptip_df[ptip_df$cor == 0 & ptip_df$sigma == 0.15,]$ptip, color = 'grey55', linetype = 'dashed') + 
  ggtitle(expression(sigma == 0.15)) +
  coord_cartesian(ylim = c(ylims[ylims$sigma == 0.15,]$ymin, ylims[ylims$sigma == 0.15,]$ymax))

plt4 <- ptip_df %>% filter(sigma == 0.2) %>% ggplot() + aes(x = cor, y = ptip) + geom_point(color = "#d7191c") + theme_classic() + 
  xlab(expression("Correlation ("*rho*")")) + ylab("Probability of tipping (%)") + 
  geom_hline(yintercept = ptip_df[ptip_df$cor == 0 & ptip_df$sigma == 0.2,]$ptip, color = 'grey55', linetype = 'dashed') + 
  ggtitle(expression(sigma == 0.2))+
  coord_cartesian(ylim = c(ylims[ylims$sigma == 0.2,]$ymin, ylims[ylims$sigma == 0.2,]$ymax))

case2_ptip_df <- case2_ptip_df %>% mutate(perc_change = (ptip-baseline)/baseline*100) %>% arrange(sigma,cor)
case2_ptip_df$pred <- 0

lm_0.05 <- lm(perc_change ~ cor, data = case2_ptip_df[case2_ptip_df$sigma == 0.05,])
lm_0.1 <- lm(perc_change ~ cor, data = case2_ptip_df[case2_ptip_df$sigma == 0.1,])
lm_0.15 <- lm(perc_change ~ cor, data = case2_ptip_df[case2_ptip_df$sigma == 0.15,])

case2_ptip_df[case2_ptip_df$sigma == 0.05,]$pred <- predict(lm_0.05, case2_ptip_df[case2_ptip_df$sigma == 0.05,])
case2_ptip_df[case2_ptip_df$sigma == 0.1,]$pred <- predict(lm_0.1, case2_ptip_df[case2_ptip_df$sigma == 0.1,])
case2_ptip_df[case2_ptip_df$sigma == 0.15,]$pred <- predict(lm_0.15, case2_ptip_df[case2_ptip_df$sigma == 0.15,])

plt5 <- case2_ptip_df %>% ggplot() + aes(x = cor, y = (ptip - baseline)/baseline*100, color = as.factor(sigma)) +
  geom_point() + theme_classic() + 
  xlab(expression("Correlation ("*rho*")")) + ylab("% change in tipping risk") + 
  theme(legend.position = "inside",
        legend.position.inside = c(0.9,0.8)) + 
  geom_line(aes(x = cor, y = pred, group = sigma, color = as.factor(sigma)), size = 1.2) + 
  scale_color_manual(expression(sigma),
                     values = c(`0.05` = "#2b83ba", `0.1` = "#abdda4", `0.15` = "#fdae61", `0.2` = "#d7191c"))

pdf("Code/figures/kde/prob_rtip_case2.pdf", height = 8, width = 8)
plot_grid(plot_grid(plt1,plt2,plt3,nrow = 1), plt5, nrow = 2, rel_heights = c(0.5, 1)) 
dev.off()
