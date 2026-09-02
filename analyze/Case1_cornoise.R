library(tidyverse)
library(gridExtra)
library(cowplot)

`%ni%` <- Negate(`%in%`)

case1_ics <- read.table("data/ICgrid_case1_alpha2theta20rF2p5to8p5in25.txt", sep = ",")

source("functions.R")

theta = 20
alpha = 2.0
dC = 2.5
dF = 2.5
rC = 1.0
rF = 2.2
F_eq = (theta*(1- alpha) - 1 + sqrt((1 + theta*(alpha - 1))^2 - 4*theta*(alpha - 1 + dF/rF - alpha*dC/rC)))/(2*theta)
C_eq = 1 - dC/(rC*(1 + theta*F_eq))

case1_ptip_df <- read_csv("data/kde/summary/case1_ptip.csv")
case1_ptip_df <- read_csv("data/kde/case1_ptip.csv")

case1_0.05 <- case1_0.05_wide %>% pivot_longer(cols = colnames(case1_0.05_wide)[5:25], names_to = "cor", values_to = "value")
case1_0.1 <- case1_0.1_wide %>% pivot_longer(cols = colnames(case1_0.05_wide)[5:25], names_to = "cor", values_to = "value")
case1_0.15 <- case1_0.15_wide %>% pivot_longer(cols = colnames(case1_0.05_wide)[5:25], names_to = "cor", values_to = "value")
case1_0.2 <- case1_0.2_wide %>% pivot_longer(cols = colnames(case1_0.05_wide)[5:25], names_to = "cor", values_to = "value")

case1_noise <- rbind(case1_0.05, case1_0.1, case1_0.15, case1_0.2)

grid_size <- diff(sort(unique(case1_0.05$Competitor)))[1]

case1_ptip_baseline <- case1_ptip_df %>% filter(cor == 0) %>% rename(baseline = ptip) %>% select(-cor)

case1_ptip_df <- merge(case1_ptip_df, case1_ptip_baseline)

ylims <- case1_ptip_df %>% mutate(diff = abs(ptip - baseline)) %>% group_by(sigma,baseline) %>% 
  summarize(max_diff = max(diff)) %>% mutate(ymin = baseline - max_diff, ymax = baseline + max_diff)

plt1 <- case1_ptip_df %>% filter(sigma == 0.05) %>% ggplot() + aes(x = cor, y = ptip) + geom_point(color = "#2b83ba") + theme_classic() + 
  xlab(expression("Correlation ("*rho*")")) + ylab("Probability of tipping (%)") + 
  geom_hline(yintercept = case1_ptip_df[case1_ptip_df$cor == 0 & case1_ptip_df$sigma == 0.05,]$ptip, color = 'grey55', linetype = 'dashed') + 
  ggtitle(expression(sigma == 0.05)) +
  coord_cartesian(ylim = c(ylims[ylims$sigma == 0.05,]$ymin, ylims[ylims$sigma == 0.05,]$ymax))

plt2 <- case1_ptip_df %>% filter(sigma == 0.1) %>% ggplot() + aes(x = cor, y = ptip) + geom_point(color = "#abdda4") + theme_classic() + 
  xlab(expression("Correlation ("*rho*")")) + ylab("Probability of tipping (%)") + 
  geom_hline(yintercept = case1_ptip_df[case1_ptip_df$cor == 0 & case1_ptip_df$sigma == 0.1,]$ptip, color = 'grey55', linetype = 'dashed') + 
  ggtitle(expression(sigma == 0.1)) +
  coord_cartesian(ylim = c(ylims[ylims$sigma == 0.1,]$ymin, ylims[ylims$sigma == 0.1,]$ymax))

plt3 <- case1_ptip_df %>% filter(sigma == 0.15) %>% ggplot() + aes(x = cor, y = ptip) + geom_point(color = "#fdae61") + theme_classic() + 
  xlab(expression("Correlation ("*rho*")")) + ylab("Probability of tipping (%)") + 
  geom_hline(yintercept = case1_ptip_df[case1_ptip_df$cor == 0 & case1_ptip_df$sigma == 0.15,]$ptip, color = 'grey55', linetype = 'dashed') + 
  ggtitle(expression(sigma == 0.15)) +
  coord_cartesian(ylim = c(ylims[ylims$sigma == 0.15,]$ymin, ylims[ylims$sigma == 0.15,]$ymax))

plt4 <- case1_ptip_df %>% filter(sigma == 0.2) %>% ggplot() + aes(x = cor, y = ptip) + geom_point(color = "#d7191c") + theme_classic() + 
  xlab(expression("Correlation ("*rho*")")) + ylab("Probability of tipping (%)") + 
  geom_hline(yintercept = case1_ptip_df[case1_ptip_df$cor == 0 & case1_ptip_df$sigma == 0.2,]$ptip, color = 'grey55', linetype = 'dashed') + 
  ggtitle(expression(sigma == 0.2))+
  coord_cartesian(ylim = c(ylims[ylims$sigma == 0.2,]$ymin, ylims[ylims$sigma == 0.2,]$ymax))

case1_ptip_df <- case1_ptip_df %>% mutate(perc_change = (ptip-baseline)/baseline*100) %>% arrange(sigma,cor)
case1_ptip_df$pred <- 0

lm_0.05 <- lm(perc_change ~ cor, data = case1_ptip_df[case1_ptip_df$sigma == 0.05,])
lm_0.1 <- lm(perc_change ~ cor, data = case1_ptip_df[case1_ptip_df$sigma == 0.1,])
lm_0.15 <- lm(perc_change ~ cor, data = case1_ptip_df[case1_ptip_df$sigma == 0.15,])
lm_0.2 <- lm(perc_change ~ cor, data = case1_ptip_df[case1_ptip_df$sigma == 0.2,])

case1_ptip_df[case1_ptip_df$sigma == 0.05,]$pred <- predict(lm_0.05, case1_ptip_df[case1_ptip_df$sigma == 0.05,])
case1_ptip_df[case1_ptip_df$sigma == 0.1,]$pred <- predict(lm_0.1, case1_ptip_df[case1_ptip_df$sigma == 0.1,])
case1_ptip_df[case1_ptip_df$sigma == 0.15,]$pred <- predict(lm_0.15, case1_ptip_df[case1_ptip_df$sigma == 0.15,])
case1_ptip_df[case1_ptip_df$sigma == 0.2,]$pred <- predict(lm_0.2, case1_ptip_df[case1_ptip_df$sigma == 0.2,])

plt5 <- case1_ptip_df %>% ggplot() + aes(x = cor, y = (ptip - baseline)/baseline*100, color = as.factor(sigma)) +
  geom_point() + theme_classic() + 
  xlab(expression("Correlation ("*rho*")")) + ylab("% change in tipping risk") + 
  theme(legend.position = "inside",
        legend.position.inside = c(0.9,0.8)) + 
  geom_line(aes(x = cor, y = pred, color = as.factor(sigma)), size = 1.2) + 
  scale_color_manual(expression(sigma),
                     values = c(`0.05` = "#2b83ba", `0.1` = "#abdda4", `0.15` = "#fdae61", `0.2` = "#d7191c"))

pdf("Code/figures/kde/prob_rtip_case1.pdf", height = 8, width = 8)
plot_grid(plot_grid(plt1,plt2,plt3,plt4,nrow = 1), plt5, nrow = 2, rel_heights = c(0.5, 1),
          labels = c("(a)", "(b)"), label_size = 12, label_fontface = "plain") 
dev.off()