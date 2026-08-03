library(tidyverse)
library(gridExtra)
library(cowplot)

`%ni%` <- Negate(`%in%`)

case2_ics <- read.table("Code/data/ICgrid/ICgrid_case2_alpha2theta40rF2p5to8p5in50.txt", sep = ",")

source("functions.R")

theta = 40
alpha = 2.0
dC = 2.5
dF = 2.5
rC = 1.0
rF = 2.2
F_eq = (theta*(1- alpha) - 1 + sqrt((1 + theta*(alpha - 1))^2 - 4*theta*(alpha - 1 + dF/rF - alpha*dC/rC)))/(2*theta)
C_eq = 1 - dC/(rC*(1 + theta*F_eq))

dir <- "Code/data/kde/case2_0.05/"

case2_files <- c(list.files(dir,pattern = ".csv"))
cor_vals <- c(-0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

case2_0.05 <- c()
ptip_0.05 <- c()
for(i in 1:length(case2_files)){
  temp <- read_csv(paste0(dir, case2_files[i]))
  
  temp2 <- clean_merge(temp, case2_ics)
  temp2$sigma <- 0.05
  temp2$cor <- cor_vals[i]
  
  ptip_temp <- data.frame(sigma = 0.05, cor = cor_vals[i], ptip = calc_ptip(temp2)*100)
  
  case2_0.05 <- rbind(case2_0.05, temp2)
  ptip_0.05 <- rbind(ptip_0.05, ptip_temp)
}

dir <- "Code/data/kde/case2_0.1/"

case2_files <- c(list.files(dir,pattern = ".csv"))
cor_vals <- c(-0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

case2_0.1 <- c()
ptip_0.1 <- c()
for(i in 1:length(case2_files)){
  temp <- read_csv(paste0(dir, case2_files[i]))
  
  temp2 <- clean_merge(temp, case2_ics)
  temp2$sigma <- 0.1
  temp2$cor <- cor_vals[i]
  
  ptip_temp <- data.frame(sigma = 0.1, cor = cor_vals[i], ptip = calc_ptip(temp2)*100)
  
  case2_0.1 <- rbind(case2_0.1, temp2)
  ptip_0.1 <- rbind(ptip_0.1, ptip_temp)
}


case2_noise <- rbind(case2_0.05, case2_0.1)

ptip_df <- rbind(ptip_0.05, ptip_0.1)

ptip_baseline <- ptip_df %>% filter(cor == 0) %>% rename(baseline = ptip) %>% select(-cor)

ptip_df <- merge(ptip_df, ptip_baseline)

pdf("Code/figures/kde/prob_rtip_c2.pdf", height = 4, width = 12)
plt1 <- ptip_df %>% filter(sigma == 0.05) %>% ggplot() + aes(x = cor, y = ptip) + geom_point() + theme_classic() + 
  xlab(expression("Correlation ("*rho*")")) + ylab("Probability of tipping (%)") + 
  geom_hline(yintercept = ptip_df[ptip_df$cor == 0 & ptip_df$sigma == 0.05,]$ptip, color = 'red', linetype = 'dashed') + 
  ggtitle(expression(sigma == 0.05)) 
plt2 <- ptip_df %>% filter(sigma == 0.1) %>% ggplot() + aes(x = cor, y = ptip) + geom_point() + theme_classic() + 
  xlab(expression("Correlation ("*rho*")")) + ylab("Probability of tipping (%)") + 
  geom_hline(yintercept = ptip_df[ptip_df$cor == 0 & ptip_df$sigma == 0.1,]$ptip, color = 'red', linetype = 'dashed') + 
  ggtitle(expression(sigma == 0.1)) 
plt3 <- ptip_df %>% ggplot() + aes(x = cor, y = (ptip - baseline)/baseline*100, color = as.factor(sigma), shape = as.factor(sigma)) +
  geom_point() + theme_classic() + 
  xlab(expression("Correlation ("*rho*")")) + ylab("% change in tipping risk") + 
  scale_color_brewer(expression(sigma), palette = "Dark2") + 
  scale_shape_manual(expression(sigma), values = c(16,17)) + 
  theme(legend.position = "inside",
        legend.position.inside = c(0.1,0.9)) + 
  ggtitle("")
grid.arrange(plt1,plt2,plt3,nrow = 1) 
dev.off()

write_csv(case2_noise, "Code/data/kde/summary/case2_kdes.csv")
write_csv(ptip_df, "Code/data/kde/summary/case2_ptip.csv")

case2_wide_prob = case2_kde %>% mutate(value = value*grid_size^2) %>%
  pivot_wider(names_from = cor, values_from = value) %>% 
  rename(cor0.0 = `0`, cor0.5 = `0.5`, `corn0.5` = `-0.5`, cor0.9 = `0.9`, corn0.9 = `-0.9`)

case2_err = case2_wide_prob %>% mutate(diff_0.0_0.5 = (cor0.0 - cor0.5), 
                                       diff_0.0_n0.5 = (cor0.0 - corn0.5),
                                       diff_0.0_0.9 = (cor0.0 - cor0.9), 
                                       diff_0.0_n0.9 = (cor0.0 - corn0.9),
                                       diff_0.5_n0.5 = cor0.5 - corn0.5,
                                       diff_0.9_n0.9 = cor0.9 - corn0.9) 

case2_err %>% group_by(sigma) %>% summarize(mse0_0.5 = sum(diff_0.0_0.5^2),
                                            mse0_n0.5 = sum(diff_0.0_n0.5^2),
                                            mse0_0.9 = sum(diff_0.0_0.9^2),
                                            mse0_n0.9 = sum(diff_0.0_n0.9^2),
                                            mse0.5_n0.5 = sum(diff_0.5_n0.5^2),
                                            mse0.9_n0.9 = sum(diff_0.9_n0.9^2))

calc_ptip(c2_11)*100
calc_ptip(c2_12)*100
calc_ptip(c2_13)*100
calc_ptip(c2_14)*100
calc_ptip(c2_15)*100

scale_max = case2_kde %>% group_by(sigma) %>% summarize(max = max(value*tip*grid_size^2))

c2_0.0_plts <- plot_dyn_cor(c2_11, scale_max[scale_max$sigma == 0.05,]$max, 0.05, cor = 0, digits = 1)
c2_0.5_plts <- plot_dyn_cor(c2_12, scale_max[scale_max$sigma == 0.05,]$max, 0.05, cor = 0.5, digits = 1)
c2_n0.5_plts <- plot_dyn_cor(c2_13, scale_max[scale_max$sigma == 0.05,]$max, 0.05, cor = -0.5, digits = 1)
c2_0.9_plts <- plot_dyn_cor(c2_14, scale_max[scale_max$sigma == 0.05,]$max, 0.05, cor = 0.9, digits = 1)
c2_n0.9_plts <- plot_dyn_cor(c2_15, scale_max[scale_max$sigma == 0.05,]$max, 0.05, cor = 0.9, digits = 1)

pdf("Code/figures/case2_0.05_cor2.pdf",height = 20, width = 15)
grid.arrange(c2_0.0_plts, c2_0.5_plts, c2_n0.5_plts, c2_0.9_plts, c2_n0.9_plts, nrow = 5)
dev.off()

calc_ptip(c2_21)*100
calc_ptip(c2_22)*100
calc_ptip(c2_23)*100
calc_ptip(c2_24)*100
calc_ptip(c2_25)*100

c2_0.0_plts <- plot_dyn_cor(c2_21, scale_max[scale_max$sigma == 0.1,]$max, 0.1, cor = 0, digits = 1)
c2_0.5_plts <- plot_dyn_cor(c2_22, scale_max[scale_max$sigma == 0.1,]$max, 0.1, cor = 0.5, digits = 1)
c2_n0.5_plts <- plot_dyn_cor(c2_23, scale_max[scale_max$sigma == 0.1,]$max, 0.1, cor = -0.5, digits = 1)
c2_0.9_plts <- plot_dyn_cor(c2_24, scale_max[scale_max$sigma == 0.1,]$max, 0.1, cor = 0.9, digits = 1)
c2_n0.9_plts <- plot_dyn_cor(c2_25, scale_max[scale_max$sigma == 0.1,]$max, 0.1, cor = -0.9, digits = 1)

pdf("Code/figures/case2_0.1_cor2.pdf",height = 20, width = 15)
grid.arrange(c2_0.0_plts, c2_0.5_plts, c2_n0.5_plts, c2_0.9_plts, c2_n0.9_plts, nrow = 5)
dev.off()

upper = 6e-4
lower = -6e-4

plt1 <- case2_err %>% 
  ggplot() + aes(x = Competitor,y = Facilitator, fill = diff_0.0_n0.9) + geom_tile(color = NA) + 
  theme_classic() +
  annotate(geom = "point", x = C_eq, y = F_eq, color = "black", size = 0.5) + 
  facet_wrap(~sigma, nrow = 2) + 
  scale_fill_gradient2(expression(Delta*"p"), high = 'blue', low = 'red', mid = 'white', midpoint = 0, limits = c(lower, upper)) + 
  ggtitle("p[0] - p[-0.9]") +
  theme(legend.position = "left")

plt2 <- case2_err %>% 
  ggplot() + aes(x = Competitor,y = Facilitator, fill = diff_0.0_n0.5) + geom_tile(color = NA) + 
  theme_classic() +
  annotate(geom = "point", x = C_eq, y = F_eq, color = "black", size = 0.5) + 
  facet_wrap(~sigma, nrow = 2) + 
  scale_fill_gradient2(expression(Delta*"p"), high = 'blue', low = 'red', mid = 'white', midpoint = 0, limits = c(lower, upper)) + 
  ggtitle("p[0] - p[-0.5]") +
  theme(legend.position = "none")

plt3 <- case2_err %>% 
  ggplot() + aes(x = Competitor,y = Facilitator, fill = diff_0.0_0.5) + geom_tile(color = NA) + 
  theme_classic() +
  annotate(geom = "point", x = C_eq, y = F_eq, color = "black", size = 0.5) + 
  facet_wrap(~sigma, nrow = 2) + 
  scale_fill_gradient2(expression(Delta*"p"), high = 'blue', low = 'red', mid = 'white', midpoint = 0, limits = c(lower, upper)) + 
  ggtitle("p[0] - p[0.5]") +
  theme(legend.position = "none")

plt4 <- case2_err %>% 
  ggplot() + aes(x = Competitor,y = Facilitator, fill = diff_0.0_0.9) + geom_tile(color = NA) + 
  theme_classic() +
  annotate(geom = "point", x = C_eq, y = F_eq, color = "black", size = 0.5) + 
  facet_wrap(~sigma, nrow = 2) + 
  scale_fill_gradient2(expression(Delta*"p"), high = 'blue', low = 'red', mid = 'white', midpoint = 0, limits = c(lower, upper)) + 
  ggtitle("p[0] - p[0.9]") +
  theme(legend.position = "none")

pdf("Code/figures/case2_kde_differences.pdf",height = 6, width = 10)
grid.arrange(plt1,plt2,plt3,plt4,nrow = 1, widths = c(1,0.8,0.8,0.8))
dev.off()
