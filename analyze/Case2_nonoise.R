library(tidyverse)
library(gridExtra)
library(cowplot)
library(scico)

`%ni%` <- Negate(`%in%`)

source('functions.R')

case2_ics <- read.table("Code/data/ICgrid_case2_alpha2theta40rF2p5to8p5in50.txt", sep = ",")

case2_kde_g2.5_1 <- read_csv("Code/data/nostoch/kde_case2_0.0_g2.5_10k_r1.csv")
case2_kde_g2.5_2 <- read_csv("Code/data/nostoch/kde_case2_0.0_g2.5_10k_r2.csv")
case2_kde_g2.5_3 <- read_csv("Code/data/nostoch/kde_case2_0.0_g2.5_10k_r3.csv")
case2_kde_g2.5_4 <- read_csv("Code/data/nostoch/kde_case2_0.0_g2.5_10k_r4.csv")
case2_kde_g2.5_5 <- read_csv("Code/data/nostoch/kde_case2_0.0_g2.5_10k_r5.csv")

case2_kde_g2.5_6 <- read_csv("Code/data/nostoch/kde_case2_0.5_g2.5_10k_r1.csv")
case2_kde_g2.5_7 <- read_csv("Code/data/nostoch/kde_case2_0.5_g2.5_10k_r2.csv")
case2_kde_g2.5_8 <- read_csv("Code/data/nostoch/kde_case2_0.5_g2.5_10k_r3.csv")
case2_kde_g2.5_9 <- read_csv("Code/data/nostoch/kde_case2_0.5_g2.5_10k_r4.csv")
case2_kde_g2.5_10 <- read_csv("Code/data/nostoch/kde_case2_0.5_g2.5_10k_r5.csv")

case2_kde_g2.5_11 <- read_csv("Code/data/nostoch/kde_case2_-0.5_g2.5_10k_r1.csv")
case2_kde_g2.5_12 <- read_csv("Code/data/nostoch/kde_case2_-0.5_g2.5_10k_r2.csv")
case2_kde_g2.5_13 <- read_csv("Code/data/nostoch/kde_case2_-0.5_g2.5_10k_r3.csv")
case2_kde_g2.5_14 <- read_csv("Code/data/nostoch/kde_case2_-0.5_g2.5_10k_r4.csv")
case2_kde_g2.5_15 <- read_csv("Code/data/nostoch/kde_case2_-0.5_g2.5_10k_r5.csv")

theta = 40
alpha = 2.0
dC = 2.5
dF = 2.5
rC = 1.0
rF = 2.2
F_eq = (theta*(1- alpha) - 1 + sqrt((1 + theta*(alpha - 1))^2 - 4*theta*(alpha - 1 + dF/rF - alpha*dC/rC)))/(2*theta)
C_eq = 1 - dC/(rC*(1 + theta*F_eq))

c2_1 <- clean_merge(case2_kde_g2.5_1, case2_ics)
c2_2 <- clean_merge(case2_kde_g2.5_2, case2_ics)
c2_3 <- clean_merge(case2_kde_g2.5_3, case2_ics)
c2_4 <- clean_merge(case2_kde_g2.5_4, case2_ics)
c2_5 <- clean_merge(case2_kde_g2.5_5, case2_ics)

c2_6 <- clean_merge(case2_kde_g2.5_6, case2_ics)
c2_7 <- clean_merge(case2_kde_g2.5_7, case2_ics)
c2_8 <- clean_merge(case2_kde_g2.5_8, case2_ics)
c2_9 <- clean_merge(case2_kde_g2.5_9, case2_ics)
c2_10 <- clean_merge(case2_kde_g2.5_10, case2_ics)

c2_11 <- clean_merge(case2_kde_g2.5_11, case2_ics)
c2_12 <- clean_merge(case2_kde_g2.5_12, case2_ics)
c2_13 <- clean_merge(case2_kde_g2.5_13, case2_ics)
c2_14 <- clean_merge(case2_kde_g2.5_14, case2_ics)
c2_15 <- clean_merge(case2_kde_g2.5_15, case2_ics)

c2_1$round <- 1
c2_2$round <- 2
c2_3$round <- 3
c2_4$round <- 4
c2_5$round <- 5

c2_6$round <- 6
c2_7$round <- 7
c2_8$round <- 8
c2_9$round <- 9
c2_10$round <- 10

c2_11$round <- 11
c2_12$round <- 12
c2_13$round <- 13
c2_14$round <- 14
c2_15$round <- 15

case2_kde <- rbind(c2_1, c2_2, c2_3, c2_4, c2_5,
                   c2_6, c2_7, c2_8, c2_9, c2_10,
                   c2_11, c2_12, c2_13, c2_14, c2_15)

case2_kde %>% group_by(Competitor, Facilitator, tip) %>%
  ggplot() + aes(x = Competitor,y = Facilitator, fill = value) + geom_tile(color = NA) + 
  theme_classic() + facet_wrap(~round) +
  scale_fill_viridis_c(option = 'turbo') + 
  annotate(geom = "point", x = C_eq, y = F_eq, color = "red", size = 0.5)

case2_kde %>% group_by(Competitor, Facilitator, tip) %>%
  summarize(mean_value = mean(value)) %>% 
  ggplot() + aes(x = Competitor,y = Facilitator, fill = mean_value) + geom_tile(color = NA) + 
  theme_classic() +
  scale_fill_viridis_c(option = 'turbo')

case2_kde_avg <- case2_kde %>% group_by(Competitor, Facilitator, tip) %>%
  summarize(mean_value = mean(value)) %>% rename(value = mean_value)

write_csv(case2_kde_avg, "Code/data/case2_no_noise_150avg.csv")

c2_g2_plts <- plot_dyn(case2_kde_avg)

pdf("Code/figures/case2_results_10k_50sim.pdf",height = 6, width = 15)
c2_g2_plts
dev.off()

calc_ptip(case2_kde_avg)

case2_sim_small <- case2_sim %>% filter(time >= 475*10000, time <= 610*10000)

case2_sim_small <- case2_sim_small %>% arrange(time)

case2_sim_small_subsamp = case2_sim_small[seq(1, nrow(case2_sim_small), 20),]

comp_col <- '#2596be'
fac_col <- "#e59f02"

change_df <- data.frame(x1 = 510, x2 = 560, y1 = -0.025, y2 = 0)

png("Code/figures/case2_ts_10k_ex6.png", height = 1200, width = 1500,  units = "px", res = 300)
case2_sim_small_subsamp %>% ggplot() + geom_line(aes(x = time/10000, y = Competitor, color = "Competitor")) + 
  geom_line(aes(x = time/10000, y = Facilitator, color = "Facilitator")) + 
  scale_color_manual("", values = c("Facilitator" = fac_col, "Competitor" = comp_col)) + 
  theme_bw(base_size = 15) + 
  ylab("Scaled density") + xlab("Time (years)") +
  theme(legend.position = c(0.175,0.9),
        legend.background = element_blank()) + 
  coord_cartesian(xlim = c(490,600)) + 
  geom_rect(data = change_df, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), 
            fill = "grey55", alpha = 0.4, color = "grey55")
dev.off()

png("Code/figures/case2_sim_10k_ex3.png", height = 1200, width = 1500,  units = "px", res = 300)
case2_sim_small_subsamp %>% filter(time >= 490*10000, time <= 600*10000) %>% 
  ggplot() + aes(x = Competitor, y = Facilitator, color = gamma) +
  geom_path() + theme_bw(base_size = 15) + 
  scale_color_scico("d", palette = "roma", limits = c(2.5, 8.5)) + 
  xlab("Scaled competitor density") + ylab("Scaled facilitator density")
dev.off()


