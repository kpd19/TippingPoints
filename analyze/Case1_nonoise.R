library(tidyverse)
library(gridExtra)
library(cowplot)

`%ni%` <- Negate(`%in%`)

source("functions.R")

case1_ics <- read.table("Code/data/ICgrid_case1_alpha2theta20rF2p5to8p5in25.txt", sep = ",")

case1_kde_g2.5_1 <- read_csv("Code/data/nostoch/kde_case1_0.0_g2.5_10k_r1.csv")
case1_kde_g2.5_2 <- read_csv("Code/data/nostoch/kde_case1_0.0_g2.5_10k_r2.csv")
case1_kde_g2.5_3 <- read_csv("Code/data/nostoch/kde_case1_0.0_g2.5_10k_r3.csv")
case1_kde_g2.5_4 <- read_csv("Code/data/nostoch/kde_case1_0.0_g2.5_10k_r4.csv")
case1_kde_g2.5_5 <- read_csv("Code/data/nostoch/kde_case1_0.0_g2.5_10k_r5.csv")

case1_kde_g2.5_6 <- read_csv("Code/data/nostoch/kde_case1_0.5_g2.5_10k_r1.csv")
case1_kde_g2.5_7 <- read_csv("Code/data/nostoch/kde_case1_0.5_g2.5_10k_r2.csv")
case1_kde_g2.5_8 <- read_csv("Code/data/nostoch/kde_case1_0.5_g2.5_10k_r3.csv")
case1_kde_g2.5_9 <- read_csv("Code/data/nostoch/kde_case1_0.5_g2.5_10k_r4.csv")
case1_kde_g2.5_10 <- read_csv("Code/data/nostoch/kde_case1_0.5_g2.5_10k_r5.csv")

case1_kde_g2.5_11 <- read_csv("Code/data/nostoch/kde_case1_-0.5_g2.5_10k_r1.csv")
case1_kde_g2.5_12 <- read_csv("Code/data/nostoch/kde_case1_-0.5_g2.5_10k_r2.csv")
case1_kde_g2.5_13 <- read_csv("Code/data/nostoch/kde_case1_-0.5_g2.5_10k_r3.csv")
case1_kde_g2.5_14 <- read_csv("Code/data/nostoch/kde_case1_-0.5_g2.5_10k_r4.csv")
case1_kde_g2.5_15 <- read_csv("Code/data/nostoch/kde_case1_-0.5_g2.5_10k_r5.csv")

theta = 20
alpha = 2.0
dC = 2.5
dF = 2.5
rC = 1.0
rF = 2.2
F_eq = (theta*(1- alpha) - 1 + sqrt((1 + theta*(alpha - 1))^2 - 4*theta*(alpha - 1 + dF/rF - alpha*dC/rC)))/(2*theta)
C_eq = 1 - dC/(rC*(1 + theta*F_eq))

c1_1 <- clean_merge(case1_kde_g2.5_1, case1_ics)
c1_2 <- clean_merge(case1_kde_g2.5_2, case1_ics)
c1_3 <- clean_merge(case1_kde_g2.5_3, case1_ics)
c1_4 <- clean_merge(case1_kde_g2.5_4, case1_ics)
c1_5 <- clean_merge(case1_kde_g2.5_5, case1_ics)

c1_6 <- clean_merge(case1_kde_g2.5_6, case1_ics)
c1_7 <- clean_merge(case1_kde_g2.5_7, case1_ics)
c1_8 <- clean_merge(case1_kde_g2.5_8, case1_ics)
c1_9 <- clean_merge(case1_kde_g2.5_9, case1_ics)
c1_10 <- clean_merge(case1_kde_g2.5_10, case1_ics)

c1_11 <- clean_merge(case1_kde_g2.5_11, case1_ics)
c1_12 <- clean_merge(case1_kde_g2.5_12, case1_ics)
c1_13 <- clean_merge(case1_kde_g2.5_13, case1_ics)
c1_14 <- clean_merge(case1_kde_g2.5_14, case1_ics)
c1_15 <- clean_merge(case1_kde_g2.5_15, case1_ics)

c1_1$round <- 1
c1_2$round <- 2
c1_3$round <- 3
c1_4$round <- 4
c1_5$round <- 5

c1_6$round <- 6
c1_7$round <- 7
c1_8$round <- 8
c1_9$round <- 9
c1_10$round <- 10

c1_11$round <- 11
c1_12$round <- 12
c1_13$round <- 13
c1_14$round <- 14
c1_15$round <- 15

case1_kde <- rbind(c1_1, c1_2, c1_3, c1_4, c1_5, 
                   c1_6, c1_7, c1_8, c1_9, c1_10, 
                   c1_11, c1_12, c1_13, c1_14, c1_15)

case1_kde %>% 
  ggplot() + aes(x = Competitor,y = Facilitator, fill = value) + geom_tile(color = NA) + 
  theme_classic() +
  scale_fill_viridis_c(option = 'plasma') +
  annotate(geom = "point", x = C_eq, y = F_eq, color = "black", size = 0.5) + 
  facet_wrap(~round) + 
  coord_cartesian(xlim = c(0,0.5), ylim = c(0,0.25))

case1_kde %>% group_by(Competitor, Facilitator, tip) %>%
  summarize(mean_value = mean(value)) %>% 
  ggplot() + aes(x = Competitor,y = Facilitator, fill = mean_value) + geom_tile(color = NA) + 
  theme_classic() +
  scale_fill_viridis_c(option = 'plasma') +
  annotate(geom = "point", x = C_eq, y = F_eq, color = "black", size = 0.5)

case1_kde_avg <- case1_kde %>% group_by(Competitor, Facilitator, tip) %>%
  summarize(mean_value = mean(value)) %>% rename(value = mean_value)

write_csv(case1_kde_avg, "Code/data/case1_no_noise_150avg.csv")

grid_size <- diff(sort(unique(case1_kde$Competitor)))[1]

scale_max = case1_kde_avg %>% ungroup() %>% summarize(max = max(value*tip*grid_size^2)) %>% pull(max)

c1_g2_plts <- plot_dyn(case1_kde_avg, scale_max)

c1_g2_p2 <- plot_dyn_p2(case1_kde)

pdf("Code/figures/case1_results_10k_2.5_150sim.pdf",height = 6, width = 15)
c1_g2_plts
dev.off()

calc_ptip(case1_kde_avg)*100

case1_kde_avg %>% filter(tip != -1) %>% ungroup() %>% 
  mutate(val_scale = value*grid_size*grid_size*tip) %>% 
  summarize(ptip = sum(val_scale)) %>% pull(ptip)

pdf("Code/figures/case1_10k_sim.pdf",height = 6, width = 15)
case1_sim %>% filter(time >= 250*10000, time <= 750*10000) %>% 
  ggplot() + aes(x = Competitor, y = Facilitator) + geom_path() + theme_classic()
dev.off()

comp_col <- '#2596be'
fac_col <- "#e59f02"

case1_sim_small <- case1_sim %>% filter(time >= 475*10000, time <= 610*10000)

case1_sim_small <- case1_sim_small %>% arrange(time)

case1_sim_small_subsamp = case1_sim_small[seq(1, nrow(case1_sim_small), 20),]

comp_col <- '#2596be'
fac_col <- "#e59f02"

change_df <- data.frame(x1 = 510, x2 = 535, y1 = -0.025, y2 = 0)

png("Code/figures/case1_ts_10k_coex2.png", height = 1200, width = 1500,  units = "px", res = 300)
case1_sim_small_subsamp %>% ggplot() + geom_line(aes(x = time/10000, y = Competitor, color = "Competitor")) + 
  geom_line(aes(x = time/10000, y = Facilitator, color = "Facilitator")) + 
  scale_color_manual("", values = c("Facilitator" = fac_col, "Competitor" = comp_col)) + 
  theme_bw(base_size = 15) + 
  ylab("Scaled density") + xlab("Time (years)") +
  theme(legend.position = c(0.175,0.9),
        legend.background = element_blank())+ 
  coord_cartesian(xlim = c(490,600)) + 
  geom_rect(data = change_df, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), 
            fill = "grey55", alpha = 0.4, color = "grey55")
dev.off()

png("Code/figures/case1_sim_10k_coex2.png", height = 1200, width = 1500,  units = "px", res = 300)
case1_sim_small_subsamp %>% filter(time >= 490*10000, time <= 600*10000) %>% 
  ggplot() + aes(x = Competitor, y = Facilitator, color = gamma) +
  geom_path() + theme_bw(base_size = 15) + 
  scale_color_scico("d", palette = "roma", limits = c(2.5, 8.5)) + 
  xlab("Scaled Competitor Density") + 
  ylab("Scaled Facilitator Density")
dev.off()

case1_0.5_g2.5_1 <- read_csv("Code/data/kde_case1_0.5_g2.5_10k_r1.csv")
case1_0.5_g2.5_2 <- read_csv("Code/data/kde_case1_0.5_g2.5_10k_r2.csv")

c1_0.5_1 <- clean_merge(case1_0.5_g2.5_1, case1_ics)
c1_0.5_2 <- clean_merge(case1_0.5_g2.5_2, case1_ics)
#c1_g3 <- clean_merge(case1_kde_g3.5, case1_ics)
#c1_g4 <- clean_merge(case1_kde_g4.5, case1_ics)

c1_0.5_1$round <- 1
c1_0.5_2$round <- 2

c1_0.5_1$cor <- 0.5
c1_0.5_2$cor <- 0.5

c1_1$cor <- 0.0
c1_2$cor <- 0.0

case1_cor <- rbind(c1_1, c1_2, c1_0.5_1, c1_0.5_2)

case1_cor %>% group_by(Competitor, Facilitator, tip, cor) %>%
  summarize(mean_value = mean(value)) %>% 
  ggplot() + aes(x = Competitor,y = Facilitator, z = mean_value, group = cor, color = as.factor(cor)) + geom_contour() + 
  theme_classic() +
  scale_color_brewer(palette = 'Dark2') +
  annotate(geom = "point", x = C_eq, y = F_eq, color = "black", size = 0.5) 




# Integrated Squared Difference (ISD): 
#v Calculate the sum of the squared differences over the grid of the two KDEs.

# Kullback-Leibler (KL) Divergence: Measures how much the probability distribution 
# \(P\) diverges from an expected probability distribution \(Q\).