library(tidyverse)
library(gridExtra)
library(cowplot)
library(latex2exp)
library(scico)

source("functions.R")

`%ni%` <- Negate(`%in%`)

case1_ics <- read.table("data/ICgrid/ICgrid_case1_alpha2theta20rF2p5to8p5in25.txt", sep = ",")

case1_kde_g2.5_1 <- read_csv("data/kde/tau_case1_0.0_0.0_g2.5_1k.csv")

case1_sim <- read_csv("data/simulations/kde_case1_0.25_2.5to8.5in25_2.csv")
case1_time <- read_csv("data/simulations/kde_case1_0.25_2.5to8.5in25_time_2.csv")
case1_gamma <- read_csv("data/simulations/kde_case1_0.25_2.5to8.5in25_gamma_2.csv")

colnames(case1_sim) <- c("Competitor","Facilitator")
case1_sim$time <- case1_time$Column1
case1_sim$gamma <- case1_gamma$Column1

theta = 20
alpha = 2.0
dC = 2.5
dF = 2.5
rC = 1.0
rF = 2.2
F_eq = (theta*(1- alpha) - 1 + sqrt((1 + theta*(alpha - 1))^2 - 4*theta*(alpha - 1 + dF/rF - alpha*dC/rC)))/(2*theta)
C_eq = 1 - dC/(rC*(1 + theta*F_eq))

case1_kde <- clean_merge(case1_kde_g2.5_1, case1_ics)

case1_kde %>% group_by(Competitor, Facilitator, tip) %>%
  summarize(mean_value = mean(value)) %>% 
  ggplot() + aes(x = Competitor,y = Facilitator, fill = mean_value) + geom_tile(color = NA) + 
  theme_classic() +
  scale_fill_viridis_c(option = 'plasma') +
  annotate(geom = "point", x = C_eq, y = F_eq, color = "black", size = 0.5)
  
case1_kde_avg <- case1_kde %>% group_by(Competitor, Facilitator, tip) %>%
  summarize(mean_value = mean(value)) %>% rename(value = mean_value) %>% mutate(cor = 0.0)

write_csv(case1_kde_avg, "Code/data/kde/summary/case1_0.0.csv")

grid_size <- diff(sort(unique(case1_kde$Competitor)))[1]

scale_max <- case1_kde_avg %>% ungroup() %>% summarize(max = max(value*tip*grid_size^2)) %>% pull(max)

c1_g2_plts <- plot_dyn(case1_kde_avg, scale_max)

c1_g2_p2 <- plot_dyn_p2(case1_kde)

pdf("figures/kde/case1_sig0.pdf",height = 6, width = 15)
c1_g2_plts
dev.off()

calc_ptip(case1_kde_avg)


pdf("figures/case1_10k_sim.pdf",height = 6, width = 15)
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

pdf("figures/kde/case1_ts_10k_coex2.pdf", height = 6, width = 10)
plt1 <- case1_sim_small_subsamp %>% ggplot() + geom_line(aes(x = time/10000, y = Competitor, color = "Competitor")) + 
  geom_line(aes(x = time/10000, y = Facilitator, color = "Facilitator")) + 
  scale_color_manual("", values = c("Facilitator" = fac_col, "Competitor" = comp_col)) + 
  theme_bw(base_size = 15) + 
  ylab("Scaled density") + xlab("Time (years)") +
  theme(legend.position = c(0.175,0.9),
        legend.background = element_blank())+ 
  coord_cartesian(xlim = c(480,600)) + 
  geom_rect(data = change_df, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), 
            fill = "grey55", alpha = 0.4, color = "grey55")

plt2 <- case1_sim_small_subsamp %>% filter(time >= 480*10000, time <= 600*10000) %>% 
  ggplot() + aes(x = Competitor, y = Facilitator, color = gamma) +
  geom_path() + theme_bw(base_size = 15) + 
  scale_color_scico("d", palette = "roma", limits = c(2.5, 8.5)) + 
  xlab("Scaled Competitor Density") + 
  ylab("Scaled Facilitator Density")
plot_grid(plt1,plt2, align = "h", rel_widths = c(1,1.2))
dev.off()