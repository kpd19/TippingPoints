library(tidyverse)
library(gridExtra)
library(cowplot)
library(scico)

source('functions.R')

`%ni%` <- Negate(`%in%`)

case2_ics <- read.table("data/ICgrid/ICgrid_case2_alpha2theta40rF2p5to8p5in50.txt", sep = ",")

case2_kde_g2.5_1 <- read_csv("data/tau_case2_0.0_0.0_g2.5_1k.csv")

case2_ext_sim <- read_csv("data/simulations/kde_case2_sim_0.0_N2_2.csv")
case2_ext_time <- read_csv("data/simulations/kde_case2_sim_0.0_t_2.csv")
case2_ext_d <- read_csv("data/simulations/kde_case2_sim_0.0_d_2.csv")

case2_coex_sim <- read_csv("data/simulations/kde_case2_sim_0.0_N2_3.csv")
case2_coex_time <- read_csv("data/simulations/kde_case2_sim_0.0_t_3.csv")
case2_coex_d <- read_csv("data/simulations/kde_case2_sim_0.0_d_3.csv")

colnames(case2_ext_sim) <- c("Competitor","Facilitator")
case2_ext_sim$time <- case2_ext_time$Column1
case2_ext_sim$gamma <- case2_ext_d$Column1

colnames(case2_coex_sim) <- c("Competitor","Facilitator")
case2_coex_sim$time <- case2_coex_time$Column1
case2_coex_sim$gamma <- case2_coex_d$Column1

theta = 40
alpha = 2.0
dC = 2.5
dF = 2.5
rC = 1.0
rF = 2.2
F_eq = (theta*(1- alpha) - 1 + sqrt((1 + theta*(alpha - 1))^2 - 4*theta*(alpha - 1 + dF/rF - alpha*dC/rC)))/(2*theta)
C_eq = 1 - dC/(rC*(1 + theta*F_eq))
  
#####################
#
# Calculating r-tipping risk
#
#####################

case2_kde <- clean_merge(case2_kde_g2.5_1, case2_ics)

case2_kde %>% group_by(Competitor, Facilitator, tip) %>%
  ggplot() + aes(x = Competitor,y = Facilitator, fill = value) + geom_tile(color = NA) + 
  theme_classic() + 
  scale_fill_viridis_c(option = 'turbo') + 
  annotate(geom = "point", x = C_eq, y = F_eq, color = "red", size = 0.5)

grid_size <- diff(sort(unique(case2_kde$Competitor)))[1]

scale_max <- case2_kde %>% ungroup() %>% summarize(max = max(value*tip*grid_size^2)) %>% pull(max)

c2_g2_plts <- plot_dyn(case2_kde, scale_max)

pdf("Code/figures/kde/case2_sig0.pdf",height = 6, width = 15)
c2_g2_plts
dev.off()

calc_ptip(case2_kde)

case2_sim_small <- case2_sim %>% filter(time >= 475*10000, time <= 610*10000)

case2_sim_small <- case2_sim_small %>% arrange(time)

case2_sim_small_subsamp = case2_sim_small[seq(1, nrow(case2_sim_small), 20),]

#####################
#
# Plotting simulation
#
#####################

comp_col <- '#2596be'
fac_col <- "#e59f02"

change_df <- data.frame(x1 = 510, x2 = 560, y1 = -0.025, y2 = 0)

pdf("figures/kde/case2_ts_10k_coex2.pdf", height = 6, width = 12)
plt1 <- case2_coex_sim %>% ggplot() + geom_line(aes(x = time/10000, y = Competitor, color = "Competitor")) + 
  geom_line(aes(x = time/10000, y = Facilitator, color = "Facilitator")) + 
  scale_color_manual("", values = c("Facilitator" = fac_col, "Competitor" = comp_col)) + 
  theme_bw(base_size = 15) + 
  ylab("Scaled density") + xlab("Time (years)") +
  theme(legend.position = c(0.8,0.94),
        legend.background = element_blank())+ 
  coord_cartesian(xlim = c(480,700)) + 
  geom_rect(data = change_df, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), 
            fill = "grey55", alpha = 0.4, color = "grey55")

plt2 <- case2_coex_sim %>% filter(time >= 480*10000, time <= 700*10000) %>% 
  ggplot() + aes(x = Competitor, y = Facilitator, color = gamma) +
  geom_path() + theme_bw(base_size = 15) + 
  scale_color_scico("d", palette = "roma", limits = c(2.5, 8.5)) + 
  xlab("Scaled Competitor Density") + 
  ylab("Scaled Facilitator Density")
plot_grid(plt1,plt2, align = "h", rel_widths = c(1,1.2))
dev.off()

pdf("figures/kde/case2_ts_10k_extinct.pdf", height = 6, width = 12)
plt1 <- case2_ext_sim %>% ggplot() + geom_line(aes(x = time/10000, y = Competitor, color = "Competitor")) + 
  geom_line(aes(x = time/10000, y = Facilitator, color = "Facilitator")) + 
  scale_color_manual("", values = c("Facilitator" = fac_col, "Competitor" = comp_col)) + 
  theme_bw(base_size = 15) + 
  ylab("Scaled density") + xlab("Time (years)") +
  theme(legend.position = c(0.8,0.94),
        legend.background = element_blank())+ 
  coord_cartesian(xlim = c(480,700)) + 
  geom_rect(data = change_df, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), 
            fill = "grey55", alpha = 0.4, color = "grey55")

plt2 <- case2_ext_sim %>% filter(time >= 480*10000, time <= 700*10000) %>% 
  ggplot() + aes(x = Competitor, y = Facilitator, color = gamma) +
  geom_path() + theme_bw(base_size = 15) + 
  scale_color_scico("d", palette = "roma", limits = c(2.5, 8.5)) + 
  xlab("Scaled Competitor Density") + 
  ylab("Scaled Facilitator Density")
plot_grid(plt1,plt2, align = "h", rel_widths = c(1,1.2))
dev.off()

