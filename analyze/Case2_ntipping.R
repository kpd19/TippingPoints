library(tidyverse)
library(ggh4x)

`%ni%` <- Negate(`%in%`)

KC = 10000

dir <- "Code/data/extinctions/case2/"

case2_files <- c(list.files(dir,pattern = ".csv"))

case2 <- c()
for(i in 1:length(case2_files)){
  temp <- read_csv(paste0(dir, case2_files[i]))
  
  case2 <- rbind(case2, temp)
}

case2_zero <- case2 %>% filter(cor == 0, sigma ==0)
case2_zero <- case2_zero %>% select(-cor)
unique(case2_zero$d)

case2 <- case2 %>% filter(sigma !=0)

cor_vals <- data.frame(cor = seq(-1,1,0.25))

case2_zero <- merge(case2_zero,cor_vals)

case2 <- rbind(case2,case2_zero)

case2 <- case2 %>% mutate(type = case_when(t_ext == 999.99 ~ "Coexistence",
                                           t_ext < 999.99 & F_end >= 0.8 ~ "Facilitator only",
                                           t_ext < 999.99 & F_end <= 0.8 ~ "Mutual extinction")) 



write_csv(case2, "Code/data/extinctions/summary/case2.csv")

case2 <- read_csv("Code/data/extinctions/summary/case2.csv")

case2 <- case2 %>% filter(d == 2.5)

case2_category <- case2 %>% group_by(sigma,d,cor) %>% count(type) 

case2_ext <- case2_category %>% filter(type != "Coexistence") %>% group_by(sigma,d,cor) %>% summarize(sum_n = sum(n))

ext_col <- "#797979"
coex_col <- "#0f8563"
fac_col <- "#e59f02"

pdf("Code/figures/extinction/case2_sigma_d2.5_2.pdf",height = 4, width = 12)
plt1 <- case2_category %>% filter(sigma >0, sigma <= 0.3) %>% 
  ggplot() + aes(x = cor, y = n/1000, group = type, fill = type) + geom_bar(stat = 'identity') + 
  theme_bw() + 
  #scale_fill_brewer(palette = 'Dark2') +
  facet_grid(~sigma, labeller = label_bquote(cols = sigma == .(sigma))) + 
  xlab(expression("Correlation ("*rho*")")) + 
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'red') + 
  ylab("Proportion") +
  theme(legend.position = 'top') + 
  scale_fill_manual("", values = c("Coexistence" = coex_col, "Facilitator only" = fac_col, "Mutual extinction" = ext_col))
plt1
dev.off()

case2_category2 <- case2_category %>% filter(cor != 0)

case2_cat_zero <- case2_category %>% filter(cor == 0) %>% rename(baseline = n) %>% ungroup() %>% select(-cor)

case2_category2 <- merge(case2_category2, case2_cat_zero, all = TRUE)
case2_category2 <- case2_category2 %>% mutate(n = ifelse(is.na(n), 0, n),
                                              baseline = ifelse(is.na(baseline), 0, baseline))

grid.arrange(plt1,plt2)

case2_category2 %>% filter(sigma == 0.1, d == 2.5) %>% arrange(cor)

pdf("Code/figures/extinction/case2_1k_allrho.pdf",height = 8, width = 8)
uni_d <- unique(case2_category2$d)

for (i in 1:length(uni_d)){
  plt <- case2_category2 %>% filter(sigma >= 0.05, d == uni_d[i]) %>% 
    ggplot() + aes(x = cor, y = n - baseline, fill = type) +
    geom_bar(stat = 'identity') + 
    theme_bw() +
    facet_wrap(~sigma, labeller = label_bquote(sigma == .(sigma))) + 
    xlab(expression("Correlation ("*rho*")")) + 
    ylab("Difference with positive correlation") + 
    theme(legend.position = "none") +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey55') + 
    scale_fill_manual(values = c("Coexistence" = coex_col, "Facilitator only" = fac_col, "Mutual extinction" = ext_col)) + 
    ggtitle(paste0("d = ", uni_d[i]))
  print(plt)
}
dev.off()

pdf("Code/figures/extinction/case2_t_ext1.pdf",height = 4, width = 8)
case2 %>% filter(type != "Coexistence") %>%
  ggplot() + aes(x = t_ext/KC, group = sigma, fill = as.factor(sigma)) + geom_histogram() + 
  scale_fill_brewer(expression(sigma), palette = 'Spectral') +
  theme_bw() + 
  facet_grid(cor~d, scales = 'free_y') + 
  xlab("Time to 1+ species extinction (years)")
dev.off()

pdf("Code/figures/extinction/case2_text_d2.5.pdf",height = 4, width = 12)
case2 %>% #filter(type %in% c("Facilitator only")) %>% 
  filter(type %ni% c("Coexistence"), sigma >0, d != 7, sigma <= 0.3) %>%
  ggplot() + aes(x = cor, y = t_ext, group = cor, fill = factor(cor)) +
  geom_boxplot(outliers = FALSE) + 
  scale_fill_brewer(expression(rho), palette = 'RdYlBu') +
  theme_bw() + 
  facet_grid(~sigma, labeller = label_bquote(cols = sigma == .(sigma)), scales = 'free_y') + 
  #ggh4x::facet_grid2(d~sigma, scales = 'free', independent = "y") + 
  ylab("Time to 1+ species extinction (years)") + 
  xlab(expression("Correlation ("*rho*")")) + 
  geom_text(data = case2_ext[case2_ext$sigma <= 0.3 & case2_ext$sigma >=0.05,], aes(x = cor, y = 1100, label = sum_n), vjust = 1, size = 2, alpha = 1) + 
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1)) + 
  scale_y_continuous(breaks = c(0,250, 500, 750, 1000)) +
  theme(legend.position = 'top') +
  guides(fill = guide_legend(nrow = 1)) 
dev.off()

df1 <- read_csv("/Volumes/My Book/CriticalTransitions/tau_case2_0.0_2.5_extinctions_ptip_N.csv")
df2 <- read_csv("/Volumes/My Book/CriticalTransitions/tau_case2_0.0_2.5_extinctions_ptip_R.csv")

df1 <- merge(df1,df2)

df1 <- df1 %>% mutate(type = case_when(t_ext == 999.99 ~ "Coexistence",
                                       t_ext < 999.99 & F_end >= 0.8 ~ "Facilitator only",
                                       t_ext < 999.99 & F_end <= 0.8 ~ "Mutual extinction")) 


pdf("Code/figures/trajectories.pdf", height = 5, width = 5)
df1 %>% filter(type != "Coexistence") %>%
  ggplot() + aes(x = C, y = F, color = type, group = rep) + geom_path() + theme_bw() +
  scale_color_brewer("Tipping result", palette = "Set1") +
  theme(legend.position = 'top') + 
  facet_wrap(~type)
dev.off()
