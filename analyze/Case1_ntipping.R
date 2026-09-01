library(tidyverse)

`%ni%` <- Negate(`%in%`)

KC = 10000

case1 <- read_csv("data/extinctions/case1_d2.5.csv")

case1_category <- case1 %>% group_by(sigma,d,cor) %>% count(type) 

case1_ext <- case1_category %>% filter(type != "Coexistence") %>% group_by(sigma,d,cor) %>% summarize(sum_n = sum(n))

ext_col <- "#797979"
coex_col <- "#0f8563"
fac_col <- "#e59f02"

plt1 <- case1_category %>% filter(sigma >0, sigma <= 0.3) %>% 
  ggplot() + aes(x = cor, y = n/1000, group = type, fill = type) + geom_bar(stat = 'identity') + 
  theme_bw() + 
  #scale_fill_brewer(palette = 'Dark2') +
  facet_grid(~sigma, labeller = label_bquote(cols = sigma == .(sigma))) + 
  xlab(expression("Correlation ("*rho*")")) + 
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'red') + 
  ylab("Proportion") +
  theme(legend.position = 'top') + 
  scale_fill_manual("", values = c("Coexistence" = coex_col, "Facilitator only" = fac_col, "Mutual extinction" = ext_col))

case1_ext2 <- data.frame(sigma = 0.05, d = 2.5, sum_n = 0)

case1_ext2 <- merge(case1_ext2, data.frame(cor = seq(-1,1,0.25)))

case1_ext <- rbind(case1_ext, case1_ext2)

plt2 <- case1 %>% #filter(type %in% c("Facilitator only")) %>% 
  filter(type %ni% c("Coexistence"), sigma >0, sigma <= 0.3) %>%
  ggplot() + aes(x = cor, y = t_ext, group = cor, fill = factor(cor)) +
  geom_boxplot(outliers = FALSE) + 
  scale_fill_brewer(expression(rho), palette = 'RdYlBu') +
  theme_bw() + 
  facet_grid(~sigma, labeller = label_bquote(cols = sigma == .(sigma)), scales = 'free_y') + 
  ylab("Time to 1+ species extinction (years)") + 
  xlab(expression("Correlation ("*rho*")")) + 
  geom_text(data = case1_ext[case1_ext$sigma <= 0.3,], aes(x = cor, y = 1100, label = sum_n), vjust = 1, size = 2.5, alpha = 1) + 
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1)) + 
  scale_y_continuous(breaks = c(0,250, 500, 750, 1000)) +
  theme(legend.position = 'top') +
  guides(fill = guide_legend(nrow = 1)) 
plt2

pdf("figures/extinction/case1_extinctions.pdf",height = 8, width = 14)
plot_grid(plt1,plt2,nrow = 2, labels = c("(a)", "(b)"), label_size = 16, label_fontface = "plain")
dev.off()