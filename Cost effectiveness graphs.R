rm(list=ls())
library(tidyverse)

CASES_AVERTED_THRESHOLD_ARR <- 2e5/c(1500, 2500, 3500)
cases_averted_res <- read_csv("cases_averted_res.csv")

###############################
# Code needed for processing
###############################
MECHANISM_TABLE <- tribble(
  ~mech_num, ~"Mechanism of control", ~N_factor, ~R_factor, ~rho_factor,
  1, "Depopulation\n(freq-dep)", 1, 0,0,
  2, "Vaccination\n(freq-dep)", 0, 0, 1,
  3, "Density-dep\ncontrol", 1, 1, 0,
  4, "Reduced\ntransmission only", 0, 1, 0 
)

MECHANISM_TABLE2 <- tribble(
  ~mech_num, ~"Mechanism of control", ~N_factor, ~R_factor, ~rho_factor,
  3, "Density-dep\ncontrol", 1, 1, 0,
  4, "Reduced\ntransmission only", 0, 1, 0 
)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

man_color <- gg_color_hue(4)
man_line <- MECHANISM_TABLE$mech_num[c(2, 1, 3, 4)]

###############################
# PLOTS
#
# Note that for now I've filtered our R = 1.5 (purely for visualization)
#
###############################

g_expected_cases <- ggplot(cases_averted_res %>% filter (R > 1.5)) +
  facet_wrap(~R) +
  geom_line(aes(x=control, y = mean_inf, col = as.factor(mech_num), linetype = as.factor(mech_num)), size = 1) +
  scale_linetype_manual(name = NULL, values = man_line,  labels = MECHANISM_TABLE$`Mechanism of control`) +
  #  scale_shape_manual(name = 'Baseline R', values = man_shape) +
  xlab(expression('Degree of control,' ~ gamma)) + ylab('Expected cases\nper 100 days') + expand_limits(y=0) +
  guides(color = "none", shape = "none") + theme(legend.position = "top") +
  theme(legend.text = element_text(size = 10),legend.title = element_text(size = 12),text = element_text(size=16))
(g_expected_cases)
ggsave('expected_cases.jpg',plot = g_expected_cases, scale = 1.25)

g_cases_averted <- ggplot(cases_averted_res %>% filter (R > 1.5, mech_num > 2)) +
  geom_line(aes(x=control, y = cases_averted, col = as.factor(R), linetype = as.factor(mech_num)), size = 1) +
  scale_linetype_manual(name = NULL, values = man_line,  labels = MECHANISM_TABLE2$`Mechanism of control`) +
  #  scale_shape_manual(name = 'Baseline R', values = man_shape) +
  xlab(expression('Degree of control,' ~ gamma)) + ylab('Expected cases\naverted per 100 days') + expand_limits(y=0) +
  theme(legend.position = "top") + labs(colour="R") +
  geom_hline(yintercept = CASES_AVERTED_THRESHOLD_ARR) +
  theme(legend.text = element_text(size = 10),legend.title = element_text(size = 12),text = element_text(size=16))
(g_cases_averted)
ggsave('cases_averted.jpg',plot = g_cases_averted, scale = 1.25)
