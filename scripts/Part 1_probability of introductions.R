# ---------------------------------------
# Get distribution of outbreak sizes among facilities
# Chris Hoover Jan 2021
# ---------------------------------------

rm(list=ls())
library(janitor)
library(dplyr)
library(tidyverse)
library(zoo)
library(ggplot2)
library(ggtext)
library(knitr)
library(kableExtra)
library(readr)
library(extrafont)
library(ggprism)
library(rlang)

#setwd("C:/Users/liliwes/Desktop/mmed2023/covid19 in prisons/prison introductions/Chris")

# Load data (from 0) -------------------
dat <- readRDS("state_prisons_pop_cases_fin2021-02-02.rds") %>% 
  filter(!is.na(Facility) & 
           !grepl("CHCF", Facility) & 
           !grepl("SATF", Facility)) # Both of these facilities had 0 cases. Both seem to be specilized for heatlhcare/treatment, so makes sense

head(dat)

# Identify outbreaks using 10 days wash-out period ----
outbreaks_df10day <- dat %>% 
  dplyr::select(c(Facility:Residents_Active | matches("County"))) %>% 
  group_by(Facility) %>% 
  mutate( # Identify outbreaks as new cases emerging following 10 days with no cases
    new_cases_10day = zoo::rollsum(New_Residents_Confirmed_rmv_neg, k = 10, 
                                   na.pad = T, align = "right"),
    new_cases_10day_lead1 = lead(new_cases_10day),
    outbreak_start = if_else(new_cases_10day == 0 & new_cases_10day_lead1 > 0, 1, 0),
    # Give each outbreak a unique identifier
    outbreak_num = cumsum(if_else(is.na(outbreak_start), 0, outbreak_start)) + outbreak_start*0,
    Facility_Outbreak = paste0(Facility, " Outbreak ", outbreak_num),
    plot7day_cases = if_else(new_cases_10day == 0, NA_real_, New_Residents_Confirmed_7day)) %>% 
  ungroup() %>% 
  filter(!is.na(outbreak_num) & outbreak_num > 0)


# Plot time series with outbreaks delineated ----
outbreaks_df10day %>% 
  mutate(Facility2 = str_replace(Facility, " State Prison", "")) %>% 
  ggplot() +
  geom_line(aes(x = Date, y = plot7day_cases, col = as.factor(outbreak_num))) +
  theme_classic() +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0,10,100,1000)) +
  theme(strip.text = element_text(size = 8)) +
  facet_wrap(facets = "Facility2",
             nrow = 4, ncol = 8, 
             labeller = label_wrap_gen(20)) +
  labs(y = "7-day average of incident resident cases",
       col = "Facility\noutbreak\nnumber")


# Calculate the number of cases for each introduction or outbreak ----
outbreak_size10day <- outbreaks_df10day %>% 
  group_by(Facility_Outbreak) %>% 
  summarise(outbreak_size = sum(New_Residents_Confirmed_rmv_neg))

outbreaks_df10day <- left_join(outbreak_size10day, outbreaks_df10day, by = c('Facility_Outbreak'))


# Community prevalence ----
# Get counties represented in the facility data
FIPS <- unique(outbreaks_df10day %>%
                 group_by(Facility, Facility_Outbreak) %>%
                 mutate(outbreak_sizes = sum(New_Residents_Confirmed_rmv_neg, na.rm = TRUE)) %>%
                 ungroup() %>%
                 filter(outbreak_sizes>0 & outbreak_start==1))

FIPS1 <- unique(dat$County.FIPS)
FIPS1 <- FIPS1[!is.na(FIPS1)]
FIPS1

# Import county (i.e., community) daily COVID cases and population size
# Source: https://usafacts.org/visualizations/coronavirus-covid-19-spread-map/

# COVID-19 cases
county_daily_cases <- read_csv("covid_confirmed_usafacts.csv")%>%
  pivot_longer(5:1227, names_to="date", values_to="cases") %>%
  mutate(`County Name`= str_remove(`County Name`,"County")) %>%
  mutate(`County Name`=trimws(`County Name`)) %>%
  filter(countyFIPS%in%FIPS1)
county_daily_cases

# County population
county_pop <- read_csv("covid_county_population_usafacts.csv") %>%
  filter(countyFIPS%in%FIPS1)
county_pop

# Daily county prevalence
prev <- county_daily_cases %>%
  full_join(county_pop%>%select(countyFIPS, population), by="countyFIPS") %>%
  mutate(community_prev=round(cases/population,3)) %>%
  mutate(date=ymd(date)) %>%
  select(date, countyFIPS, community_prev)
prev

# Merge county prevalence with facility data ----
mdata <-  outbreaks_df10day %>%
  select(Facility_Outbreak, Facility, Date, `Introduction no.` = outbreak_num, `Outbreak size` = outbreak_size, 
         County.FIPS, N0, outbreak_start, new_cases_10day_lead1) %>%
  mutate(County.FIPS=ifelse(is.na(County.FIPS ), 0, County.FIPS ) )%>%
  group_by(Facility) %>%
  mutate(County.FIPS=replace_na( max(County.FIPS))) %>%
  ungroup() %>%
  mutate(County.FIPS=ifelse(Facility%in%"California Institution for Men (CIM)", 6071,
                            ifelse(Facility%in%"Deuel Vocational Institution (DVI)", 6077,
                                   ifelse(Facility%in%"Pelican Bay State Prison (PBSP)", 6015, County.FIPS)))) %>%
  left_join(prev, by=c("Date"="date", "County.FIPS"="countyFIPS" )) %>% 
  mutate(infectious_contacts = N0 * community_prev) %>% 
  group_by(Facility) %>% 
  mutate(no_outbreaks = n_distinct(Facility_Outbreak)) %>% 
  mutate(infectious_contacts = N0 * community_prev) %>% 
  mutate(avg_county_prevalence = mean(community_prev)) %>% 
  mutate(avg_infectious_contacts = N0 * avg_county_prevalence) %>% 
  mutate(no_cases_at_outbreak_start = ifelse(outbreak_start == 1, new_cases_10day_lead1, NA_real_)) %>% 
  select(Facility, County = County.FIPS, `Residents no.` = N0, Date, `County prevalence` = community_prev,
         `Avg. County prevalence` = avg_county_prevalence, `Infectious contacts` = infectious_contacts,
          Facility_Outbreak, `Introduction no.`, `Outbreak size`, `No. introductions` = no_outbreaks,
         `No. cases at introduction` = no_cases_at_outbreak_start, `Avg infectious contacts` = avg_infectious_contacts)

mdata
tabyl(mdata$`No. introductions`)

# Regression for the relationship between COVID-19 prevalence in staff (using community prevalence as proxy) ----
# and introductions or outbreaks

data2 <- distinct(mdata, Facility, .keep_all = TRUE)
kable(data2, format = "html") %>%
  kable_styling(font_size = 10, bootstrap_options = c("striped", "hover"))

# Function for the plots
lmplot <- function (data, x, y, title){
  ggplot(data, aes({{x}}, {{y}})) +
    theme_prism() +
    geom_point(size=4) +
    labs(title = {{title}}) +
    stat_smooth(method = lm) + 
    theme(axis.text.x = element_markdown(color = "grey20", size = 12, family="Calibri"),
          axis.text.y = element_markdown(color = "grey20", size = 12, angle = 0, hjust = 0, vjust = 0, family="Calibri"),
          axis.title.x = element_markdown(size = 12),
          axis.title.y = element_markdown(size = 12),
          legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(plot.title = element_markdown(family = "Calibri", face = "bold", colour = "black", size = 20, hjust =0)) 
}

# Community prevalence versus number of introductions
a <- lmplot(data2, `No. introductions`,  `Avg. County prevalence`,
            "A") +
  labs(y = "Average county prevalence", x = "No. introductions per facility") +
  scale_y_continuous(breaks=seq(0, 0.1, 0.02)) + 
  scale_x_continuous(breaks=seq(0,8, 1))

#testing correlation
cor.test(data2$`No. introductions`, data2$`Avg. County prevalence`)

model <- lm(data2$`No. introductions` ~ data2$`Avg. County prevalence`)
summary(model)

# Infectious contacts versus number of introductions
b <- lmplot(data2, `No. introductions`,  `Avg infectious contacts`,
            "B") +
  labs(y = "Average infectious contacts", x = "No. introductions per facility") +
  scale_y_continuous(breaks=seq(0, 300, 50)) +
  scale_x_continuous(breaks=seq(0,8, 1))

#testing correlation
cor.test(data2$`No. introductions`, data2$`Avg infectious contacts`)

model <- lm(data2$`No. introductions` ~ data2$`Avg infectious contacts`)
summary(model) 

plot(model$residuals, pch = 16, col = "red")

# Explore the relationship between community prevalence (during the outbreaks) and size of outbreaks 
data3 <- mdata %>% ungroup() %>% 
  group_by(Facility_Outbreak) %>% 
  mutate(large_outbreak = ifelse(`Outbreak size` >5, 1, 0)) %>% 
  mutate(avg_outbreak_prevalence = mean(`County prevalence`)) %>% 
  distinct(Facility_Outbreak, .keep_all = T)

kable(data3, format = "html") %>%
  kable_styling(font_size = 10, bootstrap_options = c("striped", "hover"))

# Community prevalence 
c <- lmplot(data3, `Outbreak size`,  `avg_outbreak_prevalence`,
       "C") +
  labs(y = "Average county prevalence", x = "Total cases in an outbreak" ) +
  scale_x_continuous(trans='log2')

cor.test(data3$`Outbreak size`, data3$`avg_outbreak_prevalence`)
model <- lm(data3$`Outbreak size` ~ data3$`avg_outbreak_prevalence`)
summary(model)

# infectious contacts (community prevalence * number of residents)
data5 <- data3 %>% 
  mutate(avg_outbreak_inf_contacts = `Residents no.` * avg_outbreak_prevalence)

d <- lmplot(data5, `Outbreak size`,  avg_outbreak_inf_contacts,
       "D") +
  labs(y = "Average infectious contacts)", x = "Total cases in an outbreak" ) +
  scale_x_continuous(trans='log2') +
  scale_y_continuous(breaks = seq(0, 600, 100))

# combine plots
grid.arrange(a, b, c, d, nrow = 2)

grid.arrange(a, b, nrow = 1)
