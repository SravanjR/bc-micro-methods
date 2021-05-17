#Clear Environment 
rm(list = ls())

library(glue)
library(R.matlab)
library(tidyverse)
library(EnvStats)
library(priceR)
library(quantreg)
library(dplyr)
library(doBy)

path <- "~/Code"
data_1 = read.csv(glue("{path}/diffindiff_data.csv"))

# Perform Naive Difference-In-Difference to capture the change in price, pre and post treatment
# The treatment group consists of the products operated by US Airways and American Airlines
# The Product Definition consists of a destination and origin combination and ticketing carrier.
# I focus only on the top 20 busiest airports in the US in 2016. 

#Create Product ID for Carrier and Market and 
data_1<- data_1 %>% group_by(MkID,TICKET_CARRIER) %>% mutate(prodCarrID = cur_group_id())

#Post Descriptive Statistics of Data, Mean and Standard Deviation

library(TAM)
# Premerger - Weighted Average Price
pre_df <- subset(data_1, merger == 0)
pre_df_treat <- subset(pre_df, treatment == 1)
pre_df_control <- subset(pre_df, treatment == 0)
m_pre_control = weighted_mean(pre_df_control$BinFare, pre_df_control$weight)
m_pre_treat = weighted_mean(pre_df_treat$BinFare, pre_df_treat$weight)

# Standard Deviation 

sd_pre_control = weighted_sd(pre_df_control$BinFare, pre_df_control$weight)
sd_pre_treat = weighted_sd(pre_df_treat$BinFare, pre_df_treat$weight)

# Quantiles

q_pre_control <- weighted_quantile(pre_df_control$BinFare, pre_df_control$weight, probs=seq(0,1,.25) )
q_pre_treat <- weighted_quantile(pre_df_treat$BinFare, pre_df_treat$weight, probs=seq(0,1,.25) )


# Number of Consumers
cons_pre_control <- sum(pre_df_control$total_passengers_by_product)
cons_pre_treat <- sum(pre_df_treat$total_passengers_by_product)

# Postmerger - Weight Average Price
post_df <- subset(data_1, merger == 1)
post_df_treat <- subset(post_df, treatment == 1)
post_df_control <- subset(post_df, treatment == 0)
m_post_control = weighted_mean(post_df_control$BinFare, post_df_control$weight)
m_post_treat = weighted_mean(post_df_treat$BinFare, post_df_treat$weight)

# SD

sd_post_control = weighted_sd(post_df_control$BinFare, post_df_control$weight)
sd_post_treat = weighted_sd(post_df_treat$BinFare, post_df_treat$weight)

# Quantiles

q_post_control <- weighted_quantile(post_df_control$BinFare, post_df_control$weight, probs=seq(0,1,.25) )
q_post_treat <- weighted_quantile(post_df_treat$BinFare, post_df_treat$weight, probs=seq(0,1,.25) )

# Number of Consumers
cons_post_control <- sum(post_df_control$total_passengers_by_product)
cons_post_treat <- sum(post_df_treat$total_passengers_by_product)


#Create Diff-in-Diff Regression
Model  <- lm(BinFare ~ treatment + as.factor(Year) + I(treatment*merger), data = data_1, weights = weight)


#Display Results of Regression
summary(Model)

library(stargazer)

stargazer(Model)


