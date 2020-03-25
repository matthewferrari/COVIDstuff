library(lubridate)
library(circular)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(ggplot2)

#
#grab data from https://covidtracking.com/api/states/daily.csv

dt <- read.csv("https://covidtracking.com/api/states/daily.csv")
dt$date <-as.Date(trunc(strptime(dt$dateChecked,format="%Y-%m-%d", tz="EST"),"day"))

dt <- tbl_df(dt)

dt <- dt %>% 
  group_by(date,state) %>%
  mutate(test_pos = positive / (positive + negative), test_neg = negative / (positive + negative), tests = (positive + negative)) 
pa <- dt %>%
  ungroup() %>%
  dplyr::filter(state == "PA")  # set a horizontal line for your state of reference

# plot the number of tests per state over time
# size of point indicates the test-positive rate
# dashed line is number of tests in PA for reference
  gg <- ggplot(dt, aes(x=date,y=tests,size=test_pos,color=state)) +
    geom_point() + 
    scale_y_continuous(trans='log10') +
    geom_hline(yintercept=(pa$positive[1] + pa$negative[1]), linetype="dashed", color = "red") +
    facet_wrap(~state,ncol=10) + 
    scale_size_continuous(name = "% test positive") +
    ylab("cumulative tests conducted")
  
  gg

# plot the number of tests per state over time and the number of postive cases over time
# if these lines are diverging, that implies that testing is increasing faster than cases are increasing
# if the lines are converging, that implies that cases are increasing faster than testing
# this could be used to provide an index of expected trend over time in cases
# dashed line is number of tests in PA for reference
ff <- ggplot(dt, aes(x=date,y=tests,color=state)) +
  geom_point() + 
  scale_y_continuous(trans='log10') +
  geom_hline(yintercept=(pa$positive[1] + pa$negative[1]), linetype="dashed", color = "red") +
  facet_wrap(~state,ncol=10) + 
  scale_size_continuous(name = "% test positive") +
  ylab("cumulative tests conducted")+
  geom_point(data=dt,aes(x=date,y=positive,color=state),shape=15)

ff
