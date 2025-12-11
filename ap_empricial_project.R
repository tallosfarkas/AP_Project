### Setup ###

library(tidyverse)
library(quantmod)

### Data Download ###

# Fama-French 3 Factor
ff3 <- read.csv("F-F_Research_Data_Factors_daily.csv") %>%
  mutate(
    date = as.Date(date, format = "%Y%m%d"),
    Mkt.RF = Mkt.RF / 100,
    SMB = SMB / 100,
    HML = HML / 100,
  )

# Stock Market Data
tickers <- c("AAPL","MSFT","AMZN","GOOGL","META","JPM")

x <- getSymbols(tickers, from = "1900-01-01")
adj <-
