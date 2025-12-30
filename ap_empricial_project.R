### Setup ###

library(tidyverse)
library(tidyr)
library(purrr)
library(dplyr)
library(lubridate)
library(quantmod)

### Data Download ###

# Fama-French 3 Factor
ff3 <- read.csv("F-F_Research_Data_Factors.csv") %>%
  mutate(
    date = as.Date(paste0(date, "01"), format = "%Y%m%d"),
    Mkt.RF = Mkt.RF / 100,
    SMB = SMB / 100,
    HML = HML / 100,
    RF = RF / 100
  )

# Stock Market Data
tickers <- c("AAPL", "MSFT", "AMZN", "GOOGL", "META", "JPM")
stock_returns <- map_df(
  tickers,
  \(ticker) {
    prices <- getSymbols(
      ticker,
      src = "yahoo",
      from = "1900-01-01",
      auto.assign = FALSE
    )
    returns <- monthlyReturn(Ad(prices))
    tibble(
      date = as.Date(index(returns)),
      ticker = ticker,
      return = as.numeric(returns)
    )
  }
)

# Merge FF3 with returns (by month)
returns_with_factors <- stock_returns %>%
  mutate(date = floor_date(date, "month")) %>%
  inner_join(ff3, by = "date") %>%
  mutate(excess_ret = return - RF)



### Step 1: Estimate factor betas ###

estimate_betas <- function(returns_with_factors) {
  tickers <- sort(unique(returns_with_factors$ticker))

  map_dfr(
    tickers,
    \(t) {
      model_data <- returns_with_factors %>% filter(ticker == t)

      fit <- lm(excess_ret ~ Mkt.RF + SMB + HML, data = model_data)

      tibble(
        ticker = t,
        alpha = coef(fit)[["(Intercept)"]],
        beta_mkt = coef(fit)[["Mkt.RF"]],
        beta_smb = coef(fit)[["SMB"]],
        beta_hml = coef(fit)[["HML"]],
        n = nobs(fit)
      )
    }
  )
}

betas <- estimate_betas(returns_with_factors)

### Step 2: Estimate lambdas for each date x factor ###

estimate_lambdas <- function(returns_with_factors, betas_df) {
  returns_with_factors %>%
    select(date, ticker, excess_ret) %>%
    left_join(betas_df, by = "ticker") %>%
    group_by(date) %>%
    group_modify(
      \(df, key) {
        df <- df %>% drop_na(excess_ret, beta_mkt, beta_smb, beta_hml)
        fit <- lm(excess_ret ~ beta_mkt + beta_smb + beta_hml, data = df)

        tibble(
          alpha_t = coef(fit)[["(Intercept)"]],
          lambda_mkt = coef(fit)[["beta_mkt"]],
          lambda_smb = coef(fit)[["beta_smb"]],
          lambda_hml = coef(fit)[["beta_hml"]],
          r2 = summary(fit)$r.squared,
          n = nobs(fit)
        )
      }
    ) %>%
    ungroup()
}

lambdas <- estimate_lambdas(returns_with_factors, betas)

### Step 3: Expected risk premia and pricing errors ###

# Expected (avg) lambdas (risk premia) and their variances across time

estimate_expected_lambda <- function(lambdas) {
  lambdas %>%
    summarize(
      alpha_bar = mean(alpha_t, na.rm = TRUE),
      lambda_mkt_bar = mean(lambda_mkt, na.rm = TRUE),
      lambda_smb_bar = mean(lambda_smb, na.rm = TRUE),
      lambda_hml_bar = mean(lambda_hml, na.rm = TRUE),
      var_alpha = var(alpha_t, na.rm = TRUE),
      var_lambda_mkt = var(lambda_mkt, na.rm = TRUE),
      var_lambda_smb = var(lambda_smb, na.rm = TRUE),
      var_lambda_hml = var(lambda_hml, na.rm = TRUE)
    )
}

lambda_bar <- estimate_expected_lambda(lambdas)

# Pricing errors: average residuals when pricing with average lambdas

estimate_pricing_errors <- function(returns_with_factors, lambda_bar) {
  returns_with_factors %>%
    left_join(betas, by = "ticker") %>%
    mutate(
      predicted = lambda_bar$alpha_bar +
                  beta_mkt * lambda_bar$lambda_mkt_bar +
                  beta_smb * lambda_bar$lambda_smb_bar +
                  beta_hml * lambda_bar$lambda_hml_bar,
      alpha_hat = excess_ret - predicted
    ) %>%
    group_by(ticker) %>%
    summarize(
      pe = mean(alpha_hat, na.rm = TRUE),
      var_pe = var(alpha_hat, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
}

pricing_errors <- estimate_pricing_errors(returns_with_factors, lambda_bar)
