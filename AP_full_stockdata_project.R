# ==============================================================================
# PHASE A: CONSTRUCT S&P 500
# ==============================================================================
library(RPostgres)
library(dbplyr)
library(lubridate)
library(tidyverse)

# # 1. Connect
# wrds <- dbConnect(
#   Postgres(),
#   host = "wrds-pgdata.wharton.upenn.edu",
#   dbname = "wrds",
#   port = 9737,
#   sslmode = "require",
#   user = Sys.getenv("WRDS_USER"),
#   password = Sys.getenv("WRDS_PASSWORD")
# )
# 
# # 2. Access Monthly Stock File 
# msf_db <- tbl(wrds, I("crsp.msf")) 
# mse_db <- tbl(wrds, I("crsp.msenames"))
# 
# # 3. Filter and Join
# universe_query <- msf_db |>
#   select(permno, date, ret, prc, shrout) |>
#   filter(date >= "1960-01-01" & date <= "2024-12-31") |>
#   # Join with Event Names to get Share Codes
#   inner_join(
#     mse_db |> select(permno, namestart = namedt, nameend = nameendt, shrcd, exchcd),
#     by = "permno"
#   ) |>
#   # Ensure the date matches the valid name/share code range
#   filter(date >= namestart & date <= nameend) |>
#   # KEEP: Common Stocks (10,11) ONLY
#   filter(shrcd %in% c(10, 11)) |>
#   # KEEP: NYSE (1), AMEX (2), NASDAQ (3) ONLY
#   filter(exchcd %in% c(1, 2, 3))
# 
# # 4. Download the Data
# raw_data <- universe_query |>
#   select(date, permno, ret, prc, shrout) |>
#   collect()

#save(raw_data, file = "sp500_universe_rawdata.RData")

load("sp500_universe_rawdata.RData")

# 5. Construct the "Top 500" Universe Locally

stock_returns <- raw_data |>
  mutate(date = as.Date(date)) |>
  # 1. Calculate Market Cap (Handling bid/ask negatives)
  mutate(mktcap = abs(prc) * shrout) |>
  drop_na(mktcap, ret) |>
  
  # 2. Calculate LAGGED Market Cap
  arrange(permno, date) |>
  group_by(permno) |>
  mutate(mktcap_lag = lag(mktcap)) |>
  ungroup() |>
  
  # 3. Rank and Filter
  group_by(date) |>
  mutate(rank = min_rank(desc(mktcap))) |> 
  filter(rank <= 500) |>
  ungroup() |>
  
  # 4. Final Polish
  mutate(ticker = as.character(permno)) |>
  select(date, ticker, ret, mktcap, mktcap_lag) |>
  arrange(ticker, date)

print(paste("Average stocks per month:", round(mean(table(stock_returns$date)))))


# ==============================================================================
# VALIDATION - CHART COMPARISON
# ==============================================================================
library(quantmod)
library(scales)

# 1. Calculate Synthetic Index LAGGED Market Cap

synthetic_index <- stock_returns |>
  # sort by ticker and date to lag
  arrange(ticker, date) |>
  group_by(ticker) |>
  # Create Lagged Market Cap (Weight at t-1)
  mutate(mktcap_lag = lag(mktcap)) |>
  # Remove the first month for each stock
  drop_na(mktcap_lag, ret) |>
  group_by(date) |>
  # Weight by the MARKET CAP AT START OF MONTH
  summarise(
    synthetic_ret = weighted.mean(ret, mktcap_lag, na.rm = TRUE), 
    .groups = "drop"
  ) |>
  mutate(date = floor_date(date, "month"))

# 2. Download Official S&P 500 (Price Index)
# Note: We compare against ^GSPC (Price)
getSymbols("^GSPC", src = "yahoo", from = "1960-01-01", to = "2024-12-31")
official_index <- monthlyReturn(Ad(GSPC))
official_df <- data.frame(date = index(official_index), official_ret = as.numeric(official_index)) |>
  mutate(date = floor_date(date, "month"))

# 3. Merge and Normalize
validation_data <- synthetic_index |>
  inner_join(official_df, by = "date") |>
  arrange(date) |>
  mutate(
    # Cumulative Wealth (Log Scale)
    Wealth_Synthetic = 100 * cumprod(1 + synthetic_ret),
    Wealth_Official  = 100 * cumprod(1 + official_ret)
  ) |>
  select(date, Wealth_Synthetic, Wealth_Official) |>
  pivot_longer(cols = -date, names_to = "Index", values_to = "Value")

# 4. Plot
ggplot(validation_data, aes(x = date, y = Value, color = Index)) +
  geom_line(linewidth = 0.8) +
  scale_y_log10(labels = comma) + 
  scale_color_manual(values = c("Wealth_Official" = "black", "Wealth_Synthetic" = "red")) +
  labs(
    title = "Validation: Synthetic (Total Return) vs Official (Price Return)",
    subtitle = "Red > Black due to Dividends, but same shape.",
    y = "Index Value (Log Scale)"
  ) +
  theme_minimal()

# 5. Correlation Check
cor_check <- synthetic_index |>
  inner_join(official_df, by = "date") |>
  summarise(correlation = cor(synthetic_ret, official_ret))

print(paste("Correlation:", round(cor_check$correlation, 4)))

# ==============================================================================
# PHASE B: BETA ESTIMATION (TIME SERIES REGRESSION)
# ==============================================================================
library(frenchdata)
library(broom)

# 1. Download Fama-French 3 Factors (if not already loaded)
# We need to match the dates with your stock_returns
message("Fetching FF3 Factors...")
ff_factors <- download_french_data("Fama/French 3 Factors")$subsets$data[[1]] |>
  mutate(
    date = floor_date(ymd(paste0(date, "01")), "month"),
    across(c(`Mkt-RF`, SMB, HML, RF), ~as.numeric(.) / 100)
  ) |>
  rename(mkt_excess = `Mkt-RF`, smb = SMB, hml = HML, rf = RF) |>
  select(date, mkt_excess, smb, hml, rf)

# 2. Join Returns with Factors
data_for_betas <- stock_returns |>
  # Align dates to be safe (ensure both are first of month)
  mutate(date = floor_date(date, "month")) |>
  inner_join(ff_factors, by = "date") |>
  mutate(excess_ret = ret - rf) |>
  drop_na(excess_ret, mkt_excess, smb, hml)

# 3. Estimate Betas Efficiently
# We filter for stocks that have at least 24 months of data to ensure stability
message("Estimating Betas for ", length(unique(data_for_betas$ticker)), " stocks...")

stock_betas <- data_for_betas |>
  group_by(ticker) |>
  filter(n() >= 24) |> 
  summarise(
    # Run regression for each stock
    model = list(lm(excess_ret ~ mkt_excess + smb + hml)), 
    .groups = "drop"
  ) |>
  # Extract coefficients cleanly
  mutate(coefs = map(model, tidy)) |>
  unnest(coefs) |>
  select(ticker, term, estimate) |>
  pivot_wider(names_from = term, values_from = estimate) |>
  rename(alpha = `(Intercept)`, beta_mkt = mkt_excess, beta_smb = smb, beta_hml = hml)

print(head(stock_betas))







# ==============================================================================
# PHASE C: RISK PREMIA ESTIMATION (CROSS-SECTIONAL REGRESSION)
# ==============================================================================

# 1. Join Betas back to the monthly data
# Note: We are using "Full Sample Betas" here
# Advanced students might use "Rolling Betas"
fmb_data <- data_for_betas |>
  inner_join(stock_betas, by = "ticker")

# 2. Run Cross-Sectional Regression for EACH Month
message("Running Cross-Sectional Regressions...")

fmb_lambdas <- fmb_data |>
  group_by(date) |>
  # We need enough stocks in a single month to run a regression (e.g., > 10)
  filter(n() > 10) |>
  summarise(
    model = list(lm(excess_ret ~ beta_mkt + beta_smb + beta_hml)),
    .groups = "drop"
  ) |>
  mutate(coefs = map(model, tidy)) |>
  unnest(coefs) |>
  select(date, term, estimate) |>
  pivot_wider(names_from = term, values_from = estimate) |>
  rename(lambda_0 = `(Intercept)`, lambda_mkt = beta_mkt, lambda_smb = beta_smb, lambda_hml = beta_hml)

print(head(fmb_lambdas))







# ==============================================================================
# PHASE D: FINAL STATISTICS
# ==============================================================================

final_stats <- fmb_lambdas |>
  summarise(
    # 1. Average Risk Premia (Lambda)
    mean_lambda_0   = mean(lambda_0),
    mean_lambda_mkt = mean(lambda_mkt),
    mean_lambda_smb = mean(lambda_smb),
    mean_lambda_hml = mean(lambda_hml),
    
    # 2. T-Statistics (Mean / Standard Error)
    # SE = SD / sqrt(T)
    t_lambda_0   = mean(lambda_0) / (sd(lambda_0) / sqrt(n())),
    t_lambda_mkt = mean(lambda_mkt) / (sd(lambda_mkt) / sqrt(n())),
    t_lambda_smb = mean(lambda_smb) / (sd(lambda_smb) / sqrt(n())),
    t_lambda_hml = mean(lambda_hml) / (sd(lambda_hml) / sqrt(n()))
  ) |>
  pivot_longer(everything(), names_to = "stat", values_to = "value")

print(final_stats)




# ==============================================================================
# PHASE E: STEP 4 - SUB-PERIOD ANALYSIS
# ==============================================================================