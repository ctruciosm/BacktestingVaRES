################################
##  Size of calibration tests ##
################################
library(rugarch)
library(GAS)
library(MSGARCH)
library(esback)


# Setup
risk_level <- 0.025
alpha <- 0.05
model_type <- "GARCH" # "GAS" or "MSGARCH"
MC <- 10000
distri <- "std"

for (i in 1:MC) {
  set.seed(i + 1234)
  simulated_data <- DGP(2000, model_type, distribution = distri)
  r <- simulated_data$r
  var <- simulated_data$var
  sigma <- simulated_data$sigma
  es <- simulated_data$es
  
  test_01_var_250 <- BacktestVaR(r[1:250], var[1:250], alpha = risk_level, Lags = 4)
  test_01_var_500 <- BacktestVaR(r[1:500], var[1:500], alpha = risk_level, Lags = 4)
  test_01_var_1000 <- BacktestVaR(r[1:1000], var[1:1000], alpha = risk_level, Lags = 4)
  test_01_var_2000 <- BacktestVaR(r[1:2000], var[1:2000], alpha = risk_level, Lags = 4)
  
  test_02_var_250 <- VaR_VQR(r[1:250], var[1:250], risk_level)
  test_02_var_500 <- VaR_VQR(r[1:500], var[1:500], risk_level)
  test_02_var_1000 <- VaR_VQR(r[1:1000], var[1:1000], risk_level)
  test_02_var_2000 <- VaR_VQR(r[1:2000], var[1:2000], risk_level)
  
  test_01_es_250 <- ESTest(alpha = risk_level, r[1:250], es[1:250], var[1:250], conf.level = 0.95,  boot = TRUE, n.boot = 5000)
  test_01_es_500 <- ESTest(alpha = risk_level, r[1:500], es[1:500], var[1:500], conf.level = 0.95,  boot = TRUE, n.boot = 5000)
  test_01_es_1000 <- ESTest(alpha = risk_level, r[1:1000], es[1:1000], var[1:1000], conf.level = 0.95,  boot = TRUE, n.boot = 5000)
  test_01_es_2000 <- ESTest(alpha = risk_level, r[1:2000], es[1:2000], var[1:2000], conf.level = 0.95,  boot = TRUE, n.boot = 5000)
  
  test_02_es_250 <- cc_backtest(r[1:250], var[1:250], es[1:250], alpha = risk_level)
  test_02_es_500 <- cc_backtest(r[1:500], var[1:500], es[1:500], alpha = risk_level)
  test_02_es_1000 <- cc_backtest(r[1:1000], var[1:1000], es[1:1000], alpha = risk_level)
  test_02_es_2000 <- cc_backtest(r[1:2000], var[1:2000], es[1:2000], alpha = risk_level)
  
  test_03_es_250 <- cc_backtest(r[1:250], var[1:250], es[1:250], s = sigma[1:250], alpha = risk_level)
  test_03_es_500 <- cc_backtest(r[1:500], var[1:500], es[1:500], s = sigma[1:500], alpha = risk_level)
  test_03_es_1000 <- cc_backtest(r[1:1000], var[1:1000], es[1:1000], s = sigma[1:1000], alpha = risk_level)
  test_03_es_2000 <- cc_backtest(r[1:2000], var[1:2000], es[1:2000], s = sigma[1:2000], alpha = risk_level)
  
  test_04_es_250 <- esr_backtest(r[1:250], var[1:250], es[1:250], alpha = risk_level, B = 5000, version = 1)
  test_04_es_500 <- esr_backtest(r[1:500], var[1:500], es[1:500], alpha = risk_level, B = 5000, version = 1)
  test_04_es_1000 <- esr_backtest(r[1:1000], var[1:1000], es[1:1000], alpha = risk_level, B = 5000, version = 1)
  test_04_es_2000 <- esr_backtest(r[1:2000], var[1:2000], es[1:2000], alpha = risk_level, B = 5000, version = 1)
  
  test_05_es_250 <- esr_backtest(r[1:250], var[1:250], es[1:250], alpha = risk_level, B = 5000, version = 2)
  test_05_es_500 <- esr_backtest(r[1:500], var[1:500], es[1:500], alpha = risk_level, B = 5000, version = 2)
  test_05_es_1000 <- esr_backtest(r[1:1000], var[1:1000], es[1:1000], alpha = risk_level, B = 5000, version = 2)
  test_05_es_2000 <- esr_backtest(r[1:2000], var[1:2000], es[1:2000], alpha = risk_level, B = 5000, version = 2)
  
  test_06_es_250 <- esr_backtest(r[1:250], var[1:250], es[1:250], alpha = risk_level, B = 5000, version = 3)
  test_06_es_500 <- esr_backtest(r[1:500], var[1:500], es[1:500], alpha = risk_level, B = 5000, version = 3)
  test_06_es_1000 <- esr_backtest(r[1:1000], var[1:1000], es[1:1000], alpha = risk_level, B = 5000, version = 3)
  test_06_es_2000 <- esr_backtest(r[1:2000], var[1:2000], es[1:2000], alpha = risk_level, B = 5000, version = 3)
  
}



