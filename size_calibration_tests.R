################################
##  Size of calibration tests ##
################################
library(rugarch)
library(GAS)
library(MSGARCH)
library(esback)
source("DGPs.R")


# Setup
risk_level <- 0.025
alpha <- 0.05
model_type <- "GARCH" # "GAS" or "MSGARCH"
MC <- 10000
distri <- "std"
n <- 2500

Table_500 <- matrix(NA, ncol = 10, nrow = MC)
Table_1000 <- matrix(NA, ncol = 10, nrow = MC)
Table_2500 <- matrix(NA, ncol = 10, nrow = MC)

for (i in 1:MC) {
  print(sprintf("Replication %i of %i", i, MC))
  set.seed(i + 1234)
  simulated_data <- simulate_garch(n, 0.01, 0.1, 0.8, 7, risk_level)
  r <- simulated_data$r
  var <- simulated_data$var
  h <- simulated_data$h
  es <- simulated_data$es
  
  test_01_var_500 <- BacktestVaR(r[1:500], var[1:500], alpha = risk_level, Lags = 4)
  test_01_var_1000 <- BacktestVaR(r[1:1000], var[1:1000], alpha = risk_level, Lags = 4)
  test_01_var_2500 <- BacktestVaR(r[1:n], var[1:n], alpha = risk_level, Lags = 4)
  
  test_02_var_500 <- VaR_VQR(r[1:500], var[1:500], risk_level)
  test_02_var_1000 <- VaR_VQR(r[1:1000], var[1:1000], risk_level)
  test_02_var_2500 <- VaR_VQR(r[1:n], var[1:n], risk_level)
  
  test_01_es_500 <- ESTest(alpha = risk_level, r[1:500], es[1:500], var[1:500], conf.level = 0.95,  boot = FALSE)
  test_01_es_1000 <- ESTest(alpha = risk_level, r[1:1000], es[1:1000], var[1:1000], conf.level = 0.95,  boot = FALSE)
  test_01_es_2500 <- ESTest(alpha = risk_level, r[1:n], es[1:n], var[1:n], conf.level = 0.95,  boot = FALSE)
  
  test_02_es_500 <- cc_backtest(r[1:500], var[1:500], es[1:500], sqrt(h[1:500]), alpha = risk_level)
  test_02_es_1000 <- cc_backtest(r[1:1000], var[1:1000], es[1:1000], sqrt(h[1:1000]),alpha = risk_level)
  test_02_es_2500 <- cc_backtest(r[1:n], var[1:n], es[1:n], sqrt(h[1:n]), alpha = risk_level)
  
  test_03_es_500 <- esr_backtest(r[1:500], var[1:500], es[1:500], alpha = risk_level, B = 0, version = 1)
  test_03_es_1000 <- esr_backtest(r[1:1000], var[1:1000], es[1:1000], alpha = risk_level, B = 0, version = 1)
  test_03_es_2500 <- esr_backtest(r[1:n], var[1:n], es[1:n], alpha = risk_level, B = 0, version = 1)
  
  test_04_es_500 <- esr_backtest(r[1:500], var[1:500], es[1:500], alpha = risk_level, B = 0, version = 2)
  test_04_es_1000 <- esr_backtest(r[1:1000], var[1:1000], es[1:1000], alpha = risk_level, B = 0, version = 2)
  test_04_es_2500 <- esr_backtest(r[1:n], var[1:n], es[1:n], alpha = risk_level, B = 0, version = 2)
  
  test_05_es_500 <- esr_backtest(r[1:500], var[1:500], es[1:500], alpha = risk_level, B = 0, version = 3)
  test_05_es_1000 <- esr_backtest(r[1:1000], var[1:1000], es[1:1000], alpha = risk_level, B = 0, version = 3)
  test_05_es_2500 <- esr_backtest(r[1:n], var[1:n], es[1:n], alpha = risk_level, B = 0, version = 3)
  
  Table_500[i, ] <- c(ifelse(test_01_var_500$LRuc[2] < 0.05, 1, 0),
                      ifelse(test_01_var_500$LRcc[2] < 0.05, 1, 0),
                      ifelse(test_02_var_500 < 0.05, 1, 0),
                      ifelse(test_01_es_500$p.value < 0.05, 1, 0),
                      ifelse(test_02_es_500$pvalue_twosided_simple < 0.05, 1, 0),
                      ifelse(test_02_es_500$pvalue_twosided_general < 0.05, 1, 0),
                      ifelse(test_03_es_500$pvalue_twosided_asymptotic < 0.05, 1, 0),
                      ifelse(test_04_es_500$pvalue_twosided_asymptotic < 0.05, 1, 0),
                      ifelse(test_05_es_500$pvalue_twosided_asymptotic < 0.05, 1, 0),
                      ifelse(test_05_es_500$pvalue_onesided_asymptotic < 0.05, 1, 0))
  
  Table_1000[i, ] <- c(ifelse(test_01_var_1000$LRuc[2] < 0.05, 1, 0),
                      ifelse(test_01_var_1000$LRcc[2] < 0.05, 1, 0),
                      ifelse(test_02_var_1000 < 0.05, 1, 0),
                      ifelse(test_01_es_1000$p.value < 0.05, 1, 0),
                      ifelse(test_02_es_1000$pvalue_twosided_simple < 0.05, 1, 0),
                      ifelse(test_02_es_1000$pvalue_twosided_general < 0.05, 1, 0),
                      ifelse(test_03_es_1000$pvalue_twosided_asymptotic < 0.05, 1, 0),
                      ifelse(test_04_es_1000$pvalue_twosided_asymptotic < 0.05, 1, 0),
                      ifelse(test_05_es_1000$pvalue_twosided_asymptotic < 0.05, 1, 0),
                      ifelse(test_05_es_1000$pvalue_onesided_asymptotic < 0.05, 1, 0))
  
  Table_2500[i, ] <- c(ifelse(test_01_var_2500$LRuc[2] < 0.05, 1, 0),
                      ifelse(test_01_var_2500$LRcc[2] < 0.05, 1, 0),
                      ifelse(test_02_var_2500 < 0.05, 1, 0),
                      ifelse(test_01_es_2500$p.value < 0.05, 1, 0),
                      ifelse(test_02_es_2500$pvalue_twosided_simple < 0.05, 1, 0),
                      ifelse(test_02_es_2500$pvalue_twosided_general < 0.05, 1, 0),
                      ifelse(test_03_es_2500$pvalue_twosided_asymptotic < 0.05, 1, 0),
                      ifelse(test_04_es_2500$pvalue_twosided_asymptotic < 0.05, 1, 0),
                      ifelse(test_05_es_2500$pvalue_twosided_asymptotic < 0.05, 1, 0),
                      ifelse(test_05_es_2500$pvalue_onesided_asymptotic < 0.05, 1, 0))
  
}

Table <- rbind(apply(Table_500, 2, mean),
               apply(Table_1000, 2, mean),
               apply(Table_2500, 2, mean))
row.names(Table) <- c(500, 1000, 2500)
colnames(Table) <- c("UC", "CC", "DQ", "VDQ", "ES", "CoC", "ESR1", "ESR2", "ESR3", "ESR3_1")
xtable::xtable(Table, digits = 3)
