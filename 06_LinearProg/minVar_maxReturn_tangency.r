### Assignment 7: Introduction to Portfolio Theory
# To clean up the memory of your current R session run the following line
rm(list=ls(all=TRUE))

# Load relevant packages
library("PerformanceAnalytics")
library("zoo")


## Part 1: Loading in your data set
# Load relevant packages
library("PerformanceAnalytics")
library("zoo")

# Load the data
data <- url("http://s3.amazonaws.com/assets.datacamp.com/course/compfin/lab8.RData")
load(data)
# Explore the data set
head(returns_df)
tail(returns_df)

## Part 2: The CER model
# The returns_df data is preloaded in your workspace

# Estimate the parameters: multivariate
mu_hat_annual <- apply(returns_df,2,mean)*12   
sigma2_annual <- apply(returns_df,2,var)*12 
sigma_annual <- sqrt(sigma2_annual) 
#sigma_annual <- apply(returns_df,2,sd)*12 
cov_mat_annual <- cov(returns_df)*12 
cov_hat_annual <- cov(returns_df)[1,2]*12    
rho_hat_annual <- cor(returns_df)[1,2]

# The annual estimates of the CER model parameters for Boeing and Microsoft
mu_boeing <- mu_hat_annual["rboeing"]
mu_msft <- mu_hat_annual["rmsft"]
sigma2_boeing <-  sigma2_annual["rboeing"]
sigma2_msft <- sigma2_annual["rmsft"]
sigma_boeing <- sigma_annual["rboeing"]
sigma_msft <- sigma_annual["rmsft"]
sigma_boeing_msft <- cov_hat_annual
rho_boeing_msft <- rho_hat_annual

## Part3: A portfolio of Boeing and Microsoft stock
# All data and CER parameters are preloaded in your workspace.
# Type "ls()" in the console to see them.

# The ratio Boeing stock vs Microsoft stock (adds up to 1)
boeing_weights <- seq(from=-1, to=2, by=0.1)
msft_weights <- 1- boeing_weights
  
# Portfolio parameters
mu_portfolio <- boeing_weights*mu_boeing + msft_weights*mu_msft
sigma2_portfolio <- boeing_weights^2 * sigma2_boeing + msft_weights^2 * sigma2_msft + 2 * boeing_weights * msft_weights * sigma_boeing_msft
sigma_portfolio <- sqrt(sigma2_portfolio)

# Plotting the different portfolios
plot(sigma_portfolio, mu_portfolio, type="b", pch=16, ylim=c(0, max(mu_portfolio)), xlim=c(0, max(sigma_portfolio)), xlab=expression(sigma[p]), ylab=expression(mu[p]),col=c(rep("green", 18), rep("red", 13)))
text(x=sigma_boeing, y=mu_boeing, labels="Boeing", pos=4)
text(x=sigma_msft, y=mu_msft, labels="Microsoft", pos=4)

## Part 4: Adding T-bills to your portfolios
# Annual risk-free rate of 3% per year for the T-bill
t_bill_rate <- 0.03
  
# Ratio Boeing stocks
boeing_weights <- seq(from=-1, to=2, by=0.1)

# Portfolio parameters
mu_portfolio_boeing_bill <- t_bill_rate + boeing_weights*(mu_boeing - t_bill_rate)
sigma_portfolio_boeing_bill <- boeing_weights * sigma_boeing
  
# Plot previous exercise
plot(sigma_portfolio, mu_portfolio, type="b", pch=16, ylim=c(0, max(mu_portfolio)), xlim=c(0, max(sigma_portfolio)), xlab=expression(sigma[p]), ylab=expression(mu[p]),col=c(rep("green", 18), rep("red", 13)))
text(x=sigma_boeing, y=mu_boeing, labels="Boeing", pos=4)
text(x=sigma_msft, y=mu_msft, labels="MSFT", pos=4)
# Portfolio Combination Boeing and T-bills
points(sigma_portfolio_boeing_bill, mu_portfolio_boeing_bill, type="b", col="blue", pch=16)

## Part 5: The Sharpe Slope
# Sharp ratio Boeing
sharp_ratio_boeing <- (mu_boeing - t_bill_rate) / sigma_boeing
sharp_ratio_boeing

# Plot previous exercises
plot(sigma_portfolio, mu_portfolio, type="b", pch=16, ylim=c(0, max(mu_portfolio)), xlim=c(0, max(sigma_portfolio)), xlab=expression(sigma[p]), ylab=expression(mu[p]),col=c(rep("green", 18), rep("red", 13)))
text(x=sigma_boeing, y=mu_boeing, labels="Boeing", pos=4)
text(x=sigma_msft, y=mu_msft, labels="MSFT", pos=4)

# Plot the position of the global minimum variance portfolio
text(x=global_min_var_portfolio$sd, y=global_min_var_portfolio$er, labels="Global min", pos=2)

##Part 6: Global Minimum Variance Portfolio

globalMin.portfolio <-
  function(er, cov.mat)
  {
    # Compute global minimum variance portfolio
    #
    # inputs:
    # er				N x 1 vector of expected returns
    # cov.mat		N x N return covariance matrix
    #
    # output is portfolio object with the following elements
    # call			original function call
    # er				portfolio expected return
    # sd				portfolio standard deviation
    # weights		N x 1 vector of portfolio weights
    call <- match.call()
    
    #
    # check for valid inputs
    #
    asset.names <- names(er)
    er <- as.vector(er)					# assign names if none exist
    cov.mat <- as.matrix(cov.mat)
    if(length(er) != nrow(cov.mat))
      stop("invalid inputs")
    if(any(diag(chol(cov.mat)) <= 0))
      stop("Covariance matrix not positive definite")
    # remark: could use generalized inverse if cov.mat is positive semi-definite
    
    #
    # compute global minimum portfolio
    #
    cov.mat.inv <- solve(cov.mat)
    one.vec <- rep(1,length(er))
    #  w.gmin <- cov.mat.inv %*% one.vec/as.vector(one.vec %*% cov.mat.inv %*% one.vec)
    w.gmin <- rowSums(cov.mat.inv) / sum(cov.mat.inv)
    w.gmin <- as.vector(w.gmin)
    names(w.gmin) <- asset.names
    er.gmin <- crossprod(w.gmin,er)
    sd.gmin <- sqrt(t(w.gmin) %*% cov.mat %*% w.gmin)
    gmin.port <- list("call" = call,
                      "er" = as.vector(er.gmin),
                      "sd" = as.vector(sd.gmin),
                      "weights" = w.gmin)
    class(gmin.port) <- "portfolio"
    gmin.port
  }


# The global minimum variance portfolio
# The global minimum variance portfolio
global_min_var_portfolio <- globalMin.portfolio(mu_hat_annual, cov_mat_annual)
global_min_var_portfolio

# Summary of global_min_var_portfolio that takes into account the annual risk-free rate of 3% per year
summary(global_min_var_portfolio, risk.free=0.03)

# Portfolio weights Boeing and Microsoft
# Portfolio weights Boeing and Microsoft
plot(global_min_var_portfolio)

# Plot previous exercises
plot(sigma_portfolio, mu_portfolio, type="b", pch=16, ylim=c(0, max(mu_portfolio)), xlim=c(0, max(sigma_portfolio)), xlab=expression(sigma[p]), ylab=expression(mu[p]),col=c(rep("green", 18), rep("red", 13)))
text(x=sigma_boeing, y=mu_boeing, labels="Boeing", pos=4)
text(x=sigma_msft, y=mu_msft, labels="MSFT", pos=4)

# Plot the position of the global minimum variance portfolio
text(x=global_min_var_portfolio$sd, y=global_min_var_portfolio$er, labels="Global min", pos=2)


## part 7: Question

## Part 8: Tangency Portfolio
# The tangency portfolio
tangency.portfolio <- 
  function(er,cov.mat,risk.free)
  {
    # compute tangency portfolio
    #
    # inputs:
    # er				   N x 1 vector of expected returns
    # cov.mat		   N x N return covariance matrix
    # risk.free		 scalar, risk-free rate
    #
    # output is portfolio object with the following elements
    # call			  captures function call
    # er				  portfolio expected return
    # sd				  portfolio standard deviation
    # weights		 N x 1 vector of portfolio weights
    call <- match.call()
    
    #
    # check for valid inputs
    #
    asset.names <- names(er)
    if(risk.free < 0)
      stop("Risk-free rate must be positive")
    er <- as.vector(er)
    cov.mat <- as.matrix(cov.mat)
    if(length(er) != nrow(cov.mat))
      stop("invalid inputs")
    if(any(diag(chol(cov.mat)) <= 0))
      stop("Covariance matrix not positive definite")
    # remark: could use generalized inverse if cov.mat is positive semi-definite
    
    #
    # compute global minimum variance portfolio
    #
    gmin.port <- globalMin.portfolio(er,cov.mat)
    if(gmin.port$er < risk.free)
      stop("Risk-free rate greater than avg return on global minimum variance portfolio")
    
    # 
    # compute tangency portfolio
    #
    cov.mat.inv <- solve(cov.mat)
    w.t <- cov.mat.inv %*% (er - risk.free) # tangency portfolio
    w.t <- as.vector(w.t/sum(w.t))	# normalize weights
    names(w.t) <- asset.names
    er.t <- crossprod(w.t,er)
    sd.t <- sqrt(t(w.t) %*% cov.mat %*% w.t)
    tan.port <- list("call" = call,
                     "er" = as.vector(er.t),
                     "sd" = as.vector(sd.t),
                     "weights" = w.t)
    class(tan.port) <- "portfolio"
    tan.port
  }

tangency_portfolio <-  tangency.portfolio(mu_hat_annual, cov_mat_annual, risk.free=0.03)

tangency_portfolio

# Summary of tangency_portfolio with annual risk free rate of 3%
summary(tangency_portfolio, risk.free = 0.03)

# Portfolio weights Boeing and Microsoft
plot(tangency_portfolio)

## Part 9: Tangency portfolio and T-bills
# Annual risk-free rate of 3% per year for the T-bill
t_bill_rate <- 0.03

# Set of tangency portfolio weights
tangency_weights <- seq(from=0, to=2, by=0.1)

# Portfolio parameters
# Tangent Portfolio
# Annual risk-free rate of 3% per year for the T-bill
t_bill_rate <- 0.03

# Set of tangency portfolio weights
tangency_weights <- seq(from=0, to=2, by=0.1)

# Portfolio parameters
tangency_portfolio <-  tangency.portfolio(mu_hat_annual, cov_mat_annual, risk.free=0.03)
tangency_portfolio
tangency_portfolio$er
tangency_portfolio$sd

mu_portfolio_tangency_bill <- (1 -tangency_weights) * t_bill_rate + tangency_weights * tangency_portfolio$er
sigma_portfolio_tangency_bill <- tangency_weights * tangency_portfolio$sd

# Plot previous exercises
plot(sigma_portfolio, mu_portfolio, type="b", pch=16, ylim=c(0, max(mu_portfolio)), xlim=c(0, max(sigma_portfolio)), xlab=expression(sigma[p]), ylab=expression(mu[p]),col=c(rep("green", 18), rep("red", 13)))
points(sigma_portfolio_tangency_bill,mu_portfolio_tangency_bill, col="blue", type="b", pch=16)
text(x=sigma_boeing, y=mu_boeing, labels="Boeing", pos=4)
text(x=sigma_msft, y=mu_msft, labels="MSFT", pos=4)

# Plot portfolio combinations of tangency portfolio and T-bills
text(x=tangency_portfolio$sd, y=tangency_portfolio$er, labels="Tangency", pos=2)

## Part 10: An Efficient Portfolio with 30% Tangency
# Define the portfolio ratio's
tangency_weight <- 0.3
t_bill_weight <- 1 - 0.3

# Define the portfolio parameters
mu_portfolio_efficient <- t_bill_rate * t_bill_weight + tangency_weight * tangency_portfolio$er

sd_portfolio_efficient <- tangency_weight * tangency_portfolio$sd

# Plot previous exercises
plot(sigma_portfolio, mu_portfolio, type="b", pch=16, ylim=c(0, max(mu_portfolio)), xlim=c(0, max(sigma_portfolio)), xlab=expression(sigma[p]), ylab=expression(mu[p]),col=c(rep("green", 18), rep("red", 13)))
text(x=sigma_boeing, y=mu_boeing, labels="Boeing", pos=4)
text(x=sigma_msft, y=mu_msft, labels="MSFT", pos=4)
text(x=tangency_portfolio$sd, y=tangency_portfolio$er, labels="Tangency", pos=2)
points(sigma_portfolio_tangency_bill, mu_portfolio_tangency_bill, type="b", col="blue", pch=16)

# Plot Efficient Portfolio with 30% Tangency
points(x=sd_portfolio_efficient, y=mu_portfolio_efficient, type="b", col="orange", pch=16, cex=2)
text(x=sd_portfolio_efficient, y=mu_portfolio_efficient, labels="Efficient Portfolio with 30% Tangency", pos=4, cex=0.75)

## Part 11: An Efficient Portfolio with the SD of Boeing
tangency_weight = sigma_boeing / tangency_portfolio$sd
t_bill_weight = 1 - tangency_weight

# Calculate the portfolio parameters
mu_portfolio_efficient <- t_bill_rate * t_bill_weight + tangency_weight * tangency_portfolio$er
sd_portfolio_efficient <- tangency_weight * tangency_portfolio$sd    #sigma_boeing

# Plot previous exercises
plot(sigma_portfolio, mu_portfolio,bg="NA", type="b", pch=16, ylim=c(0, max(mu_portfolio)), xlim=c(0, max(sigma_portfolio)), xlab=expression(sigma[p]), ylab=expression(mu[p]),col=c(rep("green", 18), rep("red", 13)))
text(x=sigma_boeing, y=mu_boeing, labels="Boeing", pos=4)
text(x=sigma_msft, y=mu_msft, labels="MSFT", pos=4)
text(x=tangency_portfolio$sd, y=tangency_portfolio$er, labels="Tangency", pos=2)
points(sigma_portfolio_tangency_bill, mu_portfolio_tangency_bill, type="b", col="blue", pch=16)

# Plot Efficient Portfolio with the same risk as Boeing
points(sd_portfolio_efficient, mu_portfolio_efficient, type="p", col="orange", pch=16, cex=2)
text(x=sd_portfolio_efficient, y=mu_portfolio_efficient, labels="Efficient Portfolio with same risk as Boeing", pos=2, cex=0.75)

### QUIZ QUESTIONS
#Q1: Sharpe Slope of Boeing
(mu_boeing - 0.03) / sigma_boeing

# Q2: Sharpe Slope of Microsoft
(mu_msft - 0.03) / sigma_msft

# Q3:  What is the Sharpe slope of the global minimum variance portfolio? 
(global_min_var_portfolio$er - 0.03) / global_min_var_portfolio$sd

# Q4:  What is the Sharpe slope of the tangency portfolio? 
(tangency_portfolio$er - 0.03) / tangency_portfolio$sd

# Q5: What is the Sharpe slope of a portfolio that has 10% in the tangency portfolio and 90% in T-bills? 
# Define the portfolio ratio's
tangency_weight <- 0.1
t_bill_weight <- 1 - tangency_weight

# Define the portfolio parameters
mu_portfolio_efficient <- t_bill_rate * t_bill_weight + tangency_weight * tangency_portfolio$er
sd_portfolio_efficient <- tangency_weight * tangency_portfolio$sd

(mu_portfolio_efficient - 0.03) / sd_portfolio_efficient

# Q6:  What is the Sharpe slope of the efficient portfolio (combination of T-bills and tangency portfolio) that has the same risk (SD) as Microsoft?
tangency_weight = sigma_msft / tangency_portfolio$sd
t_bill_weight = 1 - tangency_weight

mu_portfolio_efficient <- t_bill_rate * t_bill_weight + tangency_weight * tangency_portfolio$er
sd_portfolio_efficient <- tangency_weight * tangency_portfolio$sd 

(mu_portfolio_efficient - 0.03) / sd_portfolio_efficient
sigma_msft
sd_portfolio_efficient

# Q7:  What is the portfolio weight of Microsoft in the global minimum variance portfolio?
global_min_var_portfolio

# Q8:  What is the portfolio weight of Microsoft in the tangency portfolio?
tangency_portfolio

# Q9: What is the expected return of the efficient portfolio (combination of T-bills and tangency portfolio) that has the same risk (SD) as Microsoft? 
mu_portfolio_efficient

# Q10:  For the efficient portfolio (combination of T-bills and tangency portfolio) that has the same risk (SD) as Microsoft, 
# what is the percentage of wealth invested into T-bills? 
t_bill_weight

