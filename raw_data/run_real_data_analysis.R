get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  path <- sub(file_arg, "", args[grep(file_arg, args)])
  if (length(path) > 0) return(dirname(path[1]))
  return(getwd())
}

data_dir <- file.path("data")
raw_data_dir <- file.path("raw_data")

evalue_dir <- file.path(data_dir, "evalues")

dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(evalue_dir, showWarnings = FALSE, recursive = TRUE)

# Install missing packages if necessary
list.of.packages <- c("dplyr", "stlplus", "ks", "tseries", "urca", "tidyquant")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

# Load necessary libraries
library(stats)
library(dplyr)
library(stlplus)
library(ks)
library(tseries)
library(urca)
library(tidyquant)

# New KDE function
adaptive_bivariate_kernel_density <- function(x, y, eval_points_x, eval_points_y, alpha = 0.5, h_x_adaptive = 15, h_y_adaptive = 15, adaptive = TRUE) {
  n <- length(x)
  
  if (adaptive) {
    h_global <- Hpi(cbind(x, y))
    h_x <- sqrt(h_global[1, 1])
    h_y <- sqrt(h_global[2, 2])
    
    h_x_adaptive <- h_x
    h_y_adaptive <- h_y
  }
  
  weights <- dnorm((eval_points_x - x) / h_x) * dnorm((eval_points_y - y) / h_y) / (h_x * h_y)
  weights_t <- dnorm((eval_points_x - x) / h_x) * h_x
  density_matrix <- sum(weights) / sum(weights_t)
  return(density_matrix)
}

# Real Data 5.2 Processing
message("--- Processing Real Data 5.2 (NYC Taxi) ---")

# Define the path for the nyc_taxi.csv file
nyc_taxi_file <- file.path(data_dir, "nyc_taxi.csv")
nyc_taxi_url <- "https://raw.githubusercontent.com/numenta/NAB/refs/heads/master/data/realKnownCause/nyc_taxi.csv"

# Download the file if it doesn't exist
if (!file.exists(nyc_taxi_file)) {
  message("Downloading nyc_taxi.csv...")
  download.file(nyc_taxi_url, destfile = nyc_taxi_file, method = "auto")
  message("Download complete.")
} else {
  message("nyc_taxi.csv already exists. Skipping download.")
}

data_test <- read.csv(nyc_taxi_file)
data <- data_test

data_test <- stlplus(
  x = data$value,
  t = data$timestamp,
  n.p = 48*14,
  s.window = "periodic",
  t.window = 50,
  t.degree = 1,
  inner = 20,
  outer = 5
)

X_calib<-  data_test[[1]][[4]][1:2000]
X_test <- data_test[[1]][[4]][-(1:2000)]

Z <- c(X_calib,X_test)
Z0 <- X_calib

n_cal <- length(X_calib)
n_test <- length(X_test)
init <- n_cal
m <- length(Z)
h <- 5
N <- 20#floor(5*h)
f_t <- rep(0,n_test)

# # Initialize progress bar
# pb <- txtProgressBar(min = init + 1, max = m, style = 3)

# for (i in (init+1):m){
#   tmp_t <- max(1,i-N):(i-1)
#   tmp_x <- Z[max(1,i-N):(i-1)]
#   f_t[i-init]<- adaptive_bivariate_kernel_density(tmp_t, tmp_x, i, Z[i], alpha = 0.5)
#   # Update progress bar
#   setTxtProgressBar(pb, i)
# }

# # Close progress bar
# close(pb)

# f_0 <- rep(0,n_test)


# f_0_func <- density(Z0)


# for (i in (init+1):m){
#   f_0[i-init]<- approx(f_0_func$x, f_0_func$y, xout = Z[i])$y
# }

# e_LR <- f_t/f_0

# input_e_LR <- data.frame(evalue = e_LR)
# # drop na rows
# input_e_LR <- na.omit(input_e_LR)

# # Output e-values to CSV
# write.csv(input_e_LR, file.path(evalue_dir, "nyc_taxi_evalues.csv"), row.names = FALSE)
# message("NYC Taxi e-values saved to data/nyc_taxi_evalues.csv")

# Real Data 5.3 Processing
message("--- Processing Real Data 5.3 (NASDAQ) ---")


# Cache data - save locally if not exists, load if exists

nasdaq_csv_path <- file.path(raw_data_dir, "nasdaq_data.csv")
if (!file.exists(nasdaq_csv_path)) {
  message("Downloading NASDAQ data...")
  fred_url <- "https://fred.stlouisfed.org/graph/fredgraph.csv?id=NASDAQCOM"
  download.file(fred_url, destfile = nasdaq_csv_path, mode = "wb")
}
nasdaq_data <- read.csv(nasdaq_csv_path)
colnames(nasdaq_data) <- c("date", "price")
nasdaq_data$date <- as.Date(nasdaq_data$date)
# select dates between "1970-01-01", to = "2008-12-31"
nasdaq_data <- nasdaq_data %>%
  filter(date >= as.Date("1970-01-01") & date <= as.Date("2008-12-31"))


data_daily<- nasdaq_data[!is.na(nasdaq_data$price),]

data2 <- data_daily$price

data_weekly <- nasdaq_data %>%
  tq_transmute(select = price,
               mutate_fun = to.period,
               period = "weeks",
               col_rename = "weekly_price")
colnames(data_weekly)[2] <- "price"

ncal <- round(1/3*length(data2))
X_calib2<-  data2[1:ncal]
X_test2 <- data2[-(1:ncal)]

Z2 <- data2
Z02 <- X_calib2

ntest <- length(X_test2)
init2 <- ncal

m2 <- length(Z2)
h <- 10
N <- floor(10*h)
f_t2 <- rep(0,ntest)

alpha = 0.01
n_min <- 100
p_values2 <- numeric(ntest)

# Initialize progress bar
pb2 <- txtProgressBar(min = init2 + 1, max = m2, style = 3)

for(i in (init2+1):m2){
  adf_result <- ur.df(Z2[max(1, (i - n_min)):i],type = "drift", lags = 1)
  p_values2[i-init2] <- summary(adf_result)@testreg$coefficients["z.lag.1", "Pr(>|t|)"]
  # Update progress bar
  setTxtProgressBar(pb2, i)
}

# Close progress bar
close(pb2)

e_LR2 <- numeric(ntest)
for(i in (init2+1):m2){
  Y_tmp2 <- Z2[(i-n_min+2):i]
  Y_tmp <- Z2[(i-n_min+1):(i-1)]
  X_tmp <- t(rbind(rep(1,(n_min-1)),(seq(1,(n_min-1))-n_min/2),Y_tmp))
  Se2 <- 1/(n_min-4)*t(Y_tmp)%*%(diag(n_min-1)-X_tmp%*%solve(t(X_tmp)%*%X_tmp)%*%t(X_tmp))%*%Y_tmp
  sigma02 <- 1/(n_min-1)*sum(Y_tmp2-Y_tmp)^2
  phi2 <- 1/(3*Se2)*((n_min-1)*sigma02)-(n_min-4)/3
  phi3 <- 1/(2*Se2)*((n_min-1)*(sigma02-mean(Y_tmp2-Y_tmp)^2))-(n_min-4)/2
  e_LR2[i-init2] <- (pmax(0,1+3/(n_min-4)*phi2))^(1/2)
}


input_e_LR2 <- data.frame(evalue = e_LR2)
input_e_LR2 <- na.omit(input_e_LR2)

# Output e-values to CSV
write.csv(input_e_LR2, file.path(evalue_dir, "nasdaq_evalues.csv"), row.names = FALSE)
message("NASDAQ e-values saved to data/nasdaq_evalues.csv")
