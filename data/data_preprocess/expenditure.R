# Total Economy: General Government, Final Consumption Expenditures

# Load necessary packages
library(tidyverse)
library(lubridate)
library(zoo)  # for as.yearqtr

# Step 1: Import the data
# df <- read_csv("expenditure.csv")
df <- read_csv("expenditures.csv")  # Total Economy: Household and non-profit, Final Consumption Expenditures

# Step 2: Pivot to p x T format (Reference area x TIME_PERIOD)
time_series_df <- df %>%
  dplyr::select(`Reference area`, TIME_PERIOD, OBS_VALUE) %>%
  pivot_wider(
    names_from = TIME_PERIOD,
    values_from = OBS_VALUE
  )

# Step 3: Remove rows with any NA values
time_series_df_clean <- time_series_df %>%
  drop_na()

# Step 4: Sort time columns chronologically
# Extract column names excluding the first ("Reference area")
time_cols <- names(time_series_df_clean)[-1]

# Sort using proper quarterly order
sorted_time_cols <- time_cols[order(as.yearqtr(time_cols, format = "%Y-Q%q"))]

# Reassemble the data with time columns in correct order
time_series_df_clean <- time_series_df_clean %>%
  dplyr::select(`Reference area`, all_of(sorted_time_cols))

# View cleaned and sorted multivariate time series
# View(time_series_df_clean)

df1 = time_series_df_clean[-c(7,11,17,18,23,31,41,44), ]
# View(df1)
df1 <- df1 %>%
  drop_na()

countries = c(df1[, 1])
df1 = as.data.frame(df1[, -1])
df1 = t(df1)

indx1 = which(rownames(df1) == "2007-Q3")
indx2 = which(rownames(df1) == "2009-Q4")
indx3 = which(rownames(df1) == "2019-Q4")

df2 = apply(log(df1), 2, diff)
# df2 = df2[1:86, ]
# View(df2)
# indx = which(rownames(df2) == "2009-Q2")
# rownames(df2) = NULL
exp1 = df2[1:(indx1-1), ] # T x p format
exp2 = df2[(indx2-1):(indx3-1), ]


results = rep(0, 2)
Y = t(rbind(exp1, exp2))
cat("\n", dim(Y))
tm1 = tic()
model = OUTAR(Y = Y[, 1:84], K = 2, n1 = 44, n2 = 40, n.iter = 10000, verbose = TRUE, adaptA = FALSE, adaptAR = FALSE)
tm2 = toc()

library(foreach)
library(doParallel)

# Set up parallel backend
num_cores <- 8  # Use all but one core
cl <- makePSOCKcluster(num_cores)
registerDoParallel(cl)

# ****************** Parallelized loop *******************
foreach(k = 2:8,
        .packages = c("BDgraph", "MASS", "pracma", "INLA", "Matrix", "astsa", "igraph", "combinat", "parallel"),
        .export = c("OUTAR", "Y")) %dopar% {
          # Assuming OUTAR is from a package; load it in .packages if needed
          model <- OUTAR(
            Y = 1000 * Y[, 1:78],
            K = k,
            n1 = 38,
            n2 = 40,
            n.iter = 10000,
            verbose = TRUE,
            adaptA = FALSE,
            adaptAR = FALSE
          )
          # Save the model
          saveRDS(model, paste0("exp1_results/model", k, ".rds"))
        }

# Stop the cluster
stopCluster(cl)
# ***********************************************************
