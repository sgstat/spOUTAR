# Total Economy: Exports

# Load necessary packages
library(tidyverse)
library(lubridate)
library(zoo)  # for as.yearqtr

# Step 1: Import the data
df <- read_csv("exports.csv")

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

df1 = time_series_df_clean[-c(5,6,10,15,19,34,36,44), ]
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

export1 = df2[1:(indx1-1), ] # T x p format
export2 = df2[(indx2-1):(indx3-1), ]

results = rep(0, 2)
Y = t(rbind(export1, export2))
cat("\n", dim(Y))
# tm1 = tic()
# model = OUTAR(Y = 100*Y[, 1:78], K = 2, n1 = 38, n2 = 40, n.iter = 10000, verbose = TRUE, adaptA = FALSE, adaptAR = FALSE)
# tm2 = toc()

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
            Y = 10000 * Y[, 1:78],
            K = k,
            n1 = 38,
            n2 = 40,
            n.iter = 10000,
            verbose = TRUE,
            adaptA = FALSE,
            adaptAR = FALSE
          )
          # Save the model
          dir.create(paste0("OECD results new/export_results"), showWarnings = FALSE, recursive = TRUE)
          saveRDS(model, paste0("OECD results new/export_results/model", k, ".rds"))
        }

# Stop the cluster
stopCluster(cl)


omegas1 = model$Omega1.samples[5001:10000,,] # before
omegas2 = model$Omega2.samples[5001:10000,,] # after

omegas_diff = omegas2 - omegas1

Omega1.pred = apply(model$Omega1.samples[5001:10000,,], c(2, 3), mean) # before
Omega2.pred = apply(model$Omega2.samples[5001:10000,,], c(2, 3), mean) # after

Omega_diff.pred = Omega2.pred - Omega1.pred  # posterior mean

lower = apply(omegas_diff, c(2,3), quantile, probs=0.05)
upper = apply(omegas_diff, c(2,3), quantile, probs=0.975)

adjmat = apply((lower > 0 | upper < 0), c(1,2), as.integer)  # flag 1 if 95% credible interval doesn't contain 0.
W.pred = graph_from_adjacency_matrix(adjmat, mode = "undirected", diag=FALSE)
E(W.pred)$width = 1.5
E(W.pred)$length = 4
V(W.pred)$name = countries$`Reference area`
# layout <- layout_with_fr(W.pred)     # Fruchterman-Reingold (force-directed)
# layout <- layout_with_kk(W.pred)  # Kamada-Kawai (spring layout)
# layout <- layout_in_circle(W.pred)  # Circular layout
# layout <- layout_with_drl(W.pred)   # Large graphs (DRL layout)
# plot(W.pred,
#      layout = layout,
#      vertex.size = 15,          # Increase node size
#      vertex.label.cex = 0.7,    # Decrease font size of labels
#      vertex.label.color = "black",  # Optional: set label color
#      vertex.color = "skyblue")

# Determine the sign of significant change
edge_signs = matrix(0, nrow = nrow(adjmat), ncol = ncol(adjmat))
edge_signs[lower > 0] = 1   # Positive change
edge_signs[upper < 0] = -1  # Negative change

# Extract edge list to match color with edge
edge_list = as_edgelist(W.pred, names = FALSE)

# Initialize color vector
edge_colors = character(ecount(W.pred))

# Assign color based on sign of change
for (i in seq_len(nrow(edge_list))) {
  v1 <- edge_list[i, 1]
  v2 <- edge_list[i, 2]
  sign_val <- edge_signs[v1, v2]
  
  edge_colors[i] <- if (sign_val == 1) {
    "blue"     # positive effect
  } else if (sign_val == -1) {
    "red"      # negative effect
  } else {
    "black"    # fallback (should not happen with filtered adjmat)
  }
}

# Assign to graph
E(W.pred)$color <- edge_colors


# Total number of nodes
n <- vcount(W.pred)

# Example: split nodes into inner (1/3) and outer (2/3) circles
n_inner <- ceiling(n / 3)
n_outer <- n - n_inner

# Radii for inner and outer circles
r_inner <- 1
r_outer <- 2

# Create empty layout matrix
layout <- matrix(NA, nrow = n, ncol = 2)

# Inner circle coordinates
angles_inner <- seq(0, 2 * pi, length.out = n_inner + 1)[- (n_inner + 1)]
layout[1:n_inner, ] <- cbind(
  r_inner * cos(angles_inner),
  r_inner * sin(angles_inner)
)

# Outer circle coordinates
angles_outer <- seq(0, 2 * pi, length.out = n_outer + 1)[- (n_outer + 1)]
layout[(n_inner + 1):n, ] <- cbind(
  r_outer * cos(angles_outer),
  r_outer * sin(angles_outer)
)

# Plot with custom layout
plot(
  W.pred,
  layout = layout,
  vertex.size = 17,
  vertex.label.cex = 0.7,
  vertex.label.color = "black",  # Optional: set label color
  vertex.color = "oldlace", # ivory, ghostwhite, whiteSmoke, oldlace
  vertex.frame.width = 2
)

# Add legend
legend(
  "topright",
  legend = c("Increase after Recession", "Decrease after Recession"),
  col = c("blue", "red"),
  lty = 1,
  lwd = 2,
  title = "Shifts in Export Expenditure Conditional \nDependencies After the 2009 Recession",
  box.lwd = 0,
  inset = 0.001
)
