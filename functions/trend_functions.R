#######################################################################################################################################

# Trend functions

# Linear trend
linear_trend <- function(j,lambda,sample_size){lambda*(j-1)/(sample_size-1)}

# Linear trend only second period
# sample size: vector of dim 2, with sample size per periods
linear_trend2 <- function(j,lambda,sample_size){ifelse(j<=sample_size[1],0,lambda*(j-1)/(sum(sample_size)-1))}

# Stepwise trend
sw_trend <- function(cj,lambda){lambda*cj}

#######################################################################################################################################