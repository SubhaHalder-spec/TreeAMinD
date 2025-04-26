#' Asymptotic Min-D test for testing equality means against tree-ordered alternatives in one-way ANOVA
#' @export
#' @param sample_data list
#' @param significance_level numeric
#' @return Critical value numeric
#' @return Test statistic value numeric
#' @return Result Character
#' @details Testing of H_0:mu_0 = mu_1 = ... = mu_k vs H_1:mu_0 <= mu_i (at least one strict inequality), for all i >= 1, where mu_i represents the population means of the i-th treatment. The input consists of two variables: sample_data and significance_level. The output consists of the critical value, the AMin-D test statistic value, and the result, which indicates whether to reject or not reject the null hypothesis.
#' @importFrom MASS mvrnorm
#' @import stats
#' @author Subha Halder
TreeAMinD <- function(sample_data, significance_level){
  set.seed(456)
  sample_data <- lapply(sample_data, function(x) x[!is.na(x)])
  num_samples = 100000
  num_datasets <- length(sample_data)
  n <- sapply(sample_data, length)
  proportions <- n / sum(n)
  var_data <- sapply(1:num_datasets, function(j) var(sample_data[[j]]))
  b_sq <- proportions[1] * var_data[-1] + proportions[-1] * var_data[1]
  s <- num_datasets - 1
  u_t <- matrix(0, s, s)
  l_t <- matrix(0, s, s)
  B <- diag(sqrt(b_sq))

  for (i in 1:(s - 1)) {
    for (j in (i + 1):s) {
      u_t[i, j] <- sqrt(proportions[i + 1] * proportions[j + 1]) * var_data[1]
      l_t[j, i] <- u_t[i, j]
    }
  }
  cov_matrix <- u_t + l_t + B %*% t(B)
  D <- solve(B) %*% cov_matrix %*% solve(B)
  D_star_amin_values <- numeric(num_samples)
  for (k in 1:num_samples) {
    bootstrap_s <- MASS::mvrnorm(n = 1, mu = rep(0, s), Sigma = D)
    D_star_amin_values[k] <- min(bootstrap_s)
  }
  quantile_value <- quantile(sort(D_star_amin_values), probs = 1-significance_level)
  D_Amin <- min(sapply(2:num_datasets, function(j) (mean(sample_data[[j]]) - mean(sample_data[[1]])) /
                         sqrt((var(sample_data[[j]]) / length(sample_data[[j]])) +
                                (var(sample_data[[1]]) / length(sample_data[[1]])))))
  if (D_Amin > quantile_value) {
    result <- "Reject null hypothesis"
  } else {
    result <- "Do not reject null hypothesis"
  }
  return(paste("Critical value:", quantile_value, "; D_AMin Test statistic:", D_Amin, "; Result:", result))
}


