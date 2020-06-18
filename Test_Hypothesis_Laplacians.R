# Test the equality of the Fréchet mean and variance of graph laplacians using the test from 
# "Fréchet Analysis of Variance for Random Objects" by Paromita Dubey and Hans-Georg Müller. 
# The hypothesis is written as
#     H_0: \mu_F1 = \mu_F2 and V_F1 = V_F2 vs H_a: ~ H_0 
# and the test for H_0 is based on the test statistic T_n from (11) in Dubey et al.
# Code Author: Shane Lubold (sl223@uw.edu)
# Last updated: June 17, 2020

# ----------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------

TestH0 = function(n_1, n_2, graph_size, Laplacian_Array_1, Laplacian_Array_2, alpha){

# Compute Fréchet means for groups 1 and 2.
FrechetMean_G1 = apply(simplify2array(Laplacian_Array_1), 1:2, mean) # this is the sample Fréchet mean for group 1
FrechetMean_G2 = apply(simplify2array(Laplacian_Array_2), 1:2, mean) # this is the sample Fréchet mean for group 2

# Compute the distances to Frechet mean for group 1.
dist_G1 = rep(0, n_1)
for(i in 1:n_1){
  dist_G1[i] = norm(Laplacian_Array_1[,,i] - FrechetMean_G1, type = 'F')
}
V_F1 = mean(dist_G1^2)

# Compute the distances to Frechet mean for group 2.
dist_G2 = rep(0, n_2)
for(i in 1:n_2){
  dist_G2[i] =  norm(Laplacian_Array_2[,,i] - FrechetMean_G2, type = 'F')
}
V_F2 = mean(dist_G2^2)


# Compute the sample pooled Frechet variance
n = n_1 + n_2  # pooled sample size
pooled_data = abind::abind(Laplacian_Array_1, Laplacian_Array_2)
FrechetMean_Pooled = apply(simplify2array(pooled_data), 1:2, mean) # compute pooled Fréchet mean.
dist_Pooled = rep(0, n)
for(i in 1:n){
  dist_Pooled[i] = norm(pooled_data[,,i] - FrechetMean_Pooled, type = 'F')
}

V_F = mean(dist_Pooled^2)

# Now compute the test statistic
lambda_1n = n_1/n
lambda_2n = n_2/n
F_n = V_F -  (lambda_1n * V_F1 + lambda_2n * V_F2)
hat_Sigma_1 = mean(dist_G1^4) - (mean(dist_G1^2))^2
hat_Sigma_2 = mean(dist_G2^4) - (mean(dist_G2^2))^2
U_n = lambda_1n * lambda_2n * (V_F1 - V_F2)^2 / ( hat_Sigma_1 * hat_Sigma_2)

T_n = n * U_n / (lambda_1n/hat_Sigma_1 + lambda_2n/hat_Sigma_2) + n * F_n^2 / (lambda_1n^2 * hat_Sigma_1 +  lambda_2n^2 * hat_Sigma_2)
return(as.numeric(T_n > cutoff))
}



