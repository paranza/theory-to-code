# Simulation: Var(beta_hat) ∝ Var(u)  and  Var(beta_hat) ∝ 1/(n * Var(X))
# Model: Y = beta0 + beta1 * X + u
# Under standard OLS assumptions (homoskedastic u), approx:
#   Var(beta1_hat | X) = sigma_u^2 / Sxx
# and for random X with Var(X)=vX, Sxx ≈ n * vX,
# so Var(beta1_hat) ≈ sigma_u^2 / (n * Var(X)).

set.seed(123)

sim_beta1_var <- function(n, varX, varU, R = 5000, beta0 = 0, beta1 = 2) {
  b1hat <- numeric(R)
  
  for (r in 1:R) {
    # Generate regressor with desired variance:
    # If X ~ N(0, varX), then Var(X) = varX
    X <- rnorm(n, mean = 0, sd = sqrt(varX))
    
    # Noise with Var(u) = varU
    u <- rnorm(n, mean = 0, sd = sqrt(varU))
    
    Y <- beta0 + beta1 * X + u
    fit <- lm(Y ~ X)
    b1hat[r] <- coef(fit)[["X"]]
  }
  
  # Return empirical variance and the theoretical approximation using n*Var(X)
  empirical_var <- var(b1hat)
  theory_approx <- varU / (n * varX)
  
  c(empirical_var = empirical_var, theory_approx = theory_approx)
}


# 1) Inverse in n (hold VarX and VarU fixed)
n_grid <- c(25, 50, 100, 200, 400, 800)
varX_fixed <- 4
varU_fixed <- 9

res_n <- t(sapply(n_grid, function(n) sim_beta1_var(n, varX_fixed, varU_fixed)))
res_n <- data.frame(n = n_grid,
                    empirical = res_n[, "empirical_var"],
                    theory = res_n[, "theory_approx"],
                    ratio_emp_to_theory = res_n[, "empirical_var"] / res_n[, "theory_approx"])
print(res_n)

# Optional plot: Var(beta1_hat) vs 1/n should look linear
plot(1 / res_n$n, res_n$empirical,
     xlab = "1/n", ylab = "Empirical Var(beta1_hat)",
     main = "Var(beta1_hat) falls ~ linearly with 1/n")
abline(lm(empirical ~ I(1/n), data = res_n), lwd = 2)

# 2) Inverse in Var(X) (hold n and VarU fixed)
varX_grid <- c(0.5, 1, 2, 4, 8, 16)
n_fixed <- 200
varU_fixed <- 9

res_x <- t(sapply(varX_grid, function(vx) sim_beta1_var(n_fixed, vx, varU_fixed)))
res_x <- data.frame(varX = varX_grid,
                    empirical = res_x[, "empirical_var"],
                    theory = res_x[, "theory_approx"],
                    ratio_emp_to_theory = res_x[, "empirical_var"] / res_x[, "theory_approx"])
print(res_x)

#Var(beta1_hat) vs 1/Var(X) should look linear
plot(1 / res_x$varX, res_x$empirical,
     xlab = "1/Var(X)", ylab = "Empirical Var(beta1_hat)",
     main = "Var(beta1_hat) falls as Var(X) grows")
abline(lm(empirical ~ I(1/varX), data = res_x), lwd = 2)


# 3) Direct in Var(U) (hold n and VarX fixed)

varU_grid <- c(0.25, 1, 4, 9, 16, 25)
n_fixed <- 200
varX_fixed <- 4

res_u <- t(sapply(varU_grid, function(vu) sim_beta1_var(n_fixed, varX_fixed, vu)))
res_u <- data.frame(varU = varU_grid,
                    empirical = res_u[, "empirical_var"],
                    theory = res_u[, "theory_approx"],
                    ratio_emp_to_theory = res_u[, "empirical_var"] / res_u[, "theory_approx"])
print(res_u)

# Var(beta1_hat) vs Var(U) should look linear
plot(res_u$varU, res_u$empirical,
     xlab = "Var(U)", ylab = "Empirical Var(beta1_hat)",
     main = "Var(beta1_hat) rises ~ linearly with Var(U)")
abline(lm(empirical ~ varU, data = res_u), lwd = 2)

# ratio_emp_to_theory should be near a constant ~ 1
# (Monte Carlo noise shrinks as R increases).
