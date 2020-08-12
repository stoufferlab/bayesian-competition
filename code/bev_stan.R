
// generated with brms 2.10.0
functions {
}
data {
  int<lower=1> N;  // number of observations
  int Y[N];  // response variable
  int<lower=1> K_lambda;  // number of population-level effects
  matrix[N, K_lambda] X_lambda;  // population-level design matrix
  int<lower=1> K_alphaii;  // number of population-level effects
  matrix[N, K_alphaii] X_alphaii;  // population-level design matrix
  int<lower=1> K_alphaij;  // number of population-level effects
  matrix[N, K_alphaij] X_alphaij;  // population-level design matrix
  // covariate vectors
  vector[N] C_1;
  vector[N] C_2;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  vector[K_lambda] b_lambda;  // population-level effects
  vector[K_alphaii] b_alphaii;  // population-level effects
  vector[K_alphaij] b_alphaij;  // population-level effects
}
transformed parameters {
}
model {
  // initialize linear predictor term
  vector[N] nlp_lambda = X_lambda * b_lambda;
  // initialize linear predictor term
  vector[N] nlp_alphaii = X_alphaii * b_alphaii;
  // initialize linear predictor term
  vector[N] nlp_alphaij = X_alphaij * b_alphaij;
  // initialize non-linear predictor term
  vector[N] mu;
  for (n in 1:N) {
    // compute non-linear predictor values
    mu[n] = (nlp_lambda[n] / (1 + (nlp_alphaii[n] * C_1[n]) + (nlp_alphaij[n] * C_2[n])));
  }
  // priors including all constants
  target += normal_lpdf(b_lambda | 0, 10);
  target += normal_lpdf(b_alphaii | 0, 10);
  target += normal_lpdf(b_alphaij | 0, 10);
  // likelihood including all constants
  if (!prior_only) {
    target += poisson_lpmf(Y | mu);
  }
}
generated quantities {
}