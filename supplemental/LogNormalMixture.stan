data {
# Data dimensions
int<lower=1> n_sites;     # Number of surveys
int<lower=1> n_ints;      # Number of censored detection intervals

# Covariates
int<lower=0> n_ba;        # Number of abundance fixed effects
int<lower=0> n_bd;        # Number of detection fixed effects
matrix[n_sites,n_ba] Xa;  # Explanatory variables for abundance
matrix[n_sites,n_bd] Xd;  # Explanatory variables for detection

# Random effect ids
int<lower=1> n_ra;         # Number of abundance random effects (e.g., Site or Year)
int<lower=1> n_ras[n_ra];  # Number of levels for each random effect
int<lower=1> n_rds;        # Number of levels for the lone detection random effect (Observer)
int<lower=1> cn_ras[n_ra]; # Cumulative number of abundance random effect levels... used for indexing
int<lower=1> cn_rds;       # Same as n_rds, because there is only one detection random effect
int<lower=1> va_id[cn_ras[n_ra]]; # Effect-category ID for each effect level --- this is a vector of indices.
     # length(va_id) = total number of abundance random effect levels across all random effects = sum(n_ras).
     # There is no vd_id, because there is only one detection random effect.
int ia[n_sites,n_ra];      # IDs what abundance random effect levels go with each survey
int id[n_sites];           # IDs what detection random effect level goes with each survey  

# Data
real<lower=0> tau[n_ints];       # Endpoints of observation intervals
int<lower=0> y[n_sites,n_ints];  # Counts by survey and observation interval
}



transformed data {
int<lower=0> ii[n_sites];         # An index identifying where in the vectorized data each survey begins (minus 1)
int<lower=0> yv[n_sites*n_ints];  # Vectorized data
int<lower=0> obsN[n_sites];       # Total count at each survey across all intervals

# Calculates each of the above
for (s in 1:n_sites) {
  ii[s] = (s-1)*n_ints;
  obsN[s] = 0;
  for (i in 1:n_ints) {
    yv[ii[s]+i] = y[s,i];
    obsN[s] = obsN[s] + y[s,i];
  }
}
}



parameters {
real intcpt_a;     # Abundance intercept
real intcpt_d;     # Detection intercept
vector[n_ba] ba;   # Abundance fixed effect estimates
vector[n_bd] bd;   # Detection fixed effect estimates
real<lower=0> sigma_det;        # Shape parameter for the Lognormal distribution
real<lower=0,upper=1> gamma;    # Heterogeneity mixing parameter

# Random effects
vector<lower=0>[n_ra] sigma_a;  # Abundance random effect variances
real<lower=0> sigma_d;          # Detection random effect variance
vector[cn_ras[n_ra]] ra;        # Abundance random effect estimates
vector[cn_rds] rd;              # Detection random effect estimates
}



transformed parameters {
vector[n_sites] log_lambda;             # Survey-level log(Expected abundance)
vector[n_sites] mu_det;                 # Survey-level mu from the LogNormal distribution (= -log(phi))
vector<upper=0>[n_ints] log_p[n_sites]; # Interval-specific log(detection probability)
vector[n_sites*n_ints] log_mu_ab;       # Interval-specific log(Expected count = lambda * Pr(detection in interval))
vector<lower=0>[n_sites] lambda;        # Survey-level Expected abundance

# Calculating Survey-level expected abundance
log_lambda = Xa*ba;
for (s in 1:n_sites) {
  log_lambda[s] = log_lambda[s] + intcpt_a;
  for (i in 1:n_ra) log_lambda[s] = log_lambda[s] + ra[ia[s,i]];
}
lambda = exp(log_lambda);

# Calculating Survey-level rate of detection
mu_det = -Xd*bd;  # Note: the negative sign is in keeping with the other distributions, since E[T] \propto 1/rate
for (s in 1:n_sites) {
  mu_det[s] = mu_det[s] - intcpt_d - rd[id[s]]; # WE COULD ADD 0.5*(sigma_det)^2 TO MODEL THE TRUE MEAN
}

# Calculating interval-specific log-scale detection probabilities (log_p) and log-scale Poisson means (log_mu_ab)
for (s in 1:n_sites) {
  log_p[s,1] = log(gamma * lognormal_cdf(tau[1],mu_det[s],sigma_det) + 1 - gamma);
  log_mu_ab[ii[s]+1] = log_lambda[s] + log_p[s,1];
  for (i in 2:n_ints) {
    log_p[s,i] = log(gamma) + log_diff_exp(lognormal_lccdf(tau[i-1]|mu_det[s],sigma_det),
                                            lognormal_lccdf(tau[i]|mu_det[s],sigma_det));
    log_mu_ab[ii[s]+i] = log_lambda[s]+log_p[s,i];
  }
}
}



model {
# Fixed effect priors
ba ~ normal(0,1); 
bd ~ normal(0,1.2); 
intcpt_a ~ normal(1,1);
intcpt_d ~ normal(-2.3,1.2);
sigma_det ~ cauchy(0,1);
gamma ~ beta(1,1);

# Random effect priors
for (i in 1:cn_ras[n_ra]) ra[i] ~ normal(0,sigma_a[va_id[i]]);
for (i in 1:cn_rds) rd[i] ~ normal(0,sigma_d);
sigma_a ~ cauchy(0,1);
sigma_d ~ cauchy(0,1);

# Data model
yv ~ poisson(log_mu_ab);
}


generated quantities {
int<lower=0> unobserved[n_sites];    # Estimated uncounted birds
int<lower=0> totN[n_sites];          # Estimated total birds
int<lower=0> yrep[n_sites,n_ints];   # Replicate data based on estimates of expected abundance and detection rates
real dev1;

real lpn_BK[n_sites,n_ints];         # log(p(n|\beta's,\ksi's))
vector [n_ints+1] p_int;             # Holding vector for interval-specific p_det at current site

dev1 = 0;
for (s in 1:n_sites) {
  unobserved[s] = poisson_rng(exp(log_lambda[s] + log(fmax(1-exp(log_sum_exp(log_p[s])),1e-8)))); # A computational precision issue sometimes leads to p_det very close to 1
  totN[s] = obsN[s] + unobserved[s];
  for (i in 1:n_ints) {
    yrep[s,i] = poisson_rng(exp(log_mu_ab[ii[s]+i]));
    p_int[i] = exp(log_p[s,i]);     # Estimated p[s,i]
    lpn_BK[s,i] = poisson_lpmf(y[s,i] | exp(log_mu_ab[ii[s]+i]));
    dev1 = dev1 - 2 * lpn_BK[s,i];
  }
}
}
