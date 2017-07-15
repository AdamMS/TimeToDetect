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
int<lower=1> n_ra; # Number of sets of random effects
int<lower=1> n_ras[n_ra]; # Number of levels for each random effect --- this is an (n_ra)-length vector
int<lower=1> n_rds;
int<lower=1> cn_ras[n_ra]; # Cumulative number of random effects
int<lower=1> cn_rds;
int<lower=1> va_id[cn_ras[n_ra]]; # Variance id for each random effect --- this is a vector of indices.  Its length = total number of abundance random effect levels across all random effects.  Each index specifies to which \psi_j each random effect level belongs.  So, for a 2-psi model, va_id is composed of all 1's and 2's.  There is no vd_id, because there is only one detection random effect.
int ia[n_sites,n_ra];
int id[n_sites]; # Note the different dimension versus 'ia' because there's only one random detection effect.

# Data
real<lower=0> tau[n_ints];
int<lower=0> y[n_sites,n_ints];
}



transformed data {
int ii[n_sites];          # 'ii' is an index of sites relative to the data vector
int yv[n_sites*n_ints];   # The data vector
int obsN[n_sites];        # Total count at each site



# Vectorize data matrix
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
# Fixed effects
real intcpt_a;     # Abundance intercept
real intcpt_d;     # Detection intercept
vector[n_ba] ba;
vector[n_bd] bd;
real<lower=0,upper=1> gamma;

# Random effects
vector<lower=0>[n_ra] sigma_a;
real<lower=0> sigma_d;
vector[cn_ras[n_ra]] ra;
vector[cn_rds] rd;
}



transformed parameters {
vector[n_sites] log_lambda;
vector[n_sites] log_rho;
vector[n_sites] rho;
real<upper=0> log_p[n_sites,n_ints];
vector[n_sites*n_ints] log_mu_ab;
vector[n_sites] lambda;

# Abundance
log_lambda = Xa*ba;
for (s in 1:n_sites) {
  log_lambda[s] = log_lambda[s] + intcpt_a;
  for (i in 1:n_ra) log_lambda[s] = log_lambda[s] + ra[ia[s,i]];
}
lambda = exp(log_lambda);

# Detection
log_rho = Xd*bd;
for (s in 1:n_sites) {
  log_rho[s] = intcpt_d + log_rho[s] + rd[id[s]];
}
rho = exp(log_rho);

# Calculate site-specific interval detection probabilities
# and Poisson means
for (s in 1:n_sites) {
  log_p[s,1] = log(gamma * exponential_cdf(tau[1], rho[s]) + 1 - gamma);
  log_mu_ab[ii[s]+1] = log_lambda[s] + log_p[s,1];
  for (i in 2:n_ints) {
    log_p[s,i] = log(gamma) + log_diff_exp(exponential_lccdf(tau[i-1] | rho[s]), 
                                           exponential_lccdf(tau[i] | rho[s]));
    log_mu_ab[ii[s]+i] = log_lambda[s]+log_p[s,i];
  }
}
}



model {
# Fixed effects
ba ~ normal(0,1); 
bd ~ normal(0,2.25); 
intcpt_a ~ normal(1,1);
intcpt_d ~ normal(-2.7,2.25);
gamma ~ beta(1,1);

# Random effects
for (i in 1:cn_ras[n_ra]) ra[i] ~ normal(0,sigma_a[va_id[i]]);
for (i in 1:cn_rds) rd[i] ~ normal(0,sigma_d);
sigma_a ~ cauchy(0,1);
sigma_d ~ cauchy(0,1);

# Data models
yv ~ poisson_log(log_mu_ab);
}



generated quantities {
int<lower=0> uncounted;              # Total uncounted birds
int<lower=0> abundance;              # Total uncounted birds
real<lower=0> p_global;              # Estimated overall detection probability
int<lower=0> unobserved[n_sites];    # Estimated uncounted birds
int<lower=0> totN[n_sites];          # Estimated total birds
# int<lower=0> yrep[n_sites,n_ints];   # Replicate data based on estimates of expected abundance and detection rates
real dev1;                           # Total deviance

real lpn_BK[n_sites,n_ints];         # log(p(n|\beta's,\ksi's))
# vector [n_ints+1] p_int;             # Holding vector for interval-specific p_det at current site

dev1 = 0;
for (s in 1:n_sites) {
  unobserved[s] = poisson_log_rng(log_lambda[s] + log(gamma) + exponential_lccdf(tau[n_ints] | rho[s]));
  totN[s] = obsN[s] + unobserved[s];
  for (i in 1:n_ints) {
    # yrep[s,i] = poisson_rng(exp(log_mu_ab[ii[s]+i]));
    # p_int[i] = exp(log_p[s,i]);     # Estimated p[s,i]
    lpn_BK[s,i] = poisson_log_lpmf(y[s,i] | log_mu_ab[ii[s]+i]);
    dev1 = dev1 - 2 * lpn_BK[s,i];
  }
}
uncounted = sum(unobserved);
p_global  = sum(obsN); 
p_global  = p_global / (p_global + uncounted);
abundance = sum(totN); 
}
