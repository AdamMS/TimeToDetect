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
int<lower=1> cn_ras[n_ra]; # Cumulative number of abundance random effect levels
int<lower=1> cn_rds;       # Same as n_rds, because there is only one detection random effect
int<lower=1> va_id[cn_ras[n_ra]];   # Variance ID for each abundance random effect level -- a vector of indices.  For each random effect level, it tells me what random effect it belongs to (e.g., 2011 belongs to Year random effects).  So, for a 2-effect model, va_id is composed of all 1's and 2's
int<lower=1> ia[n_sites,n_ra];      # IDs what abundance random effect levels go with each site survey
int<lower=1> id[n_sites];           # IDs what detection random effect level goes with each site survey  

# Data
real<lower=0> tau[n_ints];       # Endpoints of observation intervals
int<lower=0> y[n_sites,n_ints];  # Counts by site and observation interval
}



transformed data {
int<lower=0> ii[n_sites];         # An index identifying where in the vectorized data each site begins (minus 1)
int<lower=0> yv[n_sites*n_ints];  # Vectorized data
int<lower=0> obsN[n_sites];       # Total count at each site across all intervals

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
real<upper=15> intcpt_a;     # Abundance intercept
real intcpt_d;     # Detection intercept
vector[n_ba] ba;   # Abundance fixed effect estimates
vector[n_bd] bd;   # Detection fixed effect estimates
real<lower=0> alpha;  # Universal alpha from the standard Gamma distribution

# Random effects
vector<lower=0>[n_ra] sigma_a;       # Abundance random effect variance
real<lower=0> sigma_d;               # Detection random effect variance
vector[cn_ras[n_ra]] ra;             # Abundance random effect estimates
vector[cn_rds] rd; # Detection random effect variance
}



transformed parameters {
vector[n_sites] log_lambda;             # Site-level log(Expected abundance) --- 'upper' is a constraint on early iterations
vector[n_sites] log_beta;               # Site-level log(Expected detection rate) -- this is Beta in the standard Gamma distribution
vector<lower=0>[n_sites] beta;          # Site-level detection rate
vector<upper=0>[n_ints] log_p[n_sites]; # Interval-specific log(detection probability)
vector[n_sites*n_ints] log_mu_ab;       # Interval-specific log(Expected count = lambda * Pr(detection in interval))
vector<lower=0>[n_sites] lambda;        # Site-level Expected abundance
vector[n_ints] cutoffs;                 # Temporary vector to hold cdf at cutoffs (saves computation of gamma_cdf)

# Calculating site-level expected abundance
log_lambda = Xa*ba;
for (s in 1:n_sites) {
  log_lambda[s] = log_lambda[s] + intcpt_a;
  for (i in 1:n_ra) log_lambda[s] = log_lambda[s] + ra[ia[s,i]];
}
lambda = exp(log_lambda);

# In this version of the Gamma model, we model log(beta), not log(mu)
# So instead of the mean param: log(beta) = Xd*bd + ObsEff + log(alpha)
# We now use the simpler      : log(beta) = Xd*bd + ObsEff
# We can move back and forth between the two.  In the mean parameterization, the interpretation of intcpt_d was -log(mu).  Now, the interpretation is intcpt_d = log(beta) = log(alpha)-log(mu)
log_beta = Xd*bd;
for (s in 1:n_sites) {
  log_beta[s] = intcpt_d + log_beta[s] + rd[id[s]];
}
beta = exp(log_beta);

# Calculating interval-specific log-scale detection probabilities (log_p) and log-scale Poisson means (log_mu)
for (s in 1:n_sites) {
  for (i in 1:n_ints) {
    cutoffs[i] = gamma_lcdf(tau[i] | alpha,beta[s]);
  }
  log_p[s,1] = cutoffs[1];
  log_mu_ab[ii[s]+1] = log_lambda[s] + log_p[s,1];
  for (i in 2:n_ints) {
    log_p[s,i] = log_diff_exp(cutoffs[i], cutoffs[i-1]);
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
alpha ~ cauchy(0,1);

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
  unobserved[s] = poisson_log_rng(log_lambda[s] + gamma_lccdf(tau[n_ints] | alpha, beta[s]));
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
