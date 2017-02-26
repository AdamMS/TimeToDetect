data {
# Data dimensions
int<lower=1> n_sites;
int<lower=1> n_ints;

# Data
real<lower=0> tau[n_ints];       # Endpoints of observation intervals
int<lower=0> y[n_sites,n_ints];  # Counts by site and observation interval
}



transformed data {
int<lower=0> ii[n_sites];         # An index identifying where in the vectorized data each site begins (minus 1)
int<lower=0> yv[n_sites*n_ints];  # Vectorized data
int<lower=0> obsN[n_sites];       # Total count at each site across all intervals

# Vectorize the data and tally counts by survey
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
real intcpt_a;                       # Abundance intercept
real intcpt_d;                       # Detection intercept
real<lower=0> sigma_det;             # Variance of the lognormal TTD
}



transformed parameters {
real mu_det;                   # Site-level mu from the LogNormal distribution
real<upper=20> log_lambda;     # Site-level log(Expected abundance) --- 'upper' is a constraint on early iterations
vector<upper=0>[n_ints] log_p; # Interval-specific log(detection probability)
real log_mu_ab[n_ints];        # Interval-specific log(Expected count = lambda * Pr(detection in interval))
vector[n_sites*n_ints] long_log_mu; # Site- and interval-specific log_mu_ab (i.e., replicated over all sites)
real<lower=0> lambda;          # Site-level Expected abundance

# Calculating site-level expected abundance
log_lambda = intcpt_a;
lambda     = exp(log_lambda);

# Calculating site-level rate of detection
# mu_det = -Xd*bd;  # Note: the negative sign is in keeping with the other distributions, since E[T] \propto 1/rate, more or less
mu_det = - intcpt_d;

# Calculating interval-specific log-scale detection probabilities (log_p) and log-scale Poisson means (log_mu_ab)
log_p[1] = lognormal_lcdf(tau[1] | mu_det,sigma_det);
log_mu_ab[1] = log_lambda + log_p[1];
for (i in 2:n_ints) {
  log_p[i] = log_diff_exp(lognormal_lccdf(tau[i-1] | mu_det,sigma_det),
                          lognormal_lccdf(tau[i] | mu_det,sigma_det));
  log_mu_ab[i] = log_lambda + log_p[i];
}
# Replicate log_mu_ab to match yv in dimension
for (s in 1:n_sites) for (i in 1:n_ints) long_log_mu[ii[s]+i] = log_mu_ab[i];
}



model {
intcpt_a ~ normal(1,1);
intcpt_d ~ normal(-2.3,1.2);
sigma_det ~ cauchy(0,1);

# Data models
yv ~ poisson_log(long_log_mu);
}


generated quantities {
int<lower=0> unobserved[n_sites];   # Estimated uncounted birds
int<lower=0> uncounted;             # Total uncounted birds
real dev1;
real lpn_BK[n_sites,n_ints];         # log(p(n|\beta's,\ksi's))
int<lower=0> totN[n_sites];          # Estimated total birds by survey
real<lower=0> p_global;              # Estimated overall detection probability

dev1 = 0;
for (s in 1:n_sites) {
  unobserved[s] = poisson_rng(exp(log_lambda + log(fmax(1-exp(log_sum_exp(log_p)),1e-8)))); # A computational precision issue sometimes leads to p_det very close to 1
  totN[s] = obsN[s] + unobserved[s];
  for (i in 1:n_ints) {
    lpn_BK[s,i] = poisson_lpmf(y[s,i] | exp(long_log_mu[(s-1)*n_ints + i]));
    dev1 = dev1 - 2 * lpn_BK[s,i];
  }
}
uncounted = sum(unobserved);
# This calculation is split into 2 pieces, because of integers, reals, and the C++ oddity that: int/int = int
p_global = sum(obsN); 
p_global = p_global / sum(totN);
}

