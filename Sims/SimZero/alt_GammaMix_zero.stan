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
real<upper=15> intcpt_a;     # Abundance intercept
real intcpt_d;     # Detection intercept
real<lower=0> alpha;  # Universal alpha from the standard Gamma distribution
real<lower=0,upper=1> gamma; # Mixing parameters
}



transformed parameters {
real log_beta;                 # Site-level log(Expected detection rate) -- this is Beta in the standard Gamma distribution
real<lower=0> beta;            # Site-level detection rate
real<upper=20> log_lambda;     # Site-level log(Expected abundance) --- 'upper' is a constraint on early iterations
vector<upper=0>[n_ints] log_p; # Interval-specific log(detection probability)
real log_mu_ab[n_ints];        # Interval-specific log(Expected count = lambda * Pr(detection in interval))
vector[n_sites*n_ints] long_log_mu; # Site- and interval-specific log_mu_ab (i.e., replicated over all sites)
real<lower=0> lambda;          # Site-level Expected abundance

# Calculating site-level expected abundance
log_lambda = intcpt_a;
lambda     = exp(log_lambda);

# Calculating site-level rate of detection -- log(beta) = Xd*bd + ObsEff + log(alpha)
log_beta = intcpt_d;
beta     = exp(log_beta);

# Calculating interval-specific log-scale detection probabilities (log_p) and log-scale Poisson means (log_mu_ab)
log_p[1] = log(gamma * gamma_cdf(tau[1], alpha, beta) + 1 - gamma);
log_mu_ab[1] = log_lambda + log_p[1];
for (i in 2:n_ints) {
  log_p[i] = log(gamma) + log_diff_exp(gamma_lcdf(tau[i] | alpha, beta),
                                       gamma_lcdf(tau[i-1] | alpha, beta));
  log_mu_ab[i] = log_lambda+log_p[i];
}
# Replicate log_mu_ab to match yv in dimension
for (s in 1:n_sites) for (i in 1:n_ints) long_log_mu[ii[s]+i] = log_mu_ab[i];
}



model {
intcpt_a ~ normal(1.6,0.5);
intcpt_d ~ normal(-2.7,2.25);
alpha ~ cauchy(0,1);
gamma ~ beta(1,1);

# Data models
yv ~ poisson_log(long_log_mu);
}



generated quantities {
int<lower=0> uncounted;             # Total uncounted birds
real dev1;
real lpn_BK[n_ints];                # log(p(n|\beta's,\ksi's))
real<lower=0> p_global;             # Estimated overall detection probability

# Because of homogeneity and the properties of Poissons, 
# I can calculate likelihood over all survey en masse rather than one survey at a time
dev1 = 0;
for (i in 1:n_ints) {
  lpn_BK[i] = poisson_lpmf(sum(y[,i]) | n_sites*exp(log_mu_ab[i]));
  dev1 = dev1 - 2 * lpn_BK[i];
}
# Draw total unobserved abundance from Poisson(lambda * (1-pdet))
uncounted = poisson_rng(n_sites*exp(log_lambda + log(fmax(1-exp(log_sum_exp(log_p)),1e-8))));

# This calculation is split into 2 pieces, because of integers, reals, and the C++ oddity that: int/int = int
p_global = sum(obsN); 
p_global = p_global / (p_global + uncounted);
}

