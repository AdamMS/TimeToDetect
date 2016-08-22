# As a prerequisite, point count survey data should be organized in a dataframe saved as 'surveydata'.
surveydata <- read.csv("surveydata.csv")

########## Code to package data as a list for the model

tau <- c(2,3,4,5,6,7,8,9,10)    # Observation interval cutoffs (does not include zero)

### Construct covariate matrices
# Abundance covariates: site age, year, stock density, logging indicator
Xa <- with(surveydata, data.frame(siteage=siteage, 
                                year=year,
                                density=1*(stockdens==9),
                                logged=1*(logged==1)))
# Detection covariates: julian date, time of day, temperature, observer's first year indicator, jd_fy interaction
Xd <- with(surveydata, data.frame(jd=jd,
                                time=time,
                                temp=temp,
                                firstyr=firstyr, 
                                jd_fy=jd*firstyr))

### Construct indices for random effect covariates
# For abundance: site, year
ia <- with(surveydata, data.frame(yearf=as.numeric(factor(year.unscale)), 
                                  stand=as.numeric(factor(standunique))))
# For detection: observer
id <- data.frame(obs=as.numeric(factor(surveydata$obs)))

### Combine multiple effect types into a single vector

# Abundance random effects
lu <- function(x) length(unique(x))  # Helper function for apply statements
n_ras <- apply(ia, 2, lu)            # Number of levels for each abundance random effect type
va_id <- rep(1:ncol(ia), n_ras)      # Index to associate random effect levels to their random effect type
# Renumber indices so that levels of different random effects do not share the same number
for (j in 2:ncol(ia)) ia[,j] <- ia[,j] + cumsum(n_ras)[j-1]

# Detection random effects
n_rds <- apply(id, 2, lu)
vd_id <- rep(1:ncol(id), n_rds)
# Note: because we have only one detection random effect,
#       there is no need to renumber indices cumulatively

### Counts by interval
# Here, W1...W9 are counts by observation interval
Ival <- subset(surveydata, select=c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9"))



### Compile data into a list for direct use by Stan
dat = list(n_sites = nrow(Ival),
           n_ints = ncol(Ival),
           n_ba = ncol(Xa),
           n_bd = ncol(Xd),
           n_ra = dim(ia)[2],
           n_rd = dim(id)[2],
           n_ras = n_ras,
           n_rds = n_rds,
           cn_ras = cumsum(n_ras),
           cn_rds = cumsum(n_rds),
           va_id = va_id,
           vd_id = vd_id,
           Xa = Xa,
           Xd = Xd,
           ia = ia,
           id = id[,1],
           tau = tau,
           y = Ival)

save(dat, file="Data for rStan.Rdata")



########## Run Stan code
load("Data for rStan.Rdata")   # Load Data

# Build list of variables to use.
params <- c("intcpt_a", "ba", "intcpt_d", "bd", "sigma_a", "sigma_d", "gamma", "ra", "rd", "shape_k")

# MCMC variables for which to store output
# - omitting lambda and varphi, since they can be reconstructed
# - 'yrep' are simulated counts from the posterior distribution
# - 'unobserved' are MCMC samples of uncounted individuals = N_s-n_s^{(obs)}
record.list <- c(params, "unobserved", "yrep")

library(rstan)
# Compile model -- takes a few minutes
m = stan_model("WeibullMix.stan", auto_write=rstan_options("auto_write" = TRUE))
# Perform a short 1000-iteration MCMC sampling, thinned by a factor of 20 -- takes 7.5 minutes
fit <- sampling(m, data=dat, chains=1, iter=1000, thin=20, pars=record.list)
