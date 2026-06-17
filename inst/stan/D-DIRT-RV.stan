functions {
  real resp_sign(int r) {
    return (r == 1) ? 1.0 : -1.0;
  }
}

data {
  int<lower=1> nObs;       // total number of observations
  int<lower=1> nPerson;    // number of persons
  int<lower=1> nItem;      // number of items

  vector<lower=1e-6>[nPerson] tauUpper;

  array[nObs] int<lower=1, upper=nPerson> person; // Subject identifier
  array[nObs] int<lower=1, upper=nItem> item; // item identifier
  vector<lower=1e-6>[nObs] rt; // reaction times vector
  array[nObs] int<lower=0, upper=1> resp; // 0 = rejected, 1 = accepted

  int<lower=1, upper=3> omega_theta_prior_family; // prior on omega_theta
  real omega_theta_prior_par1;
  real<lower=1e-6> omega_theta_prior_par2;

  int<lower=1, upper=3> omega_gamma_prior_family;
  real omega_gamma_prior_par1;
  real<lower=1e-6> omega_gamma_prior_par2;

  int<lower=1, upper=3> nu_prior_family;
  real nu_prior_par1;
  real<lower=1e-6> nu_prior_par2;

  int<lower=1, upper=3> a_prior_family;
  real a_prior_par1;
  real<lower=1e-6> a_prior_par2;

  int<lower=1, upper=3> tnd_prior_family;
  real tnd_prior_par1;
  real<lower=1e-6> tnd_prior_par2;

  int<lower=1, upper=3> s_delta_prior_family;
  real s_delta_prior_par1;
  real<lower=1e-6> s_delta_prior_par2;

  int<lower=1, upper=2> s_beta_prior_family;
  real s_beta_prior_par1;
  real<lower=1e-6> s_beta_prior_par2;

}

parameters {
  // Drift parameterization
  vector[nPerson] z_theta; // z-standardized person drift
  vector[nItem] nu;          // item drift effects
  real<lower=1e-6> omega_theta; // standard deviation of person drift effects (on the log scale)
  real<lower=1e-6> s_delta; // standard deviation of the Gaussian drift rate variability, constant across subjects and items

  // Boundary separation parameterization on log scale
  vector[nPerson] z_gamma; // z-standardized person response caution
  vector<lower=0.01>[nItem] a;     // item boundary effects
  real<lower=1e-6> omega_gamma; // standard deviation of person boundary effects (on the log scale)
  real<lower = .01, upper = .99> s_beta; // range of the unifirm starting value variability

  // Nondecision time
  vector<lower=1e-6, upper=tauUpper>[nPerson] tnd;// non-decision time is a person parameter

}

transformed parameters{
  vector[nPerson] theta;
  theta = omega_theta * z_theta; //  person drift
  vector[nPerson] log_gamma;
  vector<lower=1e-6>[nPerson] gamma;
  log_gamma = omega_gamma * z_gamma;
  gamma = exp(log_gamma); // boundary effects, assumed to be lognormal, logmean set to 0
}

model {
  // Prior for person drift variability
  if (omega_theta_prior_family == 1) {
    omega_theta ~ lognormal(omega_theta_prior_par1, omega_theta_prior_par2);
  } else if (omega_theta_prior_family == 2) {
    omega_theta ~ normal(omega_theta_prior_par1, omega_theta_prior_par2);
  }else if (omega_theta_prior_family == 3) {
    omega_theta ~ uniform(omega_theta_prior_par1, omega_theta_prior_par2);
  }
  z_theta ~ normal(0, 1);

  // Priors for person boundary
  if (omega_gamma_prior_family == 1) {
    omega_gamma ~ lognormal(omega_gamma_prior_par1, omega_gamma_prior_par2);
  } else if (omega_gamma_prior_family == 2) {
    omega_gamma ~ normal(omega_gamma_prior_par1, omega_gamma_prior_par2);
  }else if (omega_gamma_prior_family == 3) {
    omega_gamma ~ uniform(omega_gamma_prior_par1, omega_gamma_prior_par2);
  }
  z_gamma ~ normal(0, 1);  // Vielleicht doch lieber Gauss?

  // Prior for nondecision time
   for (p in 1:nPerson) {
    if (tnd_prior_family == 1) {
      tnd[p] ~ lognormal(tnd_prior_par1, tnd_prior_par2); //T[0, tnd_upper[p]];
    } else if (tnd_prior_family == 2) {
      tnd[p] ~ normal(tnd_prior_par1, tnd_prior_par2); //T[0, tnd_upper[p]];
    } else if (tnd_prior_family == 3) {
      tnd[p] ~ uniform(tnd_prior_par1, tnd_prior_par2); //T[0, tnd_upper[p]];
    }
  }

  // Priors for item drift
  if (nu_prior_family == 1) {
    nu ~ lognormal(nu_prior_par1, nu_prior_par2);
  }else if (nu_prior_family == 2) {
    nu ~ normal(nu_prior_par1, nu_prior_par2);
  }else if (nu_prior_family == 3) {
    nu ~ uniform(nu_prior_par1, nu_prior_par2);
  }


  // Priors for item boundary aka item time pressure
  if (a_prior_family == 1) {
    a ~ lognormal(a_prior_par1, a_prior_par2);
  } else if (a_prior_family == 2) {
    a ~ normal(a_prior_par1, a_prior_par2);
  }else if (a_prior_family == 3) {
    a ~ uniform(a_prior_par1, a_prior_par2);
  }

  // priors trial-to-trial drift rate variability
  if (s_delta_prior_family == 1) {
    s_delta ~ lognormal(s_delta_prior_par1, s_delta_prior_par2);
  } else if (s_delta_prior_family == 2) {
    s_delta ~ normal(s_delta_prior_par1, s_delta_prior_par2);
  }else if (s_delta_prior_family == 3) {
    s_delta ~ uniform(s_delta_prior_par1, s_delta_prior_par2);
  }

  // priors trial-to-trial stating value variability
  if (s_beta_prior_family == 1) {
    s_beta ~ beta(s_beta_prior_par1, s_beta_prior_par2);
  }else if (s_beta_prior_family == 2) {
    s_beta ~ uniform(s_beta_prior_par1, s_beta_prior_par2);
  }

  // Likelihood
  for (n in 1:nObs) {
    int p = person[n];
    int i = item[n];

    real delta = theta[p] - nu[i];
    real alpha = gamma[p] / a[i];
    real delta_eff = resp_sign(resp[n]) * delta;

    rt[n] ~ wiener(alpha, tnd[p], 0.5, delta_eff, s_delta, s_beta, 0);
  }
}
