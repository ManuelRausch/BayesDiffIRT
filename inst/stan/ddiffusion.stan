functions {
  real resp_sign(int r) {
    return (r == 1) ? 1.0 : -1.0;
  }
}

data {
  int<lower=1> nObs;       // total number of observations
  int<lower=1> nPerson;    // number of persons
  int<lower=1> nItem;      // number of items

  vector<lower=0>[nPerson] tau_upper;

  array[nObs] int<lower=1, upper=nPerson> person; // Subject identifier
  array[nObs] int<lower=1, upper=nItem> item; // item identifier
  vector<lower=0>[nObs] rt; // reaction times vector
  array[nObs] int<lower=0, upper=1> resp; // 0 = rejected, 1 = accepted

  int<lower=1, upper=3> omega_theta_prior_family; // prior on omega_theta
  real omega_theta_prior_par1;
  real<lower=0> omega_theta_prior_par2;

  int<lower=1, upper=3> omega_gamma_prior_family;
  real omega_gamma_prior_par1;
  real<lower=0> omega_gamma_prior_par2;

  int<lower=1, upper=2> nu_prior_family;
  real nu_prior_par1;
  real<lower=0> nu_prior_par2;

  int<lower=1, upper=3> a_prior_family;
  real a_prior_par1;
  real<lower=0> a_prior_par2;

  int<lower=1, upper=3> tnd_prior_family;
  real tnd_prior_par1;
  real<lower=0> tnd_prior_par2;

}

parameters {
  // Drift parameterization
  vector[nPerson] theta;     // person drift effects, assumed to be Gaussian, mean set to 0.
  vector[nItem] nu;          // item drift effects
  real<lower=0> omega_theta; // standard deviation of person drift effects (on the log scale)

  // Boundary separation parameterization on log scale
  vector<lower=0>[nPerson] gamma; // person boundary effects, assumed to be lognormal, logmean set to 0
  vector<lower=0.01>[nItem] a;     // item boundary effects
  real<lower=0> omega_gamma; // standard deviation of person boundary effects (on the log scale)

  // Nondecision time
 //  vector<lower=0>[nPerson] tau; // non-decision time is a person parameter
  vector<lower=0, upper=tau_upper>[nPerson] tnd;

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
  theta ~ normal(0, omega_theta);

  // Priors for person boundary
  if (omega_gamma_prior_family == 1) {
    omega_gamma ~ lognormal(omega_gamma_prior_par1, omega_gamma_prior_par2);
  } else if (omega_gamma_prior_family == 2) {
    omega_gamma ~ normal(omega_gamma_prior_par1, omega_gamma_prior_par2);
  }else if (omega_gamma_prior_family == 3) {
    omega_gamma ~ uniform(omega_gamma_prior_par1, omega_gamma_prior_par2);
  }
  gamma ~ lognormal(0, omega_gamma);

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
    nu ~ normal(nu_prior_par1, nu_prior_par2);
  }else if (nu_prior_family == 2) {
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
  // Likelihood
  for (n in 1:nObs) {
    int p = person[n];
    int i = item[n];

    real delta = theta[p] - nu[i];
    real alpha = gamma[p] / a[i];
    real delta_eff = resp_sign(resp[n]) * delta;

    rt[n] ~ wiener(alpha, tnd[p], 0.5, delta_eff);
  }
}

// generated quantities {
//   vector[nObs] log_lik;
//
//   for (n in 1:nObs) {
//     int p = person[n];
//     int i = item[n];
//
//     real delta = theta[p] - nu[i];
//     real alpha = gamma[p] / a[i];
//     real delta_eff = resp_sign(resp[n]) * delta;
//     log_lik[n] = wiener_lpdf(rt[n] | alpha, tnd[p], .5, delta_eff);
//   }}

