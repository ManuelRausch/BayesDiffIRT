functions {
  real resp_sign(int r) {
    return (r == 1) ? 1.0 : -1.0;
  }

  real resp_beta(int r, real beta) {
    return (r == 1) ? beta : 1.0 - beta;
  }
}

data {
  int<lower=1> nObs;       // total number of observations
  int<lower=1> nPerson;    // number of persons
  int<lower=1> nItem;      // number of items

  array[nObs] int<lower=1, upper=nPerson> person;
  array[nObs] int<lower=1, upper=nItem> item;

  vector<lower=0>[nObs] rt;               // response times
  array[nObs] int<lower=0, upper=1> resp; // 1 = upper boundary, 0 = lower boundary

  real<lower=0> tau_lower;                // lower bound for nondecision time
}

transformed data {
  real<lower=0, upper=1> beta = 0.5;      // fixed starting point
}

parameters {
  // Drift parameterization
  vector[nPerson] theta;                  // person drift effects, mean anchored at 0 via prior
  vector[nItem] nu;                       // item drift effects, free mean via mu_nu
  real<lower=0> sigma_theta;
  real mu_nu;
  real<lower=0> sigma_nu;

  // Boundary separation parameterization on log scale
  vector[nPerson] log_gamma;              // person boundary effects, mean anchored at 0 via prior
  vector[nItem] log_a;                    // item boundary effects, free mean via mu_log_a
  real<lower=0> sigma_log_gamma;
  real mu_log_a;
  real<lower=0> sigma_log_a;

  // Nondecision time
  real<lower=tau_lower> tau;
}

model {
  // Priors for drift part
  sigma_theta ~ normal(0, 1);
  sigma_nu    ~ normal(0, 1);
  mu_nu       ~ normal(0, 2);

  theta ~ normal(0, sigma_theta);
  nu    ~ normal(mu_nu, sigma_nu);

  // Priors for boundary part
  sigma_log_gamma ~ normal(0, 1);
  sigma_log_a     ~ normal(0, 1);
  mu_log_a        ~ normal(0, 1);

  log_gamma ~ normal(0, sigma_log_gamma);
  log_a     ~ normal(mu_log_a, sigma_log_a);

  // Prior for nondecision time
  tau ~ normal(0.3, 0.2);

  // Likelihood
  for (n in 1:nObs) {
    int p = person[n];
    int i = item[n];

    real delta = theta[p] - nu[i];
    real alpha = exp(log_gamma[p] - log_a[i]);
    real delta_eff = resp_sign(resp[n]) * delta;
    real beta_eff = resp_beta(resp[n], beta);

    rt[n] ~ wiener(alpha, tau, beta_eff, delta_eff);
  }
}

generated quantities {
  vector[nObs] log_lik;
  vector<lower=0>[nPerson] gamma;
  vector<lower=0>[nItem] a;

  for (p in 1:nPerson) {
    gamma[p] = exp(log_gamma[p]);
  }

  for (i in 1:nItem) {
    a[i] = exp(log_a[i]);
  }

  for (n in 1:nObs) {
    int p = person[n];
    int i = item[n];

    real delta = theta[p] - nu[i];
    real alpha = exp(log_gamma[p] - log_a[i]);
    real delta_eff = resp_sign(resp[n]) * delta;
    real beta_eff = resp_beta(resp[n], beta);

    log_lik[n] = wiener_lpdf(rt[n] | alpha, tau, beta_eff, delta_eff);
  }
}
