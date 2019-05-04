#include <vector>

#include <include/constants.h>
#include <include/cholesky.h>
#include <lib/stats/stats.hpp>

using namespace std;
using namespace stats;

/*
double sample_nln(
    const double n1,
    const double n2,
    const double sigma,
    const double rho)
{
  const double phi = n1;
  const double eta = sigma * (rho * n1 + sqrt(1 - pow(rho, 2)) * n2);

  return exp(0.5 * eta) * phi;
}

double standardize_nln(
    const double u,
    const double sigma,
    const double rho)
{
  const double mu = 0.5 * rho * sigma * exp(0.125 * pow(sigma, 2.0));
  const double var = exp(0.5 * pow(sigma, 2.0)) * (1.0 + pow(rho, 2.0) * pow(sigma, 2.0) * (1.0 - 0.25 * exp(-0.25 * pow(sigma, 2.0))));
  return (u - mu) / sqrt(var);
}

vector<double>
sample_uniform_randoms(
    const int n_randoms)
{
  vector<double> r(n_randoms);
  for (int i = 0; i < n_randoms; ++i)
  {
    r[i] = runif(0.0, 1.0);
  }
  return r;
}

vector<double> normalize_uniform_randoms(
    vector<double> &uniforms)
{
  vector<double> normalized(uniforms.size());
  for (int i = 0; i < uniforms.size(); ++i)
  {
    normalized[i] = qnorm(uniforms[i]);
  }
  return normalized;
}

vector<double>
correlate_normal_randoms(
    vector<double> &normals,
    vector<double> &correlations,
    const int n_market_risks)
{
  const int n_correlations = n_market_risks + 1;
  const int n_randoms = (2 * n_correlations);

  vector<double> cholesky = compute_cholesky(correlations, n_correlations);
  vector<double> correlated(n_randoms, 0.0);

  for (int i = 0; i < n_correlations; ++i)
  {
    for (int j = 0; j < n_correlations; ++j)
    {
      const int jx = j * n_correlations;
      correlated[i] += normals[j] * cholesky[jx + i];
      correlated[n_correlations + i] += normals[n_correlations + j] * cholesky[jx + i];
    }
  }

  return correlated;
}

vector<double>
sample_nlns(
    vector<double> &sigmas,
    vector<double> &rhos,
    vector<double> &normals,
    const int n_market_risks)
{
  const int n_nlns_per_interest_rate_risk = 3;
  const int n_nlns_per_market_risk = 1;

  vector<double> nlns(n_nlns_per_market_risk + n_nlns_per_interest_rate_risk);

  // Sample nlns for market risks
  for (int i = 0; i < n_market_risks; ++i)
  {
    const int n1 = normals[i];
    const int n2 = normals[(normals.size() / 2) + i];
    const double u = sample_nln(n1, n2, sigmas[i], rhos[i]);
    nlns[i] = standardize_nln(u, sigmas[i], rhos[i]);
  }

  // Sample three nlns for interest rate risks
  for (int i = 0; i < n_nlns_per_interest_rate_risk; ++i)
  {
    const int n1 = normals[n_market_risks + i];
    const int n2 = normals[(normals.size() / 2) + n_market_risks + i];
    const double u = sample_nln(n1, n2, sigmas[n_market_risks + 1], rhos[n_market_risks + 1]);
    nlns[n_nlns_per_market_risk + i] = standardize_nln(u, sigmas[n_market_risks + 1], rhos[n_market_risks + 1]);
  }

  return nlns;
}

vector<double>
generate_shocks(
    vector<double> &sigmas,
    vector<double> &rhos,
    vector<double> &correlations,
    vector<double> &pca_variances,
    vector<double> &pca_loadings,
    const int n_scenarios,
    const int n_market_risks)
{
  vector<double> shocks(n_scenarios * (n_market_risks + 19));
  for (int i = 0; i < n_scenarios; ++i)
  {
    const int ix = i * n_scenarios;

    vector<double> uniforms = sample_uniform_randoms(2 * n_market_risks + 6);
    vector<double> normals = normalize_uniform_randoms(uniforms);
    vector<double> normals_correlated = correlate_normal_randoms(normals, correlations, n_market_risks);
    vector<double> nlns = sample_nlns(sigmas, rhos, normals, n_market_risks);

    for (int j = 0; j < n_market_risks; ++j) {
      shocks[ix + j] = sample_market_risk_change();
    }

    vector<double> interest_rate_nlns = {
      nlns[n_market_risks],
      nlns[n_market_risks + 1],
      nlns[n_market_risks + 2]
    };

  }
}

/*
vector<double>
generate_risk_changes(
    vector<double> &market_means,
    vector<double> &market_stds,
    vector<double> &forward_rate_means,
    vector<double> &forward_rate_pca_variances,
    vector<double> &forward_rate_pca_loadings,
    vector<double> &correlations,
    vector<double> &sigmas,
    vector<double> &rhos,
    vector<double> &gammas,
    const int n_risks,
    const int n_correlations)
{
  vector<double> risk_changes(n_risks);

  vector<double> cholesky = compute_cholesky(correlations, n_risks);

  // Generate correlated normal numbers
  for (int i = 0; i < n_correlations; ++i)
  {
    const int ix = i * n_correlations;

    if (i == INTEREST_RATE_CORRELATION_INDEX)
    {
      double r1 = 0.0;
      double r2 = 0.0;
      double r3 = rnorm(0.0, 1.0);
      double r4 = rnorm(0.0, 1.0);
      double r5 = rnorm(0.0, 1.0);
      double r6 = rnorm(0.0, 1.0);

      for (int j = 0; j < n_risks; ++j)
      {
        r1 += rnorm(0.0, 1.0) * cholesky[ix + j];
        r2 += rnorm(0.0, 1.0) * cholesky[ix + j];
      }

      const double u1 = sample_nln(r1, r2, sigmas[i], rhos[i]);
      const double u2 = sample_nln(r3, r4, sigmas[i], rhos[i]);
      const double u3 = sample_nln(r5, r6, sigmas[i], rhos[i]);

      const double u1_norm = standardize_nln(u1, sigmas[i], rhos[i]);
      const double u2_norm = standardize_nln(u2, sigmas[i], rhos[i]);
      const double u3_norm = standardize_nln(u3, sigmas[i], rhos[i]);

      const int n_forward_rates = FORWARD_RATE_RISK_END_INDEX - FORWARD_RATE_RISK_START_INDEX;
      for (int j = 0; j < n_forward_rates; ++j)
      {
        const int jx = 3 * j;

        const double gamma = gammas[FORWARD_RATE_RISK_START_INDEX + j];

        double epsilon = 0.0;
        epsilon += u1_norm * sqrt(forward_rate_pca_variances[jx]) * forward_rate_pca_loadings[jx];
        epsilon += u2_norm * sqrt(forward_rate_pca_variances[jx + 1]) * forward_rate_pca_loadings[jx + 1];
        epsilon += u3_norm * sqrt(forward_rate_pca_variances[jx + 2]) * forward_rate_pca_loadings[jx + 2];

        const double forward_rate_change = forward_rate_means[i] * exp(-gamma + epsilon) - NEGATIVE_FORWARD_RATE_ADJUSTMENT;

        risk_changes[FORWARD_RATE_RISK_START_INDEX + j] = forward_rate_change;
      }
    }
    else
    {
      double r1 = 0.0;
      double r2 = 0.0;
      for (int j = 0; j < n_risks; ++j)
      {
        r1 += rnorm(0.0, 1.0) * cholesky[ix + j];
        r2 += rnorm(0.0, 1.0) * cholesky[ix + j];
      }
      const double u = sample_nln(r1, r2, sigmas[i], rhos[i]);
      const double u_norm = standardize_nln(u, sigmas[i], rhos[i]);

      const double gamma = gammas[i];
      const double epsilon = market_stds[i] * u_norm;
      const double market_change =
    }
  }

  return risk_changes;
}
*/