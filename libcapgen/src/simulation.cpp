#include <vector>

#include <lib/stats/stats.hpp>

#include <include/constants.h>
#include <include/cholesky.h>
#include <include/simulation.h>

using namespace std;
using namespace stats;

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

vector<double>
normalize_uniform_randoms(
    vector<double> &uniforms)
{
  vector<double> normalized(uniforms.size());
  for (int i = 0; i < uniforms.size(); ++i)
  {
    normalized[i] = qnorm(uniforms[i]);
  }
  return normalized;
}

double
sample_nln(
    const double n1,
    const double n2,
    const double sigma,
    const double rho)
{
  const double phi = n1;
  const double eta = sigma * (rho * n1 + sqrt(1 - pow(rho, 2)) * n2);

  return exp(0.5 * eta) * phi;
}

double
standardize_nln(
    const double u,
    const double sigma,
    const double rho)
{
  const double mu = 0.5 * rho * sigma * exp(0.125 * pow(sigma, 2.0));
  const double var = exp(0.5 * pow(sigma, 2.0)) * (1.0 + pow(rho, 2.0) * pow(sigma, 2.0) * (1.0 - 0.25 * exp(-0.25 * pow(sigma, 2.0))));
  return (u - mu) / sqrt(var);
}

vector<double>
correlate_risks(
    vector<double> &normals,
    vector<double> &correlations,
    const int n_generic_risks)
{
  const int n_correlations = n_generic_risks + 1;
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
sample_standardized_nlns(
    vector<double> &sigmas,
    vector<double> &rhos,
    vector<double> &normals,
    const int n_generic_risks)
{
  const int n_nlns_per_forward_rate_risk = 3;
  const int n_nlns_per_generic_risk = 1;

  vector<double> nlns(n_nlns_per_generic_risk + n_nlns_per_forward_rate_risk);

  // Sample nlns for market risks
  for (int i = 0; i < n_generic_risks; ++i)
  {
    const int n1 = normals[i];
    const int n2 = normals[(normals.size() / 2) + i];
    const double u = sample_nln(n1, n2, sigmas[i], rhos[i]);
    nlns[i] = standardize_nln(u, sigmas[i], rhos[i]);
  }

  // Sample three nlns for interest rate risks
  for (int i = 0; i < n_nlns_per_forward_rate_risk; ++i)
  {
    const int n1 = normals[n_generic_risks + i];
    const int n2 = normals[(normals.size() / 2) + n_generic_risks + i];
    const double u = sample_nln(n1, n2, sigmas[n_generic_risks + 1], rhos[n_generic_risks + 1]);
    nlns[n_nlns_per_generic_risk + i] = standardize_nln(u, sigmas[n_generic_risks + 1], rhos[n_generic_risks + 1]);
  }

  return nlns;
}

vector<double>
generate_epsilons(
    vector<double> &stds,
    vector<double> &sigmas,
    vector<double> &rhos,
    vector<double> &eigenvalues,
    vector<double> &eigenvectors,
    vector<double> &correlations,
    const int n_pca_components,
    const int n_scenarios,
    const int n_generic_risks,
    const int n_forward_rate_risks)
{
  const int n_risks = n_generic_risks + n_forward_rate_risks;
  vector<double> epsilons(n_scenarios * n_risks);
  for (int i = 0; i < n_scenarios; ++i)
  {
    const int ix = i * n_scenarios;

    vector<double> uniforms = sample_uniform_randoms(2 * (n_generic_risks * n_pca_components));
    vector<double> normals = normalize_uniform_randoms(uniforms);
    vector<double> normals_correlated = correlate_risks(normals, correlations, n_generic_risks);
    vector<double> nlns = sample_standardized_nlns(sigmas, rhos, normals, n_generic_risks);

    const double tau = floor(log2(i + 1));

    vector<double> generic_risk_nlns = vector<double>(nlns.begin(), nlns.begin() + n_generic_risks);
    vector<double> forward_rate_risk_nlns = vector<double>(nlns.begin() + n_generic_risks, nlns.end());

    vector<double> generic_risk_epsilons = compute_generic_risk_epsilons(
        stds,
        generic_risk_nlns,
        tau,
        n_generic_risks);

    vector<double> forward_rate_risk_epsilons = compute_forward_rate_risk_epsilons(
        eigenvalues,
        eigenvectors,
        forward_rate_risk_nlns,
        tau,
        n_pca_components,
        n_forward_rate_risks);

    for (int j = 0; j < n_generic_risks; ++j)
    {
      epsilons[j] = generic_risk_epsilons[j];
    }
    for (int j = 0; j < n_forward_rate_risks; ++j)
    {
      epsilons[j] = forward_rate_risk_epsilons[j];
    }

    return epsilons;
  }
}

vector<double>
compute_generic_risk_epsilons(
    vector<double> &stds,
    vector<double> &nlns,
    const double tau,
    const int n_generic_risks)
{
  vector<double> epsilons(n_generic_risks);
  for (int i = 0; i < n_generic_risks; ++i)
  {
    epsilons[i] = stds[i] * sqrt(tau) * nlns[i];
  }
  return epsilons;
}

vector<double>
compute_forward_rate_risk_epsilons(
    vector<double> &eigenvalues,
    vector<double> &eigenvectors,
    vector<double> &nlns,
    const double tau,
    const int n_pca_components,
    const int n_forward_rate_risks)
{
  vector<double> epsilons(n_forward_rate_risks);
  for (int i = 0; i < n_forward_rate_risks; ++i)
  {
    for (int j = 0; j < n_pca_components; ++j)
    {
      const int jx = (i * eigenvalues.size()) + j;
      epsilons[i] += nlns[j] * sqrt(eigenvalues[j]) * sqrt(tau) * eigenvectors[jx];
    }
  }
  return epsilons;
}

vector<double>
compute_gammas(
    vector<double> &epsilons,
    const int n_risks,
    const int n_scenarios)
{
  vector<double> gammas(n_risks);

  for (int i = 0; i < n_scenarios; ++i)
  {
    const int ix = i * n_scenarios;
    for (int j = 0; j < n_risks; ++j)
    {
      gammas[j] += epsilons[ix + j];
    }
  }

  for (int j = 0; j < n_risks; ++j)
  {
    gammas[j] /= (double)n_scenarios;
  }

  return gammas;
}

double
evaluate_black_76(
    const double mean,
    const double gamma,
    const double epsilon)
{
  return mean * exp(-gamma + epsilon);
}

vector<double>
evaluate_risk_processes(
    vector<double> &generic_means,
    vector<double> &forward_rate_means,
    vector<double> &gammas,
    vector<double> &epsilons,
    const int n_generic_risks,
    const int n_forward_rate_risks,
    const int n_scenarios)
{
  const int n_risks = n_generic_risks + n_forward_rate_risks;
  vector<double> risk_values(n_risks);
  for (int i = 0; i < n_scenarios; ++i)
  {
    int ix = i * n_scenarios;
    for (int j = 0; j < n_generic_risks; ++j)
    {
      risk_values[ix + j] = evaluate_black_76(
          generic_means[j],
          gammas[j],
          epsilons[ix + j]);
    }

    ix += n_generic_risks;
    for (int j = 0; j < n_forward_rate_risks; ++j)
    {
      risk_values[ix + j] = evaluate_black_76(
          forward_rate_means[j],
          gammas[j],
          epsilons[ix + j]);
    }
  }

  return risk_values;
}

vector<double>
compute_risk_changes(
    vector<double> &generic_means,
    vector<double> &forward_rate_means,
    vector<double> &gammas,
    vector<double> &initial_epsilons,
    vector<double> &epsilons,
    const int n_generic_risks,
    const int n_forward_rate_risks,
    const int n_scenarios)
{
  vector<double> risk_values = evaluate_risk_processes(
      generic_means,
      forward_rate_means,
      gammas,
      initial_epsilons,
      n_generic_risks,
      n_forward_rate_risks,
      n_scenarios);
  vector<double> initial_risk_values = evaluate_risk_processes(
      generic_means,
      forward_rate_means,
      gammas,
      initial_epsilons,
      n_generic_risks,
      n_forward_rate_risks,
      1);

  const int n_risks = n_generic_risks + n_forward_rate_risks;
  vector<double> risk_changes(n_scenarios * n_risks);

  for (int j = 0; j < n_risks; ++j)
  {
    risk_changes[j] = risk_values[j] / initial_risk_values[j];
  }

  for (int i = 1; i < n_scenarios; ++i)
  {
    const int p = ((i - 1) / 2);

    const int px = p * n_risks; // Parent
    const int ix = i * n_risks; // Child

    for (int j = 0; j < n_risks; ++j)
    {
      risk_changes[ix + j] = risk_values[ix + j] / risk_values[px + j];
    }
  }

  return risk_changes;
}

vector<double>
generate_risk_changes(
    vector<double> &generic_means,
    vector<double> &generic_stds,
    vector<double> &forward_rate_means,
    vector<double> &forward_rate_eigenvalues,
    vector<double> &forward_rate_eigenvectors,
    vector<double> &correlations,
    vector<double> &sigmas,
    vector<double> &rhos,
    const int n_generic_risks,
    const int n_forward_rate_risks,
    const int n_pca_components,
    const int n_scenarios,
    const int n_correlations)
{
  // Compute epsilons
  vector<double> epsilons = generate_epsilons(
      generic_stds,
      sigmas,
      rhos,
      forward_rate_eigenvalues,
      forward_rate_eigenvectors,
      correlations,
      n_pca_components,
      n_scenarios,
      n_generic_risks,
      n_forward_rate_risks);
  vector<double> initial_epsilons = generate_epsilons(
      generic_stds,
      sigmas,
      rhos,
      forward_rate_eigenvalues,
      forward_rate_eigenvectors,
      correlations,
      n_pca_components,
      1,
      n_generic_risks,
      n_forward_rate_risks);

  // Compute gammas
  const int n_risks = n_generic_risks + n_forward_rate_risks;
  vector<double> gammas = compute_gammas(epsilons, n_risks, n_scenarios);

  // Compute risk changes
  vector<double> risk_changes = compute_risk_changes(
      generic_means,
      forward_rate_means,
      gammas,
      initial_epsilons,
      epsilons,
      n_generic_risks,
      n_forward_rate_risks,
      n_scenarios);

  return risk_changes;
}