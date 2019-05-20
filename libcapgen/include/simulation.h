#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <tuple>

using namespace std;

/**
 * Sample a number from 
 * the normal-log-normal mixture
 * distribution.
 */
double sample_nln(
    const double n1,
    const double n2,
    const double sigma,
    const double rho);

/**
 * Standardize a normal-log-normally
 * distributed number:
 * u_ = (u - mean) / std
 */
double standardize_nln(
    const double u,
    const double sigma,
    const double rho);

/**
 * Sample a sequence of uniformly
 * distributed random numbers.
 */
vector<double>
sample_uniform_randoms(
    const int n_randoms);

/**
 * Convert uniformly distributed
 * random numbers to the normal standard
 * distribution.
 */
vector<double> normalize_uniform_randoms(
    vector<double> &uniforms);

/**
 * Correlated a series
 * of normally distributed random
 * numbers
 */
vector<double>
correlate_risks(
    vector<double> &normals,
    vector<double> &correlations,
    const int n_generic_risks,
    const int n_pca_components);

/**
 * Sample a series of standardized
 * normal-log-normally distributed
 * random numbers
 */
vector<double>
sample_standardized_nlns(
    vector<double> &sigmas,
    vector<double> &rhos,
    vector<double> &normals,
    const int n_generic_risks,
    const int n_pca_components);

/**
 * Generate epsilons (shocks)
 * to be used when simulating new
 * risk process changes.
 */
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
    const int n_forward_rate_risks);

/**
 * Generate epsilons (shocks)
 * for generic type risks.
 */
vector<double>
compute_generic_risk_epsilons(
    vector<double> &stds,
    vector<double> &nlns,
    const double tau,
    const int n_generic_risks);

/**
 * Generate epsilons (shocks)
 * for interest rate risks.
 */
vector<double>
compute_forward_rate_risk_epsilons(
    vector<double> &eigenvalues,
    vector<double> &eigenvectors,
    vector<double> &nlns,
    const double tau,
    const int n_pca_components,
    const int n_forward_rate_risks);

tuple<vector<double>, vector<double>>
compute_gammas(
    vector<double> &generic_risk_stds,
    vector<double> &pca_forward_rate_risk_eigenvalues,
    vector<double> &pca_forward_rate_risk_eigenvectors,
    vector<double> &sigmas,
    vector<double> &rhos,
    vector<double> &correlations,
    const int n_trials,
    const int n_pca_components,
    const int n_generic_risks,
    const int n_forward_rate_risks,
    const int n_states,
    const int n_steps);

double
evaluate_black_76(
    const double mean,
    const double gamma,
    const double epsilon);

vector<double>
evaluate_risk_processes(
    vector<double> &generic_means,
    vector<double> &forward_rate_means,
    vector<double> &gammas,
    vector<double> &epsilons,
    const int n_generic_risks,
    const int n_forward_rate_risks,
    const int n_scenarios);

vector<double>
compute_risk_changes(
    vector<double> &generic_means,
    vector<double> &forward_rate_means,
    vector<double> &gammas,
    vector<double> &initial_epsilons,
    vector<double> &epsilons,
    const int n_generic_risks,
    const int n_forward_rate_risks,
    const int n_scenarios);

tuple<vector<double>, vector<double>>
generate_risks(
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
    const int n_scenarios);

vector<double> compute_yield_curve(
    vector<double> &forward_rates);

vector<double> generate_state_changes(
    vector<double> &initial_generic_risk_values,
    vector<double> &initial_forward_rate_risk_values,
    vector<double> &generic_risk_means,
    vector<double> &generic_risk_stds,
    vector<double> &pca_forward_rate_risk_eigenvalues,
    vector<double> &pca_forward_rate_risk_eigenvectors,
    vector<double> &sigmas,
    vector<double> &rhos,
    vector<double> &correlations,
    vector<double> &zero_coupon_tenors,
    const int n_instruments,
    const int n_generic_risks,
    const int n_forward_rate_risks,
    const int n_pca_components,
    const int n_states,
    const int n_steps);

#endif
