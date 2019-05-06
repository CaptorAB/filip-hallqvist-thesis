#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>

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
correlate_normal_randoms(
    vector<double> &normals,
    vector<double> &correlations,
    const int n_generic_risks);

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
    const int n_generic_risks);

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

#endif