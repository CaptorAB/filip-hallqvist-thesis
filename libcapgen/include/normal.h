#ifndef NORMAL_H
#define NORMAL_H

#include <tuple>
#include <vector>

std::vector<double> compute_cholesky(
    std::vector<double> &A,
    const int n);

std::vector<double> generate_risk_changes(
    std::vector<double> &means,
    std::vector<double> &standard_deviations,
    std::vector<double> &correlations,
    const int n_risks);

std::vector<double> generate_normal_scenario(
    std::vector<double> &means,
    std::vector<double> &standard_deviations,
    std::vector<double> &correlations,
    const int n_risks,
    const int n_instruments,
    const int n_scenarios);

std::vector<double> generate_normal_scenarios(
    std::vector<double> &means,
    std::vector<double> &standard_deviations,
    std::vector<double> &correlations,
    const int n_risks,
    const int n_instruments,
    const int n_scenarios);

std::tuple<std::vector<double>, std::vector<double>> generate_normal_goals(
    std::vector<double> &instrument_changes,
    const double initial_funding_ratio,
    const double target_funding_ratio,
    const int n_scenarios,
    const int n_instruments);

std::vector<double> generate_instrument_changes(
    std::vector<double> &risk_changes,
    const int n_instruments);

#endif