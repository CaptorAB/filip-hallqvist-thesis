#ifndef NORMAL_H
#define NORMAL_H

#include <tuple>
#include <vector>

std::vector<double>
generate_normal_risk_changes(std::vector<double> &means,
                             std::vector<double> &standard_deviations,
                             std::vector<double> &correlations,
                             const int n_risks);

std::vector<double>
generate_normal_scenario(std::vector<double> &means,
                         std::vector<double> &standard_deviations,
                         std::vector<double> &correlations,
                         const int n_risks,
                         const int n_instruments,
                         const int n_scenarios);

std::vector<double>
generate_normal_scenarios(std::vector<double> &means,
                          std::vector<double> &standard_deviations,
                          std::vector<double> &correlations,
                          const int n_risks,
                          const int n_instruments,
                          const int n_scenarios);

#endif