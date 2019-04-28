#ifndef GOALS_H
#define GOALS_H

#include <tuple>
#include <vector>

std::tuple<std::vector<double>, std::vector<double>>
generate_goals(std::vector<double> &instrument_changes,
               const double initial_funding_ratio,
               const double target_funding_ratio, const int n_scenarios,
               const int n_instruments);

#endif