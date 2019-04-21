#ifndef SCENARIO_H
#define SCENARIO_H

std::tuple<std::vector<double>, std::vector<double>> generate_scenarios(int n_steps);
std::vector<double> generate_goals(std::vector<double> &price_changes, int n_steps, int n_scenarios, int n_instruments, double initial_funding_ratio, double target_funding_ratio);

#endif