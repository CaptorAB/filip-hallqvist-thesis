#ifndef BOOTSTRAPPING_H
#define BOOTSTRAPPING_H

#include <vector>

void adjust_credit_par_rates(std::vector<double> &par_rates);

std::vector<double> compute_dfs_from_pars(std::vector<double> &par_rates);

#endif