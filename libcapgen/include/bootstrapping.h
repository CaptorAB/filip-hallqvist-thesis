#ifndef BOOTSTRAPPING_H
#define BOOTSTRAPPING_H

#include <vector>

void adjust_credit_par_rates(std::vector<double> &par_rates);

std::vector<double> bootstrap(std::vector<double> &par_rates);

#endif