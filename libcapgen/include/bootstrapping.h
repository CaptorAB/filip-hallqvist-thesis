#ifndef BOOTSTRAPPING_H
#define BOOTSTRAPPING_H

#include <vector>

using namespace std;

void adjust_credit_par_rates(vector<double> &par_rates);

vector<double> bootstrap_discount_factors(vector<double> &par_rates);

vector<double> compute_forward_rates(vector<double> &discount_factors);

void adjust_forward_rates_ufr(vector<double> &forward_rates);

#endif