#ifndef INSTRUMENTS_H
#define INSTRUMENTS_H

#include <vector>

double sample_domestic_equity(std::vector<double> &risk_changes);

double sample_global_equity(std::vector<double> &risk_changes);

double sample_real_estate(std::vector<double> &risk_changes);

double sample_alternative(std::vector<double> &risk_changes);

double sample_credit(std::vector<double> &risk_changes);

double sample_bonds_2y(std::vector<double> &risk_changes);

double sample_bonds_5y(std::vector<double> &risk_changes);

double sample_cash(std::vector<double> &risk_changes);

double sample_fta(std::vector<double> &risk_changes);

double sample_domestic_equity_future(std::vector<double> &risk_changes);

double sample_interest_rate_swap_2y(std::vector<double> &risk_changes);

double sample_interest_rate_swap_5y(std::vector<double> &risk_changes);

double sample_interest_rate_swap_20y(std::vector<double> &risk_changes);

double generate_instrument_change(
    const int instrument_index,
    std::vector<double> &risk_changes);

std::vector<double> generate_instrument_changes(
    std::vector<double> &risk_changes,
    const int n_instruments);

#endif