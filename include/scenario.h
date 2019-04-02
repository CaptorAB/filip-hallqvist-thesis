#ifndef SCENARIO_H
#define SCENARIO_H

// Instrument indices
const int DOMESTIC_EQUITY_INDEX = 0;
const int GLOBAL_EQUITY_INDEX = 1;
const int REAL_ESTATE_INDEX = 2;
const int ALTERNATIVE_INDEX = 3;
const int CREDIT_INDEX = 4;
const int BONDS_2Y_INDEX = 5;
const int BONDS_5Y_INDEX = 6;
const int CASH_INDEX = 7;
const int FTA_INDEX = 8;
const int DOMESTIC_EQUITY_FUTURE_INDEX = 9;
const int INTEREST_RATE_SWAP_2Y_INDEX = 10;
const int INTEREST_RATE_SWAP_5Y_INDEX = 11;
const int INTEREST_RATE_SWAP_20Y_INDEX = 12;

// Risk indices
const int DOMESTIC_MARKET_RISK_INDEX = 0;
const int GLOBAL_MARKET_RISK_INDEX = 1;
const int ALTERNATIVE_RISK_INDEX = 2;
const int INTEREST_RATE_2Y_RISK_INDEX = 3;
const int INTEREST_RATE_5Y_RISK_INDEX = 4;
const int INTEREST_RATE_20Y_RISK_INDEX = 5;
const int CREDIT_RISK_INDEX = 6;
const int CASH_RISK_INDEX = 7;

const std::vector<double> RISK_CORRELATIONS = {
    1.0, 0.0,
    0.0, 1.0};

std::tuple<std::vector<double>, std::vector<double>> generate_scenarios(int n_steps);
std::vector<double> generate_goals(std::vector<double> price_changes, int n_steps, int n_scenarios, int n_instruments, double surplus);

#endif