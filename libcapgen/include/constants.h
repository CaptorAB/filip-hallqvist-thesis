#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <vector>
#include <string>

const double NEGATIVE_FORWARD_RATE_ADJUSTMENT = 0.01;

// Instrument names
const std::vector<std::string> INSTRUMENT_NAMES = {
    "Domestic Equity",
    "Global Equity",
    "Real Estate",
    "Alternative",
    "Credit",
    "Bonds 2Y",
    "Bonds 5Y",
    "Bonds 20Y",
    "Cash",
    "Domestic Equity Future",
    "Interest Rate Swap 2Y",
    "Interest Rate Swap 5Y",
    "Interest Rate Swap 20Y"};

// Risk names
const std::vector<std::string> RISK_NAMES = {
    "Domestic Market Risk",
    "Global Market Risk",
    "Real Estate Risk",
    "Alternative Risk",
    "Credit Risk",
    "Interest Rate Risk"};

// Instrument indices
const int N_RISKS = 6;
const int N_INSTRUMENTS = 13;
const int N_DERIVATIVES = 4;
const int N_TREES = 100;

const int DOMESTIC_EQUITY_INDEX = 0;
const int GLOBAL_EQUITY_INDEX = 1;
const int REAL_ESTATE_INDEX = 2;
const int ALTERNATIVE_INDEX = 3;
const int CREDIT_INDEX = 4;
const int BONDS_2Y_INDEX = 5;
const int BONDS_5Y_INDEX = 6;
const int BONDS_20Y_INDEX = 7;
const int CASH_INDEX = 8;
const int DOMESTIC_EQUITY_FUTURE_INDEX = 9;
const int INTEREST_RATE_SWAP_2Y_INDEX = 10;
const int INTEREST_RATE_SWAP_5Y_INDEX = 11;
const int INTEREST_RATE_SWAP_20Y_INDEX = 12;

// Correlation indices
const int DOMESTIC_MARKET_CORRELATION_INDEX = 0;
const int GLOBAL_MARKET_CORRELATION_INDEX = 1;
const int REAL_ESTATE_CORRELATION_INDEX = 2;
const int ALTERNATIVE_CORRELATION_INDEX = 3;
const int CREDIT_CORRELATION_INDEX = 4;
const int INTEREST_RATE_CORRELATION_INDEX = 5;

// Risk indices
const int DOMESTIC_MARKET_RISK_INDEX = 0;
const int GLOBAL_MARKET_RISK_INDEX = 1;
const int REAL_ESTATE_RISK_INDEX = 2;
const int ALTERNATIVE_RISK_INDEX = 3;
const int CREDIT_RISK_INDEX = 4;
const int FORWARD_RATE_RISK_START_INDEX = 5;
const int FORWARD_RATE_RISK_END_INDEX = 25;

// Normal scenarios
const std::vector<double> NORMAL_DEFAULT_MEANS = {
    0.0,
    0.0,
    0.1,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0};

const std::vector<double> NORMAL_DEFAULT_STANDARD_DEVIATIONS = {
    0.1,
    0.1,
    0.1,
    0.1,
    0.1,
    0.1,
    0.1,
    0.1};

const std::vector<double> NORMAL_DEFAULT_CORRELATIONS = {
    1.0, 0.0004, 0.0, -0.0172, 0.0034, 0.0176, 0.0, 0.0155,
    0.0004, 1.0, 0.0, -0.0036, 0.0340, 0.0202, 0.0, 0.02201,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0172, -0.0036, 0.0, 1.0, 0.0032, 0.0151, 0.0, -0.0185,
    0.0034, 0.0340, 0.0, 0.0032, 1.0, -0.0185, 0.0, 0.0343,
    0.0176, 0.0202, 0.0, 0.0151, -0.0185, 1.0, 0.0, 0.0100,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0155, 0.201, 0.0, -0.0185, 0.0343, 0.01, 0.0, 1.0};

const std::vector<double> DEFAULT_PAR_RATES = {
    0.00970,
    0.01388,
    0.01918,
    0.02402,
    0.02055,
    0.02848,
    0.02975,
    0.03047,
    0.03100,
    0.03091,
    0.03192,
    0.03293,
    0.03375,
    0.03478,
    0.03581,
    0.03699,
    0.03799,
    0.03899,
    0.03999,
    0.04100};

const double PAR_CREDIT_ADJUSTMENT = -0.0035;

#endif