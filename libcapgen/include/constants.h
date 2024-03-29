#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <vector>
#include <string>

// Instrument names
const std::vector<std::string> INSTRUMENT_NAMES = {
    "Domestic Equity",
    "Global Equity",
    "Real Estate",
    "Alternative",
    "Credit",
    "Bonds 2Y",
    "Bonds 5Y",
    "Cash",
    "FTA",
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
    "Interest Rate 2Y Risk",
    "Interest Rate 5Y Risk",
    "Interest Rate 20Y Risk"};

// Instrument indices
const int N_RISKS = 8;
const int N_INSTRUMENTS = 13;
const int N_DERIVATIVES = 4;

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
const int FOREX_RISK_INDEX = 7;

// Normal scenarios
const std::vector<double> NORMAL_DEFAULT_MEANS = {
    0.000459,
    0.000414,
    0.0, // TODO
    -0.022266,
    0.0074981,
    -0.000126,
    0.0, // TODO
    0.000042};

const std::vector<double> NORMAL_DEFAULT_STANDARD_DEVIATIONS = {
    0.012387,
    0.009088,
    0.0,       // TODO
    0.0554545, // TODO: Matlab yields 0.554545 but it feels unreasonable
    0.0469321, // TODO: Matlab yields 0.469321 but it feels unreasonable
    0.016684,
    0.0, // TODO
    0.007673};

const std::vector<double> NORMAL_DEFAULT_CORRELATIONS = {
    1.0, 0.0004, 0.0, -0.0172, 0.0034, 0.0176, 0.0, 0.0155,
    0.0004, 1.0, 0.0, -0.0036, 0.0340, 0.0202, 0.0, 0.02201,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0172, -0.0036, 0.0, 1.0, 0.0032, 0.0151, 0.0, -0.0185,
    0.0034, 0.0340, 0.0, 0.0032, 1.0, -0.0185, 0.0, 0.0343,
    0.0176, 0.0202, 0.0, 0.0151, -0.0185, 1.0, 0.0, 0.0100,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0155, 0.201, 0.0, -0.0185, 0.0343, 0.01, 0.0, 1.0};

#endif