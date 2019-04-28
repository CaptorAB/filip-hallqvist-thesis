#ifndef __EMSCRIPTEN__
#define CATCH_CONFIG_RUNNER

#include <lib/catch.h>
#include <lib/random.h>

#endif

#include <include/constants.h>
#include <include/genetic.h>

#ifdef __EMSCRIPTEN__

#include <string>
#include <emscripten/bind.h>

using namespace emscripten;

EMSCRIPTEN_BINDINGS(libcapgen)
{
    constant("N_INSTRUMENTS", N_INSTRUMENTS);
    constant("N_DERIVATIVES", N_DERIVATIVES);

    constant("DOMESTIC_EQUITY_INDEX", DOMESTIC_EQUITY_INDEX);
    constant("GLOBAL_EQUITY_INDEX", GLOBAL_EQUITY_INDEX);
    constant("REAL_ESTATE_INDEX", REAL_ESTATE_INDEX);
    constant("ALTERNATIVE_INDEX", ALTERNATIVE_INDEX);
    constant("CREDIT_INDEX", CREDIT_INDEX);
    constant("BONDS_2Y_INDEX", BONDS_2Y_INDEX);
    constant("BONDS_5Y_INDEX", BONDS_5Y_INDEX);
    constant("CASH_INDEX", CASH_INDEX);
    constant("BONDS_20Y_INDEX", BONDS_20Y_INDEX);
    constant("DOMESTIC_EQUITY_FUTURE_INDEX", DOMESTIC_EQUITY_FUTURE_INDEX);
    constant("INTEREST_RATE_SWAP_2Y_INDEX", INTEREST_RATE_SWAP_2Y_INDEX);
    constant("INTEREST_RATE_SWAP_5Y_INDEX", INTEREST_RATE_SWAP_5Y_INDEX);
    constant("INTEREST_RATE_SWAP_20Y_INDEX", INTEREST_RATE_SWAP_20Y_INDEX);

    constant("N_RISKS", N_RISKS);
    constant("DOMESTIC_MARKET_RISK_INDEX", DOMESTIC_MARKET_RISK_INDEX);
    constant("GLOBAL_MARKET_RISK_INDEX", GLOBAL_MARKET_RISK_INDEX);
    constant("REAL_ESTATE_RISK_INDEX", REAL_ESTATE_RISK_INDEX);
    constant("ALTERNATIVE_RISK_INDEX", ALTERNATIVE_RISK_INDEX);
    constant("INTEREST_RATE_2Y_RISK_INDEX", INTEREST_RATE_2Y_RISK_INDEX);
    constant("INTEREST_RATE_5Y_RISK_INDEX", INTEREST_RATE_5Y_RISK_INDEX);
    constant("INTEREST_RATE_20Y_RISK_INDEX", INTEREST_RATE_20Y_RISK_INDEX);
    constant("CREDIT_RISK_INDEX", CREDIT_RISK_INDEX);

    constant("INSTRUMENT_NAMES", INSTRUMENT_NAMES);
    constant("RISK_NAMES", RISK_NAMES);

    value_array<TransactionCosts>("TransactionCosts")
        .element(&TransactionCosts::domestic_equity)
        .element(&TransactionCosts::global_equity)
        .element(&TransactionCosts::real_estate)
        .element(&TransactionCosts::alternative)
        .element(&TransactionCosts::credit)
        .element(&TransactionCosts::bonds_2y)
        .element(&TransactionCosts::bonds_5y)
        .element(&TransactionCosts::bonds_20y)
        .element(&TransactionCosts::cash)
        .element(&TransactionCosts::domestic_equity_future)
        .element(&TransactionCosts::interest_rate_swap_2y)
        .element(&TransactionCosts::interest_rate_swap_5y)
        .element(&TransactionCosts::interest_rate_swap_20y);

    value_array<InstrumentConstraints>("InstrumentConstraints")
        .element(&InstrumentConstraints::domestic_equity_min)
        .element(&InstrumentConstraints::global_equity_min)
        .element(&InstrumentConstraints::real_estate_min)
        .element(&InstrumentConstraints::alternative_min)
        .element(&InstrumentConstraints::credit_min)
        .element(&InstrumentConstraints::bonds_2y_min)
        .element(&InstrumentConstraints::bonds_5y_min)
        .element(&InstrumentConstraints::bonds_20y_min)
        .element(&InstrumentConstraints::cash_min)
        .element(&InstrumentConstraints::domestic_equity_future_min)
        .element(&InstrumentConstraints::interest_rate_swap_2y_min)
        .element(&InstrumentConstraints::interest_rate_swap_5y_min)
        .element(&InstrumentConstraints::interest_rate_swap_20y_min)
        .element(&InstrumentConstraints::domestic_equity_max)
        .element(&InstrumentConstraints::global_equity_max)
        .element(&InstrumentConstraints::real_estate_max)
        .element(&InstrumentConstraints::alternative_max)
        .element(&InstrumentConstraints::credit_max)
        .element(&InstrumentConstraints::bonds_2y_max)
        .element(&InstrumentConstraints::bonds_5y_max)
        .element(&InstrumentConstraints::bonds_20y_max)
        .element(&InstrumentConstraints::cash_max)
        .element(&InstrumentConstraints::domestic_equity_future_max)
        .element(&InstrumentConstraints::interest_rate_swap_2y_max)
        .element(&InstrumentConstraints::interest_rate_swap_5y_max)
        .element(&InstrumentConstraints::interest_rate_swap_20y_max);

    value_array<MarginConstraints>("MarginConstraints")
        .element(&MarginConstraints::domestic_equity_future)
        .element(&MarginConstraints::interest_rate_swap_2y)
        .element(&MarginConstraints::interest_rate_swap_5y)
        .element(&MarginConstraints::interest_rate_swap_20y);

    value_object<OptimizeOptions>("OptimizeOptions")
        .field("populationSize", &OptimizeOptions::population_size)
        .field("elitismCopies", &OptimizeOptions::elitism_copies)
        .field("generations", &OptimizeOptions::generations)
        .field("steps", &OptimizeOptions::steps)
        .field("mutationRate", &OptimizeOptions::mutation_rate)
        .field("crossoverRate", &OptimizeOptions::crossover_rate)
        .field("initialFundingRatio", &OptimizeOptions::initial_funding_ratio)
        .field("targetFundingRatio", &OptimizeOptions::target_funding_ratio)
        .field("transactionCosts", &OptimizeOptions::transaction_costs)
        .field("instrumentConstraints", &OptimizeOptions::instrument_constraints)
        .field("marginConstraints", &OptimizeOptions::margin_constraints);

    value_object<Result>("Result")
        .field("fitness", &Result::fitness)
        .field("individual", &Result::individual)
        .field("expectedReturn", &Result::expected_return)
        .field("expectedRisk", &Result::expected_risk)
        .field("intermediateWealths", &Result::intermediate_wealths)
        .field("finalWealths", &Result::final_wealths)
        .field("instrumentChanges", &Result::instrument_changes)
        .field("goals", &Result::goals);

    emscripten::register_vector<double>("vector<double>");
    emscripten::register_vector<std::string>("vector<string>");

    function("optimize", &optimize);
}

#else

using Random = effolkronium::random_static;
int main(int argc, char *argv[])
{
    // Random::seed(42);
    int result = Catch::Session().run(argc, argv);
    return result;
}

#endif