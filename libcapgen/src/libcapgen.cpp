#include <string>
#include <iostream>

#ifndef __EMSCRIPTEN__
#define CATCH_CONFIG_RUNNER
#include <lib/catch.h>
#endif

#include <vector>
#include <lib/random.h>

#include <include/genetic.h>
#include <include/scenario.h>

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
using namespace emscripten;

EMSCRIPTEN_BINDINGS(libcapgen)
{
    value_array<TransactionCosts>("TransactionCosts")
        .element(&TransactionCosts::domestic_equity)
        .element(&TransactionCosts::global_equity)
        .element(&TransactionCosts::real_estate)
        .element(&TransactionCosts::alternative)
        .element(&TransactionCosts::credit)
        .element(&TransactionCosts::bonds_2y)
        .element(&TransactionCosts::bonds_5y)
        .element(&TransactionCosts::cash)
        .element(&TransactionCosts::fta)
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
        .element(&InstrumentConstraints::cash_min)
        .element(&InstrumentConstraints::fta_min)
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
        .element(&InstrumentConstraints::cash_max)
        .element(&InstrumentConstraints::fta_max)
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
        .field("riskAversion", &OptimizeOptions::risk_aversion)
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
        .field("incomingWealths", &Result::incoming_wealths)
        .field("finalWealths", &Result::final_wealths)
        .field("priceChanges", &Result::price_changes)
        .field("goals", &Result::goals);

    emscripten::register_vector<double>("VectorDouble");

    function("optimize", &optimize);
}

#else

using Random = effolkronium::random_static;
int main(int argc, char *argv[])
{
    Random::seed(42);
    int result = Catch::Session().run(argc, argv);
    return result;
}

#endif