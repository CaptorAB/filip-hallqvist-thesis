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
  value_object<OptimizeOptions>("OptimizeOptions")
      .field("populationSize", &OptimizeOptions::population_size)
      .field("elitismCopies", &OptimizeOptions::elitism_copies)
      .field("generations", &OptimizeOptions::generations)
      .field("steps", &OptimizeOptions::steps)
      .field("mutationRate", &OptimizeOptions::mutation_rate)
      .field("crossoverRate", &OptimizeOptions::crossover_rate)
      .field("riskAversion", &OptimizeOptions::risk_aversion)
      .field("initialFundingRatio", &OptimizeOptions::initial_funding_ratio)
      .field("targetFundingRatio", &OptimizeOptions::target_funding_ratio);

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