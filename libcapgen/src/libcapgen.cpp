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
      .field("individual", &Result::individual);

  emscripten::register_vector<double>("VectorDouble");

  function("optimize", &optimize);
}

#else

/*
TEST_CASE("sample_black_process calculates correct values", "[sample_black_process]")
{
  double n1 = 0.12;
  double n2 = 0.33;
  double forward_price = 0.91;
  double gamma = 0.42;
  double rho = 0.9;
  double sigma = 0.78;
  double mean = 1.0;
  double variance = 0.99;

  double result = sample_black_process(n1, n2, forward_price, gamma, rho, sigma, mean, variance);
  REQUIRE(result - 0.267704205 <= 0.0001); // TODO: Verify this using pen and paper
}

TEST_CASE("DomesticMarketRiskProcess calculates correct values", "[risk]")
{
  double n1 = 0.25;
  double n2 = 0.9;
  RiskProcess *risk_process = new DomesticMarketRiskProcess();

  double current = risk_process->get_current();

  REQUIRE(current == 0.0);

  risk_process->update(n1, n2);

  current = risk_process->get_current();

  REQUIRE(current - 0.3456077342 <= 0.0001); // TODO: Verify this using pen and paper
}

TEST_CASE("DomesticEquity moves correctly", "[risk]")
{
  double n1 = 0.25;
  double n2 = 0.9;
  risks_t risks = create_default_risks();
  Instrument *instrument = new DomesticEquityInstrument();

  double current = instrument->get_current();

  REQUIRE(current == 1.0);

  instrument->update(risks);

  current = instrument->get_current();

  REQUIRE(current - 0.3456077342 <= 0.0001); // TODO: Verify this using pen and paper
}
*/
/*
TEST_CASE("initialize_chromosomes generates a population with correct capacity and values", "[genetic]")
{
  int population = 2;
  int genes = 10;
  int size = population * genes;

  std::vector<int> chromosomes = initialize_chromosomes(population, genes);

  REQUIRE(chromosomes.capacity() == size);
  for (size_t i = 0; i < size; i++)
  {
    REQUIRE((chromosomes[i] == 1 || chromosomes[i] == 0));
  }
}

TEST_CASE("select_roulette selects proper indices", "[genetic]")
{
  std::vector<double> fitnesses({0.01, 0.02, 1.25, 0.02, 0.3, 1.1, 0.9});
  std::tuple<size_t, size_t> selected = select_roulette(fitnesses);

  REQUIRE(std::get<0>(selected) == 2);
  REQUIRE(std::get<1>(selected) == 5);
}

TEST_CASE("crossover_chromosomes correctly crossovers chromosomes", "[genetic]")
{
  std::vector<int> chromosomes({0, 0, 0, 0, 1, 1, 1, 0});
  std::vector<double> evaluated({0.0, 2.0});
  std::vector<int> crossovered = crossover_chromosomes(chromosomes, evaluated, 2, 4, 1.0);

  REQUIRE(crossovered[0] == 1);
  REQUIRE(crossovered[1] == 1);
  REQUIRE(crossovered[2] == 1);
  REQUIRE(crossovered[3] == 0);
  REQUIRE(crossovered[4] == 1);
  REQUIRE(crossovered[5] == 1);
  REQUIRE(crossovered[6] == 1);
  REQUIRE(crossovered[7] == 0);
}
*/

/*
TEST_CASE("evaluated individuals", "[genetic]")
{
  double result = 0.0;
  std::vector<double> individuals = {1.0, 0.0, 0.52, 0.48, 0.79, 0.21, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0};
  std::vector<double> price_changes = {0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.5};
  std::vector<double> probabilities = {1.0, 0.5, 0.5, 0.25, 0.25, 0.25, 0.25};
  std::vector<double> goals = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  double risk_aversion = 1.0;

  std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> evaluated = evaluate_individuals(individuals, price_changes, probabilities, goals, risk_aversion, 1, 3, 7, 2);

  std::vector<double> fitnesses = std::get<0>(evaluated);
  std::vector<double> total_returns = std::get<1>(evaluated);
  std::vector<double> risks = std::get<2>(evaluated);

  REQUIRE(fitnesses.size() == 1);
  REQUIRE(total_returns.size() == 1);
  REQUIRE(risks.size() == 1);
  REQUIRE(total_returns[0] == 1.73);

  printf("Fitnesses: %f\n", fitnesses[0]);
  printf("Total returns: %f\n", total_returns[0]);
  printf("Risks: %f\n", risks[0]);
  std::cout << std::endl;
}
*/

TEST_CASE("optimization runs without crashing", "[genetic]")
{
  OptimizeOptions options;
  options.population_size = 20;
  options.elitism_copies = 2;
  options.generations = 400;
  options.steps = 8;
  options.mutation_rate = 0.05;
  options.crossover_rate = 0.05;
  options.risk_aversion = 0.0;
  options.initial_funding_ratio = 1.0;
  options.target_funding_ratio = 1.0;

  Result r = optimize(options);

  std::cout << "\nFitness:           ";
  std::printf("%.4f \n", r.fitness);

  std::cout << "Tot. Exp. Return:  ";
  std::printf("%.4f \n", r.total_return);

  std::cout << "Acc. Semivariance: ";
  std::printf("%.4f \n", r.risk);

  std::cout << "Individual: \n"
            << std::endl;

  int i = 1;
  double s = 0.0;
  for (auto const &c : r.individual)
  {
    s += c;
    std::printf("%.2f ", c);
    if (i % N_INSTRUMENTS == 0)
    {
      std::cout << " == " << s;
      s = 0;
      std::cout << "\n";
    }
    i++;
  }
  std::cout << "\n";
}

using Random = effolkronium::random_static;
int main(int argc, char *argv[])
{
  Random::seed(42);
  int result = Catch::Session().run(argc, argv);
  return result;
}

#endif