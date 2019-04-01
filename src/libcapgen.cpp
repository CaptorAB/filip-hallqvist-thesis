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
      .field("population", &OptimizeOptions::population)
      .field("elitism", &OptimizeOptions::elitism)
      .field("generations", &OptimizeOptions::generations)
      .field("bits", &OptimizeOptions::bits)
      .field("steps", &OptimizeOptions::steps);

  value_object<Result>("Result")
      .field("fitness", &Result::fitness)
      .field("chromosome", &Result::chromosome)
      .field("individual", &Result::individual);

  emscripten::register_vector<int>("VectorInt");
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

TEST_CASE("optimization runs without crashing", "[genetic]")
{
  OptimizeOptions options;
  options.population = 100;
  options.elitism = 2;
  options.generations = 100;
  options.bits = 7;
  options.steps = 4;
  options.mutation = 0.4;
  options.crossover = 0.02;

  Result r = optimize(options);

  /*
  std::cout << "\nFitness: \n";
  std::cout << r.fitness;
  std::cout << "\n Chromosome: \n";
  for (auto const &c : r.chromosome)
    std::cout << c << ' ';
  std::cout << "\n Individual: \n"
            << std::endl;
  for (auto const &c : r.individual)
    std::cout << c << ' ';
  */
}

using Random = effolkronium::random_static;
int main(int argc, char *argv[])
{
  Random::seed(42);
  int result = Catch::Session().run(argc, argv);
  return result;
}

#endif