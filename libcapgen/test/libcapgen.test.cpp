#include <ctime>
#include <iostream>

#include <lib/catch.h>

#include <include/scenario.h>
#include <include/genetic.h>

TEST_CASE("optimization runs without crashing", "[genetic]")
{
  OptimizeOptions options;
  options.population_size = 100;
  options.elitism_copies = 5;
  options.generations = 10000;
  options.steps = 3;
  options.mutation_rate = 0.20;
  options.crossover_rate = 0.5;
  options.risk_aversion = 0.0;
  options.initial_funding_ratio = 1.0;
  options.target_funding_ratio = 1.0;

  clock_t begin = std::clock();

  Result r = optimize(options);

  clock_t end = std::clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

  std::cout << "Elapsed time: " << elapsed_secs << "\n";

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