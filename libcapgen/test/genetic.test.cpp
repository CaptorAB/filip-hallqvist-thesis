#include <iostream>
#include <lib/catch.h>
#include <lib/random.h>

#include <include/genetic.h>

using Random = effolkronium::random_static;

TEST_CASE("normalize_individuals correctly normalizes individuals", "[genetic]")
{
  const int n_individuals = 3;
  const int n_instruments = 2;
  const int n_scenarios = 2;

  std::vector<double> individuals = {
      0.1, 0.3, 0.2, 0.0,
      0.0, 0.0, 1.0, 1.0,
      0.9, 0.1, 0.625, 0.975};

  normalize_individuals(individuals, n_individuals, n_instruments, n_scenarios);

  const std::vector<double> expected = {
      0.25, 0.75, 1.0, 0.0,
      0.5, 0.5, 0.5, 0.5,
      0.9, 0.1, 0.390625, 0.609375};

  for (int i = 0; i < individuals.size(); ++i)
  {
    REQUIRE(individuals[i] == Approx(expected[i]).epsilon(0.0001));
  }
}

TEST_CASE("initialize_individuals generates a population with correct capacity and values", "[genetic]")
{
  const int n_individuals = 1000;
  const int n_instruments = 10;
  const int n_scenarios = 20;

  std::vector<double> individuals = initialize_individuals(n_individuals, n_instruments, n_scenarios);

  int expected_size = n_individuals * n_instruments * n_scenarios;
  REQUIRE(individuals.size() == expected_size);

  for (int i = 0; i < individuals.size(); ++i)
  {
    REQUIRE((individuals[i] >= 0.0 && individuals[i] <= 1.0));
  }
}

TEST_CASE("select_roulette selects proper indices", "[genetic]")
{
  const int experiments = 10000;
  const int n_individuals = 10;

  double total = 0.0;
  for (int i = 0; i < experiments; ++i)
  {
    std::vector<double> fitnesses(n_individuals);
    for (int j = 0; j < n_individuals; ++j)
    {
      fitnesses[j] = Random::get<std::uniform_real_distribution<>>();
    }
    std::tuple<int, int> selected = select_roulette(fitnesses);
    total += fitnesses[std::get<0>(selected)];
    total += fitnesses[std::get<1>(selected)];
  }

  const double average = total / (2.0 * (double)experiments);

  REQUIRE(average == Approx(0.66).epsilon(0.01));
}

TEST_CASE("mutate_individuals correctly mutates individuals", "[genetic]")
{
  const double mutation_rate = 1.0;

  std::vector<double> selected = {
      0.1, 0.1, 0.1, 0.1,
      0.9, 0.9, 0.9, 0.9};

  std::vector<double> mutated(selected);

  mutate_individuals(mutated, mutation_rate);

  for (int i = 0; i < selected.size(); ++i)
  {
    REQUIRE(selected[i] != mutated[i]);
    REQUIRE(mutated[i] == Approx(selected[i]).margin(0.3));
  }
}

TEST_CASE("mutate_individuals does not mutate if mutation_rate is 0", "[genetic]")
{
  const double mutation_rate = 0.0;

  std::vector<double> selected = {
      0.1, 0.1, 0.1, 0.1,
      0.9, 0.9, 0.9, 0.9};

  std::vector<double> mutated(selected);

  mutate_individuals(mutated, mutation_rate);

  for (int i = 0; i < selected.size(); ++i)
  {
    REQUIRE(selected[i] == mutated[i]);
  }
}

TEST_CASE("crossover_chromosomes correctly crossovers chromosomes", "[genetic]")
{
  const double crossover_rate = 1.0;

  std::vector<double> selected = {
      0.1, 0.1, 0.1, 0.1,
      0.9, 0.9, 0.9, 0.9};

  std::vector<double> crossovered(selected);

  crossover_individuals(crossovered, crossover_rate);

  int crossovers = 0;
  for (int i = 0; i < selected.size(); ++i)
  {
    if (crossovered[i] != selected[i])
      crossovers++;
  }

  REQUIRE(crossovers > 0);
}

TEST_CASE("crossover_chromosomes does not perform crossover if crossover_rate is 0", "[genetic]")
{
  const double crossover_rate = 0.0;

  std::vector<double> selected = {
      0.1, 0.2, 0.3, 0.4,
      0.5, 0.6, 0.7, 0.8};

  std::vector<double> crossovered(selected);

  crossover_individuals(crossovered, crossover_rate);

  for (int i = 0; i < selected.size(); ++i)
  {
    REQUIRE(crossovered[i] == selected[i]);
  }
}

/*
TEST_CASE("evaluated individuals", "[genetic]")
{
  double result = 0.0;

  const int n_individuals = 1;
  const int n_instruments = 2;
  const int n_steps = 3;
  const int n_scenarios = 7;

  std::vector<double> individuals = {1.0, 0.0, 0.52, 0.48, 0.79, 0.21, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0};
  std::vector<double> price_changes = {0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.5};
  std::vector<double> probabilities = {1.0, 0.5, 0.5, 0.25, 0.25, 0.25, 0.25};
  std::vector<double> goals = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  double risk_aversion = 1.0;

  std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> evaluated = evaluate_individuals(individuals, price_changes, probabilities, goals, risk_aversion, n_individuals, n_steps, n_scenarios, n_instruments);

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

TEST_CASE("compute_wealth correctly computes wealth", "[genetic]")
{
  std::vector<double> current_weights = {0.1, 0.8, 0.1};
  std::vector<double> next_weights = {0.2, 0.8, 0.0};
  std::vector<double> price_changes = {0.5, -0.4, 0.3};
  std::vector<double> transaction_costs = {0.03, 0.01, 0.02};
  const double initial_wealth = 1.0;

  const double wealth = compute_wealth(current_weights, next_weights, price_changes, transaction_costs, initial_wealth);

  REQUIRE(wealth == Approx(0.7564).margin(0.0002));
}