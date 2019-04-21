#include <iostream>
#include <lib/catch.h>
#include <lib/random.h>

#include <include/constants.h>
#include <include/genetic.h>

using Random = effolkronium::random_static;

InstrumentConstraints create_default_instrument_constraints()
{
  InstrumentConstraints instrument_constraints;
  instrument_constraints.domestic_equity_min = 0.0;
  instrument_constraints.global_equity_min = 0.0;
  instrument_constraints.real_estate_min = 0.0;
  instrument_constraints.alternative_min = 0.0;
  instrument_constraints.credit_min = 0.0;
  instrument_constraints.bonds_2y_min = 0.0;
  instrument_constraints.bonds_5y_min = 0.0;
  instrument_constraints.cash_min = 0.0;
  instrument_constraints.fta_min = 0.0;
  instrument_constraints.domestic_equity_future_min = 0.0;
  instrument_constraints.interest_rate_swap_2y_min = 0.0;
  instrument_constraints.interest_rate_swap_5y_min = 0.0;
  instrument_constraints.interest_rate_swap_20y_min = 0.0;
  instrument_constraints.domestic_equity_max = 1.0;
  instrument_constraints.global_equity_max = 1.0;
  instrument_constraints.real_estate_max = 1.0;
  instrument_constraints.alternative_max = 1.0;
  instrument_constraints.credit_max = 1.0;
  instrument_constraints.bonds_2y_max = 1.0;
  instrument_constraints.bonds_5y_max = 1.0;
  instrument_constraints.cash_max = 1.0;
  instrument_constraints.fta_max = 1.0;
  instrument_constraints.domestic_equity_future_max = 1.0;
  instrument_constraints.interest_rate_swap_2y_max = 1.0;
  instrument_constraints.interest_rate_swap_5y_max = 1.0;
  instrument_constraints.interest_rate_swap_20y_max = 1.0;
  return instrument_constraints;
}

MarginConstraints create_default_margin_constraints()
{
  MarginConstraints margin_constraints;
  margin_constraints.domestic_equity_future = 0.0;
  margin_constraints.interest_rate_swap_2y = 0.0;
  margin_constraints.interest_rate_swap_5y = 0.0;
  margin_constraints.interest_rate_swap_20y = 0.0;
  return margin_constraints;
}

TEST_CASE("normalize_individuals correctly normalizes individuals", "[genetic]")
{
  const int n_individuals = 3;
  const int n_instruments = 2;
  const int n_derivatives = 0;
  const int n_scenarios = 2;

  std::vector<double> individuals = {
      0.1, 0.3, 0.2, 0.0,
      0.0, 0.0, 1.0, 1.0,
      0.9, 0.1, 0.625, 0.975};

  normalize_individuals(individuals, n_individuals, n_instruments, n_derivatives, n_scenarios);

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

TEST_CASE("compute_wealths correctly computes wealths over several steps", "[genetic]")
{
  const int n_instruments = 2;
  const int n_scenarios = 7;

  std::vector<double> individual = {0.5, 0.5, 0.55, 0.45, 0.6, 0.4, 0.55, 0.45, 0.5, 0.5, 0.45, 0.55, 0.4, 0.6};
  std::vector<double> price_changes = {0.2, -0.3, 0.15, -0.25, -0.15, 0.25, 0.5, -0.5, 0.35, -0.45, 0.05, 0.5, -0.5, -0.5};
  std::vector<double> transaction_costs = {0.03, 0.02};

  std::tuple<std::vector<double>, std::vector<double>> wealths = compute_wealths(individual, price_changes, transaction_costs, n_instruments, n_scenarios);

  std::vector<double> incoming_wealths = std::get<0>(wealths);
  std::vector<double> final_wealths = std::get<1>(wealths);

  std::vector<double> expected_incoming_wealths = {
      1.000000,
      0.947625, 0.945250,
      0.919196, 0.916898, 0.947542, 0.945155};

  std::vector<double> expected_final_wealths = {
      0.965156, 0.871053, 1.229436, 0.472578};

  for (int i = 0; i < expected_incoming_wealths.size(); ++i)
  {
    REQUIRE(incoming_wealths[i] == Approx(expected_incoming_wealths[i]).epsilon(0.0001));
  }

  for (int i = 0; i < expected_final_wealths.size(); ++i)
  {
    REQUIRE(final_wealths[i] == Approx(expected_final_wealths[i]).epsilon(0.0001));
  }
}

TEST_CASE("compute_expected_wealth correctly computes the expected wealth of an individual", "[genetic]")
{
  std::vector<double> final_wealths = {0.11, 0.21, 0.37, 0.49};
  const double wealth = compute_expected_wealth(final_wealths);
  const double expected = 0.295;

  REQUIRE(wealth == Approx(expected).epsilon(0.000001));
}

TEST_CASE("compute_expected_risk correctly computes the risk of an individual", "[genetic]")
{
  std::vector<double> incoming_wealths = {1.0, 0.7, 1.3};
  std::vector<double> goals = {1.0, 1.0, 1.0};
  const double risk = compute_expected_risk(incoming_wealths, goals);
  const double expected = 0.045;

  REQUIRE(risk == Approx(expected).epsilon(0.000001));
}

TEST_CASE("compute_fitnesses correctly computes the fitness of all individuals", "[genetic]")
{
  const int n_individuals = 3;
  const int n_instruments = 2;
  const int n_derivatives = 0;
  const int n_scenarios = 3;
  const double risk_aversion = 0;

  std::vector<double> individuals = {
      1.0, 0.0, 0.5, 0.5, 0.5, 0.5,
      0.0, 1.0, 0.5, 0.5, 0.5, 0.5,
      0.5, 0.5, 1.0, 0.0, 0.9, 0.1};
  std::vector<double> price_changes = {
      0.0, 0.5, 0.5, -0.5, -0.5, 0.5};
  std::vector<double> transaction_costs = {
      0.0, 0.0};
  std::vector<double> goals = {
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  std::vector<double> instrument_constraints = {
      0.0, 0.0, 1.0, 1.0};
  std::vector<double> margin_constraints = {
      0.0, 0.0, 0.0, 0.0};

  std::vector<double> fitnesses = compute_fitnesses(individuals, price_changes, transaction_costs, goals, instrument_constraints, margin_constraints, risk_aversion, n_individuals, n_instruments, n_derivatives, n_scenarios);

  std::vector<double> expected_fitnesses = {
      1.0, 1.5, 1.3125};

  for (int i = 0; i < expected_fitnesses.size(); ++i)
  {
    REQUIRE(fitnesses[i] == Approx(expected_fitnesses[i]).epsilon(0.0001));
  }
}

TEST_CASE("compute_penalty gives 0 penalty to valid individuals", "[genetic]")
{
  const int n_instruments = 2;
  const int n_derivatives = 0;
  const int n_scenarios = 3;
  std::vector<double> individual = {
      0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
  std::vector<double> instrument_constraints = {
      0.0, 0.0, 1.0, 1.0};
  std::vector<double> margin_constraints = {
      0.0, 0.0, 0.0, 0.0};
  const double penalty = compute_penalty(individual, instrument_constraints, margin_constraints, n_instruments, n_derivatives, n_scenarios);
  const double expected = 0.0;
  REQUIRE(penalty == Approx(expected).epsilon(0.000001));
}

TEST_CASE("compute_penalty penalizes individuals with an instrument weight below minimum", "[genetic]")
{
  const int n_instruments = 2;
  const int n_derivatives = 0;
  const int n_scenarios = 3;
  std::vector<double> individual = {
      0.8, 0.2, 0.8, 0.2, 0.1, 0.9};
  std::vector<double> instrument_constraints = {
      0.5, 0.0, 1.0, 1.0};
  std::vector<double> margin_constraints = {
      0.0, 0.0, 0.0, 0.0};
  const double penalty = compute_penalty(individual, instrument_constraints, margin_constraints, n_instruments, n_derivatives, n_scenarios);
  const double expected = 0.16;
  REQUIRE(penalty == Approx(expected).epsilon(0.000001));
}

TEST_CASE("compute_penalty penalizes individuals with an instrument weight above maximum", "[genetic]")
{
  const int n_instruments = 2;
  const int n_derivatives = 0;
  const int n_scenarios = 3;
  std::vector<double> individual = {
      0.4, 0.6, 0.7, 0.3, 0.1, 0.9};
  std::vector<double> instrument_constraints = {
      0.0, 0.0, 1.0, 0.8};
  std::vector<double> margin_constraints = {
      0.0, 0.0, 0.0, 0.0};
  const double penalty = compute_penalty(individual, instrument_constraints, margin_constraints, n_instruments, n_derivatives, n_scenarios);
  const double expected = 0.01;
  REQUIRE(penalty == Approx(expected).epsilon(0.000001));
}

TEST_CASE("optimization runs without crashing", "[genetic]")
{
  Random::seed(42);

  TransactionCosts transaction_costs;
  InstrumentConstraints instrument_constraints = create_default_instrument_constraints();
  MarginConstraints margin_constraints = create_default_margin_constraints();

  OptimizeOptions options;
  options.population_size = 10;
  options.elitism_copies = 2;
  options.generations = 10;
  options.steps = 3;
  options.mutation_rate = 0.02;
  options.crossover_rate = 0.01;
  options.risk_aversion = 0.0;
  options.initial_funding_ratio = 1.0;
  options.target_funding_ratio = 1.0;
  options.transaction_costs = transaction_costs;
  options.margin_constraints = margin_constraints;
  options.instrument_constraints = instrument_constraints;

  clock_t begin = std::clock();

  Result r = optimize(options);

  clock_t end = std::clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

  /*
  std::cout << "Elapsed time: " << elapsed_secs << "\n";

  std::cout << "\nFitness:         ";
  std::printf("%.4f \n", r.fitness);

  std::cout << "Expected return: ";
  std::printf("%.4f \n", r.expected_return);

  std::cout << "Expected risk:   ";
  std::printf("%.4f \n", r.expected_risk);

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
      std::cout << "== " << s;
      s = 0;
      std::cout << "\n";
    }
    i++;
  }
  std::cout << "\n";
  */
}