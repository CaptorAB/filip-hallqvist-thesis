#include <lib/catch.h>
#include <lib/random.h>

#include <include/constants.h>
#include <include/genetic.h>

using Random = effolkronium::random_static;
using namespace std;

TransactionCosts create_default_transaction_costs()
{
  TransactionCosts transaction_costs;
  transaction_costs.domestic_equity = 0.0;
  transaction_costs.global_equity = 0.0;
  transaction_costs.real_estate = 0.0;
  transaction_costs.alternative = 0.0;
  transaction_costs.credit = 0.0;
  transaction_costs.bonds_2y = 0.0;
  transaction_costs.bonds_5y = 0.0;
  transaction_costs.bonds_20y = 0.0;
  transaction_costs.cash = 0.0;
  transaction_costs.domestic_equity_future = 0.0;
  transaction_costs.interest_rate_swap_2y = 0.0;
  transaction_costs.interest_rate_swap_5y = 0.0;
  transaction_costs.interest_rate_swap_20y = 0.0;
  return transaction_costs;
}

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
  instrument_constraints.bonds_20y_min = 0.0;
  instrument_constraints.cash_min = 0.0;
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
  instrument_constraints.bonds_20y_max = 1.0;
  instrument_constraints.cash_max = 1.0;
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

/*
TEST_CASE("normalize_individuals correctly normalizes individuals", "[genetic]")
{
  const int n_individuals = 3;
  const int n_instruments = 2;
  const int n_derivatives = 0;
  const int n_scenarios = 2;

  std::vector<double> individuals = {0.1, 0.3, 0.2, 0.0, 0.0, 0.0,
                                     1.0, 1.0, 0.9, 0.1, 0.625, 0.975};

  normalize_individuals(individuals, n_individuals, n_instruments,
                        n_derivatives, n_scenarios);

  const std::vector<double> expected = {
      0.25, 0.75, 1.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.9, 0.1, 0.390625, 0.609375};

  for (int i = 0; i < individuals.size(); ++i)
  {
    REQUIRE(individuals[i] == Approx(expected[i]).epsilon(0.0001));
  }
}

TEST_CASE("initialize_individuals generates a population with correct capacity "
          "and values",
          "[genetic]")
{
  const int n_individuals = 1000;
  const int n_instruments = 10;
  const int n_scenarios = 20;

  std::vector<double> individuals =
      initialize_individuals(n_individuals, n_instruments, n_scenarios);

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

  std::vector<double> selected = {0.1, 0.1, 0.1, 0.1, 0.9, 0.9, 0.9, 0.9};

  std::vector<double> mutated(selected);

  mutate_individuals(mutated, mutation_rate);

  for (int i = 0; i < selected.size(); ++i)
  {
    REQUIRE(selected[i] != mutated[i]);
  }
}

TEST_CASE("mutate_individuals does not mutate if mutation_rate is 0",
          "[genetic]")
{
  const double mutation_rate = 0.0;

  std::vector<double> selected = {0.1, 0.1, 0.1, 0.1, 0.9, 0.9, 0.9, 0.9};

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

  std::vector<double> selected = {0.1, 0.1, 0.1, 0.1, 0.9, 0.9, 0.9, 0.9};

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

TEST_CASE(
    "crossover_chromosomes does not perform crossover if crossover_rate is 0",
    "[genetic]")
{
  const double crossover_rate = 0.0;

  std::vector<double> selected = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

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
  std::vector<double> instrument_changes = {0.5, -0.4, 0.3};
  std::vector<double> transaction_costs = {0.03, 0.01, 0.02};
  const int n_instruments = 3;
  const int n_derivatives = 0;
  const double initial_wealth = 1.0;

  const double wealth = compute_wealth(
      current_weights, next_weights, instrument_changes, transaction_costs,
      initial_wealth, n_instruments, n_derivatives);

  REQUIRE(wealth == Approx(0.7564).margin(0.0002));
}

TEST_CASE("compute_wealths correctly computes wealths over several steps",
          "[genetic]")
{
  const int n_instruments = 2;
  const int n_derivatives = 0;
  const int n_scenarios = 7;

  std::vector<double> individual = {0.5, 0.5, 0.55, 0.45, 0.6, 0.4, 0.55,
                                    0.45, 0.5, 0.5, 0.45, 0.55, 0.4, 0.6};
  std::vector<double> instrument_changes = {0.2, -0.3, 0.15, -0.25, -0.15,
                                            0.25, 0.5, -0.5, 0.35, -0.45,
                                            0.05, 0.5, -0.5, -0.5};
  std::vector<double> transaction_costs = {0.03, 0.02};

  std::tuple<std::vector<double>, std::vector<double>> wealths =
      compute_wealths(individual, instrument_changes, transaction_costs,
                      n_instruments, n_derivatives, n_scenarios);

  std::vector<double> intermediate_wealths = std::get<0>(wealths);
  std::vector<double> final_wealths = std::get<1>(wealths);

  std::vector<double> expected_intermediate_wealths = {
      1.000000, 0.947625, 0.945250, 0.919196, 0.916898, 0.947542, 0.945155};

  std::vector<double> expected_final_wealths = {0.965156, 0.871053, 1.229436,
                                                0.472578};

  for (int i = 0; i < expected_intermediate_wealths.size(); ++i)
  {
    REQUIRE(intermediate_wealths[i] ==
            Approx(expected_intermediate_wealths[i]).epsilon(0.0001));
  }

  for (int i = 0; i < expected_final_wealths.size(); ++i)
  {
    REQUIRE(final_wealths[i] ==
            Approx(expected_final_wealths[i]).epsilon(0.0001));
  }
}

TEST_CASE("compute_expected_wealth correctly computes the expected wealth of "
          "an individual",
          "[genetic]")
{
  std::vector<double> final_wealths = {0.11, 0.21, 0.37, 0.49};
  const double wealth = compute_expected_wealth(final_wealths);
  const double expected = 0.295;

  REQUIRE(wealth == Approx(expected).epsilon(0.000001));
}

TEST_CASE("compute_fitnesses correctly computes the fitness of all individuals",
          "[genetic]")
{
  const int n_individuals = 3;
  const int n_instruments = 2;
  const int n_derivatives = 0;
  const int n_scenarios = 3;
  const int generation = 1;

  std::vector<double> individuals = {1.0, 0.0, 0.5, 0.5, 0.5, 0.5,
                                     0.0, 1.0, 0.5, 0.5, 0.5, 0.5,
                                     0.5, 0.5, 1.0, 0.0, 0.9, 0.1};
  std::vector<double> instrument_changes = {0.0, 0.5, 0.5, -0.5, -0.5, 0.5};
  std::vector<double> transaction_costs = {0.0, 0.0};
  std::vector<double> intermediate_goals = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  std::vector<double> final_goals = {0.0, 0.0};
  std::vector<double> instrument_constraints = {0.0, 0.0, 1.0, 1.0};
  std::vector<double> margin_constraints = {0.0, 0.0, 0.0, 0.0};

  std::vector<double> fitnesses = compute_fitnesses(
      individuals, instrument_changes, transaction_costs, intermediate_goals,
      final_goals, instrument_constraints, margin_constraints, n_individuals,
      n_instruments, n_derivatives, n_scenarios, generation);

  std::vector<double> expected_fitnesses = {1.0, 1.5, 1.3125};

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
  const int generation = 1;

  std::vector<double> individual = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
  std::vector<double> instrument_constraints = {0.0, 0.0, 1.0, 1.0};
  std::vector<double> margin_constraints = {0.0, 0.0, 0.0, 0.0};
  std::vector<double> intermediate_wealths = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  std::vector<double> final_wealths = {1.0, 1.0};
  std::vector<double> intermediate_goals = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  std::vector<double> final_goals = {0.0, 0.0};

  const double penalty =
      compute_penalty(individual, instrument_constraints, margin_constraints,
                      intermediate_wealths, final_wealths, intermediate_goals,
                      final_goals, n_instruments, n_derivatives, n_scenarios, generation);

  const double expected = 0.0;
  REQUIRE(penalty == Approx(expected).epsilon(0.000001));
}

TEST_CASE("compute_penalty penalizes individuals with an instrument weight "
          "below minimum",
          "[genetic]")
{
  const int n_instruments = 2;
  const int n_derivatives = 0;
  const int n_scenarios = 3;
  const int generation = 1;

  std::vector<double> individual = {0.8, 0.2, 0.8, 0.2, 0.1, 0.9};
  std::vector<double> instrument_constraints = {0.5, 0.0, 1.0, 1.0};
  std::vector<double> margin_constraints = {0.0, 0.0, 0.0, 0.0};
  std::vector<double> intermediate_wealths = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  std::vector<double> final_wealths = {1.0, 1.0};
  std::vector<double> intermediate_goals = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  std::vector<double> final_goals = {0.0, 0.0};
  const double penalty =
      compute_penalty(individual, instrument_constraints, margin_constraints,
                      intermediate_wealths, final_wealths, intermediate_goals,
                      final_goals, n_instruments, n_derivatives, n_scenarios, generation);

  const double expected = 0.16;
  REQUIRE(penalty == Approx(expected).epsilon(0.000001));
}

TEST_CASE("compute_penalty penalizes individuals with an instrument weight "
          "above maximum",
          "[genetic]")
{
  const int n_instruments = 2;
  const int n_derivatives = 0;
  const int n_scenarios = 3;
  const int generation = 1;

  std::vector<double> individual = {0.4, 0.6, 0.7, 0.3, 0.1, 0.9};
  std::vector<double> instrument_constraints = {0.0, 0.0, 1.0, 0.8};
  std::vector<double> margin_constraints = {0.0, 0.0, 0.0, 0.0};
  std::vector<double> intermediate_wealths = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  std::vector<double> final_wealths = {1.0, 1.0};
  std::vector<double> intermediate_goals = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  std::vector<double> final_goals = {0.0, 0.0};
  const double penalty =
      compute_penalty(individual, instrument_constraints, margin_constraints,
                      intermediate_wealths, final_wealths, intermediate_goals,
                      final_goals, n_instruments, n_derivatives, n_scenarios, generation);

  const double expected = 0.01;
  REQUIRE(penalty == Approx(expected).epsilon(0.000001));
}

*/
TEST_CASE("genetic golden master", "[genetic]")
{
  Random::seed(42);

  TransactionCosts transaction_costs =
      create_default_transaction_costs();
  InstrumentConstraints instrument_constraints =
      create_default_instrument_constraints();

  MarginConstraints margin_constraints;
  margin_constraints.domestic_equity_future = 0.0;
  margin_constraints.interest_rate_swap_2y = 0.0;
  margin_constraints.interest_rate_swap_5y = 0.0;
  margin_constraints.interest_rate_swap_20y = 0.0;

  OptimizeOptions options;
  options.population_size = 2;
  options.elitism_copies = 1;
  options.generations = 10000;
  options.steps = 3;
  options.mutation_rate = 0.5;
  options.crossover_rate = 0.0;
  options.initial_funding_ratio = 1.0;
  options.target_funding_ratio = 1.0;
  options.transaction_costs = transaction_costs;
  options.margin_constraints = margin_constraints;
  options.instrument_constraints = instrument_constraints;

  printf("Running...\n");
  const int n_scenarios = 7;
  const int n_trials = 10;
  vector<double> results(10 * N_INSTRUMENTS * n_scenarios);

  for (int i = 0; i < n_trials; ++i)
  {
    Result r = optimize(options);
    printf("(%i) Return: %.2f%%\n", i, 100 * r.expected_return);

    for (int j = 0; j < r.individual.size(); ++j)
    {
      results[i * r.individual.size() + j] = r.individual[j];
    }

    for (int j = 0; j < n_scenarios; ++j) // Scenarios
    {
      for (int k = 0; k < N_INSTRUMENTS; ++k) // Instruments
      {
        printf("%.2f ", r.individual[j * N_INSTRUMENTS + k]);
      }
      printf("\n");
    }
  }

  // Compute mean and variance
  vector<double> summed(N_INSTRUMENTS * n_scenarios);
  vector<double> squared(N_INSTRUMENTS * n_scenarios);
  vector<double> means(N_INSTRUMENTS * n_scenarios);
  vector<double> variances(N_INSTRUMENTS * n_scenarios);

  for (int i = 0; i < n_scenarios; ++i)
  {
    for (int j = 0; j < N_INSTRUMENTS; ++j)
    {
      for (int k = 0; k < n_trials; ++k)
      {
        const double x = results[k * (N_INSTRUMENTS * n_scenarios) + (i * N_INSTRUMENTS) + j];
        summed[i * N_INSTRUMENTS + j] += x;
        squared[i * N_INSTRUMENTS + j] += x * x;
      }
    }
  }

  printf("\n");
  for (int i = 0; i < summed.size(); ++i)
  {
    means[i] = summed[i] / (double)n_trials;
    variances[i] = (squared[i] - ((summed[i] * summed[i]) / (double)n_trials)) / ((double)n_trials - 1.0);
  }

  for (int j = 0; j < n_scenarios; ++j) // Scenarios
  {
    for (int k = 0; k < N_INSTRUMENTS; ++k) // Instruments
    {
      printf("%.2f (±%.2f)  ", means[j * N_INSTRUMENTS + k], sqrt(variances[j * N_INSTRUMENTS + k]));
    }
    printf("\n");
  }
}

/*
TEST_CASE("intermediate_wealths")
{
  const int n_instruments = 13;
  const int n_derivatives = 4;
  const int n_scenarios = 3;
  const int n_risks = N_RISKS;

  std::vector<double> individual = {
      0.02, 0.32, 0.00, 0.00, 0.09, 0.18, 0.20, 0.18, 0.00, 1.00, 1.00, 1.00, 1.00, 0.10, 0.19, 0.00, 0.00, 0.25, 0.38, 0.07, 0.00, 0.00, 1.00, 1.00, 1.00, 1.00, 0.02, 0.10, 0.00, 0.00, 0.35, 0.07, 0.15, 0.32, 0.00, 1.00, 1.00, 1.00, 1.00};

  std::vector<double> means = NORMAL_DEFAULT_MEANS;
  std::vector<double> standard_deviations = NORMAL_DEFAULT_STANDARD_DEVIATIONS;
  std::vector<double> correlations = NORMAL_DEFAULT_CORRELATIONS;

  std::vector<double> instrument_changes =
      generate_normal_scenarios(means, standard_deviations, correlations,
                                n_risks, n_instruments, n_scenarios);
  std::vector<double> transaction_costs = std::vector<double>(n_instruments, 0.0);

  std::tuple<std::vector<double>, std::vector<double>> wealths =
      compute_wealths(individual, instrument_changes, transaction_costs,
                      n_instruments, n_derivatives, n_scenarios);

  std::vector<double> intermediate_wealths = std::get<0>(wealths);
  std::vector<double> final_wealths = std::get<1>(wealths);

  printf("Intermediate\n");
  for (auto w : intermediate_wealths)
    printf("%.2f ", w);
  printf("\n");

  printf("Final\n");
  for (auto w : final_wealths)
    printf("%.2f ", w);
  printf("\n");
}
*/