#include <vector>
#include <tuple>
#include <cmath>

#include <lib/random.h>

#include <include/constants.h>
#include <include/util.h>
#include <include/genetic.h>
#include <include/scenario.h>

using Random = effolkronium::random_static;

std::vector<double> parse_instrument_constraints(
    InstrumentConstraints instrument_constraints)
{
  return std::vector<double>({instrument_constraints.domestic_equity_min,
                              instrument_constraints.global_equity_min,
                              instrument_constraints.real_estate_min,
                              instrument_constraints.alternative_min,
                              instrument_constraints.credit_min,
                              instrument_constraints.bonds_2y_min,
                              instrument_constraints.bonds_5y_min,
                              instrument_constraints.cash_min,
                              instrument_constraints.fta_min,
                              instrument_constraints.domestic_equity_future_min,
                              instrument_constraints.interest_rate_swap_2y_min,
                              instrument_constraints.interest_rate_swap_5y_min,
                              instrument_constraints.interest_rate_swap_20y_min,
                              instrument_constraints.domestic_equity_max,
                              instrument_constraints.global_equity_max,
                              instrument_constraints.real_estate_max,
                              instrument_constraints.alternative_max,
                              instrument_constraints.credit_max,
                              instrument_constraints.bonds_2y_max,
                              instrument_constraints.bonds_5y_max,
                              instrument_constraints.cash_max,
                              instrument_constraints.fta_max,
                              instrument_constraints.domestic_equity_future_max,
                              instrument_constraints.interest_rate_swap_2y_max,
                              instrument_constraints.interest_rate_swap_5y_max,
                              instrument_constraints.interest_rate_swap_20y_max});
}

std::vector<double> parse_margin_constraints(
    MarginConstraints margin_constraints)
{
  return std::vector<double>({margin_constraints.domestic_equity_future,
                              margin_constraints.interest_rate_swap_2y,
                              margin_constraints.interest_rate_swap_5y,
                              margin_constraints.interest_rate_swap_20y});
}

void normalize_individuals(
    std::vector<double> &individuals,
    const int n_individuals,
    const int n_instruments,
    const int n_derivatives,
    const int n_scenarios)
{
  const int n_non_derivatives = n_instruments - n_derivatives;
  const int n_genes = n_instruments * n_scenarios;

  for (int i = 0; i < n_individuals; ++i)
  {
    int ix = i * n_genes;
    for (int s = 0; s < n_scenarios; ++s)
    {
      double total = 0.0;
      int sx = ix + (s * n_instruments);

      for (int j = 0; j < n_non_derivatives; ++j)
      {
        int jx = sx + j;
        total += individuals[jx];
      }

      if (total == 0.0)
      {
        for (int j = 0; j < n_non_derivatives; ++j)
        {
          int jx = sx + j;
          individuals[jx] = 1.0 / n_instruments;
        }
      }
      else
      {
        for (int j = 0; j < n_non_derivatives; ++j)
        {
          int jx = sx + j;
          individuals[jx] = individuals[jx] / total;
        }
      }

    } // Scenarios

  } // Individuals
}

std::vector<double> initialize_individuals(
    const int n_individuals,
    const int n_instruments,
    const int n_scenarios)
{
  const int n_genes = n_instruments * n_scenarios;
  const int size = n_individuals * n_genes;

  std::vector<double> individuals(size);

  // Generate genes
  int ix = 0;
  for (int i = 0; i < n_individuals; i++)
  {
    for (int j = 0; j < (n_instruments * n_scenarios); j++)
    {
      ix = (i * n_genes) + j;
      individuals[ix] = Random::get<std::uniform_real_distribution<>>();
    }
  }
  return individuals;
}

std::tuple<int, int> select_roulette(std::vector<double> &fitnesses)
{
  {
    int i1 = Random::get<std::discrete_distribution<>>(
        fitnesses.begin(), fitnesses.end());
    int i2 = Random::get<std::discrete_distribution<>>(
        fitnesses.begin(), fitnesses.end());
    return std::make_tuple(i1, i2);
  }
}

double compute_wealth(
    std::vector<double> &current_weights,
    std::vector<double> &next_weights,
    std::vector<double> &price_changes,
    std::vector<double> &transaction_costs,
    const double initial_wealth)
{
  double holdings = 0.0;
  double reallocations = 0.0;

  // Compute wealth from price changes
  for (int i = 0; i < current_weights.size(); ++i)
  {
    const double result =
        initial_wealth * current_weights[i] * (1.0 + price_changes[i]);
    holdings += result;
  }

  // Deduct transation costs
  for (int i = 0; i < current_weights.size(); ++i)
  {
    const double diff = fabs(current_weights[i] - next_weights[i]);
    reallocations += holdings * diff * transaction_costs[i];
  }

  const double wealth = holdings - reallocations;

  return wealth;
}

std::tuple<std::vector<double>, std::vector<double>> compute_wealths(
    std::vector<double> &individual,
    std::vector<double> &price_changes,
    std::vector<double> &transaction_costs,
    const int n_instruments,
    const int n_scenarios)
{
  std::vector<double> incoming_wealths(n_scenarios);
  incoming_wealths[0] = 1.0;

  std::vector<double> final_wealths(n_scenarios / 2 + 1);
  int final_index = 0;

  // Iterate through scenarios
  for (int i = 0; i < n_scenarios; ++i)
  {
    int current = i;

    int left = 2 * current + 1;
    int right = 2 * current + 2;

    // Index of current node
    const int cx = current * n_instruments;

    const double current_wealth = incoming_wealths[current];

    std::vector<double> current_changes = std::vector<double>(
        price_changes.begin() + cx,
        price_changes.begin() + cx + n_instruments);

    std::vector<double> current_weights = std::vector<double>(
        individual.begin() + cx,
        individual.begin() + cx + n_instruments);

    if (left < n_scenarios && right < n_scenarios)
    {
      const int lx = left * n_instruments;
      const int rx = right * n_instruments;

      std::vector<double> left_weights = std::vector<double>(
          individual.begin() + lx,
          individual.begin() + lx + n_instruments);

      std::vector<double> right_weights = std::vector<double>(
          individual.begin() + rx,
          individual.begin() + rx + n_instruments);

      // Evaluate left child
      incoming_wealths[left] = compute_wealth(
          current_weights,
          left_weights,
          current_changes,
          transaction_costs,
          current_wealth);

      // Evaluate right child
      incoming_wealths[right] = compute_wealth(
          current_weights,
          right_weights,
          current_changes,
          transaction_costs,
          current_wealth);
    }
    else
    {
      // We are at the leaf nodes of the scenario tree,
      // so compute the final wealth.

      final_wealths[final_index] = compute_wealth(
          current_weights,
          current_weights,
          current_changes,
          transaction_costs,
          current_wealth);

      final_index++;
    }
  }

  return std::make_tuple(incoming_wealths, final_wealths);
}

double compute_expected_wealth(std::vector<double> &final_wealths)
{
  double expected_wealth = 0.0;
  for (double wealth : final_wealths)
  {
    expected_wealth += wealth;
  }
  expected_wealth = expected_wealth / final_wealths.size();

  return expected_wealth;
}

double compute_expected_risk(
    std::vector<double> &incoming_wealths,
    std::vector<double> &goals)
{
  double risk = 0.0;
  for (int i = 0; i < incoming_wealths.size(); ++i)
  {
    const double step = floor(std::log2(i + 1));
    const double nodes_in_step = pow(2.0, step);

    risk += pow(std::min(0.0, incoming_wealths[i] - goals[i]), 2.0) / nodes_in_step;
  }
  return risk;
}

double compute_fitness(
    std::vector<double> &individual,
    std::vector<double> &incoming_wealths,
    std::vector<double> &final_wealths,
    std::vector<double> &goals,
    std::vector<double> &instrument_constraints,
    std::vector<double> &margin_constraints,
    const double risk_aversion,
    const int n_instruments,
    const int n_derivatives,
    const int n_scenarios)
{
  const double wealth = compute_expected_wealth(final_wealths);
  const double risk = compute_expected_risk(incoming_wealths, goals);
  const double penalty = compute_penalty(
      individual,
      instrument_constraints,
      margin_constraints,
      n_instruments,
      n_derivatives,
      n_scenarios);

  return ((1.0 - risk_aversion) * wealth) - (risk_aversion * risk) - penalty;
}

double compute_penalty(
    std::vector<double> &individual,
    std::vector<double> &instrument_constraints,
    std::vector<double> &margin_constraints,
    const int n_instruments,
    const int n_derivatives,
    const int n_scenarios)
{
  const int n_non_derivatives = n_instruments - n_derivatives;
  double penalty = 0.0;

  // Check instrument allocation constraints
  for (int j = 0; j < n_scenarios; ++j)
  {
    const int jx = j * n_instruments;

    // Allocation constraints
    for (int k = 0; k < n_instruments; ++k)
    {
      // Minimum
      penalty += pow(
          std::max(0.0, instrument_constraints[k] - individual[jx + k]),
          2.0);

      // Maximum
      penalty += pow(
          std::max(0.0, individual[jx + k] - instrument_constraints[k + n_instruments]),
          2.0);
    }

    // Margin constraints
    double remaining_margin =
        individual[jx + CASH_INDEX] + individual[jx + FTA_INDEX];

    for (int k = 0; k < n_derivatives; ++k)
    {
      const int kx = n_non_derivatives + k;
      const double required_margin = margin_constraints[k] * individual[jx + kx];
      remaining_margin -= required_margin;
    }
    penalty += pow(std::min(0.0, remaining_margin), 2.0);
  }

  return penalty;
}

std::vector<double> compute_fitnesses(
    std::vector<double> &individuals,
    std::vector<double> &price_changes,
    std::vector<double> &transaction_costs,
    std::vector<double> &goals,
    std::vector<double> &instrument_constraints,
    std::vector<double> &margin_constraints,
    const double risk_aversion,
    const int n_individuals,
    const int n_instruments,
    const int n_derivatives,
    const int n_scenarios)
{
  const int n_genes = n_instruments * n_scenarios;
  std::vector<double> fitnesses(n_individuals);
  for (int i = 0; i < n_individuals; ++i)
  {
    const int ix = i * n_genes;

    std::vector<double> individual = std::vector<double>(
        individuals.begin() + ix,
        individuals.begin() + ix + n_genes);

    std::tuple<std::vector<double>, std::vector<double>> wealths =
        compute_wealths(
            individual,
            price_changes,
            transaction_costs,
            n_instruments,
            n_scenarios);

    std::vector<double> incoming_wealths = std::get<0>(wealths);
    std::vector<double> final_wealths = std::get<1>(wealths);

    fitnesses[i] = compute_fitness(
        individual,
        incoming_wealths,
        final_wealths,
        goals,
        instrument_constraints,
        margin_constraints,
        risk_aversion,
        n_instruments,
        n_derivatives,
        n_scenarios);
  }
  return fitnesses;
}

void mutate_individuals(std::vector<double> &selected, const double mutation_rate)
{
  for (int j = 0; j < selected.size(); ++j)
  {
    const double random = Random::get(0.0, 1.0);
    if (random < mutation_rate)
    {
      const double d = ((double)j + 1.0) / ((double)j + 3.0);
      const double mutation = d * Random::get(-1.0, 1.0);

      selected[j] = std::min(1.0, std::max(0.0, selected[j] + mutation));
    }
  }
}

void crossover_individuals_scenario(
    std::vector<double> &selected,
    const int n_instruments,
    const int n_scenarios,
    const double crossover_rate)
{
  const int n_genes = n_instruments * n_scenarios;
  const double r = Random::get(0.0, 1.0);
  if (r < crossover_rate)
  {
    const int scenario = Random::get(0, n_scenarios - 1);
    const int ix = scenario * n_instruments;
    for (int j = 0; j < n_genes; ++j)
    {
      const double temp = selected[j];
      selected[j] = selected[j + n_genes];
      selected[j + n_genes] = temp;
    }
  }
}

void crossover_individuals(
    std::vector<double> &selected,
    const double crossover_rate)
{

  const double r1 = Random::get(0.0, 1.0);
  if (r1 < crossover_rate)
  {
    const int n_genes = selected.size() / 2;
    std::vector<double> temp(selected);

    // Construct first individual
    for (int j = 0; j < n_genes; ++j)
    {
      const double r2 = Random::get(0.0, 1.0);
      if (r2 < 0.5)
        selected[j] = temp[j];
      else
        selected[j] = temp[j + n_genes];
    }

    // Construct Second individual
    for (int j = 0; j < n_genes; ++j)
    {
      const double r2 = Random::get(0.0, 1.0);
      if (r2 < 0.5)
        selected[j + n_genes] = temp[j];
      else
        selected[j + n_genes] = temp[j + n_genes];
    }
  }
}

Result optimize(OptimizeOptions options)
{
  const int n_individuals = options.population_size;
  const int n_elitism_copies = options.elitism_copies;
  const int n_generations = options.generations;
  const int n_steps = options.steps;
  const double mutation_rate = options.mutation_rate;
  const double crossover_rate = options.crossover_rate;
  const double risk_aversion = options.risk_aversion;
  const double initial_funding_ratio = options.initial_funding_ratio;
  const double target_funding_ratio = options.target_funding_ratio;

  std::vector<double> transaction_costs = {
      options.transaction_costs.domestic_equity,
      options.transaction_costs.global_equity,
      options.transaction_costs.real_estate,
      options.transaction_costs.alternative,
      options.transaction_costs.credit,
      options.transaction_costs.bonds_2y,
      options.transaction_costs.bonds_5y,
      options.transaction_costs.cash,
      options.transaction_costs.fta,
      options.transaction_costs.domestic_equity_future,
      options.transaction_costs.interest_rate_swap_2y,
      options.transaction_costs.interest_rate_swap_5y,
      options.transaction_costs.interest_rate_swap_20y};

  std::vector<double> instrument_constraints =
      parse_instrument_constraints(options.instrument_constraints);
  std::vector<double> margin_constraints =
      parse_margin_constraints(options.margin_constraints);

  const int n_instruments = N_INSTRUMENTS;
  const int n_derivatives = N_DERIVATIVES;

  const int n_scenarios = pow(2.0, n_steps) - 1;
  const int n_genes = n_scenarios * n_instruments;

  // Will contain the final result
  double best_fitness = std::numeric_limits<double>::min();
  std::vector<double> best_individual = std::vector<double>(n_genes, 0.0);

  // Generate initial individuals
  std::vector<double> individuals =
      initialize_individuals(n_individuals, n_instruments, n_scenarios);

  // ...and normalize them.
  normalize_individuals(
      individuals,
      n_individuals,
      n_instruments,
      n_derivatives,
      n_scenarios);

  // Generate scenarios
  std::tuple<std::vector<double>, std::vector<double>> scenarios =
      generate_scenarios(n_steps);

  std::vector<double> price_changes = std::get<0>(scenarios);
  std::vector<double> probabilities = std::get<1>(scenarios);

  // Generate goals
  std::vector<double> goals = generate_goals(
      price_changes,
      n_steps,
      n_scenarios,
      n_instruments,
      initial_funding_ratio,
      target_funding_ratio);

  for (int t = 0; t < n_generations; ++t)
  {
    std::vector<double> fitnesses = compute_fitnesses(
        individuals,
        price_changes,
        transaction_costs,
        goals,
        instrument_constraints,
        margin_constraints,
        risk_aversion,
        n_individuals,
        n_instruments,
        n_derivatives,
        n_scenarios);

    // Check global best fitness
    for (int i = 0; i < n_individuals; ++i)
    {
      if (fitnesses[i] > best_fitness)
      {
        const int ix = i * n_genes;
        best_fitness = fitnesses[i];
        best_individual = std::vector<double>(
            individuals.begin() + ix,
            individuals.begin() + ix + n_genes);
      }
    }

    // Clone individuals vector
    std::vector<double> offspring(individuals);

    for (int i = 0; i < n_individuals; i += 2)
    {
      const int ix = i * n_genes;

      // Selection
      std::tuple<int, int> indices = select_roulette(fitnesses);
      int ix1 = std::get<0>(indices) * n_genes;
      int ix2 = std::get<1>(indices) * n_genes;

      std::vector<double> selected(2 * n_genes);
      for (int j = 0; j < n_genes; ++j)
      {
        selected[j] = individuals[ix1 + j];
        selected[n_genes + j] = individuals[ix2 + j];
      }

      // Crossover
      crossover_individuals_scenario(
          selected, n_instruments,
          n_scenarios,
          crossover_rate);

      // Mutation
      mutate_individuals(selected, mutation_rate);

      // Add selected individuals
      for (int j = 0; j < n_genes; ++j)
      {
        offspring[i * n_genes] = selected[j];
      }
    }

    // Elitism
    for (int i = 0; i < n_elitism_copies; ++i)
    {
      int ix = i * n_genes;
      for (int j = 0; j < n_genes; ++j)
      {
        offspring[ix + j] = best_individual[j];
      }
    }

    normalize_individuals(
        offspring,
        n_individuals,
        n_instruments,
        n_derivatives,
        n_scenarios);

    // Replace generation
    individuals = offspring;
  }

  Result result;

  std::tuple<std::vector<double>, std::vector<double>> best_wealths =
      compute_wealths(
          best_individual,
          price_changes,
          transaction_costs,
          n_instruments,
          n_scenarios);

  std::vector<double> best_incoming_wealths = std::get<0>(best_wealths);
  std::vector<double> best_final_wealths = std::get<1>(best_wealths);
  double best_expected_return = compute_expected_wealth(best_final_wealths) - 1.0;
  double best_expected_risk = compute_expected_risk(best_incoming_wealths, goals);

  result.fitness = best_fitness;
  result.individual = best_individual;
  result.incoming_wealths = best_incoming_wealths;
  result.final_wealths = best_final_wealths;
  result.expected_return = best_expected_return;
  result.expected_risk = best_expected_risk;
  result.price_changes = price_changes;
  result.goals = goals;

  return result;
}