#include <vector>
#include <tuple>
#include <iostream>
#include <string>
#include <algorithm>
#include <limits>
#include <queue>

#include <lib/random.h>

#include <include/util.h>
#include <include/genetic.h>
#include <include/scenario.h>

using Random = effolkronium::random_static;

void normalize_individuals(std::vector<double> &individuals, const int n_individuals, const int n_instruments, const int n_scenarios)
{
  const int n_genes = n_instruments * n_scenarios;

  for (int i = 0; i < n_individuals; ++i)
  {
    int ix = i * n_genes;
    for (int s = 0; s < n_scenarios; ++s)
    {
      double total = 0.0;
      int sx = ix + (s * n_instruments);

      for (int j = 0; j < n_instruments; ++j)
      {
        int jx = sx + j;
        total += individuals[jx];
      }

      if (total == 0.0)
      {
        for (int j = 0; j < n_instruments; ++j)
        {
          int jx = sx + j;
          individuals[jx] = 1.0 / n_instruments;
        }
      }
      else
      {
        for (int j = 0; j < n_instruments; ++j)
        {
          int jx = sx + j;
          individuals[jx] = individuals[jx] / total;
        }
      }

    } // Scenarios

  } // Individuals
}

std::vector<double> initialize_individuals(const int n_individuals, const int n_instruments, const int n_scenarios)
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

  // Normalize individuals
  normalize_individuals(individuals, n_individuals, n_instruments, n_scenarios);

  return individuals;
}

std::tuple<int, int> select_roulette(std::vector<double> &fitnesses)
{
  {
    int i1 = Random::get<std::discrete_distribution<>>(fitnesses.begin(), fitnesses.end());
    int i2 = Random::get<std::discrete_distribution<>>(fitnesses.begin(), fitnesses.end());
    return std::make_tuple(i1, i2);
  }
}

double compute_wealth(std::vector<double> &current_weights, std::vector<double> &next_weights, std::vector<double> &price_changes, std::vector<double> &transaction_costs, const double initial_wealth)
{
  double holdings = 0.0;
  double reallocations = 0.0;

  // Compute wealth from price changes
  for (int i = 0; i < current_weights.size(); ++i)
  {
    const double result = initial_wealth * current_weights[i] * (1.0 + price_changes[i]);
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

std::tuple<std::vector<double>, std::vector<double>> compute_wealths(std::vector<double> &individual, std::vector<double> &price_changes, std::vector<double> transaction_costs, const int n_instruments, const int n_scenarios)
{
  std::vector<double> incoming_wealths(n_scenarios);
  incoming_wealths[0] = 1.0;

  std::vector<double> final_wealths(n_scenarios / 2 + 1);
  int final_index = 0;

  for (int i = 0; i < n_scenarios; ++i)
  {
    int current = i;

    int left = 2 * current + 1;
    int right = 2 * current + 2;

    const int cx = current * n_instruments;
    const double current_wealth = incoming_wealths[current];

    std::vector<double> current_changes = std::vector<double>(price_changes.begin() + cx, price_changes.begin() + cx + n_instruments);
    std::vector<double> current_weights = std::vector<double>(individual.begin() + cx, individual.begin() + cx + n_instruments);

    if (left < n_scenarios && right < n_scenarios)
    {
      const int lx = left * n_instruments;
      const int rx = right * n_instruments;

      std::vector<double> left_weights = std::vector<double>(individual.begin() + lx, individual.begin() + lx + n_instruments);
      std::vector<double> right_weights = std::vector<double>(individual.begin() + rx, individual.begin() + rx + n_instruments);

      incoming_wealths[left] = compute_wealth(current_weights, left_weights, current_changes, transaction_costs, current_wealth);
      incoming_wealths[right] = compute_wealth(current_weights, right_weights, current_changes, transaction_costs, current_wealth);
    }
    else
    {
      final_wealths[final_index] = compute_wealth(current_weights, current_weights, current_changes, transaction_costs, current_wealth);
      final_index++;
    }
  }

  return std::make_tuple(incoming_wealths, final_wealths);
}

double compute_fitness(std::vector<double> &final_wealths)
{
  double expected_wealth = 0.0;
  for (double wealth : final_wealths)
  {
    expected_wealth += wealth;
  }
  expected_wealth = expected_wealth / final_wealths.size();

  return expected_wealth;
}

std::vector<double> compute_fitnesses(std::vector<double> &individuals, std::vector<double> price_changes, std::vector<double> transaction_costs, const int n_individuals, const int n_instruments, const int n_scenarios)
{
  const int n_genes = n_instruments * n_scenarios;
  std::vector<double> fitnesses(n_individuals);
  for (int i = 0; i < n_individuals; ++i)
  {
    const int ix = i * n_genes;

    std::vector<double> individual = std::vector<double>(individuals.begin() + ix, individuals.begin() + ix + n_genes);

    std::tuple<std::vector<double>, std::vector<double>> wealths = compute_wealths(individual, price_changes, transaction_costs, n_instruments, n_scenarios);

    std::vector<double> incoming_wealths = std::get<0>(wealths);
    std::vector<double> final_wealths = std::get<1>(wealths);

    fitnesses[i] = compute_fitness(final_wealths);
  }
  return fitnesses;
}

/*
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> evaluate_individuals(std::vector<double> &X, std::vector<double> &price_changes, std::vector<double> &probabilities, std::vector<double> &goals, const double risk_aversion, const int n_individuals, const int n_steps, const int n_scenarios, const int n_instruments)
{
  int timestamps = n_steps;
  int branching = 2;
  int instruments = n_instruments;
  int n_genes = n_instruments * n_scenarios;

  std::vector<double> fitnesses(n_individuals);
  std::vector<double> total_returns(n_individuals);
  std::vector<double> risks(n_individuals);

  for (int in = 0; in < n_individuals; ++in)
  {
    double pt = 0.0;       // Penalty
    int xi = in * n_genes; // Index of current individual

    std::vector<double> wealth(n_scenarios);
    wealth[0] = 1.0;
    std::vector<double> expected_wealth(get_nodes_in_level(timestamps - 1));
    std::vector<double> joint_probabilities(n_scenarios);
    joint_probabilities[0] = probabilities[0];

    for (int t = 1; t < timestamps; ++t)
    {
      int first_node = get_first_node_in_level(t);
      int ti = first_node * instruments; // Index of first node in step
      int ss = get_nodes_in_level(t);    // Scenarios in this step

      // Update probabilities
      double total_probability = 0.0;
      for (int s = 0; s < ss; ++s)
      {                                                               // Iterate over scenarios in step
        int si = ti + (s * n_instruments);                            // Index of current scenario
        int si_ = ((si / instruments) - 1) / branching * instruments; // Index of previous scenario

        int ri = first_node + s;           // Index of current scenario risk map
        int ri_ = get_parent_index(ri, 1); // Index of previous scenario risk map

        // Join probability with previous scenario
        joint_probabilities[ri] = joint_probabilities[ri_] * probabilities[ri];
        total_probability += joint_probabilities[ri];
      }

      // Normalize probabilities
      for (int s = 0; s < ss; ++s)
      {                                                               // Iterate over scenarios in step
        int si = ti + (s * n_instruments);                            // Index of current scenario
        int si_ = ((si / instruments) - 1) / branching * instruments; // Index of previous scenario
        int ri = first_node + s;                                      // Index of current scenario risk map
        int ri_ = get_parent_index(ri, 1);                            // Index of previous scenario risk map

        // Join probability with previous scenario
        joint_probabilities[ri] = joint_probabilities[ri] / total_probability;
      }

      // TODO: Calculate wealth from assets for final step
      for (int s = 0; s < ss; ++s)
      {                                                               // Iterate over scenarios in step
        int si = ti + (s * n_instruments);                            // Index of current scenario
        int si_ = ((si / instruments) - 1) / branching * instruments; // Index of previous scenario

        int ri = first_node + s;           // Index of current scenario risk map
        int ri_ = get_parent_index(ri, 1); // Index of previous scenario risk map

        double aw = 0.0; // Change in wealth from assets
        double cw = 0.0; // Change in wealth from reallocations

        // Calculate wealth from assets
        for (int i = 0; i < n_instruments; ++i)
        {
          int k = si + i;   // Index of current instrumet in current scenario
          int k_ = si_ + i; // Index of current instrument in previous scenario
          double v = wealth[ri_] * X[xi + k_] * (1.0 + price_changes[k]);
          aw += v;
        }

        // TODO: Incorporate transaction costs
        /*
        for (int i = 0; i < n_instruments; ++i)
        {
          int k = si + i;   // Index of current instrument in current scenario
          int k_ = si_ + i; // Index of current instrument in previous scenario

          double diff = X[xi + k] - X[xi + k_];
          double c = aw * diff; // TODO: Add transaction costs (also make sure that we can afford them)
          cw += c;
        }
        */
/*

        wealth[ri] = aw + cw;

        // Update penalty as needed
        if (wealth[ri] < goals[ri])
        {
          pt += joint_probabilities[ri] * pow(wealth[ri] - goals[ri], 2);
        }

        // Update expected wealth if we are on final timestamp
        if (t == timestamps - 1)
        {
          expected_wealth[s] = joint_probabilities[ri] * wealth[ri];
        }
      }

      // TODO: Discount this value
    }

    // Sum up final wealths
    double total_wealth = 0.0;
    for (int i = 0; i < get_nodes_in_level(timestamps - 1); i++)
    {
      total_wealth += expected_wealth[i];
    }

    total_returns[in] = total_wealth;
    risks[in] = pt;
    fitnesses[in] = (1.0 - risk_aversion) * total_wealth - risk_aversion * pt;
  }

  return std::make_tuple(fitnesses, total_returns, risks);
}
*/

void mutate_individuals(std::vector<double> &selected, const double mutation_rate)
{
  const int n_genes = floor(selected.size() / 2);
  for (int j = 0; j < n_genes; ++j)
  {
    const double random = Random::get(0.0, 1.0);
    if (random < mutation_rate)
    {
      const double d = j / (j + 1);
      selected[j] = std::min(1.0, std::max(0.0, selected[j] + d * Random::get(-1.0, 1.0)));
    }
  }
  for (int j = n_genes; j < selected.size(); ++j)
  {
    const double random = Random::get(0.0, 1.0);
    if (random < mutation_rate)
    {
      const double d = (j - n_genes) / ((j - n_genes) + 1);
      selected[j] = std::min(1.0, std::max(0.0, selected[j] + d * Random::get(-1.0, 1.0)));
    }
  }
}

void crossover_individuals_scenario(std::vector<double> &selected, const int n_instruments, const int n_scenarios, const double crossover_rate)
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

void crossover_individuals(std::vector<double> &selected, const double crossover_rate)
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

  const int n_instruments = N_INSTRUMENTS;
  const int n_scenarios = (1 << n_steps) - 1;
  const int n_genes = n_scenarios * n_instruments;

  // TODO: Use "best" instead of "max"
  double global_max_fitness = std::numeric_limits<double>::min();
  double global_max_total_return = std::numeric_limits<double>::min();
  double global_max_risk = std::numeric_limits<double>::min();

  int i_global_max_individual = 0;

  std::vector<double> global_max_individual(n_genes);

  // Generate initial individuals
  std::vector<double> individuals = initialize_individuals(n_individuals, n_instruments, n_scenarios);

  std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> evaluated;
  std::vector<double> fitnesses(n_individuals);
  std::vector<double> total_returns(n_individuals);
  std::vector<double> risks(n_individuals);

  // Generate scenarios
  std::tuple<std::vector<double>, std::vector<double>> scenarios = generate_scenarios(n_steps);
  std::vector<double> price_changes = std::get<0>(scenarios);
  std::vector<double> probabilities = std::get<1>(scenarios);

  std::cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!n";
  std::cout << "Price changes: \n";
  for (auto c : price_changes)
    std::printf("%.2f ", c);
  std::cout << "\n";

  // Define transaction costs (currently 0.0 for all instruments)
  std::vector<double> transaction_costs(n_instruments, 0.1);

  // Define goals
  std::vector<double> goals = generate_goals(price_changes, n_steps, n_scenarios, n_instruments, initial_funding_ratio, target_funding_ratio);

  for (int t = 0; t < n_generations; ++t)
  {
    fitnesses = compute_fitnesses(individuals, price_changes, transaction_costs, n_individuals, n_instruments, n_scenarios);

    // Check global_max_fitness
    for (int i = 0; i < n_individuals; ++i)
    {
      if (fitnesses[i] > global_max_fitness)
      {
        global_max_fitness = fitnesses[i];
        i_global_max_individual = i * n_genes;

        int ix = i * n_genes;
        global_max_individual = std::vector<double>(individuals.begin() + ix, individuals.begin() + ix + n_genes);

        std::cout << t << " @@@@@@@@@\n";
        std::cout << "Fitness: " << global_max_fitness << "\n";
        for (auto c : global_max_individual)
          std::printf("%.2f ", c);
        std::cout << "\n";
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
      crossover_individuals_scenario(selected, n_instruments, n_scenarios, crossover_rate);

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
        offspring[ix + j] = individuals[i_global_max_individual + j];
      }
    }

    // Normalize individuals
    normalize_individuals(offspring, n_individuals, n_instruments, n_scenarios);

    // Replace generation
    individuals = offspring;
  }

  Result result;
  result.fitness = global_max_fitness;
  result.total_return = global_max_total_return;
  result.risk = global_max_risk;
  result.price_changes = price_changes;
  result.probabilities = probabilities;
  result.goals = goals;
  result.individual = global_max_individual;

  return result;
}