#include <vector>
#include <tuple>
#include <iostream>
#include <string>

#include <lib/random.h>

#include <include/util.h>
#include <include/genetic.h>
#include <include/scenario.h>

using Random = effolkronium::random_static;

std::vector<double> normalize_individuals(std::vector<double> individuals, int n_individuals, int n_instruments, int n_scenarios)
{
  int n_genes = n_instruments * n_scenarios;
  int size = n_individuals * n_genes;

  std::vector<double> normalized(size);

  for (size_t i = 0; i < n_individuals; ++i)
  {
    size_t ix = i * n_genes;
    for (size_t s = 0; s < n_scenarios; ++s)
    {
      double total = 0.0;
      size_t sx = ix + (s * n_instruments);

      for (size_t j = 0; j < n_instruments; ++j)
      {
        size_t jx = sx + j;
        total += individuals[jx];
      }

      if (total == 0.0)
      {
        for (size_t j = 0; j < n_instruments; ++j)
        {
          size_t jx = sx + j;
          normalized[jx] = 1.0 / n_instruments;
        }
      }
      else
      {
        for (size_t j = 0; j < n_instruments; ++j)
        {
          size_t jx = sx + j;
          normalized[jx] = individuals[jx] / total;
        }
      }

    } // Scenarios

  } // Individuals

  return normalized;
}

std::vector<double>
initialize_individuals(int n_individuals, int n_instruments, int n_scenarios)
{
  int n_genes = n_instruments * n_scenarios;
  int size = n_individuals * n_genes;

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
  individuals = normalize_individuals(individuals, n_individuals, n_instruments, n_scenarios);

  return individuals;
}

std::tuple<size_t, size_t>
select_roulette(std::vector<double> fitnesses)
{
  {
    int i1 = Random::get<std::discrete_distribution<>>(fitnesses.begin(), fitnesses.end());
    int i2 = Random::get<std::discrete_distribution<>>(fitnesses.begin(), fitnesses.end());
    return std::make_tuple(i1, i2);
  }
}

// TODO: Refactor this mess
std::vector<double> evaluate_individuals(std::vector<double> X, std::vector<double> price_changes, std::vector<double> probabilities, std::vector<double> goals, double risk_aversion, double penalty, int n_individuals, int n_steps, int n_scenarios, int n_instruments)
{
  int timestamps = n_steps;
  int branching = 2;
  int instruments = n_instruments;
  int n_genes = n_instruments * n_scenarios;

  std::vector<double> fitnesses(n_individuals);

  for (int in = 0; in < n_individuals; ++in)
  {
    double wt = 1.0;
    double pt = 0.0; // Penalty
    double a = penalty;
    size_t xi = in * n_genes; // Index of current individual

    std::vector<double> wealth(n_scenarios);
    wealth[0] = 1.0; // TODO: Dynamic initial wealth (assets vs liabilities)
    std::vector<double> expected_wealth(get_nodes_in_level(timestamps - 1));
    std::vector<double> joint_probabilities(n_scenarios);
    joint_probabilities[0] = probabilities[0];

    for (int t = 1; t < timestamps; ++t)
    {
      int first_node = get_first_node_in_level(t);
      int ti = first_node * instruments; // Index of first node in step
      int ss = get_nodes_in_level(t);    // Scenarios in this step
      double ws = 0.0;                   // Total expected wealth in this step

      // TODO: Calculate wealth from assets for final step
      for (int s = 0; s < ss; ++s)
      {                                                               // Iterate over scenarios in step
        int si = ti + (s * n_instruments);                            // Index of current scenario
        int si_ = ((si / instruments) - 1) / branching * instruments; // Index of previous scenario

        int ri = first_node + s;           // Index of current scenario risk map
        int ri_ = get_parent_index(ri, 1); // Index of previous scenario risk map

        double aw = 0.0; // Change in wealth from assets
        double cw = 0.0; // Change in wealth from reallocations

        // Join probability with previous scenario
        joint_probabilities[ri] = joint_probabilities[ri_] * probabilities[ri];

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

        wealth[ri] = aw + cw;

        // Update penalty as needed
        if (wealth[ri] < goals[ri])
        {
          pt += pow(wealth[ri] - goals[ri], 2);
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

    fitnesses[in] = (1.0 - risk_aversion) * total_wealth - risk_aversion * pt;
  }

  return fitnesses;
}

std::vector<double>
mutate_individuals(std::vector<double> individuals, int n_individuals, int n_instruments, int n_scenarios, double mutation_rate)
{
  int n_genes = n_instruments * n_scenarios;
  int size = n_individuals * n_genes;

  std::vector<double> mutated(size);

  for (size_t i = 0; i < size; i++)
  {
    double random = Random::get(0.l, 1.l);
    if (random < mutation_rate)
    {
      mutated[i] = individuals[i] + Random::get<std::normal_distribution<>>(0.0, 1.0);
    }
  }

  mutated = normalize_individuals(individuals, n_individuals, n_instruments, n_scenarios);

  return mutated;
}

std::vector<double>
crossover_individuals(std::vector<double> individuals, std::vector<double> evaluated, int n_individuals, int n_genes, double crossover_rate)
{
  int size = n_individuals * n_genes;
  std::vector<double> crossovered(size);
  for (size_t i = 0; i < n_individuals; i += 2)
  {
    std::tuple<size_t, size_t> selected = select_roulette(evaluated);
    size_t cx1 = std::get<0>(selected) * n_genes;
    size_t cx2 = std::get<1>(selected) * n_genes;

    double random = Random::get(0.l, 1.l);
    size_t ix = (i * n_genes);

    // Zip individuals
    if (random < crossover_rate)
    {
      // Construct first individuals
      for (size_t j = 0; j < n_genes; j += 2)
      {
        crossovered[ix + j] = individuals[cx1 + j];
        crossovered[ix + j + 1] = individuals[cx2 + j + 1];
      }

      // Construct second individuals
      ix += n_genes;
      for (size_t j = 0; j < n_genes; j += 2)
      {
        crossovered[ix + j] = individuals[cx2 + j];
        crossovered[ix + j + 1] = individuals[cx1 + j + 1];
      }
    }
    else
    {
      for (size_t j = 0; j < n_genes; ++j)
      {
        crossovered[ix + j] = individuals[cx1 + j];
      }
      ix += n_genes;
      for (size_t j = 0; j < n_genes; ++j)
      {
        crossovered[ix + j] = individuals[cx2 + j];
      }
    }
  }
  return crossovered;
}

std::vector<double>
elitism_individuals(std::vector<double> old_individuals, std::vector<double> new_individuals, int i_elite_individual, int n_elitism_copies, int n_genes)
{
  for (size_t i = 0; i < n_elitism_copies; ++i)
  {
    size_t ix = i * n_genes;
    for (size_t j = 0; j < n_genes; ++j)
    {
      new_individuals[ix + j] = old_individuals[i_elite_individual + j];
    }
  }
  return new_individuals;
}

Result optimize(OptimizeOptions options)
{
  int n_individuals = options.population_size;
  int n_elitism_copies = options.elitism_copies;
  int n_generations = options.generations;
  int n_steps = options.steps;
  double mutation_rate = options.mutation_rate;
  double crossover_rate = options.crossover_rate;
  double risk_aversion = options.risk_aversion;
  double penalty = options.penalty_exponent;
  double surplus = options.goal_surplus;

  int n_instruments = N_INSTRUMENTS;
  int n_scenarios = (1 << n_steps) - 1;
  int n_genes = n_scenarios * n_instruments;

  double global_max_fitness = 0.0;
  int i_global_max_individual = 0;

  std::vector<double> global_max_individual(n_genes);

  // Generate initial individuals
  std::vector<double> individuals = initialize_individuals(n_individuals, n_instruments, n_scenarios);

  std::vector<double> fitnesses(n_individuals);
  std::vector<double> crossovered(n_individuals);
  std::vector<double> mutated(n_individuals);
  std::vector<double> elitismed(n_individuals);

  // Generate scenarios
  std::tuple<std::vector<double>, std::vector<double>> scenarios = generate_scenarios(n_steps);
  std::vector<double> price_changes = std::get<0>(scenarios);
  std::vector<double> probabilities = std::get<1>(scenarios);

  // Define goals
  std::vector<double> goals = generate_goals(price_changes, n_steps, n_scenarios, n_instruments, surplus);

  /*
  std::cout << "Price changes: \n";
  for (auto const &c : price_changes)
    std::cout << c << ' ';
  */

  // printf("Running...\n");

  for (size_t t = 0; t < n_generations; ++t)
  {
    fitnesses = evaluate_individuals(individuals, price_changes, probabilities, goals, risk_aversion, penalty, n_individuals, n_steps, n_scenarios, n_instruments);

    // Check global_max_fitness
    for (size_t i = 0; i < n_individuals; ++i)
    {
      if (fitnesses[i] > global_max_fitness)
      {
        global_max_fitness = fitnesses[i];
        i_global_max_individual = i * n_genes;

        // printf("(%i) Max fitness: %f\n", t, global_max_fitness);

        int ix = i * n_genes;
        for (size_t j = 0; j < n_genes; ++j)
        {
          global_max_individual[j] = individuals[ix + j];
        }
      }
    }

    crossovered = crossover_individuals(individuals, fitnesses, n_individuals, n_genes, crossover_rate);
    mutated = mutate_individuals(crossovered, n_individuals, n_instruments, n_scenarios, mutation_rate);
    elitismed = elitism_individuals(individuals, mutated, i_global_max_individual, n_elitism_copies, n_genes);

    individuals = elitismed;
  }

  Result result;
  result.fitness = global_max_fitness;
  result.individual = global_max_individual;

  return result;
}