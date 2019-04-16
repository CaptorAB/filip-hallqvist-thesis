#include <vector>
#include <tuple>
#include <iostream>
#include <string>
#include <algorithm>
#include <limits>

#include <lib/random.h>

#include <include/util.h>
#include <include/genetic.h>
#include <include/scenario.h>

using Random = effolkronium::random_static;

void normalize_individuals(std::vector<double> &individuals, int n_individuals, int n_instruments, int n_scenarios)
{
  int n_genes = n_instruments * n_scenarios;

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
          individuals[jx] = 1.0 / n_instruments;
        }
      }
      else
      {
        for (size_t j = 0; j < n_instruments; ++j)
        {
          size_t jx = sx + j;
          individuals[jx] = individuals[jx] / total;
        }
      }

    } // Scenarios

  } // Individuals
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
  normalize_individuals(individuals, n_individuals, n_instruments, n_scenarios);

  return individuals;
}

std::tuple<size_t, size_t>
select_roulette(const std::vector<double> &fitnesses)
{
  {
    int i1 = Random::get<std::discrete_distribution<>>(fitnesses.begin(), fitnesses.end());
    int i2 = Random::get<std::discrete_distribution<>>(fitnesses.begin(), fitnesses.end());
    return std::make_tuple(i1, i2);
  }
}

// TODO: Refactor this mess
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> evaluate_individuals(const std::vector<double> &X, const std::vector<double> &price_changes, const std::vector<double> &probabilities, const std::vector<double> &goals, double risk_aversion, int n_individuals, int n_steps, int n_scenarios, int n_instruments)
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
    double pt = 0.0;          // Penalty
    size_t xi = in * n_genes; // Index of current individual

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

std::vector<double>
mutate_individuals(std::vector<double> &individuals, int n_individuals, int n_instruments, int n_scenarios, double mutation_rate)
{
  int n_genes = n_instruments * n_scenarios;
  int size = n_individuals * n_genes;

  std::vector<double> mutated(size);

  for (size_t i = 0; i < size; i++)
  {
    double random = Random::get(0.l, 1.l);
    if (random < mutation_rate)
    {
      mutated[i] = std::max(0.0, individuals[i] + Random::get<std::normal_distribution<>>(0.0, 0.1));
    }
    else
    {
      mutated[i] = individuals[i];
    }
  }

  normalize_individuals(mutated, n_individuals, n_instruments, n_scenarios);

  return mutated;
}

std::vector<double>
crossover_individuals(std::vector<double> &individuals, std::vector<double> &evaluated, int n_individuals, int n_genes, double crossover_rate)
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
elitism_individuals(std::vector<double> &old_individuals, std::vector<double> &new_individuals, int i_elite_individual, int n_elitism_copies, int n_genes)
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
  double initial_funding_ratio = options.initial_funding_ratio;
  double target_funding_ratio = options.target_funding_ratio;

  int n_instruments = N_INSTRUMENTS;
  int n_scenarios = (1 << n_steps) - 1;
  int n_genes = n_scenarios * n_instruments;

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
  std::vector<double> crossovered(n_individuals);
  std::vector<double> mutated(n_individuals);
  std::vector<double> elitismed(n_individuals);

  // Generate scenarios
  std::tuple<std::vector<double>, std::vector<double>> scenarios = generate_scenarios(n_steps);
  std::vector<double> price_changes = std::get<0>(scenarios);
  std::vector<double> probabilities = std::get<1>(scenarios);

  // Define goals
  std::vector<double> goals = generate_goals(price_changes, n_steps, n_scenarios, n_instruments, initial_funding_ratio, target_funding_ratio);

  /*
  std::cout << "Price changes: \n";
  for (auto const &c : price_changes)
    std::cout << c << ' ';
  std::cout << std::endl;

  std::cout << "Goals: \n";
  for (auto const &c : goals)
    std::cout << c << ' ';
  std::cout << std::endl;
  std::cout << std::endl;
  */

  // printf("Running...\n");

  for (size_t t = 0; t < n_generations; ++t)
  {
    evaluated = evaluate_individuals(individuals, price_changes, probabilities, goals, risk_aversion, n_individuals, n_steps, n_scenarios, n_instruments);

    fitnesses = std::get<0>(evaluated);
    total_returns = std::get<1>(evaluated);
    risks = std::get<2>(evaluated);

    // Check global_max_fitness
    for (size_t i = 0; i < n_individuals; ++i)
    {
      if (fitnesses[i] > global_max_fitness)
      {
        global_max_fitness = fitnesses[i];
        global_max_total_return = total_returns[i];
        global_max_risk = risks[i];
        i_global_max_individual = i * n_genes;

        /*
        printf("(%i) Fitness: %.10f, Total return: %.10f\n", t, global_max_fitness, global_max_total_return);
        std::cout << "\n";

        int a = 1;
        double s = 0.0;
        for (int u = 0; u < n_genes; u++)
        {
          double c = individuals[i_global_max_individual + u];
          s += c;
          std::printf("%.2f ", c);
          if (a % N_INSTRUMENTS == 0)
          {
            std::cout << " == " << s;
            s = 0;
            std::cout << "\n";
          }
          a++;
        }
        std::cout << "\n";
        */

        int ix = i * n_genes;
        for (size_t j = 0; j < n_genes; ++j)
        {
          global_max_individual[j] = individuals[ix + j];
        }
      }
    }

    int m;
    crossovered = individuals; // crossover_individuals(individuals, fitnesses, n_individuals, n_genes, crossover_rate);
    mutated = mutate_individuals(crossovered, n_individuals, n_instruments, n_scenarios, mutation_rate);
    elitismed = elitism_individuals(individuals, mutated, i_global_max_individual, n_elitism_copies, n_genes);

    individuals = elitismed;
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