#include <vector>
#include <tuple>
#include <iostream>
#include <string>

#include <lib/random.h>

#include <include/util.h>
#include <include/genetic.h>
#include <include/scenario.h>

using Random = effolkronium::random_static;

std::vector<int>
initialize_chromosomes(int population, int genes)
{
  int size = population * genes;
  std::vector<int> chromosomes(size);

  for (size_t i = 0; i < size; i++)
  {
    chromosomes[i] = Random::get<int>(0, 1);
  }

  return chromosomes;
}

std::vector<double>
decode_chromosomes(std::vector<int> chromosomes, int population, int variables, int n_scenarios, int instruments, int genes, int bits)
{
  int size = population * variables;
  std::vector<double> individuals(size);

  // Parse bits into real numbers
  for (size_t i = 0; i < population; ++i)
  {
    for (size_t j = 0; j < variables; ++j)
    {
      size_t ix = (i * variables) + j;
      individuals[ix] = 0.0;
      for (size_t k = 0; k < bits; ++k)
      {
        size_t cx = (i * genes) + (j * bits) + k;
        individuals[ix] += (double)(chromosomes[cx] << (bits - k - 1));
      }
    }
  }

  // Scale reals to [0, 1]
  for (size_t i = 0; i < population; ++i)
  {
    size_t ix = i * variables;
    for (size_t s = 0; s < n_scenarios; ++s)
    {

      double total = 0.0;
      size_t sx = ix + (s * instruments);

      for (size_t j = 0; j < instruments; ++j)
      {
        size_t jx = sx + j;
        total += individuals[jx];
      }

      if (total == 0.0)
      {
        for (size_t j = 0; j < instruments; ++j)
        {
          size_t jx = sx + j;
          individuals[jx] = 1.0 / (double)instruments;
        }
      }
      else
      {
        for (size_t j = 0; j < instruments; ++j)
        {
          size_t jx = sx + j;
          individuals[jx] /= total;
        }
      }

    } // Scenarios

  } // Individuals

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
std::vector<double> evaluate_individuals(std::vector<double> X, std::vector<double> price_changes, int n_individuals, int n_variables, int n_steps, int n_instruments)
{
  int timestamps = n_steps;
  int branching = 2;
  int instruments = n_instruments;

  std::vector<double> fitnesses(n_individuals);

  for (int in = 0; in < n_individuals; ++in)
  {

    double wt = 1.0;
    size_t xi = in * n_variables; // Index of current individual

    for (int t = 1; t < timestamps; ++t)
    {
      int first_node = get_first_node_in_level(t);
      int ti = first_node * instruments; // Index of first node in step
      int ss = get_nodes_in_level(t);    // Scenarios in this step
      double ws = 0.0;                   // Total expected wealth in this step

      for (int s = 0; s < ss; ++s)
      {                                                               // Iterate over scenarios in step
        int si = ti + (s * n_instruments);                            // Index of current scenario
        int si_ = ((si / instruments) - 1) / branching * instruments; // Index of previous scenario

        double p = 1.0 / ss; // Probability of scenario occurring (TODO: should actually be likelihood)
        double aw = 0.0;     // Wealth from assets
        double cw = 0.0;     // Wealth from reallocations

        // Calculate wealth from assets
        for (int i = 0; i < n_instruments; ++i)
        {
          int k_ = si_ + i; // Index of current instrument in previous scenario
          double r_ = 1.0 + price_changes[k_];
          double v = r_ * X[xi + k_];

          aw += v;
        }

        // Calculate wealth from reallocations
        for (int i = 0; i < n_instruments; ++i)
        {
          int k = si + i;   // Index of current instrument in current scenario
          int k_ = si_ + i; // Index of current instrument in previous scenario

          double c = aw * (X[xi + k] - X[xi + k_]);
          cw += c;
        }

        ws += p * (aw + cw);
      }

      wt *= ws;
    }

    fitnesses[in] = wt;
  }

  return fitnesses;
}

std::vector<int>
mutate_chromosomes(std::vector<int> chromosomes, int population, int genes, double mutation)
{
  int size = population * genes;
  std::vector<int> mutated(size);
  for (size_t i = 0; i < size; i++)
  {
    mutated[i] = chromosomes[i] != Random::get<std::bernoulli_distribution>(mutation);
  }

  return mutated;
}

std::vector<int>
crossover_chromosomes(std::vector<int> chromosomes, std::vector<double> evaluated, int population, int genes, double crossover)
{
  int size = population * genes;
  std::vector<int> crossovered(size);
  for (size_t i = 0; i < population; i += 2)
  {
    std::tuple<size_t, size_t> selected = select_roulette(evaluated);
    size_t cx1 = std::get<0>(selected) * genes;
    size_t cx2 = std::get<1>(selected) * genes;

    double random = Random::get(0, 1);
    size_t ix = (i * genes);

    // Zip chromosomes
    if (random < crossover)
    {
      // Construct first chromosome
      for (size_t j = 0; j < genes; j += 2)
      {
        crossovered[ix + j] = chromosomes[cx1 + j];
        crossovered[ix + j + 1] = chromosomes[cx2 + j + 1];
      }

      // Construct second chromosome
      ix += genes;
      for (size_t j = 0; j < genes; j += 2)
      {
        crossovered[ix + j] = chromosomes[cx2 + j];
        crossovered[ix + j + 1] = chromosomes[cx1 + j + 1];
      }
    }
    else
    {
      for (size_t j = 0; j < genes; ++j)
      {
        crossovered[ix + j] = chromosomes[cx1 + j];
      }
      ix += genes;
      for (size_t j = 0; j < genes; ++j)
      {
        crossovered[ix + j] = chromosomes[cx2 + j];
      }
    }
  }
  return crossovered;
}

std::vector<int>
elitism_chromosomes(std::vector<int> old_chromosomes, std::vector<int> new_chromosomes, int elite, int elitism, int genes)
{
  for (size_t i = 0; i < elitism; ++i)
  {
    size_t ix = i * genes;
    for (size_t j = 0; j < genes; ++j)
    {
      new_chromosomes[ix + j] = old_chromosomes[elite + j];
    }
  }
  return new_chromosomes;
}

Result optimize(OptimizeOptions options)
{
  int n_individuals = options.population;
  int n_elitism_copies = options.elitism;
  int n_generations = options.generations;
  int n_bits = options.bits;
  int n_steps = options.steps;
  double mutation_rate = options.mutation;
  double crossover_rate = options.crossover;

  int n_instruments = 2; // TODO: Dynamic
  int n_scenarios = (1 << n_steps) - 1;
  int n_variables = n_scenarios * n_instruments;
  int n_genes = n_variables * n_bits;

  int chromosomes_size = n_individuals * n_genes;
  int individuals_size = n_individuals * n_variables;

  double global_max_fitness = 0.0;
  size_t global_max_chromosome = 0;

  std::vector<int> global_max_chromosome_vector(n_genes);
  std::vector<double> global_max_individual_vector(n_variables);

  // Generate initial chromosomes
  std::vector<int> chromosomes = initialize_chromosomes(n_individuals, n_genes);

  std::vector<double> individuals(individuals_size);
  std::vector<double> fitnesses(n_individuals);

  std::vector<int> crossovered(chromosomes_size);
  std::vector<int> mutated(chromosomes_size);
  std::vector<int> elitismed(chromosomes_size);

  // Generate scenarios
  std::vector<double> price_changes = generate_price_changes(n_steps);

  /*
  std::cout << "Price changes: \n";
  for (auto const &c : price_changes)
    std::cout << c << ' ';
  */

  for (size_t t = 0; t < n_generations; ++t)
  {
    individuals = decode_chromosomes(chromosomes, n_individuals, n_variables, n_scenarios, n_instruments, n_genes, n_bits);
    fitnesses = evaluate_individuals(individuals, price_changes, n_individuals, n_variables, n_steps, n_instruments);

    // Check global_max_fitness
    for (size_t i = 0; i < n_individuals; ++i)
    {
      if (fitnesses[i] > global_max_fitness)
      {
        global_max_fitness = fitnesses[i];
        global_max_chromosome = i * n_genes;

        size_t ix = i * n_genes;
        for (size_t j = 0; j < n_genes; ++j)
        {
          global_max_chromosome_vector[j] = chromosomes[ix + j];
        }

        ix = i * n_variables;
        for (size_t j = 0; j < n_variables; ++j)
        {
          global_max_individual_vector[j] = individuals[ix + j];
        }
      }
    }

    std::vector<int> crossovered = crossover_chromosomes(chromosomes, fitnesses, n_individuals, n_genes, crossover_rate);
    std::vector<int> mutated = mutate_chromosomes(crossovered, n_individuals, n_genes, mutation_rate);
    std::vector<int> elitismed = elitism_chromosomes(chromosomes, mutated, global_max_chromosome, n_elitism_copies, n_genes);

    chromosomes = elitismed;
  }

  Result result;
  result.fitness = global_max_fitness;
  result.chromosome = global_max_chromosome_vector;
  result.individual = global_max_individual_vector;

  return result;
}