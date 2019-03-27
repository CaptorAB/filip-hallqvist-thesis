#include <vector>
#include <tuple>

#include <lib/random.h>
#include <include/genetic.h>

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
decode_chromosomes(std::vector<int> chromosomes, int population, int variables, int instruments, int genes, int bits)
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
    double total = 0.0;
    for (size_t j = 0; j < variables; ++j)
    {
      size_t ix = (i * variables) + j;
      total += individuals[ix];
    }
    if (total == 0.0)
    {
      for (size_t j = 0; j < variables; ++j)
      {
        size_t ix = (i * variables) + j;
        individuals[ix] = 1.0 / (double)variables;
      }
    }
    else
    {
      for (size_t j = 0; j < variables; ++j)
      {
        size_t ix = (i * variables) + j;
        individuals[ix] /= total;
      }
    }
  }

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

std::vector<double>
evaluate_individuals(std::vector<double> individuals, int population, int variables)
{
  std::vector<double> fitnesses(population);

  // Use dummy fitness calculation for now
  // Reward large values in first variable
  for (size_t i = 0; i < population; ++i)
  {
    double fitness = 0.0;
    for (size_t j = 0; j < variables; ++j)
    {
      size_t ix = (i * variables) + j;
      if (j == 0)
      {
        fitness += 1000 * individuals[ix];
      }
      else
      {
        fitness -= 100 * individuals[ix];
      }
    }
    fitnesses[i] = fitness;
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

Result
optimize(OptimizeOptions options)
{
  int population = options.population;
  int elitism = options.elitism;
  int generations = options.generations;
  int bits = options.bits;
  int steps = options.steps;
  int instruments = options.instruments;
  double mutation = options.mutation;
  double crossover = options.crossover;

  int scenarios = 1 << (steps - 1);
  int variables = scenarios * instruments;
  int genes = variables * bits;

  int chromosomes_size = population * genes;
  int individuals_size = population * variables;

  double global_max_fitness = 0.0;
  size_t global_max_chromosome = 0;

  std::vector<int> global_max_chromosome_vector(genes);
  std::vector<double> global_max_individual_vector(variables);

  // Generate initial chromosomes
  std::vector<int> chromosomes = initialize_chromosomes(population, genes);

  std::vector<double> individuals(individuals_size);
  std::vector<double> fitnesses(individuals_size);

  std::vector<int> crossovered(chromosomes_size);
  std::vector<int> mutated(chromosomes_size);
  std::vector<int> elitismed(chromosomes_size);

  for (size_t t = 0; t < generations; ++t)
  {
    individuals = decode_chromosomes(chromosomes, population, variables, instruments, genes, bits);
    fitnesses = evaluate_individuals(individuals, population, variables);

    // Check global_max_fitness
    for (size_t i = 0; i < population; ++i)
    {
      if (fitnesses[i] > global_max_fitness)
      {
        global_max_fitness = fitnesses[i];
        global_max_chromosome = i * genes;

        size_t ix = i * genes;
        for (size_t j = 0; j < genes; ++j)
        {
          global_max_chromosome_vector[j] = chromosomes[ix + j];
        }

        ix = i * variables;
        for (size_t j = 0; j < variables; ++j)
        {
          global_max_individual_vector[j] = individuals[ix + j];
        }
      }
    }

    std::vector<int> crossovered = crossover_chromosomes(chromosomes, fitnesses, population, genes, crossover);
    std::vector<int> mutated = mutate_chromosomes(crossovered, population, genes, mutation);
    std::vector<int> elitismed = elitism_chromosomes(chromosomes, mutated, global_max_chromosome, elitism, genes);

    chromosomes = elitismed;
  }

  Result result;
  result.fitness = global_max_fitness;
  result.chromosome = global_max_chromosome_vector;
  result.individual = global_max_individual_vector;

  return result;
}