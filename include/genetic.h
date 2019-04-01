#ifndef GENETIC_H
#define GENETIC_H

#include <vector>
#include <tuple>

#include <lib/random.h>

struct OptimizeOptions
{
  int population;
  int elitism;
  int generations;
  int bits;
  int steps;
  double mutation;
  double crossover;
};

struct Result
{
  double fitness;
  std::vector<int> chromosome;
  std::vector<double> individual;
};

std::vector<int> initialize_chromosomes(int population, int genes);

std::vector<double> decode_chromosomes(std::vector<int> chromosomes, int population, int variables, int scenarios, int instruments, int genes, int bits);

std::tuple<size_t, size_t> select_roulette(std::vector<double> fitnesses);

std::vector<double> evaluate_individuals(std::vector<double> individuals, std::vector<double> scenarios, std::vector<double> goals, int n_individuals, int n_variables, int n_steps, int n_instruments);

std::vector<int> mutate_chromosomes(std::vector<int> chromosomes, int population, int genes, double mutation);

std::vector<int> crossover_chromosomes(std::vector<int> chromosomes, std::vector<double> evaluated, int population, int genes, double crossover);

std::vector<int> elitism_chromosomes(std::vector<int> old_chromosomes, std::vector<int> new_chromosomes, int elite, int elitism, int genes);

Result optimize(OptimizeOptions options);

#endif