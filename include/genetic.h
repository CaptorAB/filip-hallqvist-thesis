#ifndef GENETIC_H
#define GENETIC_H

#include <vector>
#include <tuple>

#include <lib/random.h>

struct OptimizeOptions
{
  int population_size;
  int elitism_copies;
  int generations;
  int steps;
  double mutation_rate;
  double crossover_rate;
  double risk_aversion;
  double penalty_exponent;
  double goal_surplus;
};

struct Result
{
  double fitness;
  std::vector<double> individual;
};

Result optimize(OptimizeOptions options);

#endif