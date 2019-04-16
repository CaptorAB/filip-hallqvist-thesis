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
  double initial_funding_ratio;
  double target_funding_ratio;
};

struct Result
{
  double fitness;
  double total_return;
  double risk;
  std::vector<double> individual;
};

Result optimize(OptimizeOptions options);

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> evaluate_individuals(const std::vector<double> &X, const std::vector<double> &price_changes, const std::vector<double> &probabilities, const std::vector<double> &goals, double risk_aversion, int n_individuals, int n_steps, int n_scenarios, int n_instruments);

#endif