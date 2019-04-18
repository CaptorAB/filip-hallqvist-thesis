#ifndef GENETIC_H
#define GENETIC_H

#include <vector>
#include <tuple>

#include <lib/random.h>

struct TransactionCosts
{
  double domestic_equity;
  double global_equity;
  double real_estate;
  double alternative;
  double credit;
  double bonds_2y;
  double bonds_5y;
  double cash;
  double fta;
  double domestic_equity_future;
  double interest_rate_swap_2y;
  double interest_rate_swap_5y;
  double interest_rate_swap_20y;
};

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
  TransactionCosts transaction_costs;
};

struct Result
{
  double fitness;
  double expected_return;
  double expected_risk;
  std::vector<double> individual;
  std::vector<double> incoming_wealths;
  std::vector<double> final_wealths;
  std::vector<double> price_changes;
  std::vector<double> goals;
};

Result optimize(OptimizeOptions options);

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> evaluate_individuals(std::vector<double> &X, std::vector<double> &price_changes, std::vector<double> &probabilities, std::vector<double> &goals, const double risk_aversion, const int n_individuals, const int n_steps, const int n_scenarios, const int n_instruments);

void normalize_individuals(std::vector<double> &individuals, const int n_individuals, const int n_instruments, const int n_scenarios);

std::vector<double> initialize_individuals(const int n_individuals, const int n_instruments, const int n_scenarios);

std::tuple<int, int> select_roulette(std::vector<double> &fitnesses);

void crossover_individuals(std::vector<double> &selected, const double crossover_rate);

void mutate_individuals(std::vector<double> &selected, const double mutation_rate);

double compute_wealth(std::vector<double> &current_weights, std::vector<double> &next_weights, std::vector<double> &price_changes, std::vector<double> &transaction_costs, const double initial_wealth);

std::tuple<std::vector<double>, std::vector<double>> compute_wealths(std::vector<double> &individual, std::vector<double> &price_changes, std::vector<double> transaction_costs, const int n_instruments, const int n_scenarios);

std::vector<double> compute_fitnesses(std::vector<double> &individuals, std::vector<double> &price_changes, std::vector<double> &transaction_costs, std::vector<double> &goals, const double risk_aversion, const int n_individuals, const int n_instruments, const int n_scenarios);

double compute_fitness(std::vector<double> &incoming_wealths, std::vector<double> &final_wealths, std::vector<double> &goals, const double risk_aversion);

double compute_expected_wealth(std::vector<double> &final_wealths);

double compute_expected_risk(std::vector<double> &incoming_wealths, std::vector<double> &goals);

#endif