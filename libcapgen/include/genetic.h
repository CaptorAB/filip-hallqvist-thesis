#ifndef GENETIC_H
#define GENETIC_H

#include <vector>
#include <tuple>

using std::tuple;
using std::vector;

/**
 * Embind object for retrieving
 * an array of transaction costs
 * from javascript.
 */
struct TransactionCosts
{
    double domestic_equity;
    double global_equity;
    double real_estate;
    double alternative;
    double credit;
    double bonds_2y;
    double bonds_5y;
    double bonds_20y;
    double cash;
    double domestic_equity_future;
    double interest_rate_swap_2y;
    double interest_rate_swap_5y;
    double interest_rate_swap_20y;
};

/**
 * Embind object for retrieving
 * an array of instrument allocation
 * constraints from javascript.
 */
struct InstrumentConstraints
{
    double domestic_equity_min;
    double global_equity_min;
    double real_estate_min;
    double alternative_min;
    double credit_min;
    double bonds_2y_min;
    double bonds_5y_min;
    double bonds_20y_min;
    double cash_min;
    double domestic_equity_future_min;
    double interest_rate_swap_2y_min;
    double interest_rate_swap_5y_min;
    double interest_rate_swap_20y_min;
    double domestic_equity_max;
    double global_equity_max;
    double real_estate_max;
    double alternative_max;
    double credit_max;
    double bonds_2y_max;
    double bonds_5y_max;
    double bonds_20y_max;
    double cash_max;
    double domestic_equity_future_max;
    double interest_rate_swap_2y_max;
    double interest_rate_swap_5y_max;
    double interest_rate_swap_20y_max;
};

/**
 * Embind object for retrieving
 * an array of derivative margin
 * constraints from javascript.
 */
struct MarginConstraints
{
    double domestic_equity_future;
    double interest_rate_swap_2y;
    double interest_rate_swap_5y;
    double interest_rate_swap_20y;
};

/**
 * Embind object containing
 * all parameters used in the optimization.
 */
struct OptimizeOptions
{
    int population_size;
    int elitism_copies;
    int generations;
    int steps;
    double mutation_rate;
    double crossover_rate;
    double initial_funding_ratio;
    double target_funding_ratio;
    TransactionCosts transaction_costs;
    InstrumentConstraints instrument_constraints;
    MarginConstraints margin_constraints;
};

/**
 * Embind object for passing the
 * optimal solution back to javascript.
 */
struct Result
{
    double fitness;
    double expected_return;
    double expected_risk;
    vector<double> individual;
};

/**
 * Run the portfolio optimizer.
 */
Result optimize(OptimizeOptions options);

/**
 * Make sure that all (non-derivative) weights sum
 * up to 1.0.
 */
void normalize_individuals(
    vector<double> &individuals,
    const int n_individuals,
    const int n_instruments,
    const int n_derivatives,
    const int n_scenarios);

/**
 * Create a fresh batch of individuals (solutions).
 */
vector<double> initialize_individuals(
    const int n_individuals,
    const int n_instruments,
    const int n_scenarios);

/**
 * Perform roulette wheel selection to select two individuals
 * with probability proportional to their fitness.
 * 
 * @see https://en.wikipedia.org/wiki/Fitness_proportionate_selection
 */
tuple<int, int> select_roulette(
    vector<double> &fitnesses);

/**
 * Perform uniform crossover between two individuals.
 * 
 * @see https://en.wikipedia.org/wiki/Crossover_(genetic_algorithm)
 */
void crossover_individuals(
    vector<double> &selected,
    const double crossover_rate);

/**
 * Perform k-point crossover, making sure that only whole scenarios
 * are swapped.
 * 
 * @see https://en.wikipedia.org/wiki/Crossover_(genetic_algorithm)
 */
void crossover_individuals_scenario(
    vector<double> &selected,
    const int n_instruments,
    const int n_scenarios,
    const double crossover_rate);

/**
 * Mutate individuals based on the uniform distribution, and a scaling
 * parameter `d`, allowing greater mutations towards the end of
 * the scenario tree, and smaller mutations close to the root.
 * 
 * @see https://en.wikipedia.org/wiki/Mutation_(genetic_algorithm)
 */
void mutate_individuals(
    vector<double> &selected,
    const double mutation_rate);

/**
 * Compute the new wealth after evaluating a single
 * step in time.
 */
double compute_wealth(
    vector<double> &current_weights,
    vector<double> &next_weights,
    vector<double> &instrument_changes,
    vector<double> &transaction_costs,
    const double initial_wealth,
    const int n_instruments,
    const int n_derivatives);

/**
 * Compute the wealth for all steps in time.
 */
tuple<vector<double>, vector<double>> compute_wealths(
    vector<double> &individual,
    vector<double> &instrument_changes,
    vector<double> &transaction_costs,
    const int n_instruments,
    const int n_derivatives,
    const int n_scenarios);

/**
 * Compute the fitness of each individual in the population.
 */
vector<double> compute_fitnesses(
    vector<double> &individuals,
    vector<double> &instrument_changes,
    vector<double> &transaction_costs,
    vector<double> &intermediate_goals,
    vector<double> &final_goals,
    vector<double> &instrument_constraints,
    vector<double> &margin_constraints,
    const int n_individuals,
    const int n_instruments,
    const int n_derivatives,
    const int n_scenarios,
    const int generation);

/**
 * Compute the fitness of a single individual.
 */
double compute_fitness(
    vector<double> &individual,
    vector<double> &intermediate_wealths,
    vector<double> &final_wealths,
    vector<double> &intermediate_goals,
    vector<double> &final_goals,
    vector<double> &margin_constraints,
    vector<double> &instrument_constraints,
    const int n_instruments,
    const int n_derivatives,
    const int n_scenarios,
    const int generation);

/**
 * Compute the expected wealth of a portfolio.
 */
double compute_expected_wealth(
    vector<double> &final_wealths);

/**
 * Compute the penalty for breaking constraints.
 */
double compute_penalty(
    vector<double> &individual,
    vector<double> &instrument_constraints,
    vector<double> &margin_constraints,
    vector<double> &intermediate_wealths,
    vector<double> &final_wealths,
    vector<double> &intermediate_goals,
    vector<double> &final_goals,
    const int n_instruments,
    const int n_derivatives,
    const int n_scenarios,
    const int generation);

/**
 * Utility function for parsing an Embind
 * object into a vector of instrument allocation constraints.
 */
vector<double> parse_instrument_constraints(
    InstrumentConstraints instrument_constraints);

/**
 * Utility function for parsing an Embind
 * object into a vector of derivative margin constraints.
 */
vector<double> parse_margin_constraints(
    MarginConstraints margin_constraints);

#endif