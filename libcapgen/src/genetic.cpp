#include <iostream>
#include <cmath>
#include <tuple>
#include <vector>

#include <lib/random.h>

#include <include/constants.h>
#include <include/genetic.h>
#include <include/goals.h>
#include <include/trees.h>
#include <include/simulation.h>

using Random = effolkronium::random_static;

std::vector<double>
parse_instrument_constraints(InstrumentConstraints instrument_constraints)
{
  return std::vector<double>(
      {instrument_constraints.domestic_equity_min,
       instrument_constraints.global_equity_min,
       instrument_constraints.real_estate_min,
       instrument_constraints.alternative_min,
       instrument_constraints.credit_min,
       instrument_constraints.bonds_2y_min,
       instrument_constraints.bonds_5y_min,
       instrument_constraints.bonds_20y_min,
       instrument_constraints.cash_min,
       instrument_constraints.domestic_equity_future_min,
       instrument_constraints.interest_rate_swap_2y_min,
       instrument_constraints.interest_rate_swap_5y_min,
       instrument_constraints.interest_rate_swap_20y_min,
       instrument_constraints.domestic_equity_max,
       instrument_constraints.global_equity_max,
       instrument_constraints.real_estate_max,
       instrument_constraints.alternative_max,
       instrument_constraints.credit_max,
       instrument_constraints.bonds_2y_max,
       instrument_constraints.bonds_5y_max,
       instrument_constraints.bonds_20y_max,
       instrument_constraints.cash_max,
       instrument_constraints.domestic_equity_future_max,
       instrument_constraints.interest_rate_swap_2y_max,
       instrument_constraints.interest_rate_swap_5y_max,
       instrument_constraints.interest_rate_swap_20y_max});
}

std::vector<double>
parse_margin_constraints(MarginConstraints margin_constraints)
{
  return std::vector<double>({margin_constraints.domestic_equity_future,
                              margin_constraints.interest_rate_swap_2y,
                              margin_constraints.interest_rate_swap_5y,
                              margin_constraints.interest_rate_swap_20y});
}

void normalize_individuals(std::vector<double> &individuals,
                           const int n_individuals,
                           const int n_instruments,
                           const int n_derivatives,
                           const int n_scenarios)
{
  const int n_non_derivatives = n_instruments - n_derivatives;
  const int n_genes = n_instruments * n_scenarios;

  for (int i = 0; i < n_individuals; ++i)
  {
    int ix = i * n_genes;
    for (int s = 0; s < n_scenarios; ++s)
    {
      double total = 0.0;
      int sx = ix + (s * n_instruments);

      for (int j = 0; j < n_non_derivatives; ++j)
      {
        int jx = sx + j;
        total += individuals[jx];
      }

      if (total == 0.0)
      {
        for (int j = 0; j < n_non_derivatives; ++j)
        {
          int jx = sx + j;
          individuals[jx] = 1.0 / n_instruments;
        }
      }
      else
      {
        for (int j = 0; j < n_non_derivatives; ++j)
        {
          int jx = sx + j;
          individuals[jx] = individuals[jx] / total;
        }
      }

    } // Scenarios

  } // Individuals
}

std::vector<double> initialize_individuals(const int n_individuals,
                                           const int n_instruments,
                                           const int n_scenarios)
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
  return individuals;
}

std::tuple<int, int> select_roulette(std::vector<double> &fitnesses)
{
  {
    int i1 = Random::get<std::discrete_distribution<>>(fitnesses.begin(),
                                                       fitnesses.end());
    int i2 = Random::get<std::discrete_distribution<>>(fitnesses.begin(),
                                                       fitnesses.end());
    return std::make_tuple(i1, i2);
  }
}

double compute_wealth(std::vector<double> &current_weights,
                      std::vector<double> &next_weights,
                      std::vector<double> &instrument_changes,
                      std::vector<double> &transaction_costs,
                      const double initial_wealth,
                      const int n_instruments,
                      const int n_derivatives)
{
  const int n_non_derivatives = n_instruments - n_derivatives;
  double holdings = 0.0;
  double reallocations = 0.0;

  // Compute wealth from price changes
  for (int i = 0; i < n_instruments; ++i)
  {
    double result = 0.0;
    // printf("(%i) initial = %.2f, current = %.2f, change = %.2f \n", i, initial_wealth, current_weights[i], instrument_changes[i]);
    if (i < n_non_derivatives)
    {
      result = initial_wealth * current_weights[i] * (1.0 + instrument_changes[i]);
    }
    else
    {
      result = initial_wealth * current_weights[i] * instrument_changes[i];
    }
    holdings += result;
  }

  // Deduct transation costs
  for (int i = 0; i < current_weights.size(); ++i)
  {
    const double diff = fabs(current_weights[i] - next_weights[i]);
    reallocations += holdings * diff * transaction_costs[i];
  }

  const double wealth = holdings - reallocations;

  return std::max(0.0, wealth);
}

std::tuple<std::vector<double>, std::vector<double>>
compute_wealths(std::vector<double> &individual,
                std::vector<double> &instrument_changes,
                std::vector<double> &transaction_costs,
                const int n_instruments,
                const int n_derivatives,
                const int n_scenarios)
{
  std::vector<double> intermediate_wealths(n_scenarios);
  intermediate_wealths[0] = 1.0;

  std::vector<double> final_wealths(n_scenarios / 2 + 1);
  int final_index = 0;

  // Iterate through scenarios
  for (int i = 0; i < n_scenarios; ++i)
  {
    int current = i;

    int left = 2 * current + 1;
    int right = 2 * current + 2;

    // Index of current node
    const int cx = current * n_instruments;

    const double current_wealth = intermediate_wealths[current];

    std::vector<double> current_changes =
        std::vector<double>(instrument_changes.begin() + cx,
                            instrument_changes.begin() + cx + n_instruments);

    std::vector<double> current_weights = std::vector<double>(
        individual.begin() + cx, individual.begin() + cx + n_instruments);

    if (left < n_scenarios && right < n_scenarios)
    {
      const int lx = left * n_instruments;
      const int rx = right * n_instruments;

      std::vector<double> left_weights = std::vector<double>(
          individual.begin() + lx, individual.begin() + lx + n_instruments);

      std::vector<double> right_weights = std::vector<double>(
          individual.begin() + rx, individual.begin() + rx + n_instruments);

      // Evaluate left child
      // printf("Evaluating left child %i \n", left);
      intermediate_wealths[left] = compute_wealth(
          current_weights, left_weights, current_changes, transaction_costs,
          current_wealth, n_instruments, n_derivatives);

      // Evaluate right child
      // printf("Evaluating right child %i \n", right);
      intermediate_wealths[right] = compute_wealth(
          current_weights, right_weights, current_changes, transaction_costs,
          current_wealth, n_instruments, n_derivatives);
    }
    else
    {
      // We are at the leaf nodes of the scenario tree,
      // so compute the final wealth.

      // printf("Computing final index %i \n", final_index);

      final_wealths[final_index] = compute_wealth(
          current_weights, current_weights, current_changes, transaction_costs,
          current_wealth, n_instruments, n_derivatives);

      final_index++;
    }
  }

  return std::make_tuple(intermediate_wealths, final_wealths);
}

double compute_expected_wealth(std::vector<double> &final_wealths)
{
  double expected_wealth = 0.0;
  for (double wealth : final_wealths)
  {
    expected_wealth += wealth;
  }
  expected_wealth = expected_wealth / final_wealths.size();

  return expected_wealth;
}

double compute_expected_risk(std::vector<double> &intermediate_wealths,
                             std::vector<double> &final_wealths,
                             std::vector<double> &intermediate_goals,
                             std::vector<double> &final_goals)
{
  double risk = 0.0;

  // Intermediate goals
  for (int i = 0; i < intermediate_wealths.size(); ++i)
  {
    const double step = floor(std::log2(i + 1));
    const double nodes_in_step = pow(2.0, step);

    risk +=
        pow(std::min(0.0, intermediate_wealths[i] - intermediate_goals[i]),
            2.0) /
        nodes_in_step;
  }

  // Final goals
  for (int i = 0; i < final_wealths.size(); ++i)
  {
    risk += pow(std::min(0.0, final_wealths[i] - final_goals[i]), 2.0) /
            final_wealths.size();
  }
  return risk;
}

double compute_fitness(std::vector<double> &individual,
                       std::vector<double> &intermediate_wealths,
                       std::vector<double> &final_wealths,
                       std::vector<double> &intermediate_goals,
                       std::vector<double> &final_goals,
                       std::vector<double> &instrument_constraints,
                       std::vector<double> &margin_constraints,
                       const int n_instruments,
                       const int n_derivatives,
                       const int n_scenarios,
                       const int generation)
{
  const double wealth = compute_expected_wealth(final_wealths);
  const double penalty =
      compute_penalty(individual, instrument_constraints, margin_constraints,
                      intermediate_wealths, final_wealths, intermediate_goals,
                      final_goals, n_instruments, n_derivatives, n_scenarios, generation);

  return std::max(0.0, wealth - penalty);
}

double compute_penalty(std::vector<double> &individual,
                       std::vector<double> &instrument_constraints,
                       std::vector<double> &margin_constraints,
                       std::vector<double> &intermediate_wealths,
                       std::vector<double> &final_wealths,
                       std::vector<double> &intermediate_goals,
                       std::vector<double> &final_goals,
                       const int n_instruments,
                       const int n_derivatives,
                       const int n_scenarios,
                       const int generation)
{
  const int n_non_derivatives = n_instruments - n_derivatives;
  double penalty = 0.0;

  // Check instrument allocation constraints
  for (int j = 0; j < n_scenarios; ++j)
  {
    const int jx = j * n_instruments;

    // Allocation constraints
    for (int k = 0; k < n_instruments; ++k)
    {
      // Minimum
      penalty += pow(
          std::max(0.0, instrument_constraints[k] - individual[jx + k]), 2.0);

      // Maximum
      penalty +=
          pow(std::max(0.0, individual[jx + k] -
                                instrument_constraints[k + n_instruments]),
              2.0);
    }

    // Margin constraints
    double remaining_margin =
        individual[jx + CASH_INDEX] + individual[jx + BONDS_20Y_INDEX];

    for (int k = 0; k < n_derivatives; ++k)
    {
      const int kx = n_non_derivatives + k;
      const double required_margin =
          margin_constraints[k] * individual[jx + kx];
      remaining_margin -= required_margin;
    }
    penalty += pow(std::min(0.0, remaining_margin), 2.0);
  }

  // Intermediate goals
  for (int i = 0; i < intermediate_wealths.size(); ++i)
  {
    const double step = floor(std::log2(i + 1));
    const double nodes_in_step = pow(2.0, step);

    penalty +=
        pow(std::min(0.0, intermediate_wealths[i] - intermediate_goals[i]),
            2.0) /
        nodes_in_step;
  }

  // Final goals
  for (int i = 0; i < final_wealths.size(); ++i)
  {
    penalty += pow(std::min(0.0, final_wealths[i] - final_goals[i]), 2.0) /
               final_wealths.size();
  }

  return penalty;
}

std::vector<double>
compute_fitnesses(std::vector<double> &individuals,
                  std::vector<double> &instrument_changes,
                  std::vector<double> &transaction_costs,
                  std::vector<double> &intermediate_goals,
                  std::vector<double> &final_goals,
                  std::vector<double> &instrument_constraints,
                  std::vector<double> &margin_constraints,
                  const int n_individuals,
                  const int n_instruments,
                  const int n_derivatives,
                  const int n_scenarios,
                  const int generation)
{
  const int n_genes = n_instruments * n_scenarios;
  std::vector<double> fitnesses(n_individuals);
  for (int i = 0; i < n_individuals; ++i)
  {
    const int ix = i * n_genes;

    std::vector<double> individual = std::vector<double>(
        individuals.begin() + ix, individuals.begin() + ix + n_genes);

    std::tuple<std::vector<double>, std::vector<double>> wealths =
        compute_wealths(individual, instrument_changes, transaction_costs,
                        n_instruments, n_derivatives, n_scenarios);

    std::vector<double> intermediate_wealths = std::get<0>(wealths);
    std::vector<double> final_wealths = std::get<1>(wealths);

    // printf("Wealths for (%i): %.4f, %.4f\n", i, final_wealths[0], final_wealths[1]);

    fitnesses[i] = compute_fitness(
        individual, intermediate_wealths, final_wealths, intermediate_goals,
        final_goals, instrument_constraints, margin_constraints, n_instruments,
        n_derivatives, n_scenarios, generation);
  }
  return fitnesses;
}

void mutate_individuals(std::vector<double> &selected,
                        const double mutation_rate)
{
  for (int j = 0; j < selected.size(); ++j)
  {
    const double random = Random::get(0.0, 1.0);
    if (random < mutation_rate)
    {
      const double d = 1.0; //((double)j + 1.0) / ((double)j + 3.0);
      const double mutation = d * Random::get(-0.1, 0.1);

      selected[j] = std::min(1.0, std::max(0.0, selected[j] + mutation));
    }
  }
}

void crossover_individuals_scenario(std::vector<double> &selected,
                                    const int n_instruments,
                                    const int n_scenarios,
                                    const double crossover_rate)
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

void crossover_individuals(std::vector<double> &selected,
                           const double crossover_rate)
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
  const double initial_funding_ratio = options.initial_funding_ratio;
  const double target_funding_ratio = options.target_funding_ratio;

  std::vector<double> transaction_costs = {
      options.transaction_costs.domestic_equity,
      options.transaction_costs.global_equity,
      options.transaction_costs.real_estate,
      options.transaction_costs.alternative,
      options.transaction_costs.credit,
      options.transaction_costs.bonds_2y,
      options.transaction_costs.bonds_5y,
      options.transaction_costs.bonds_20y,
      options.transaction_costs.cash,
      options.transaction_costs.domestic_equity_future,
      options.transaction_costs.interest_rate_swap_2y,
      options.transaction_costs.interest_rate_swap_5y,
      options.transaction_costs.interest_rate_swap_20y};

  std::vector<double> instrument_constraints =
      parse_instrument_constraints(options.instrument_constraints);
  std::vector<double> margin_constraints =
      parse_margin_constraints(options.margin_constraints);

  const int n_generic_risks = N_GENERIC_RISKS;
  const int n_forward_rate_risks = N_FORWARD_RATE_RISKS;
  const int n_risks = N_RISKS;
  const int n_instruments = N_INSTRUMENTS;
  const int n_derivatives = N_DERIVATIVES;
  const int n_trees = N_TREES;
  const int n_pca_components = N_PCA_COMPONENTS;

  const int n_scenarios = pow(2.0, n_steps) - 1;
  const int n_genes = n_scenarios * n_instruments;

  vector<double> initial_generic_risk_values = {1.0, 1.0, 1.0, 1.0, 1.0};
  vector<double> initial_forward_rate_risk_values = {0.017682734852308, 0.022014527395177, 0.03055152080786, 0.034314728892416, 0.037618658501593, 0.035941958895264, 0.037028231412524, 0.036388694178369, 0.034014781540592, 0.030376299796611, 0.032005189961554, 0.029276283475554};

  vector<double> generic_risk_means = {0.5, -0.5, -0.5, -0.5, -0.5};
  vector<double> generic_risk_stds = {0.1, 0.1, 0.1, 0.1, 0.1};

  vector<double> pca_forward_rate_risk_eigenvalues = {0.223282705428122, 0.082751342396589, 0.048454437715918, 0.040635667202019, 0.02539470506845, 0.023281350019023, 0.019433965208767, 0.017605370378461, 0.016953588708945, 0.013347508862663, 0.01194883154095, 0.011169415442773};
  vector<double> pca_forward_rate_risk_eigenvectors = {
      0.17806310748978, 0.336052289069525, 0.352104530046595, 0.389724410069539, 0.349021453294063, 0.315252232125525, 0.26285833324634, 0.249608295798516, 0.256096107818295, 0.254773941504214, 0.224079066124409, 0.214929099673031,
      -0.596316760839539, -0.527975500405015, -0.251583154671794, -0.041465903559625, 0.094526967296751, 0.170116753835499, 0.184378018804841, 0.242160135472741, 0.226630505320397, 0.266462268435327, 0.180848288798239, 0.122691541379718,
      0.486340482637313, 0.028187296071813, -0.352705203001156, -0.309926291591076, -0.339610973326499, -0.119838777109285, -0.022143933180858, 0.094019501839615, 0.174841594801978, 0.354831655605287, 0.301300466007176, 0.394891277100978,
      0.555751828998025, -0.427628181205936, -0.25046231532465, 0.107727039741423, 0.191438014022832, 0.25543099259259, 0.120487281923441, 0.17635438532532, 0.046956136603755, 0.06384449363972, -0.292724115760584, -0.440973458515732,
      -0.192920659287202, 0.317837716410373, 0.11595703534322, -0.13446138242021, -0.355170758402066, 0.077751328915102, -0.225727533137616, 0.458589514159322, 0.166812881987041, 0.381386846616376, -0.431471044620772, -0.278094197420766,
      -0.000355596490075, -0.185988898171873, 0.022073459160166, 0.374034215628392, 0.070527463087904, -0.032454523259506, -0.327600397301577, -0.387403848946857, -0.200809131193615, 0.546817723715433, -0.349769215176918, 0.316094219500701,
      0.139006029300245, -0.430508228855957, 0.65231454940463, -0.231715691776853, -0.138560093344234, 0.154092078890021, -0.248061539720448, 0.182960002777793, -0.306348730270312, 0.067706766490438, 0.282215591764525, -0.010099810588702,
      0.00709041088476, -0.073091426233262, 0.32415209903802, -0.283900639827112, -0.205004584886281, 0.113736921858844, 0.447229551572455, -0.579390190608922, 0.377497309651679, 0.188318259633643, -0.13217425799885, -0.151079319248406,
      -0.062248423020668, 0.145605494122674, -0.131066410849474, 0.104650660594197, 0.095972375661584, -0.066946768093708, -0.267587511061236, -0.244911350780297, 0.013799056699233, 0.3466385918266, 0.558978675125192, -0.607215703088809,
      0.005347660322724, 0.053862754606072, 0.095286425069807, -0.352006363050202, 0.469552689329943, -0.551000579719726, 0.328797893389751, 0.153227752248187, -0.294744351269166, 0.325365087966101, -0.108894151744519, -0.061801688839778,
      0.060581508643696, -0.101683553907503, 0.127693442205872, -0.242126986785104, 0.439325113241745, -0.182042215299357, -0.514870515443296, -0.028886315075429, 0.61981228537522, -0.144573079389127, -0.087306912690354, 0.079348038577479,
      0.056103920008638, -0.252816668074894, 0.196240506077734, 0.5029683197517, -0.319720817495426, -0.638018884962708, 0.131438908544935, 0.158951119033581, 0.273024562645784, -0.051680006452083, 0.055263517107068, -0.096696332791446};

  vector<double> correlations = {
      1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};

  vector<double> sigmas = {
      1.1, 1.1, 1.1, 1.1, 1.1, 1.7};
  vector<double> rhos = {
      -0.15, -0.15, -0.15, -0.15, -0.15, -0.1};

  vector<double> zero_coupon_tenors = {
      1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20};

  // Will contain the final result
  double best_fitness = 0.0;
  std::vector<double> best_individual = std::vector<double>(n_genes, 0.0);

  // Generate initial individuals
  std::vector<double> individuals =
      initialize_individuals(n_individuals, n_instruments, n_scenarios);

  // ...and normalize them.
  normalize_individuals(individuals, n_individuals, n_instruments,
                        n_derivatives, n_scenarios);

  // Compute gammas
  const int n_gamma_trials = 100;
  tuple<vector<double>, vector<double>> gammas = compute_gammas(
      generic_risk_stds,
      pca_forward_rate_risk_eigenvalues,
      pca_forward_rate_risk_eigenvectors,
      sigmas,
      rhos,
      correlations,
      n_gamma_trials,
      n_pca_components,
      n_generic_risks,
      n_forward_rate_risks,
      n_scenarios,
      n_steps);

  printf("Generating trees... ");
  std::vector<ScenarioTree> scenario_trees(n_trees);
  for (int i = 0; i < n_trees; ++i)
  {
    ScenarioTree tree;
    tree.instrument_changes =
        generate_state_changes(
            initial_generic_risk_values,
            initial_forward_rate_risk_values,
            generic_risk_means,
            generic_risk_stds,
            pca_forward_rate_risk_eigenvalues,
            pca_forward_rate_risk_eigenvectors,
            sigmas,
            rhos,
            correlations,
            zero_coupon_tenors,
            gammas,
            n_instruments,
            n_generic_risks,
            n_forward_rate_risks,
            n_pca_components,
            n_scenarios,
            n_steps);

    std::tuple<std::vector<double>, std::vector<double>> goals =
        generate_goals(tree.instrument_changes, initial_funding_ratio,
                       target_funding_ratio, n_scenarios, n_instruments);

    tree.intermediate_goals = std::get<0>(goals);
    tree.final_goals = std::get<1>(goals);

    scenario_trees[i] = tree;
  }
  printf("Done!\n");

  for (int t = 0; t < n_generations; ++t)
  {

    std::vector<double> fitnesses(n_individuals);
    double generation_best_fitness = 0.0;

    for (int s = 0; s < n_trees; ++s)
    {
      ScenarioTree tree = scenario_trees[s];
      std::vector<double> temp_fitnesses = compute_fitnesses(
          individuals, tree.instrument_changes, transaction_costs, tree.intermediate_goals,
          tree.final_goals, instrument_constraints, margin_constraints, n_individuals,
          n_instruments, n_derivatives, n_scenarios, t);
      for (int i = 0; i < n_individuals; ++i)
      {
        fitnesses[i] += temp_fitnesses[i];
      }
    }

    // Compute average fitnesses
    for (int i = 0; i < n_individuals; ++i)
    {
      fitnesses[i] = fitnesses[i] / n_trees;
    }

    // Print individuals
    printf("Individuals in generation %i\n", t);
    for (int i = 0; i < n_individuals; ++i)
    {
      for (int j = 0; j < n_scenarios; ++j)
      {
        const int ix = i * n_genes + j * n_instruments;
        printf("(%.2f, %.2f) ", individuals[ix + DOMESTIC_EQUITY_INDEX], individuals[ix + DOMESTIC_EQUITY_FUTURE_INDEX]);
      }
      printf("\n");
    }
    printf("\n");

    // Check global best fitness
    for (int i = 0; i < n_individuals; ++i)
    {
      if (fitnesses[i] > best_fitness || best_fitness == 0.0)
      {
        // printf("%i %.4f %.4f (%.32f) \n", t, best_fitness, fitnesses[i], (fitnesses[i] / best_fitness) - 1.0);

        const int ix = i * n_genes;
        best_fitness = fitnesses[i];
        for (int j = 0; j < n_genes; ++j)
        {
          best_individual[j] = individuals[ix + j];
        }
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

      // printf("Selecting %i and %i\n", ix1 / n_genes, ix2 / n_genes);

      std::vector<double> selected(2 * n_genes);
      for (int j = 0; j < n_genes; ++j)
      {
        selected[j] = individuals[ix1 + j];
        selected[n_genes + j] = individuals[ix2 + j];
      }

      // Crossover
      crossover_individuals_scenario(
          selected,
          n_instruments,
          n_scenarios,
          crossover_rate);

      // Mutation
      mutate_individuals(selected, mutation_rate);

      // Add selected individuals
      for (int j = 0; j < n_genes; ++j)
      {
        offspring[ix + j] = selected[j];
        offspring[ix + n_genes + j] = selected[n_genes + j];
      }
    }

    // Elitism
    for (int i = 0; i < n_elitism_copies; ++i)
    {
      int ix = i * n_genes;
      for (int j = 0; j < n_genes; ++j)
      {
        offspring[ix + j] = best_individual[j];
      }
    }

    normalize_individuals(offspring, n_individuals, n_instruments,
                          n_derivatives, n_scenarios);

    // Replace generation
    individuals = offspring;
  }

  Result result;

  double average_expected_return = 0.0;
  double average_expected_risk = 0.0;

  for (int s = 0; s < n_trees; ++s)
  {
    ScenarioTree t = scenario_trees[s];
    std::tuple<std::vector<double>, std::vector<double>> best_wealths =
        compute_wealths(best_individual, t.instrument_changes, transaction_costs,
                        n_instruments, n_derivatives, n_scenarios);

    std::vector<double> best_intermediate_wealths = std::get<0>(best_wealths);
    std::vector<double> best_final_wealths = std::get<1>(best_wealths);
    average_expected_return +=
        compute_expected_wealth(best_final_wealths) - 1.0;
    average_expected_risk +=
        compute_expected_risk(best_intermediate_wealths, best_final_wealths, t.intermediate_goals, t.final_goals);
  }

  average_expected_return /= n_trees;
  average_expected_risk /= n_trees;

  result.fitness = best_fitness;
  result.individual = best_individual;
  result.expected_return = average_expected_return;
  result.expected_risk = average_expected_risk;

  /*
  int i = 0;
  printf("Best: %.2f%%\n", 100 * average_expected_return);
  for (auto w : best_individual)
  {
    i = (i + 1) % n_instruments;
    printf("%.2f ", w);
    if (i == 0)
    {
      printf("\n");
    }
  }
  printf("\n");
  */

  return result;
}