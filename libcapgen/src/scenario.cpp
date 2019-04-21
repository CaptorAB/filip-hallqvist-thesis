#include <vector>
#include <tuple>
#include <cmath>

#include <lib/random.h>

#include <include/constants.h>
#include <include/util.h>
#include <include/scenario.h>

using Random = effolkronium::random_static;

// Util

static double get_random_normal(double mu, double sigma)
{
  return Random::get<std::normal_distribution<>>(mu, sigma);
}

static double get_random_nln(double n1, double n2, double rho, double sigma)
{
  return sigma * (rho * n1 + (sqrt(1 - pow(rho, 2))) * n2);
}

static double sample_black_process(double n1, double n2, double forward_price, double gamma, double rho, double sigma, double mean, double variance)
{
  double u = get_random_nln(n1, n2, rho, sigma);
  double u_norm = (u - mean) / variance;
  double epsilon = variance * u_norm; // TODO: Why immediately cancel the variance?
  return forward_price * exp(-gamma + epsilon);
}

// Risk sampling functions
static double sample_domestic_market_risk(double previous_change, double n1, double n2)
{
  return get_random_normal(0.5, 0.0); // TODO: Use black
}

static double sample_global_market_risk(double previous_change, double n1, double n2)
{
  return get_random_normal(-0.5, 0.0); // TODO: Use black
}

static double sample_alternative_risk(double previous_change, double n1, double n2)
{
  return get_random_normal(-0.5, 0.0); // TODO: Use black
}

static double sample_interest_rate_2y_risk(double previous_change, double n1, double n2)
{
  return get_random_normal(-0.5, 0.0); // TODO: Use black
}

static double sample_interest_rate_5y_risk(double previous_change, double n1, double n2)
{
  return get_random_normal(-0.5, 0.0); // TODO: Use black
}

static double sample_interest_rate_20y_risk(double previous_change, double n1, double n2)
{
  return get_random_normal(-0.5, 0.0); // TODO: Use black
}

static double sample_credit_risk(double previous_change, double n1, double n2)
{
  return get_random_normal(-0.5, 0.0); // TODO: Use black
}

static double sample_cash_risk(double previous_change, double n1, double n2)
{
  return get_random_normal(-0.5, 0.0); // TODO: Use black
}

// Instrument sampling functions

static double sample_domestic_equity(double previous_change, std::vector<double> &risk_changes)
{
  return risk_changes[DOMESTIC_MARKET_RISK_INDEX];
}

static double sample_global_equity(double previous_change, std::vector<double> &risk_changes)
{
  return risk_changes[GLOBAL_MARKET_RISK_INDEX];
}

static double sample_real_estate(double previous_change, std::vector<double> &risk_changes)
{
  return risk_changes[REAL_ESTATE_INDEX];
}

static double sample_alternative(double previous_change, std::vector<double> &risk_changes)
{
  return risk_changes[ALTERNATIVE_RISK_INDEX];
}

static double sample_credit(double previous_change, std::vector<double> &risk_changes)
{
  return risk_changes[CREDIT_RISK_INDEX];
}

static double sample_bonds_2y(double previous_change, std::vector<double> &risk_changes)
{
  return risk_changes[INTEREST_RATE_2Y_RISK_INDEX];
}

static double sample_bonds_5y(double previous_change, std::vector<double> &risk_changes)
{
  return risk_changes[INTEREST_RATE_5Y_RISK_INDEX];
}

static double sample_cash(double previous_change, std::vector<double> &risk_changes)
{
  return risk_changes[CASH_RISK_INDEX];
}

static double sample_fta(double previous_change, std::vector<double> &risk_changes)
{
  return risk_changes[INTEREST_RATE_20Y_RISK_INDEX];
}

static double sample_domestic_equity_future(double previous_change, std::vector<double> &risk_changes)
{
  return risk_changes[DOMESTIC_MARKET_RISK_INDEX];
}

static double sample_interest_rate_swap_2y(double previous_change, std::vector<double> &risk_changes)
{
  return risk_changes[INTEREST_RATE_2Y_RISK_INDEX];
}

static double sample_interest_rate_swap_5y(double previous_change, std::vector<double> &risk_changes)
{
  return risk_changes[INTEREST_RATE_5Y_RISK_INDEX];
}

static double sample_interest_rate_swap_20y(double previous_change, std::vector<double> &risk_changes)
{
  return risk_changes[INTEREST_RATE_20Y_RISK_INDEX];
}

// Goals

std::vector<double> generate_goals(std::vector<double> &price_changes, int n_steps, int n_scenarios, int n_instruments, double initial_funding_ratio, double target_funding_ratio)
{
  std::vector<double> values(n_scenarios);
  std::vector<double> goals(n_scenarios);

  values[0] = 1.0f / initial_funding_ratio;
  goals[0] = values[0] * target_funding_ratio;

  for (int t = 1; t < n_steps; ++t)
  {
    int first_node_in_level = get_first_node_in_level(t);
    int i_first_node_in_step = first_node_in_level * n_instruments; // Index of first node in step
    int n_scenarios_in_step = get_nodes_in_level(t);                // Scenarios in this step
    for (int s = 0; s < n_scenarios_in_step; ++s)
    {
      int si = i_first_node_in_step + (s * n_instruments); // Index of current scenario
      int si_ = get_parent_index(si, n_instruments);       // Index of previous scenario

      int ri = first_node_in_level + s;  // Index of current scenario goal
      int ri_ = get_parent_index(ri, 1); // Index of previous scenario goal

      // TODO: Implement the interest rate swap
      // goals[ri] = goals[ri_] * price_changes[si + INTEREST_RATE_SWAP_20Y_INDEX] * target_funding_ratio;
      values[ri] = values[ri_] * 1.0;
      goals[ri] = values[ri] * target_funding_ratio;
    }
  }
  return goals;
}

// TODO: Refactor this mess
std::tuple<std::vector<double>, std::vector<double>>
generate_scenarios(int n_steps)
{
  int n_scenarios = (1 << n_steps) - 1;
  int n_risks = N_RISKS;
  int n_instruments = N_INSTRUMENTS;
  int n_risk_changes = n_risks * n_scenarios;
  int n_instrument_changes = n_instruments * n_scenarios;

  std::vector<double> risk_changes(n_risk_changes);
  std::vector<double> instrument_changes(n_instrument_changes);
  std::vector<double> probabilities(n_scenarios);

  // Generate new scenarios
  for (int t = 0; t < n_steps; ++t)
  {
    int first_node_in_level = get_first_node_in_level(t);
    int i_first_node_in_step = first_node_in_level * n_instruments; // Index of first node in step
    int n_scenarios_in_step = get_nodes_in_level(t);                // Scenarios in this step

    for (int s = 0; s < n_scenarios_in_step; ++s)
    {                                                      // Iterate over scenarios in step
      int si = i_first_node_in_step + (s * n_instruments); // Index of current scenario
      int si_ = get_parent_index(si, n_instruments);       // Index of previous scenario

      int ri = first_node_in_level + s;  // Index of current scenario risk map
      int ri_ = get_parent_index(ri, 1); // Index of previous scenario risk map

      probabilities[ri] = 1.0 / n_scenarios_in_step;

      // Generate new normals
      std::vector<double> normals(2 * n_risks);
      for (int i = 0; i < n_risks; ++i)
      {
        int ix = 2 * i;
        normals[ix] = get_random_normal(0.0, 1.0);
        normals[ix + 1] = get_random_normal(0.0, 1.0);
      }

      // TODO: Correlate normals

      // Update risks
      // TODO: Break out into own function
      std::vector<double> temp_risk_changes(n_risks);
      for (int i = 0; i < n_risks; ++i)
      {
        int nx = 2 * i;
        int rx = ri + i;
        double n1 = normals[nx];
        double n2 = normals[nx + 1];

        double old_change = risk_changes[ri_]; // TODO: To be used in new_change
        double new_change;

        if (i == DOMESTIC_MARKET_RISK_INDEX)
        {
          new_change = sample_domestic_market_risk(old_change, n1, n2);
        }
        else if (i == GLOBAL_MARKET_RISK_INDEX)
        {
          new_change = sample_global_market_risk(old_change, n1, n2);
        }
        else if (i == ALTERNATIVE_RISK_INDEX)
        {
          new_change = sample_alternative_risk(old_change, n1, n2);
        }
        else if (i == INTEREST_RATE_2Y_RISK_INDEX)
        {
          new_change = sample_interest_rate_2y_risk(old_change, n1, n2);
        }
        else if (i == INTEREST_RATE_5Y_RISK_INDEX)
        {
          new_change = sample_interest_rate_5y_risk(old_change, n1, n2);
        }
        else if (i == INTEREST_RATE_20Y_RISK_INDEX)
        {
          new_change = sample_interest_rate_20y_risk(old_change, n1, n2);
        }
        else if (i == CREDIT_RISK_INDEX)
        {
          new_change = sample_credit_risk(old_change, n1, n2);
        }
        else if (i == CASH_RISK_INDEX)
        {
          new_change = sample_cash_risk(old_change, n1, n2);
        }
        else
        {
          throw std::logic_error("not implemented");
        }

        risk_changes[rx] = new_change;
        temp_risk_changes[i] = new_change;
      }

      // Update instruments
      // TODO: Break out into own function
      for (int i = 0; i < n_instruments; ++i)
      {
        int k = si + i;   // Index of current instrument in current scenario
        int k_ = si_ + i; // Index of current instrument in previous scenario

        double old_change = instrument_changes[k_];
        double new_change = 0.0;

        // TODO: Break out into own function
        if (i == DOMESTIC_EQUITY_INDEX)
        {
          new_change = sample_domestic_equity(old_change, temp_risk_changes);
        }
        else if (i == GLOBAL_EQUITY_INDEX)
        {
          new_change = sample_global_equity(old_change, temp_risk_changes);
        }
        else if (i == REAL_ESTATE_INDEX)
        {
          new_change = sample_real_estate(old_change, temp_risk_changes);
        }
        else if (i == ALTERNATIVE_INDEX)
        {
          new_change = sample_alternative(old_change, temp_risk_changes);
        }
        else if (i == CREDIT_INDEX)
        {
          new_change = sample_credit(old_change, temp_risk_changes);
        }
        else if (i == BONDS_2Y_INDEX)
        {
          new_change = sample_bonds_2y(old_change, temp_risk_changes);
        }
        else if (i == BONDS_5Y_INDEX)
        {
          new_change = sample_bonds_5y(old_change, temp_risk_changes);
        }
        else if (i == CASH_INDEX)
        {
          new_change = sample_cash(old_change, temp_risk_changes);
        }
        else if (i == FTA_INDEX)
        {
          new_change = sample_fta(old_change, temp_risk_changes);
        }
        else if (i == DOMESTIC_EQUITY_FUTURE_INDEX)
        {
          new_change = sample_domestic_equity_future(old_change, temp_risk_changes);
        }
        else if (i == INTEREST_RATE_SWAP_2Y_INDEX)
        {
          new_change = sample_interest_rate_swap_2y(old_change, temp_risk_changes);
        }
        else if (i == INTEREST_RATE_SWAP_5Y_INDEX)
        {
          new_change = sample_interest_rate_swap_5y(old_change, temp_risk_changes);
        }
        else if (i == INTEREST_RATE_SWAP_20Y_INDEX)
        {
          new_change = sample_interest_rate_swap_20y(old_change, temp_risk_changes);
        }

        instrument_changes[k] = new_change;
      }
    }
  }

  return std::make_tuple(instrument_changes, probabilities);
}