#include <vector>
#include <tuple>
#include <map>
#include <string>

#include <lib/random.h>
#include <include/util.h>
#include <include/scenario.h>

using Random = effolkronium::random_static;

/*
  Equities = BlackEq
  Global equities = BlackEq * Normal
  Real Estate = BlackIr
  Alternative = BlackEq
  Credit = BlackIr * Normal
  2y Bonds = BlackIr
  5y Bonds = BlackIr
  Cash = BlackIr?
  FTA = BlackIr
*/

// Instrument indices
static const int DOMESTIC_EQUITY_INDEX = 0;
static const int GLOBAL_EQUITY_INDEX = 1;
static const int REAL_ESTATE_INDEX = 2;
static const int ALTERNATIVE_INDEX = 3;
static const int CREDIT_INDEX = 4;
static const int BONDS_2Y_INDEX = 5;
static const int BONDS_5Y_INDEX = 6;
static const int CASH_INDEX = 7;
static const int FTA_INDEX = 8;
static const int DOMESTIC_EQUITY_FUTURE_INDEX = 9;
static const int INTEREST_RATE_SWAP_2Y_INDEX = 10;
static const int INTEREST_RATE_SWAP_5Y_INDEX = 11;
static const int INTEREST_RATE_SWAP_20Y_INDEX = 12;

// Risk indices
static const int DOMESTIC_MARKET_RISK_INDEX = 0;
static const int GLOBAL_MARKET_RISK_INDEX = 1;
static const int ALTERNATIVE_RISK_INDEX = 2;
static const int INTEREST_RATE_RISK_INDEX = 3;
static const int CREDIT_RISK_INDEX = 4;
static const int CASH_RISK_INDEX = 5;

static const std::vector<double> RISK_CORRELATIONS = {
    1.0, 0.0,
    0.0, 1.0};

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

static double sample_normal_process(double mu, double sigma)
{
  return get_random_normal(mu, sigma);
}

// Risk sampling functions
static double sample_domestic_market_risk(double previous_change, double n1, double n2)
{
  return get_random_normal(0.5, 0.1);
}

static double sample_global_market_risk(double previous_change, double n1, double n2)
{
  return get_random_normal(-0.5, 0.1);
}

// Instrument sampling functions
static double sample_domestic_equity(double previous_change, std::vector<double> risk_changes)
{
  return risk_changes[DOMESTIC_MARKET_RISK_INDEX];
}

static double sample_global_equity(double previous_change, std::vector<double> risk_changes)
{
  return risk_changes[GLOBAL_MARKET_RISK_INDEX];
}

// TODO: Refactor this mess
static std::vector<double> generate_price_changes(int n_steps)
{
  int n_scenarios = (1 << n_steps) - 1;
  int n_risks = 5;
  int n_instruments = 2;
  int n_risk_changes = n_risks * n_scenarios;
  int n_instrument_changes = n_instruments * n_scenarios;

  std::vector<double> risk_changes(n_risk_changes);
  std::vector<double> instrument_changes(n_instrument_changes);

  // Set initial risk changes
  for (int i = 0; i < n_risks; ++i)
  {
    risk_changes[i] = 0.0;
  }

  // Set initial instrument changes
  for (int i = 0; i < n_instruments; ++i)
  {
    instrument_changes[i] = 0.0;
  }

  // Generate new scenarios
  for (int t = 1; t < n_steps; ++t)
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
        int rx_ = ri_ + i;
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

        if (i == DOMESTIC_EQUITY_INDEX)
        {
          instrument_changes[k] = sample_domestic_equity(old_change, temp_risk_changes);
        }
        else if (i == GLOBAL_EQUITY_INDEX)
        {
          instrument_changes[k] = sample_global_equity(old_change, temp_risk_changes);
        }
      }
    }
  }

  return instrument_changes;
}