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

// Util

double get_random_normal(double mu, double sigma)
{
  return Random::get<std::normal_distribution<>>(mu, sigma);
}

double get_random_nln(double n1, double n2, double rho, double sigma)
{
  return sigma * (rho * n1 + (sqrt(1 - pow(rho, 2))) * n2);
}

double sample_black_process(double n1, double n2, double forward_price, double gamma, double rho, double sigma, double mean, double variance)
{
  double u = get_random_nln(n1, n2, rho, sigma);
  double u_norm = (u - mean) / variance;
  double epsilon = variance * u_norm; // TODO: Why immediately cancel the variance?
  return forward_price * exp(-gamma + epsilon);
}

risks_t create_default_risks()
{
  risks_t risks;
  risks[RiskType::DomesticMarket] = new DomesticMarketRiskProcess();
  return risks;
}

instruments_t create_default_instruments()
{
  instruments_t instruments = {new DomesticEquityInstrument(), new DomesticEquityInstrument()};
  return instruments;
}

// Risks

double RiskProcess::get_current()
{
  return current;
};

void DomesticMarketRiskProcess::update(double n1, double n2)
{
  current = current * (1.0 + sample_black_process(n1, n2, forward_price, gamma, rho, sigma, mean, variance));
}

// Instruments
double Instrument::get_current()
{
  return current;
};

void DomesticEquityInstrument::update(risks_t risks)
{
  current = risks[RiskType::DomesticMarket]->get_current();
}

// Generate scenario

std::vector<double>
generate_price_change(instruments_t instruments, risks_t risks, correlations_t correlations)
{
  std::vector<double> price_changes; // TODO: Define length

  // Sample normals
  std::map<RiskType, std::tuple<double, double>> normals;
  for (auto const &risk : risks)
  {
    double n1 = get_random_normal(0.0, 1.0);
    double n2 = get_random_normal(0.0, 1.0);
    normals[risk.first] = std::make_tuple(n1, n2);
  }

  // TODO: Correlate normals

  // Evaluate risk processes
  for (auto &risk : risks)
  {
    double n1 = std::get<0>(normals[risk.first]);
    double n2 = std::get<1>(normals[risk.first]);
    risk.second->update(n1, n2); // TODO: Copies element each time?
  }

  // Evaluate instruments
  for (auto &instrument : instruments)
  {
    instrument->update(risks);
    price_changes.push_back(instrument->get_current());
  }

  return price_changes;
};

// TODO: Refactor this mess
std::vector<double> generate_price_changes(int n_steps, instruments_t instruments, risks_t risks, correlations_t correlations)
{
  int n_scenarios = (1 << n_steps) - 1;
  int n_risks = risks.size();
  int n_instruments = instruments.size();
  int n_instrument_changes = n_instruments * n_scenarios;

  std::vector<std::map<RiskType, double>> risk_changes(n_scenarios);
  std::vector<double> instrument_changes(n_instrument_changes);

  // Set initial risk values
  std::map<RiskType, double> changes_map;
  risk_changes[0] = changes_map;
  for (auto const &risk : risks)
  {
    risk_changes[0][risk.first] = 1.0;
  }

  // Set starting prices
  for (int i = 0; i < n_instruments; ++i)
  {
    instrument_changes[i] = 1.0;
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
      std::map<RiskType, std::tuple<double, double>> normals;
      for (auto const &risk : risks)
      {
        double n1 = get_random_normal(0.0, 1.0);
        double n2 = get_random_normal(0.0, 1.0);
        normals[risk.first] = std::make_tuple(n1, n2);
      }

      // TODO: Correlate normals

      // Update risks
      std::map<RiskType, double> changes_map;
      risk_changes[ri] = changes_map;
      for (auto const &risk : risks)
      {
        double n1 = std::get<0>(normals[risk.first]);
        double n2 = std::get<1>(normals[risk.first]);

        double old_change = risk_changes[ri_][risk.first]; // TODO: To be used in new_change
        double new_change = sample_black_process(n1, n2, 0.25, 1.0, 1.0, 1.0, 0.0, 1.0);

        risk_changes[ri][risk.first] = new_change;
      }

      // Update instruments
      for (int i = 0; i < n_instruments; ++i)
      {
        int k = si + i;   // Index of current instrument in current scenario
        int k_ = si_ + i; // Index of current instrument in previous scenario

        if (i == 0)
        {
          instrument_changes[k] = risk_changes[ri][RiskType::DomesticMarket];
        }
        else
        {
          instrument_changes[k] = (1.0 - risk_changes[ri][RiskType::DomesticMarket]);
        }
      }
    }
  }

  return instrument_changes;
}