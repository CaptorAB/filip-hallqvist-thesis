#include <vector>
#include <tuple>
#include <map>
#include <string>

#include <lib/random.h>
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

/*
enum class Instrument
{
  DomesticEquity,
  GlobalEquity,
  RealEstate,
  Alternative,
  InterestRate2Y,
  InterestRate5Y,
  InterestRate20Y,
  Credit,
  Cash,
  DomesticEquityFuture,
  InterestRateSwap2Y,
  InterestRateSwap5Y,
  InterestRateSwap20Y,
};
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
  risks[RiskType::DomesticMarketRisk] = new DomesticMarketRiskProcess();
  return risks;
}

// Risks

double RiskProcess::get_current()
{
  return current;
};

void DomesticMarketRiskProcess::update(double n1, double n2)
{
  current = sample_black_process(n1, n2, forward_price, gamma, rho, sigma, mean, variance);
}

// Instruments
double Instrument::get_current()
{
  return current;
};

void DomesticEquityInstrument::update(risks_t risks)
{
  current = current * risks[RiskType::DomesticMarketRisk]->get_current();
}

// Generate scenario

double
generate_scenario(instruments_t instruments, risks_t risks, correlations_t correlations)
{
  // std::map<Instrument, double> instrument_values;
  // std::map<Risk, double> risk_values;

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
  }

  return 0.0;
};