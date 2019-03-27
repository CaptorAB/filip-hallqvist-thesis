#ifndef SCENARIO_H
#define SCENARIO_H

#include <vector>
#include <tuple>
#include <map>
#include <string>

#include <lib/random.h>

double get_random_normal(double mu, double sigma);
double get_random_nln(double n1, double n2, double rho, double sigma);
double sample_black_process(double n1, double n2, double forward_price, double gamma, double rho, double sigma, double mean, double variance);

enum class RiskType
{
  DomesticMarketRisk,
  GlobalMarketRisk,
  AlternativeRisk,
  InterestRateRisk,
  CreditRisk,
  CashRisk,
};

class RiskProcess
{
protected:
  RiskProcess(const double initial) : current(initial) {}
  double current;

public:
  virtual void update(double n1, double n2) = 0;
  double get_current();
};

class DomesticMarketRiskProcess : public RiskProcess
{
private:
  double current = 0.0;
  const double forward_price = 1.0;
  const double gamma = 1.0;
  const double rho = -0.15;
  const double sigma = 1.1;
  const double mean = 1.0;
  const double variance = 1.0;

public:
  DomesticMarketRiskProcess() : RiskProcess(0.0){};
  void update(double n1, double n2);
};

typedef std::map<RiskType, RiskProcess *> risks_t;

risks_t create_default_risks();

class Instrument
{
protected:
  Instrument(const double initial) : current(initial) {}
  double current;

public:
  virtual void update(risks_t risks) = 0;
  double get_current();
};

class DomesticEquityInstrument : public Instrument
{

public:
  DomesticEquityInstrument() : Instrument(1.0){};
  void update(risks_t risks);
};

typedef std::vector<Instrument *> instruments_t;
typedef std::vector<double> correlations_t;

double generate_scenario(instruments_t instruments, risks_t risks, correlations_t correlations);

#endif