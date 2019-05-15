#include <vector>

#include <include/constants.h>

using namespace std;

double sample_domestic_equity(std::vector<double> &risk_values)
{
  return risk_values[DOMESTIC_MARKET_RISK_INDEX];
}

double sample_global_equity(std::vector<double> &risk_values)
{
  return risk_values[GLOBAL_MARKET_RISK_INDEX];
}

double sample_real_estate(std::vector<double> &risk_values)
{
  return risk_values[REAL_ESTATE_RISK_INDEX];
}

double sample_alternative(std::vector<double> &risk_values)
{
  return risk_values[ALTERNATIVE_RISK_INDEX];
}

double sample_credit(std::vector<double> &risk_values)
{
  return risk_values[CREDIT_RISK_INDEX];
}

double sample_bonds_2y(std::vector<double> &risk_values)
{
  double next_price = 0.0;
  return next_price;
}

double sample_bonds_5y(std::vector<double> &risk_values)
{
  return 0.0;
}

double sample_bonds_20y(std::vector<double> &risk_values)
{
  return 0.0;
}

double sample_cash(std::vector<double> &risk_values)
{
  return 0.0;
}

double sample_domestic_equity_future(std::vector<double> &risk_values)
{
  return risk_values[DOMESTIC_MARKET_RISK_INDEX];
}

double sample_interest_rate_swap_2y(std::vector<double> &risk_values)
{
  return 0.0;
}

double sample_interest_rate_swap_5y(std::vector<double> &risk_values)
{
  return 0.0;
}

double sample_interest_rate_swap_20y(std::vector<double> &risk_values)
{
  return 0.0;
}

double generate_instrument_change(
    const int instrument_index,
    std::vector<double> &risk_values)
{
  if (instrument_index == DOMESTIC_EQUITY_INDEX)
  {
    return sample_domestic_equity(risk_values);
  }
  else if (instrument_index == GLOBAL_EQUITY_INDEX)
  {
    return sample_global_equity(risk_values);
  }
  else if (instrument_index == REAL_ESTATE_INDEX)
  {
    return sample_real_estate(risk_values);
  }
  else if (instrument_index == ALTERNATIVE_INDEX)
  {
    return sample_alternative(risk_values);
  }
  else if (instrument_index == CREDIT_INDEX)
  {
    return sample_credit(risk_values);
  }
  else if (instrument_index == BONDS_2Y_INDEX)
  {
    return sample_bonds_2y(risk_values);
  }
  else if (instrument_index == BONDS_5Y_INDEX)
  {
    return sample_bonds_5y(risk_values);
  }
  else if (instrument_index == BONDS_20Y_INDEX)
  {
    return sample_bonds_20y(risk_values);
  }
  else if (instrument_index == CASH_INDEX)
  {
    return sample_cash(risk_values);
  }
  else if (instrument_index == DOMESTIC_EQUITY_FUTURE_INDEX)
  {
    return sample_domestic_equity_future(risk_values);
  }
  else if (instrument_index == INTEREST_RATE_SWAP_2Y_INDEX)
  {
    return sample_interest_rate_swap_2y(risk_values);
  }
  else if (instrument_index == INTEREST_RATE_SWAP_5Y_INDEX)
  {
    return sample_interest_rate_swap_5y(risk_values);
  }
  else if (instrument_index == INTEREST_RATE_SWAP_20Y_INDEX)
  {
    return sample_interest_rate_swap_20y(risk_values);
  }
  else
  {
    // This should never occur, but currently Embind does
    // not allow us to throw a proper exception.
    return 0.0;
  }
}

std::vector<double> generate_instrument_changes(
    std::vector<double> &risk_changes,
    const int n_instruments)
{
  std::vector<double> instrument_changes(n_instruments);

  for (int i = 0; i < n_instruments; ++i)
  {
    instrument_changes[i] = generate_instrument_change(i, risk_changes);
  }

  return instrument_changes;
}