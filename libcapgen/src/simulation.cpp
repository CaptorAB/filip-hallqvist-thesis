#include <vector>
#include <tuple>

#include <lib/stats/stats.hpp>

#include <include/constants.h>
#include <include/cholesky.h>
#include <include/simulation.h>
#include <include/bootstrapping.h>

using namespace std;
using namespace stats;

vector<double>
sample_uniform_randoms(
    const int n_randoms)
{
  vector<double> r(n_randoms);
  for (int i = 0; i < n_randoms; ++i)
  {
    r[i] = runif(0.0, 1.0);
  }
  return r;
}

vector<double>
normalize_uniform_randoms(
    vector<double> &uniforms)
{
  vector<double> normalized(uniforms.size());
  for (int i = 0; i < uniforms.size(); ++i)
  {
    normalized[i] = qnorm(uniforms[i]);
  }
  return normalized;
}

double
sample_nln(
    const double n1,
    const double n2,
    const double sigma,
    const double rho)
{
  const double phi = n1;
  const double eta = sigma * (rho * n1 + sqrt(1.0 - pow(rho, 2.0)) * n2);

  return exp(0.5 * eta) * phi;
}

double
standardize_nln(
    const double u,
    const double sigma,
    const double rho)
{
  const double mu = 0.5 * rho * sigma * exp(0.125 * pow(sigma, 2.0));
  const double var = exp(0.5 * pow(sigma, 2.0)) * (1.0 + pow(rho, 2.0) * pow(sigma, 2.0) * (1.0 - 0.25 * exp(-0.25 * pow(sigma, 2.0))));
  return (u - mu) / sqrt(var);
}

vector<double>
correlate_risks(
    vector<double> &normals,
    vector<double> &correlations,
    const int n_generic_risks,
    const int n_pca_components)
{
  const int n_correlations = n_generic_risks + n_pca_components;
  const int n_randoms = 2 * n_correlations;

  vector<double> cholesky = compute_cholesky(correlations, n_correlations);
  vector<double> correlated(n_randoms, 0.0);

  for (int i = 0; i < n_correlations; ++i)
  {
    for (int j = 0; j < n_correlations; ++j)
    {
      const int jx = j * n_correlations;
      correlated[i] += normals[j] * cholesky[jx + i];
      correlated[n_correlations + i] += normals[n_correlations + j] * cholesky[jx + i];
    }
  }

  return correlated;
}

vector<double> compute_yield_curve(
    vector<double> &forward_rates)
{
  vector<double> yield_curve(forward_rates.size());
  yield_curve[0] = forward_rates[0];
  double z = yield_curve[0];
  for (int i = 1; i < forward_rates.size(); ++i)
  {
    yield_curve[i] = pow(pow((1.0 + forward_rates[i]) * (1.0 + z), i - 1.0), 1.0 / i) - 1.0;
    z = yield_curve[i];
  }
  return yield_curve;
}

vector<double> generate_correlated_nlns(
    vector<double> &sigmas,
    vector<double> &rhos,
    vector<double> &correlations,
    const int n_generic_risks,
    const int n_pca_components)
{
  const int n_uniforms = 2 * correlations.size();
  vector<double> uniforms = sample_uniform_randoms(n_uniforms);
  vector<double> normals = normalize_uniform_randoms(uniforms);
  vector<double> normals_correlated = correlate_risks(normals, correlations, n_generic_risks, n_pca_components);

  vector<double> nlns(n_generic_risks + n_pca_components);

  // Sample nlns for market risks
  for (int i = 0; i < n_generic_risks; ++i)
  {
    const double n1 = normals[i];
    const double n2 = normals[(normals.size() / 2) + i];
    const double u = sample_nln(n1, n2, sigmas[i], rhos[i]);
    nlns[i] = standardize_nln(u, sigmas[i], rhos[i]);
  }

  // Sample three nlns for interest rate risks
  for (int i = 0; i < n_pca_components; ++i)
  {
    const double n1 = normals[n_generic_risks + i];
    const double n2 = normals[(normals.size() / 2) + n_generic_risks + i];
    const double u = sample_nln(n1, n2, sigmas[n_generic_risks + 1], rhos[n_generic_risks + 1]);
    nlns[n_generic_risks + i] = standardize_nln(u, sigmas[n_generic_risks + 1], rhos[n_generic_risks + 1]);
  }

  return nlns;
}

tuple<vector<double>, vector<double>> compute_gammas(
    vector<double> &generic_risk_stds,
    vector<double> &pca_forward_rate_risk_eigenvalues,
    vector<double> &pca_forward_rate_risk_eigenvectors,
    vector<double> &sigmas,
    vector<double> &rhos,
    vector<double> &correlations,
    const int n_trials,
    const int n_pca_components,
    const int n_generic_risks,
    const int n_forward_rate_risks,
    const int n_states,
    const int n_steps)
{
  const int n_risks = n_generic_risks + n_forward_rate_risks;

  vector<double> generic_risk_gammas(n_steps * n_generic_risks);
  vector<double> forward_rate_risk_gammas(n_steps * n_forward_rate_risks);

  for (int t = 0; t < n_trials; ++t)
  {
    for (int s = 0; s < n_steps; ++s)
    {
      vector<double> nlns = generate_correlated_nlns(
          sigmas, rhos, correlations,
          n_generic_risks, n_pca_components);

      const int sx1 = s * n_generic_risks;
      for (int j = 0; j < n_generic_risks; ++j)
      {
        generic_risk_gammas[sx1 + j] += exp(generic_risk_stds[j] * sqrt(1.0 + s) * nlns[j]);
      }

      const int sx2 = s * n_forward_rate_risks;
      for (int j = 0; j < n_forward_rate_risks; ++j)
      {
        const double lambda = sqrt(pca_forward_rate_risk_eigenvalues[j]);
        double epsilon = 0.0;
        for (int k = 0; k < n_pca_components; ++k)
        {
          const int row = j * n_forward_rate_risks;
          epsilon += nlns[n_generic_risks + k] * lambda * sqrt(12.0 * (1.0 + s)) * pca_forward_rate_risk_eigenvectors[row + k];
        }
        forward_rate_risk_gammas[sx2 + j] += exp(epsilon);
      }
    }
  }

  for (int i = 0; i < generic_risk_gammas.size(); ++i)
  {
    generic_risk_gammas[i] = log(generic_risk_gammas[i] / (double)n_trials);
  }

  for (int i = 0; i < forward_rate_risk_gammas.size(); ++i)
  {
    forward_rate_risk_gammas[i] = log(forward_rate_risk_gammas[i] / (double)n_trials);
  }

  return make_tuple(generic_risk_gammas, forward_rate_risk_gammas);
}

double compute_domestic_equity_price(
    vector<double> &generic_risk_values,
    vector<double> &zero_coupon_rates)
{
  return generic_risk_values[DOMESTIC_MARKET_RISK_INDEX];
}

double compute_global_equity_price(
    vector<double> &generic_risk_values,
    vector<double> &zero_coupon_rates)
{
  return generic_risk_values[GLOBAL_MARKET_RISK_INDEX];
}

double compute_real_estate_price(vector<double> &generic_risk_values, vector<double> &zero_coupon_rates)
{
  return generic_risk_values[REAL_ESTATE_RISK_INDEX];
}

double compute_alternative_price(vector<double> &generic_risk_values, vector<double> &zero_coupon_rates)
{
  return generic_risk_values[ALTERNATIVE_RISK_INDEX];
}

double compute_credit_price(vector<double> &generic_risk_values, vector<double> &zero_coupon_rates)
{
  return 1.0; // TODO: credit price is computed both from the yield curve and credit risk
}

double compute_cash_price(vector<double> &generic_risk_values, vector<double> &zero_coupon_rates)
{
  return 1.0;
}

double compute_bonds_2y_price(vector<double> &generic_risk_values, vector<double> &zero_coupon_rates)
{
  return 1.0 / pow(1.0 + zero_coupon_rates[1], 2.0);
}

double compute_bonds_5y_price(vector<double> &generic_risk_values, vector<double> &zero_coupon_rates)
{
  return 1.0 / pow(1.0 + zero_coupon_rates[4], 5.0);
}

double compute_bonds_20y_price(vector<double> &generic_risk_values, vector<double> &zero_coupon_rates)
{
  return 1.0 / pow(1.0 + zero_coupon_rates[19], 20.0);
}

double compute_domestic_equity_future_price(
    vector<double> &generic_risk_values,
    vector<double> &zero_coupon_rates)
{
  return generic_risk_values[DOMESTIC_MARKET_RISK_INDEX];
}

double compute_interest_rate_swap_2y_price(vector<double> &generic_risk_values, vector<double> &zero_coupon_rates)
{
  return 1.0 / pow(1.0 + zero_coupon_rates[1], 2.0);
}

double compute_interest_rate_swap_5y_price(vector<double> &generic_risk_values, vector<double> &zero_coupon_rates)
{
  return 1.0 / pow(1.0 + zero_coupon_rates[4], 5.0);
}

double compute_interest_rate_swap_20y_price(vector<double> &generic_risk_values, vector<double> &zero_coupon_rates)
{
  return 1.0 / pow(1.0 + zero_coupon_rates[19], 20.0);
}

double compute_instrument_price(
    const int instrument_index,
    vector<double> &generic_risk_values,
    vector<double> &zero_coupon_rates)
{
  if (instrument_index == DOMESTIC_EQUITY_INDEX)
  {
    return compute_domestic_equity_price(generic_risk_values, zero_coupon_rates);
  }
  else if (instrument_index == GLOBAL_EQUITY_INDEX)
  {
    return compute_global_equity_price(generic_risk_values, zero_coupon_rates);
  }
  else if (instrument_index == REAL_ESTATE_INDEX)
  {
    return compute_real_estate_price(generic_risk_values, zero_coupon_rates);
  }
  else if (instrument_index == ALTERNATIVE_INDEX)
  {
    return compute_alternative_price(generic_risk_values, zero_coupon_rates);
  }
  else if (instrument_index == CREDIT_INDEX)
  {
    return compute_credit_price(generic_risk_values, zero_coupon_rates);
  }
  else if (instrument_index == BONDS_2Y_INDEX)
  {
    return compute_bonds_2y_price(generic_risk_values, zero_coupon_rates);
  }
  else if (instrument_index == BONDS_5Y_INDEX)
  {
    return compute_bonds_5y_price(generic_risk_values, zero_coupon_rates);
  }
  else if (instrument_index == BONDS_20Y_INDEX)
  {
    return compute_bonds_20y_price(generic_risk_values, zero_coupon_rates);
  }
  else if (instrument_index == CASH_INDEX)
  {
    return compute_cash_price(generic_risk_values, zero_coupon_rates);
  }
  else if (instrument_index == DOMESTIC_EQUITY_FUTURE_INDEX)
  {
    return compute_domestic_equity_future_price(generic_risk_values, zero_coupon_rates);
  }
  else if (instrument_index == INTEREST_RATE_SWAP_2Y_INDEX)
  {
    return compute_interest_rate_swap_2y_price(generic_risk_values, zero_coupon_rates);
  }
  else if (instrument_index == INTEREST_RATE_SWAP_5Y_INDEX)
  {
    return compute_interest_rate_swap_5y_price(generic_risk_values, zero_coupon_rates);
  }
  else if (instrument_index == INTEREST_RATE_SWAP_20Y_INDEX)
  {
    return compute_interest_rate_swap_20y_price(generic_risk_values, zero_coupon_rates);
  }
  else
  {
    // This should never occur, but currently Embind does
    // not allow us to throw a proper exception.
    return 0.0;
  }
}

vector<double>
interpolate_zero_curve(
    vector<double> &zero_curve,
    vector<double> &tenors)
{
  vector<double> interpolated(tenors[tenors.size() - 1]); // Assume sorted
  interpolated[0] = zero_curve[0];
  int i = 1;
  int j = 1;

  while (i < zero_curve.size())
  {
    const int d = tenors[i] - tenors[i - 1];
    if (d > 1)
    {
      while ((j + 1) < tenors[i])
      {
        interpolated[j] = ((((j + 1) - tenors[i - 1]) * zero_curve[i]) + (tenors[i] - (j + 1)) * zero_curve[i - 1]) / d;
        j++;
      }
    }
    interpolated[j] = zero_curve[i];
    i++;
    j++;
  }

  return interpolated;
}

tuple<vector<double>, vector<double>> compute_instrument_prices(
    vector<double> &initial_forward_rate_risk_values,
    vector<double> &initial_generic_risk_values,
    vector<double> &intermediate_generic_risk_values,
    vector<double> &intermediate_forward_rate_risk_values,
    vector<double> &zero_curve_tenors,
    const int n_generic_risks,
    const int n_forward_rate_risks,
    const int n_instruments,
    const int n_states)
{
  // Compute initial zero curve
  vector<double> initial_discount_factors = compute_discount_factors(initial_forward_rate_risk_values);
  vector<double> initial_zero_curve = compute_zero_coupon_rates(initial_discount_factors);
  vector<double> interpolated_initial_zero_curve = interpolate_zero_curve(initial_zero_curve, zero_curve_tenors);

  vector<double> initial_prices(n_instruments);
  for (int j = 0; j < n_instruments; ++j)
  {
    initial_prices[j] = compute_instrument_price(j, initial_generic_risk_values, interpolated_initial_zero_curve);
  }

  // Compute intermediate instrument prices from intermediate risk values
  vector<double> intermediate_prices(n_states * n_instruments);
  vector<double> instrument_changes(n_states * n_instruments);
  for (int i = 0; i < n_states; ++i)
  {
    const int ix = i * n_instruments;
    const int gx = i * n_generic_risks;
    const int fx = i * n_forward_rate_risks;

    vector<double> current_generic_risk_values(
        intermediate_generic_risk_values.begin() + gx,
        intermediate_generic_risk_values.begin() + gx + n_generic_risks);
    vector<double> current_forward_rate_risk_values(
        intermediate_forward_rate_risk_values.begin() + fx,
        intermediate_forward_rate_risk_values.begin() + fx + n_forward_rate_risks);

    vector<double> current_discount_factors = compute_discount_factors(current_forward_rate_risk_values);
    vector<double> current_zero_curve = compute_zero_coupon_rates(current_discount_factors);
    vector<double> interpolated_current_zero_curve = interpolate_zero_curve(current_zero_curve, zero_curve_tenors);

    for (int j = 0; j < n_instruments; ++j)
    {
      intermediate_prices[ix + j] = compute_instrument_price(j, current_generic_risk_values, interpolated_current_zero_curve);
    }
  }

  return make_tuple(initial_prices, intermediate_prices);
}

double compute_instrument_change(
    const int instrument_index,
    vector<double> &current_instrument_prices,
    vector<double> &previous_instrument_prices,
    vector<double> &current_zero_curve)
{
  double change = 0.0;
  if (instrument_index == DOMESTIC_EQUITY_INDEX)
  {
    return (current_instrument_prices[instrument_index] / previous_instrument_prices[instrument_index]) - 1.0;
  }
  else if (instrument_index == GLOBAL_EQUITY_INDEX)
  {
    return (current_instrument_prices[instrument_index] / previous_instrument_prices[instrument_index]) - 1.0;
  }
  else if (instrument_index == REAL_ESTATE_INDEX)
  {
    return (current_instrument_prices[instrument_index] / previous_instrument_prices[instrument_index]) - 1.0;
  }
  else if (instrument_index == ALTERNATIVE_INDEX)
  {
    return (current_instrument_prices[instrument_index] / previous_instrument_prices[instrument_index]) - 1.0;
  }
  else if (instrument_index == CREDIT_INDEX)
  {
    return (current_instrument_prices[instrument_index] / previous_instrument_prices[instrument_index]) - 1.0;
  }
  else if (instrument_index == BONDS_2Y_INDEX)
  {
    const double new_price = 1.0 / pow(1.0 + current_zero_curve[0], 1.0);
    return (new_price / previous_instrument_prices[instrument_index]) - 1.0;
  }
  else if (instrument_index == BONDS_5Y_INDEX)
  {
    const double new_price = 1.0 / pow(1.0 + current_zero_curve[3], 4.0);
    return (new_price / previous_instrument_prices[instrument_index]) - 1.0;
  }
  else if (instrument_index == BONDS_20Y_INDEX)
  {
    const double new_price = 1.0 / pow(1.0 + current_zero_curve[18], 19.0);
    return (new_price / previous_instrument_prices[instrument_index]) - 1.0;
  }
  else if (instrument_index == CASH_INDEX)
  {
    return 0.0;
  }
  else if (instrument_index == DOMESTIC_EQUITY_FUTURE_INDEX)
  {
    return (current_instrument_prices[instrument_index] / previous_instrument_prices[instrument_index]) - 1.0;
  }
  else if (instrument_index == INTEREST_RATE_SWAP_2Y_INDEX)
  {
    const double new_price = 1.0 / pow(1.0 + current_zero_curve[0], 1.0);
    return (new_price / previous_instrument_prices[instrument_index]) - 1.0;
  }
  else if (instrument_index == INTEREST_RATE_SWAP_5Y_INDEX)
  {
    const double new_price = 1.0 / pow(1.0 + current_zero_curve[3], 4.0);
    return (new_price / previous_instrument_prices[instrument_index]) - 1.0;
  }
  else if (instrument_index == INTEREST_RATE_SWAP_20Y_INDEX)
  {
    const double new_price = 1.0 / pow(1.0 + current_zero_curve[18], 19.0);
    return (new_price / previous_instrument_prices[instrument_index]) - 1.0;
  }

  return change;
}

vector<double> compute_instrument_changes(
    vector<double> &initial_instrument_prices,
    vector<double> &intermediate_instrument_prices,
    vector<double> &intermediate_forward_rates,
    vector<double> &zero_curve_tenors,
    const int n_forward_rate_risks,
    const int n_instruments,
    const int n_states)
{
  vector<double> instrument_changes(n_states * n_instruments);

  // Initial changes
  vector<double> current_forward_rate_risk_values(
      intermediate_forward_rates.begin(),
      intermediate_forward_rates.begin() + n_forward_rate_risks);

  vector<double> current_discount_factors = compute_discount_factors(current_forward_rate_risk_values);
  vector<double> current_zero_curve = compute_zero_coupon_rates(current_discount_factors);
  vector<double> interpolated_zero_curve = interpolate_zero_curve(current_zero_curve, zero_curve_tenors);

  vector<double> current_instrument_prices(
      intermediate_instrument_prices.begin(),
      intermediate_instrument_prices.begin() + n_instruments);

  vector<double> initial_prices(n_instruments);
  for (int j = 0; j < n_instruments; ++j)
  {
    instrument_changes[j] = compute_instrument_change(j, current_instrument_prices, initial_instrument_prices, interpolated_zero_curve);
  }

  // Intermediate changes
  for (int i = 1; i < n_states; ++i)
  {
    const int p = ((i - 1) / 2);

    const int px = p * n_instruments; // Parent
    const int ix = i * n_instruments; // Child
    const int fx = i * n_forward_rate_risks;

    vector<double> current_instrument_prices(
        intermediate_instrument_prices.begin() + ix,
        intermediate_instrument_prices.begin() + ix + n_instruments);
    vector<double> previous_instrument_prices(
        intermediate_instrument_prices.begin() + px,
        intermediate_instrument_prices.begin() + px + n_instruments);

    vector<double> current_forward_rate_risk_values(
        intermediate_forward_rates.begin() + fx,
        intermediate_forward_rates.begin() + fx + n_forward_rate_risks);

    vector<double> current_discount_factors = compute_discount_factors(current_forward_rate_risk_values);
    vector<double> current_zero_curve = compute_zero_coupon_rates(current_discount_factors);
    vector<double> interpolated_zero_curve = interpolate_zero_curve(current_zero_curve, zero_curve_tenors);

    for (int j = 0; j < n_instruments; ++j)
    {
      instrument_changes[ix + j] = compute_instrument_change(j, current_instrument_prices, previous_instrument_prices, interpolated_zero_curve);
    }
  }

  return instrument_changes;
}

vector<double> generate_state_changes(
    vector<double> &initial_generic_risk_values,
    vector<double> &initial_forward_rate_risk_values,
    vector<double> &generic_risk_means,
    vector<double> &generic_risk_stds,
    vector<double> &pca_forward_rate_risk_eigenvalues,
    vector<double> &pca_forward_rate_risk_eigenvectors,
    vector<double> &sigmas,
    vector<double> &rhos,
    vector<double> &correlations,
    vector<double> &zero_curve_tenors,
    const int n_instruments,
    const int n_generic_risks,
    const int n_forward_rate_risks,
    const int n_pca_components,
    const int n_states,
    const int n_steps)
{
  const int n_risks = n_generic_risks + n_forward_rate_risks;

  // Compute gammas
  vector<double> generic_risk_sigmas(sigmas.begin(), sigmas.begin() + n_generic_risks);
  vector<double> generic_risk_rhos(rhos.begin(), rhos.begin() + n_generic_risks);

  vector<double> forward_rate_risk_sigmas(sigmas.begin() + n_generic_risks, sigmas.end());
  vector<double> forward_rate_risk_rhos(rhos.begin() + n_generic_risks, rhos.end());

  // TODO: Compute outside of this function
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
      n_states,
      n_steps);

  vector<double> generic_gammas = get<0>(gammas);
  vector<double> forward_rate_gammas = get<1>(gammas);

  // Generate intermediate risk values
  vector<double> intermediate_generic_risk_values(n_states * n_generic_risks);
  vector<double> intermediate_forward_rate_risk_values(n_states * n_forward_rate_risks);

  for (int i = 0; i < n_states; ++i)
  {
    vector<double> nlns = generate_correlated_nlns(
        sigmas, rhos, correlations,
        n_generic_risks, n_pca_components);

    const double t = floor(log2(i + 1)) + 1.0;

    const int ix1 = i * n_generic_risks;
    const int tx1 = t * n_generic_risks;
    for (int j = 0; j < n_generic_risks; ++j)
    {
      const double s0 = initial_generic_risk_values[j];
      const double mu = generic_risk_means[j];
      const double std = generic_risk_stds[j];
      const double gamma = generic_gammas[tx1 + j];
      const double epsilon = std * sqrt(t) * nlns[j];
      intermediate_generic_risk_values[ix1 + j] = s0 * exp(mu * t - gamma + epsilon);
    }

    const int ix2 = i * n_forward_rate_risks;
    const int tx2 = t * n_forward_rate_risks;
    for (int j = 0; j < n_forward_rate_risks; ++j)
    {
      const double s0 = initial_forward_rate_risk_values[j];
      const double lambda = sqrt(pca_forward_rate_risk_eigenvalues[j]);
      const double gamma = forward_rate_gammas[tx2 + j];
      double epsilon = 0.0;
      for (int k = 0; k < n_pca_components; ++k)
      {
        const int row = j * n_forward_rate_risks;

        epsilon += nlns[n_generic_risks + j] * lambda * sqrt(12.0 * t) * pca_forward_rate_risk_eigenvectors[row + k];
      }
      intermediate_forward_rate_risk_values[ix2 + j] = s0 * exp(epsilon - gamma);
    }
  }

  // Compute instrument prices
  tuple<vector<double>, vector<double>> instrument_prices = compute_instrument_prices(
      initial_forward_rate_risk_values,
      initial_generic_risk_values,
      intermediate_generic_risk_values,
      intermediate_forward_rate_risk_values,
      zero_curve_tenors,
      n_generic_risks,
      n_forward_rate_risks,
      n_instruments,
      n_states);

  // Compute instrument changes
  vector<double> initial_instrument_prices = get<0>(instrument_prices);
  vector<double> intermediate_instrument_prices = get<1>(instrument_prices);

  vector<double> instrument_changes = compute_instrument_changes(
      initial_instrument_prices,
      intermediate_instrument_prices,
      intermediate_forward_rate_risk_values,
      zero_curve_tenors,
      n_forward_rate_risks,
      n_instruments,
      n_states);

  return instrument_changes;
}