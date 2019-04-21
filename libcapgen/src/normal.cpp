#include <tuple>
#include <vector>

#include <lib/random.h>

#include <include/constants.h>
#include <include/instruments.h>

using Random = effolkronium::random_static;

/**
 * Compute the cholesky decomposition of a matrix.
 * 
 * @see https://rosettacode.org/wiki/Cholesky_decomposition#C
 */
std::vector<double> compute_cholesky(
    std::vector<double> &A,
    const int n)
{
  std::vector<double> L(n * n);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < (i + 1); j++)
    {
      double s = 0;
      for (int k = 0; k < j; k++)
      {
        s += L[i * n + k] * L[j * n + k];
      }

      if (i == j)
      {
        L[i * n + j] =
            sqrt(A[i * n + i] - s);
      }
      else
      {
        L[i * n + j] =
            (1.0 / L[j * n + j] * (A[i * n + j] - s));
      }
    }
  }
  return L;
}

std::vector<double> generate_risk_changes(
    std::vector<double> &means,
    std::vector<double> &standard_deviations,
    std::vector<double> &correlations,
    const int n_risks)
{
  std::vector<double> risk_changes(n_risks);

  std::vector<double> cholesky = compute_cholesky(
      correlations, n_risks);

  // Generate correlated normal numbers
  for (int i = 0; i < n_risks; ++i)
  {
    const int ix = i * n_risks;
    double correlated = 0.0;
    for (int j = 0; j < n_risks; ++j)
    {
      double random =
          Random::get<std::normal_distribution<>>(means[i], standard_deviations[i]);
      correlated += random * cholesky[ix + j];
    }

    risk_changes[i] = means[i] + correlated * standard_deviations[i];
  }

  return risk_changes;
}

std::vector<double> generate_normal_scenario(
    std::vector<double> &means,
    std::vector<double> &standard_deviations,
    std::vector<double> &correlations,
    const int n_risks,
    const int n_instruments,
    const int n_scenarios)
{
  std::vector<double> risk_changes = generate_risk_changes(
      means,
      standard_deviations,
      correlations,
      n_risks);

  std::vector<double> instrument_changes = generate_instrument_changes(
      risk_changes,
      n_instruments);

  return instrument_changes;
}

std::vector<double> generate_normal_scenarios(
    std::vector<double> &means,
    std::vector<double> &standard_deviations,
    std::vector<double> &correlations,
    const int n_risks,
    const int n_instruments,
    const int n_scenarios)
{
  std::vector<double> instrument_changes(n_scenarios);

  for (int i = 0; i < n_scenarios; ++i)
  {
    int current = i;

    int left = 2 * current + 1;
    int right = 2 * current + 2;

    // Index of current node
    const int cx = current * n_instruments;

    if (left < n_scenarios && right < n_scenarios)
    {
      const int lx = left * n_instruments;
      const int rx = right * n_instruments;

      // Evaluate left child
      std::vector<double> left_scenario = generate_normal_scenario(
          means,
          standard_deviations,
          correlations,
          n_risks,
          n_instruments,
          n_scenarios);

      // Evaluate right child
      std::vector<double> right_scenario = generate_normal_scenario(
          means,
          standard_deviations,
          correlations,
          n_risks,
          n_instruments,
          n_scenarios);

      // Assign branches
      for (int j = 0; j < n_instruments; ++j)
      {
        instrument_changes[lx + j] = left_scenario[j];
        instrument_changes[rx + j] = right_scenario[j];
      }
    }
  }

  return instrument_changes;
}

std::tuple<std::vector<double>, std::vector<double>> generate_normal_goals(
    std::vector<double> &instrument_changes,
    const double initial_funding_ratio,
    const double target_funding_ratio,
    const int n_scenarios,
    const int n_instruments)
{
  std::vector<double> intermediate_goals(n_scenarios);
  std::vector<double> final_goals(n_scenarios / 2 + 1);

  intermediate_goals[0] = 1.0 / initial_funding_ratio;
  int final_index = 0;

  // Iterate through scenarios
  for (int i = 0; i < n_scenarios; ++i)
  {
    int current = i;

    int left = 2 * current + 1;
    int right = 2 * current + 2;

    // Index of current node
    const double current_goal = intermediate_goals[current];
    const double current_change = instrument_changes[FTA_INDEX];

    if (left < n_scenarios && right < n_scenarios)
    {
      const int lx = left * n_instruments;
      const int rx = right * n_instruments;

      // Evaluate left child
      intermediate_goals[left] = target_funding_ratio * current_change * current_goal;
      intermediate_goals[right] = target_funding_ratio * current_change * current_goal;
    }
    else
    {
      final_goals[final_index] = target_funding_ratio * current_change * current_goal;
      final_index++;
    }
  }

  return std::make_tuple(intermediate_goals, final_goals);
}