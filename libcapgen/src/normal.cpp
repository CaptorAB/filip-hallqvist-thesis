#include <tuple>
#include <vector>

#include <lib/random.h>

#include <include/cholesky.h>
#include <include/constants.h>
#include <include/instruments.h>

using Random = effolkronium::random_static;

std::vector<double>
generate_normal_risk_changes(std::vector<double> &means,
                             std::vector<double> &standard_deviations,
                             std::vector<double> &correlations,
                             const int n_risks)
{
  std::vector<double> risk_changes(n_risks);

  std::vector<double> cholesky = compute_cholesky(correlations, n_risks);

  // Generate correlated normal numbers
  for (int i = 0; i < n_risks; ++i)
  {
    /*
    const int ix = i * n_risks;
    double correlated = 0.0;
    for (int j = 0; j < n_risks; ++j)
    {
      double random = Random::get<std::normal_distribution<>>(
          means[i], standard_deviations[i]);
      correlated += random * cholesky[ix + j];
    }

    risk_changes[i] = means[i] + correlated * standard_deviations[i];
    */

    risk_changes[i] = Random::get<std::normal_distribution<>>(
        means[i], standard_deviations[i]);
  }

  return risk_changes;
}

std::vector<double>
generate_normal_scenario(std::vector<double> &means,
                         std::vector<double> &standard_deviations,
                         std::vector<double> &correlations,
                         const int n_risks,
                         const int n_instruments,
                         const int n_scenarios)
{
  std::vector<double> risk_changes =
      generate_normal_risk_changes(means, standard_deviations, correlations, n_risks);

  std::vector<double> instrument_changes =
      generate_instrument_changes(risk_changes, n_instruments);

  return instrument_changes;
}

std::vector<double>
generate_normal_scenarios(std::vector<double> &means,
                          std::vector<double> &standard_deviations,
                          std::vector<double> &correlations,
                          const int n_risks,
                          const int n_instruments,
                          const int n_scenarios)
{
  std::vector<double> instrument_changes(n_instruments * n_scenarios);

  for (int i = 0; i < n_scenarios; ++i)
  {
    std::vector<double> current_scenario =
        generate_normal_scenario(means, standard_deviations, correlations,
                                 n_risks, n_instruments, n_scenarios);

    const int cx = i * n_instruments;
    for (int j = 0; j < n_instruments; ++j)
    {
      instrument_changes[cx + j] = current_scenario[j];
    }
  }

  return instrument_changes;
}
