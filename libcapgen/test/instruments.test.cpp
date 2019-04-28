#include <lib/catch.h>
#include <lib/random.h>

#include <include/constants.h>
#include <include/instruments.h>
#include <include/normal.h>

TEST_CASE("generate_instrument_changes runs without crashing", "[normal]")
{
  const int n_risks = N_RISKS;
  const int n_instruments = N_INSTRUMENTS;
  std::vector<double> means = NORMAL_DEFAULT_MEANS;
  std::vector<double> standard_deviations = NORMAL_DEFAULT_STANDARD_DEVIATIONS;
  std::vector<double> correlations = NORMAL_DEFAULT_CORRELATIONS;
  std::vector<double>
      risk_changes = generate_normal_risk_changes(
          means,
          standard_deviations,
          correlations,
          n_risks);

  std::vector<double> instrument_changes = generate_instrument_changes(
      risk_changes,
      n_instruments);
}
