#include <lib/catch.h>
#include <lib/random.h>

#include <include/constants.h>
#include <include/normal.h>

using Random = effolkronium::random_static;

TEST_CASE("generate_risk_changes runs without crashing", "[normal]")
{
  const int n_risks = N_RISKS;
  std::vector<double> means = NORMAL_DEFAULT_MEANS;
  std::vector<double> standard_deviations = NORMAL_DEFAULT_STANDARD_DEVIATIONS;
  std::vector<double> correlations = NORMAL_DEFAULT_CORRELATIONS;
  std::vector<double>
      risk_changes = generate_normal_risk_changes(
          means,
          standard_deviations,
          correlations,
          n_risks);
}

/*
TEST_CASE("generate_normal_scenarios yields reasonable changes", "[normal]")
{
  const int n_risks = N_RISKS;
  const int n_instruments = N_INSTRUMENTS;
  const int n_scenarios = 3;
  std::vector<double> means = NORMAL_DEFAULT_MEANS;
  std::vector<double> standard_deviations = NORMAL_DEFAULT_STANDARD_DEVIATIONS;
  std::vector<double> correlations = NORMAL_DEFAULT_CORRELATIONS;
  std::vector<double> instrument_changes = generate_normal_scenarios(
      means,
      standard_deviations,
      correlations,
      n_risks,
      n_instruments,
      n_scenarios);

  for (int i = 0; i < instrument_changes.size(); ++i)
  {
    REQUIRE(instrument_changes[i] < 0.1);
    REQUIRE(instrument_changes[i] > -0.1);
  }
}
*/