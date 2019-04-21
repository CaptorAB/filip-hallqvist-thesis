#include <iostream>
#include <lib/catch.h>
#include <lib/random.h>

#include <include/constants.h>
#include <include/normal.h>

using Random = effolkronium::random_static;

TEST_CASE("cholesky correctly performs cholesky decomposition", "[normal]")
{
  std::vector<double> correlations({1.00, 0.50, 0.25, -0.10,
                                    0.50, 1.00, 0.75, -0.20,
                                    0.25, 0.75, 1.00, 0.40,
                                    -0.10, -0.20, 0.40, 1.00});
  const int n_instruments = 4;

  std::vector<double> cholesky = compute_cholesky(correlations, n_instruments);

  const std::vector<double> expected = {
      1.0, 0.0, 0.0, 0.0,
      0.5, 0.86603, 0.0, 0.0,
      0.25, 0.72169, 0.6455, 0.0,
      -0.10, -0.1732, 0.85206, 0.48374};

  for (int i = 0; i < cholesky.size(); ++i)
  {
    REQUIRE(cholesky[i] == Approx(expected[i]).epsilon(0.0001));
  }
}

TEST_CASE("generate_risk_changes runs without crashing", "[normal]")
{
  const int n_risks = N_RISKS;
  std::vector<double> means = NORMAL_DEFAULT_MEANS;
  std::vector<double> standard_deviations = NORMAL_DEFAULT_STANDARD_DEVIATIONS;
  std::vector<double> correlations = NORMAL_DEFAULT_CORRELATIONS;
  std::vector<double>
      risk_changes = generate_risk_changes(
          means,
          standard_deviations,
          correlations,
          n_risks);
}

TEST_CASE("generate_instrument_changes runs without crashing", "[normal]")
{
  const int n_risks = N_RISKS;
  const int n_instruments = N_INSTRUMENTS;
  std::vector<double> means = NORMAL_DEFAULT_MEANS;
  std::vector<double> standard_deviations = NORMAL_DEFAULT_STANDARD_DEVIATIONS;
  std::vector<double> correlations = NORMAL_DEFAULT_CORRELATIONS;
  std::vector<double>
      risk_changes = generate_risk_changes(
          means,
          standard_deviations,
          correlations,
          n_risks);

  std::vector<double> instrument_changes = generate_instrument_changes(
      risk_changes,
      n_instruments);
}

TEST_CASE("generate_normal_scenarios runs without crashing", "[normal]")
{
  const int n_risks = N_RISKS;
  const int n_instruments = N_INSTRUMENTS;
  const int n_scenarios = 3;
  std::vector<double> means = NORMAL_DEFAULT_MEANS;
  std::vector<double> standard_deviations = NORMAL_DEFAULT_STANDARD_DEVIATIONS;
  std::vector<double> correlations = NORMAL_DEFAULT_CORRELATIONS;
  std::vector<double> price_changes = generate_normal_scenarios(
      means,
      standard_deviations,
      correlations,
      n_risks,
      n_instruments,
      n_scenarios);
}