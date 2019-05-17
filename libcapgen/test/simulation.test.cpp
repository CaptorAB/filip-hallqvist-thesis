#include <lib/catch.h>
#include <lib/random.h>
#include <lib/stats/stats.hpp>

#include <include/constants.h>
#include <include/simulation.h>

using namespace std;
using namespace stats;

double compute_correlation(vector<double> X, vector<double> Y, const int n)
{
  double sum_X = 0;
  double sum_Y = 0;
  double sum_XY = 0;

  double squareSum_X = 0;
  double squareSum_Y = 0;

  for (int i = 0; i < n; i++)
  {
    sum_X = sum_X + X[i];
    sum_Y = sum_Y + Y[i];
    sum_XY = sum_XY + X[i] * Y[i];
    squareSum_X = squareSum_X + X[i] * X[i];
    squareSum_Y = squareSum_Y + Y[i] * Y[i];
  }

  return ((double)n * sum_XY - sum_X * sum_Y) / sqrt(((double)n * squareSum_X - sum_X * sum_X) * ((double)n * squareSum_Y - sum_Y * sum_Y));
}

TEST_CASE("standardize_nln correctly standardizes nln distributed numbers", "[simulation]")
{
  const int n_randoms = 10000;
  const double sigma = 0.5;
  const double rho = 0.9;

  vector<double> nlns(n_randoms);
  for (int i = 0; i < n_randoms; ++i)
  {
    const double n1 = rnorm();
    const double n2 = rnorm();
    const double nln = sample_nln(n1, n2, sigma, rho);

    nlns[i] = standardize_nln(nln, sigma, rho);
  }

  double sum = accumulate(nlns.begin(), nlns.end(), 0.0);
  double mean = sum / n_randoms;

  double accumulator = 0.0;
  std::for_each(nlns.begin(), nlns.end(), [&](const double d) {
    accumulator += (d - mean) * (d - mean);
  });

  double stdev = sqrt(accumulator / (n_randoms - 1));

  REQUIRE(mean == Approx(0.0).margin(0.02));
  REQUIRE(stdev == Approx(1.0).margin(0.02));
}

/*
TEST_CASE("correlate_risks correctly correlates risks", "[simulation]")
{
  const int n_samples = 1000;
  const int n_generic_risks = 2;
  const int n_pca_components = 3;
  const int n_randoms = n_generic_risks + n_pca_components;
  const int n_randoms_per_sample = 2 * n_randoms;

  vector<double> correlations = {
      1.0, 0.1, 0.2,
      0.1, 1.0, -0.5,
      0.2, -0.5, 1.0};

  vector<vector<double>> correlated;

  for (int i = 0; i < n_randoms; ++i)
  {
    correlated.push_back(vector<double>());
  }

  for (int i = 0; i < n_samples; ++i)
  {
    const int ix = i * n_randoms_per_sample;
    vector<double> uniforms = sample_uniform_randoms(2 * (n_generic_risks * n_pca_components));
    vector<double> normals = normalize_uniform_randoms(uniforms);
    vector<double> normals_correlated = correlate_risks(normals, correlations, n_generic_risks);

    for (int j = 0; j < n_randoms; ++j)
    {
      const int jx = j * n_randoms;
      for (int k = 0; k < 2; ++k)
      {
        correlated[j].push_back(normals_correlated[jx + k]);
      }
    }
  }

  vector<vector<double>> coefficients(n_randoms);
  for (int i = 0; i < n_randoms; ++i) {
    const int ix = i * n_randoms;
    for (int j = 0; j < n_randoms; ++j) {
      coefficients[ix + j] = 
    }
  }
}
*/

TEST_CASE("simulation golden master", "[simulation]")
{
  vector<double> initial_generic_risk_values = {1.0, 1.0};
  vector<double> initial_forward_rate_risk_values = {0.01768, 0.02201, 0.03055, 0.03431, 0.03762};
  vector<double> generic_risk_means = {0.1, 0.1};
  vector<double> generic_risk_stds = {0.1, 0.1};

  vector<double> pca_forward_rate_risk_eigenvalues = {0.223282705428122, 0.082751342396589, 0.048454437715918, 0.040635667202019, 0.02539470506845};
  vector<double> pca_forward_rate_risk_eigenvectors = {
      0.17806310748978, 0.336052289069525, 0.352104530046595, 0.389724410069539, 0.349021453294063,
      -0.596316760839539, -0.527975500405015, -0.251583154671794, -0.041465903559625, 0.094526967296751,
      0.486340482637313, 0.028187296071813, -0.352705203001156, -0.309926291591076, -0.339610973326499,
      0.555751828998025, -0.427628181205936, -0.25046231532465, 0.107727039741423, 0.191438014022832,
      -0.192920659287202, 0.317837716410373, 0.11595703534322, -0.13446138242021, -0.355170758402066};

  vector<double> correlations = {
      1.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 1.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 1.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 1.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 1.0};

  vector<double> sigmas = {
      1.1, 1.1, 1.7};
  vector<double> rhos = {
      -0.15, -0.15, -0.1};

  const int n_instruments = 2;
  const int n_generic_risks = 2;
  const int n_forward_rate_risks = 5;
  const int n_pca_components = 3;
  const int n_states = 3;
  const int n_steps = 2;

  vector<double> instrument_changes = generate_state_changes(
      initial_generic_risk_values,
      initial_forward_rate_risk_values,
      generic_risk_means,
      generic_risk_stds,
      pca_forward_rate_risk_eigenvalues,
      pca_forward_rate_risk_eigenvectors,
      sigmas,
      rhos,
      correlations,
      n_instruments,
      n_generic_risks,
      n_forward_rate_risks,
      n_pca_components,
      n_states,
      n_steps);

  for (int i = 0; i < n_states; ++i)
  {
    const int ix = i * n_instruments;
    for (int j = 0; j < n_instruments; ++j)
    {
      printf("%.4f ", instrument_changes[ix + j]);
    }
    printf("\n");
  }
}
