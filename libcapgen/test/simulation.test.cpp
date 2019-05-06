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