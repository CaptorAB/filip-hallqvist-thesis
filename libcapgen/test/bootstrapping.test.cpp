#include <lib/catch.h>
#include <lib/random.h>

#include <include/constants.h>
#include <include/bootstrapping.h>

using namespace std;

TEST_CASE("adjust_credit_par_rates", "[bootstrapping]")
{
  std::vector<double> actual = {0.05, -1.0, 0.1, -1.0};
  adjust_credit_par_rates(actual);

  std::vector<double> expected = {0.05, -1.0, 0.1, -1.0};

  for (int i = 0; i < actual.size(); ++i)
  {
    if (expected[i] == -1.0)
    {
      REQUIRE(actual[i] == -1.0);
    }
    else
    {
      REQUIRE(actual[i] == expected[i] + PAR_CREDIT_ADJUSTMENT);
    }
  }
}

TEST_CASE("golden master bootstrapping test", "[bootstrapping]")
{
  vector<double> par_rates = {
      0.0,
      0.009700,
      0.011775,
      0.014200,
      0.016580,
      0.018580,
      0.020130,
      0.021400,
      0.022430,
      0.023280,
      0.023950,
      -1.0,
      0.02495,
      -1.0,
      -1.0,
      0.02590,
      -1.0,
      -1.0,
      -1.0,
      -1.0,
      0.02690};

  vector<double> expected_discount_factors = {
      1.0,
      0.9904,
      0.9768,
      0.9585,
      0.9360,
      0.9113,
      0.8861,
      0.8605,
      0.8350,
      0.8099,
      0.7856,
      0.7621,
      0.7392,
      0.7172,
      0.6959,
      0.6752,
      0.6549,
      0.6353,
      0.6162,
      0.5977,
      0.5797};

  vector<double> actual_discount_factors = bootstrap_discount_factors(par_rates);

  // Disount factors
  REQUIRE(actual_discount_factors.size() == expected_discount_factors.size());
  for (int i = 0; i < expected_discount_factors.size(); ++i)
  {
    REQUIRE(actual_discount_factors[i] == Approx(expected_discount_factors[i]).margin(0.0001));
  }

  vector<double> expected_forward_rates = {
      0.009700,
      0.013879,
      0.019177,
      0.024019,
      0.027055,
      0.028479,
      0.029752,
      0.030472,
      0.030998,
      0.030913,
      0.030911,
      0.030911,
      0.030654,
      0.030654,
      0.030654,
      0.030966,
      0.030966,
      0.030966,
      0.030966,
      0.030966};
  vector<double> actual_forward_rates = compute_forward_rates(expected_discount_factors);

  // Forward rates
  REQUIRE(actual_forward_rates.size() == expected_forward_rates.size());
  for (int i = 0; i < expected_forward_rates.size(); ++i)
  {
    REQUIRE(actual_forward_rates[i] == Approx(expected_forward_rates[i]).margin(0.001));
  }

  // UFR-weighted forward rates
  vector<double> expected_ufr_forward_rates = {
      0.009700,
      0.013879,
      0.019177,
      0.024019,
      0.027055,
      0.028479,
      0.029752,
      0.030472,
      0.030998,
      0.030913,
      0.031920,
      0.032928,
      0.033749,
      0.034780,
      0.035811,
      0.036985,
      0.037988,
      0.038991,
      0.039994,
      0.040997};

  vector<double> actual_ufr_forward_rates = expected_forward_rates;
  adjust_ufr_forward_rates(actual_ufr_forward_rates);

  // Forward rates
  REQUIRE(actual_forward_rates.size() == expected_ufr_forward_rates.size());
  for (int i = 0; i < expected_ufr_forward_rates.size(); ++i)
  {
    REQUIRE(actual_ufr_forward_rates[i] == Approx(expected_ufr_forward_rates[i]).margin(0.00001));
  }
}