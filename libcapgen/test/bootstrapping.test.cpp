#include <lib/catch.h>
#include <lib/random.h>

#include <include/constants.h>
#include <include/bootstrapping.h>

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

/*
TEST_CASE("compute_dfs_from_pars", "[bootstrapping]")
{
  std::vector<double> pars = {
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
  std::vector<double> expected_dfs = {
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

  std::vector<double> actual_dfs = compute_dfs_from_pars(pars);

  for (int i = 0; i < expected_dfs.size(); ++i)
  {
    REQUIRE(actual_dfs[i] == Approx(expected_dfs[i]).margin(0.0001));
  }
}
*/