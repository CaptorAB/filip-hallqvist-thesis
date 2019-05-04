#include <vector>
#include <math.h>

#include <include/constants.h>

using namespace std;

void adjust_credit_par_rates(std::vector<double> &par_rates)
{
  for (int i = 0; i < par_rates.size(); ++i)
  {
    if (par_rates[i] != -1.0)
    {
      par_rates[i] = par_rates[i] + PAR_CREDIT_ADJUSTMENT;
    }
  }
}

double f(
    const double r,
    const double p,
    std::vector<double> &df,
    const int a,
    const int b)
{
  double dfb = exp(-1.0 * r * b);
  double fr = -1.0 * (log(dfb) - log(df[a])) / (b - a);
  double sum = 0.0;

  for (int t = 0; t <= a; ++t)
  {
    sum += df[t];
  }
  for (int t = a + 1; t <= b; ++t)
  {
    sum += df[a] / pow(1.0 + fr, t - a);
  }

  sum *= p;
  return sum - (1 - df[b]);
}

double compute_intermediate_discount_factor(
    const int t,
    const int ta,
    const int tb,
    const double df_ta,
    const double df_tb)
{
  const double k = ((double)(t - ta)) / ((double)(tb - ta));
  return pow(df_ta, 1.0 - k) * pow(df_tb, k);
}

double
newton_raphson_minimize(
    const double p,
    const double s0,
    const double df_ta,
    const double df_tb,
    const int ta,
    const int tb)
{
  double s = s0;
  int t = ta + 1;
  while (t <= tb)
  {
    const double v = compute_intermediate_discount_factor(t, ta, tb, df_ta, df_tb);
    printf("computed discount factor at %i is %.4f\n", t, v);
    s += v;
    t++;
  }
  printf("p = %.4f, s = %.4f, df_tb = %.4f\n", p, s, df_tb);
  return p * s - (1.0 - df_tb);
}

double
newton_raphson(
    const double p,
    const double s,
    const double df_ta,
    const int ta,
    const int tb)
{
  const int max_iterations = 20;
  const double tolerance = 0.00001;

  printf("Computing for %i\n", tb);
  int i = 0;
  double df_tb = df_ta; // Initial guess
  do
  {

    printf("Computing f\n");
    const double f = newton_raphson_minimize(p, s, df_ta, df_tb, ta, tb);
    printf("Computing df\n");
    const double df = (newton_raphson_minimize(p, s, df_ta, df_tb + 0.1, ta, tb) - f) / 0.1;
    const double d = f / df;

    printf("f = %.4f, df = %.4f, d = %.4f\n", f, df, d);

    if (abs(d) < tolerance)
    {
      return df_tb;
    }

    df_tb = df_tb - d;

  } while (i++ < max_iterations);
}

vector<double>
bootstrap(
    vector<double> &par_rates)
{
  const int T = par_rates.size();
  vector<double> discount_factors(T);
  discount_factors[0] = 1.0;

  double s = 0.0;

  int t = 1;
  while (t < T)
  {
    printf("t = %i\n", t);
    double p = par_rates[t];
    if (p == -1.0)
    {
      int ta = t - 1;
      int tb = t + 1;
      while (par_rates[tb] == -1.0)
        tb++;

      p = par_rates[tb];

      const double df_ta = discount_factors[ta];
      const double df_tb = newton_raphson(p, s, df_ta, ta, tb);

      printf("df_tb = %.4f\n", df_tb);

      // Add discount factors
      int i = t;
      while (i <= tb)
      {
        discount_factors[i] = compute_intermediate_discount_factor(i, ta, tb, df_ta, df_tb);
        s += discount_factors[i];
        i++;
      }

      t = tb + 1;
    }
    else
    {
      discount_factors[t] = (1.0 - p * s) / (1.0 + p);
      s += discount_factors[t];
      t++;
    }
  }

  return discount_factors;
}