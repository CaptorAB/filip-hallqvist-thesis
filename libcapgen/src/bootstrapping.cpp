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
    s += v;
    t++;
  }
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

  int i = 0;
  double df_tb = df_ta; // Initial guess
  do
  {

    const double f = newton_raphson_minimize(p, s, df_ta, df_tb, ta, tb);
    const double df = (newton_raphson_minimize(p, s, df_ta, df_tb + 0.1, ta, tb) - f) / 0.1;
    const double d = f / df;

    if (abs(d) < tolerance)
    {
      return df_tb;
    }

    df_tb = df_tb - d;

  } while (i++ < max_iterations);
}

vector<double>
compute_forward_rates(
    vector<double> &discount_factors)
{
  vector<double> forward_rates(discount_factors.size() - 1);
  for (int i = 1; i < discount_factors.size(); ++i)
  {
    forward_rates[i - 1] = (discount_factors[i - 1] / discount_factors[i]) - 1.0;
  }
  return forward_rates;
}

void adjust_ufr_forward_rates(
    vector<double> &forward_rates)
{
  const int T1 = 10;
  const int T2 = 20;

  for (int t = T1; t < forward_rates.size(); ++t)
  {
    const int T = t + 1;
    double w = 0.0;

    if (T1 < T && T <= T2)
    {
      w = ((double)(T - T1)) / ((double)(T2 - T1 + 1.0));
    }
    else if (T2 < T)
    {
      w = 1.0;
    }

    forward_rates[t] = (1.0 - w) * forward_rates[t] + (w * ULTIMATE_FORWARD_RATE);
  }
}

vector<double>
bootstrap_discount_factors(
    vector<double> &par_rates)
{
  const int T = par_rates.size();
  vector<double> discount_factors(T);
  discount_factors[0] = 1.0;

  double s = 0.0;

  int t = 1;
  while (t < T)
  {
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