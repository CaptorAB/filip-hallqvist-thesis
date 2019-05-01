#include <vector>
#include <math.h>

#include <include/constants.h>

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

double g(
    const double z,
    const double par,
    std::vector<double> &dfs,
    const int b)
{
  double sum = 0.0;
  for (int i = 0; i < b; ++i)
  {
    sum += dfs[i];
  }
  return pow(par * sum - (1.0 - (1.0 / (1.0 + pow(z, b)))), 2.0);
}

double g_prime(
    const double z,
    const double par,
    std::vector<double> &dfs,
    const int a,
    const int b)
{
  double sum = 0.0;
  for (int t = a + 1; t < b; ++t)
  {
    sum += pow(dfs[a], ((b - t + a) / (b - a))) * pow((1.0 + pow(z, b)), (t - a) / (b - a));
  }
  return par * sum - pow(1.0 + pow(z, b), -2.0);
}

double f(
    const double z,
    const double par,
    std::vector<double> &dfs,
    const int b)
{
  return pow(g(z, par, dfs, b), 2.0);
}

double f_prime(
    const double z,
    const double par,
    std::vector<double> &dfs,
    const int a,
    const int b)
{
  return 2.0 * f(z, par, dfs, b) * g_prime(z, par, dfs, a, b);
}

std::vector<double> compute_dfs_from_pars(std::vector<double> &pars)
{
  std::vector<double> dfs(pars.size());
  dfs[0] = 1.0;
  int i = 1;
  while (i < dfs.size())
  {
    if (pars[i] == -1.0)
    {
      int j = i + 1;
      while (pars[j] == -1.0)
      {
        j++;
      }

      // We have known par rates at i - 1 and j
      const int a = i - 1;
      const int b = j;
      printf("We have known par rates at %i and %i\n", a, b);

      // Run Newton-Raphson
      int runs = 0;
      double epsilon = 0.00001;
      double d = 0.0;
      double z = 1.0;

      do
      {
        // Set temporary discount factors
        int t = a + 1;
        printf("dfs ");
        while (t <= b)
        {
          dfs[t] = pow(dfs[a], ((b - t + a) / (b - a))) * pow(1.0 + z, -b * (t - a) / (b - a));
          printf("[%i] = %.4f ", t, dfs[t]);
          t++;
        }

        printf("\n");

        d = f(z, pars[b], dfs, b) / f_prime(z, pars[b], dfs, a, b);
        printf("d = %.6f, z = %.6f\n", d, z);
        z = z - d;
        if (runs++ > 20)
        {
          break;
        }
      } while (abs(d) > epsilon);

      i = b + 1;
    }
    else
    {
      double sum = 0.0;
      int a = 1;
      int b = i;
      while (a < b)
      {
        sum += dfs[a];
        a++;
      }
      printf("Sum is %.6f\n", sum);
      dfs[i] = (1.0 - pars[i] * sum) / (1.0 + pars[i]);
      i++;
    }
  }
  return dfs;
}
