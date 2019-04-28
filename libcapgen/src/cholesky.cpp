#include <cmath>
#include <vector>

/**
 * Compute the cholesky decomposition of a matrix.
 *
 * @see https://rosettacode.org/wiki/Cholesky_decomposition#C
 */
std::vector<double> compute_cholesky(std::vector<double> &A, const int n)
{
  std::vector<double> L(n * n);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < (i + 1); j++)
    {
      double s = 0;
      for (int k = 0; k < j; k++)
      {
        s += L[i * n + k] * L[j * n + k];
      }

      if (i == j)
      {
        L[i * n + j] = std::sqrt(A[i * n + i] - s);
      }
      else
      {
        L[i * n + j] = (1.0 / L[j * n + j] * (A[i * n + j] - s));
      }
    }
  }
  return L;
}