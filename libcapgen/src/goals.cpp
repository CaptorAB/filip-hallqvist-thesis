
#include <tuple>
#include <vector>

#include <include/constants.h>

std::tuple<std::vector<double>, std::vector<double>>
generate_goals(std::vector<double> &instrument_changes,
               const double initial_funding_ratio,
               const double target_funding_ratio, const int n_scenarios,
               const int n_instruments)
{
  std::vector<double> intermediate_goals(n_scenarios);
  std::vector<double> final_goals(n_scenarios / 2 + 1);

  intermediate_goals[0] = 1.0 / initial_funding_ratio;
  int final_index = 0;

  // Iterate through scenarios
  for (int i = 0; i < n_scenarios; ++i)
  {
    int current = i;

    int left = 2 * current + 1;
    int right = 2 * current + 2;

    // Index of current node
    const double current_goal = intermediate_goals[current];
    const double current_change = instrument_changes[BONDS_20Y_INDEX];

    if (left < n_scenarios && right < n_scenarios)
    {
      const int lx = left * n_instruments;
      const int rx = right * n_instruments;

      // Evaluate left child
      intermediate_goals[left] =
          target_funding_ratio * current_change * current_goal;
      intermediate_goals[right] =
          target_funding_ratio * current_change * current_goal;
    }
    else
    {
      final_goals[final_index] =
          target_funding_ratio * current_change * current_goal;
      final_index++;
    }
  }

  return std::make_tuple(intermediate_goals, final_goals);
}