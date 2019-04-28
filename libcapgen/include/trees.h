#ifndef TREES_H
#define TREES_H

#include <vector>
#include <tuple>

struct ScenarioTree
{
  std::vector<double> instrument_changes;
  std::vector<double> intermediate_goals;
  std::vector<double> final_goals;
};

#endif