#include <cmath>

int get_first_node_in_level(int level)
{
  int l = level;
  int b = 2;
  if (l > 0)
  {
    return ((b * (pow(b, l - 1) - 1)) / (b - 1)) + 1;
  }
  return 0;
}

int get_nodes_in_level(int level)
{
  return 1 << level;
}

int get_parent_index(int child_index, int items_per_node)
{
  if (child_index == 0)
  {
    return 0;
  }
  return (((child_index / items_per_node) - 1) / 2) * items_per_node;
}