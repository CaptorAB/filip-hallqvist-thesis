import { Container } from "unstated";

import { libcapgen } from "./libcapgen";
import { format } from "date-fns";
import { N_INSTRUMENTS } from "./constants";

/*
const _preOrderTraversal = (array, current, path, paths) => {
  path.push(array[current]);
  const left = 2 * current + 1;
  const right = 2 * current + 2;
  if (left <= array.length && right <= array.length) {
    _preOrderTraversal(array, left, [...path], paths);
    _preOrderTraversal(array, right, [...path], paths);
  } else {
    paths.push(path);
  }
};

const preOrderTraversal = array => {
  const paths = [];
  _preOrderTraversal(array, 0, [], paths);
  return paths;
};
*/

export class SimulatorContainer extends Container {
  state = {
    history: [],
    loading: false,
    metrics: {
      fitness: 0,
      expectedReturn: 0,
      expectedRisk: 0
    },
    weights: []
  };
  optimize = async options => {
    const instance = await libcapgen();
    const result = instance.optimize(options);
    const metrics = {
      fitness: result.fitness,
      expectedReturn: result.expectedReturn,
      expectedRisk: result.expectedRisk
    };

    const weights = [];
    for (let i = 0; i < N_INSTRUMENTS; ++i) {
      weights.push(result.individual.get(i));
    }

    const state = {
      title: "Backtest",
      timestamp: format(new Date(), "yyyy-MM-dd HH:mm:ss"),
      metrics: { ...metrics },
      weights: [...weights]
    };

    this.setState({
      ...this.state,
      ...state,
      history: [state, ...this.state.history]
    });
  };
}
