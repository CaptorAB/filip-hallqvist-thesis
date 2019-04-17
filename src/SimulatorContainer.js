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
      totalReturn: 0,
      risk: 0
    },
    scenarios: []
  };
  optimize = async options => {
    const instance = await libcapgen();
    const result = instance.optimize(options);
    const metrics = {
      totalReturn: result.totalReturn,
      risk: result.risk
    };

    const events = [];
    const nEvents = 2 ** options.steps - 1;

    // Build events from memory linear arrays
    for (let i = 0; i < nEvents; i++) {
      const ix = i * N_INSTRUMENTS;
      const event = {
        priceChanges: [],
        goal: result.goals.get(i),
        probability: result.probabilities.get(i),
        weights: []
      };

      // Update instruments and weights
      for (let j = 0; j < N_INSTRUMENTS; j++) {
        const jx = ix + j;
        event.priceChanges.push(result.priceChanges.get(jx));
        event.weights.push(result.individual.get(jx));
      }

      events.push(event);
    }

    // Update scenarios
    // const scenarios = preOrderTraversal(events);
    const scenarios = [];

    this.setState({
      ...this.state,
      history: [
        {
          title: "Backtest",
          timestamp: format(new Date(), "yyyy-MM-dd HH:mm:ss"),
          metrics: { ...metrics }
        },
        ...this.state.history
      ],
      metrics: { ...metrics },
      scenarios
    });
  };
}
