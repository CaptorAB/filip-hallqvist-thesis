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
    weights: [],
    finalReturns: []
  };
  optimize = async options => {
    const instance = await libcapgen();
    const result = instance.optimize(options);
    const metrics = {
      fitness: result.fitness,
      expectedReturn: result.expectedReturn,
      expectedRisk: result.expectedRisk
    };

    const nFinalWealths = 2 ** (options.steps - 1);

    const finalReturns = [];
    let iBestScenario = -1;
    let iWorstScenario = -1;
    for (let i = 0; i < nFinalWealths; ++i) {
      const finalReturn = result.finalWealths.get(i) - 1;
      finalReturns.push(finalReturn);

      if (iBestScenario < 0 || finalReturn > finalReturns[iBestScenario]) {
        iBestScenario = i;
      }
      if (iWorstScenario < 0 || finalReturn < finalReturns[iWorstScenario]) {
        iWorstScenario = i;
      }
    }

    // Backtrack from best and worst scenarios
    // to find path from the root
    const bestPath = [result.finalWealths.get(iBestScenario)];
    let i = nFinalWealths + iBestScenario - 1;
    while (i >= 0) {
      bestPath.push(result.intermediateWealths.get(i));
      i = Math.floor((i - 1) / 2);
    }
    bestPath.reverse();

    // Do the same for worst path
    const worstPath = [result.finalWealths.get(iWorstScenario)];
    let j = nFinalWealths + iWorstScenario - 1;
    while (j >= 0) {
      worstPath.push(result.intermediateWealths.get(j));
      j = Math.floor((j - 1) / 2);
    }
    worstPath.reverse();

    const weights = [];
    for (let i = 0; i < N_INSTRUMENTS; ++i) {
      weights.push(result.individual.get(i));
    }

    console.log(weights.map(w => w.toFixed(2)).join(" "));

    const state = {
      title: "Backtest",
      timestamp: format(new Date(), "yyyy-MM-dd HH:mm:ss"),
      metrics: { ...metrics },
      weights: [...weights],
      finalReturns: [...finalReturns],
      bestPath: [...bestPath],
      worstPath: [...worstPath]
    };

    this.setState({
      ...this.state,
      ...state,
      history: [state, ...this.state.history]
    });
  };
}
