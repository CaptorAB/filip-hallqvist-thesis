import React from "react";
import { Pane, Heading, Button, withTheme } from "evergreen-ui";
import { Parameters } from "./Parameters";
import { Results } from "./Results";
import { Consumer as Libcapgen } from "../Libcapgen";

const runSimulation = libcapgen => {
  const result = libcapgen.optimize({
    populationSize: 10,
    elitismCopies: 2,
    generations: 10,
    mutationRate: 0.02,
    crossoverRate: 0.02,
    steps: 4,
    riskAversion: 0.5,
    initialFundingRatio: 1.3,
    targetFundingRatio: 1.3
  });

  console.log(result);
};

export const Simulator = withTheme(({ theme, ...rest }) => (
  <Pane display="flex" flexDirection="column" {...rest}>
    <Pane
      borderBottom={`1px solid ${theme.colors.border.default}`}
      display="flex"
      padding={16}
      alignItems="center"
    >
      <Pane flex={1}>
        <Heading>Backtest</Heading>
      </Pane>
      <Pane>
        <Libcapgen>
          {({ libcapgen }) => (
            <Button
              appearance="primary"
              iconBefore="play"
              onClick={() => runSimulation(libcapgen)}
            >
              Run
            </Button>
          )}
        </Libcapgen>
      </Pane>
    </Pane>
    <Parameters />
    <Results />
  </Pane>
));
