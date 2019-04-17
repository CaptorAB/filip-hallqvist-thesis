import React from "react";
import { Pane, Heading, Button, withTheme } from "evergreen-ui";
import { Subscribe } from "unstated";

import { Parameters } from "./Parameters";
import { Results } from "./Results";
import { SimulatorContainer } from "../SimulatorContainer";

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
        <Subscribe to={[SimulatorContainer]}>
          {simulator => (
            <Button
              appearance="primary"
              iconBefore="play"
              onClick={() =>
                simulator.optimize({
                  populationSize: 10,
                  elitismCopies: 2,
                  generations: 10,
                  mutationRate: 0.02,
                  crossoverRate: 0.02,
                  steps: 4,
                  riskAversion: 0.0,
                  initialFundingRatio: 2.0,
                  targetFundingRatio: 1.3
                })
              }
            >
              Run
            </Button>
          )}
        </Subscribe>
      </Pane>
    </Pane>
    <Parameters />
    <Results />
  </Pane>
));
