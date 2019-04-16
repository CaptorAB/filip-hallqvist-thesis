import React, { useState } from "react";
import { Pane, Heading, TextInputField, withTheme } from "evergreen-ui";

import { Tabs, Tab } from "../Tabs/Tabs";

export const INITIAL_PARAMETERS = {
  populationSize: 100,
  elitismCopies: 2,
  generations: 100,
  mutationRate: 0.02,
  crossoverRate: 0.02,
  steps: 4,
  riskAversion: 0.5,
  initialFundingRatio: 1.3,
  targetFundingRatio: 1.3
};

export const Parameters = withTheme(({ theme }) => {
  const [parameters, setParameters] = useState({ ...INITIAL_PARAMETERS });
  return (
    <>
      <Pane
        background="tint1"
        borderBottom={`1px solid ${theme.colors.border.default}`}
        display="flex"
        paddingX={16}
        paddingY={8}
        flexDirection="column"
      >
        <Pane flex={1}>
          <Heading size={200}>Parameters</Heading>
        </Pane>
      </Pane>
      <Pane
        borderBottom={`1px solid ${theme.colors.border.default}`}
        display="flex"
        padding={16}
        flexDirection="column"
      >
        <Tabs>
          <Tab title="Portfolio">
            <PortfolioParameters />
          </Tab>
          <Tab title="Simulation">
            <SimulationParameters />
          </Tab>
          <Tab title="Optimizer">
            <OptimizerParameters />
          </Tab>
          <Tab title="Correlations" disabled={true}>
            Correlations
          </Tab>
        </Tabs>
      </Pane>
    </>
  );
});

const Parameter = props => (
  <TextInputField flexBasis="28%" marginX={8} {...props} />
);

export const PortfolioParameters = () => {
  return (
    <Pane display="flex" flexWrap="wrap" marginX={-8}>
      <Parameter
        label="Risk aversion"
        hint="Set to 0.0 to optimize without considering risk."
        placeholder="0.5"
        value={0.5}
      />
      <Parameter
        label="Initial funding ratio"
        hint="Initial funding ratio of the portfolio."
        placeholder="1.3"
        value={1.3}
      />
      <Parameter
        label="Target funding ratio"
        hint="Penalize portfolios with funding ratio below this value."
        placeholder="1.3"
        value={1.3}
      />
    </Pane>
  );
};

export const SimulationParameters = () => {
  return (
    <Pane display="flex" flexWrap="wrap" marginX={-8}>
      <Parameter
        label="Steps"
        hint="Number of steps to simulate."
        placeholder="4"
        value={4}
      />
    </Pane>
  );
};

export const OptimizerParameters = () => {
  return (
    <Pane display="flex" flexWrap="wrap" marginX={-8}>
      <Parameter
        label="Population size"
        hint="Number of individuals in the population."
        placeholder="100"
        value={100}
      />
      <Parameter
        label="Elitism copies"
        hint="Keep a number of clones of the best individuals at iteration."
        placeholder="2"
        value={2}
      />
      <Parameter
        label="Generations"
        hint="Iterate through this many generations."
        placeholder="100"
        value={100}
      />
      <Parameter
        label="Mutation rate"
        hint="Probability of mutating a gene."
        placeholder="0.02"
        value={0.02}
      />
      <Parameter
        label="Crossover rate"
        hint="Probability of crossing two individuals."
        placeholder="0.02"
        value={0.02}
      />
    </Pane>
  );
};
