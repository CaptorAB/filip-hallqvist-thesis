import React from "react";
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

export const Parameters = withTheme(({ theme, values, handleChange }) => {
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
            <PortfolioParameters values={values} handleChange={handleChange} />
          </Tab>
          <Tab title="Simulation">
            <SimulationParameters values={values} handleChange={handleChange} />
          </Tab>
          <Tab title="Optimizer">
            <OptimizerParameters values={values} handleChange={handleChange} />
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

export const PortfolioParameters = ({ values, handleChange }) => {
  return (
    <Pane display="flex" flexWrap="wrap" marginX={-8}>
      <Parameter
        label="Risk aversion"
        hint="Set to 0.0 to optimize without considering risk."
        placeholder="0.5"
        value={values.riskAversion}
        name="riskAversion"
        type="number"
        onChange={handleChange}
      />
      <Parameter
        label="Initial funding ratio"
        hint="Initial funding ratio of the portfolio."
        placeholder="1.3"
        value={values.initialFundingRatio}
        name="initialFundingRatio"
        type="number"
        onChange={handleChange}
      />
      <Parameter
        label="Target funding ratio"
        hint="Penalize portfolios with funding ratio below this value."
        placeholder="1.3"
        value={values.targetFundingRatio}
        name="targetFundingRatio"
        type="number"
        onChange={handleChange}
      />
    </Pane>
  );
};

export const SimulationParameters = ({ values, handleChange }) => {
  return (
    <Pane display="flex" flexWrap="wrap" marginX={-8}>
      <Parameter
        label="Steps"
        hint="Number of steps to simulate."
        placeholder="4"
        value={values.steps}
        name="steps"
        type="number"
        onChange={handleChange}
      />
    </Pane>
  );
};

export const OptimizerParameters = ({ values, handleChange }) => {
  return (
    <Pane display="flex" flexWrap="wrap" marginX={-8}>
      <Parameter
        label="Population size"
        hint="Number of individuals in the population."
        placeholder="100"
        value={values.populationSize}
        name="populationSize"
        type="number"
        onChange={handleChange}
      />
      <Parameter
        label="Elitism copies"
        hint="Keep a number of clones of the best individuals at iteration."
        placeholder="2"
        value={values.elitismCopies}
        name="elitismCopies"
        type="number"
        onChange={handleChange}
      />
      <Parameter
        label="Generations"
        hint="Iterate through this many generations."
        placeholder="100"
        value={values.generations}
        name="generations"
        type="number"
        onChange={handleChange}
      />
      <Parameter
        label="Mutation rate"
        hint="Probability of mutating a gene."
        placeholder="0.02"
        value={values.mutationRate}
        name="mutationRate"
        type="number"
        onChange={handleChange}
      />
      <Parameter
        label="Crossover rate"
        hint="Probability of crossing two individuals."
        placeholder="0.02"
        value={values.crossoverRate}
        name="crossoverRate"
        type="number"
        onChange={handleChange}
      />
    </Pane>
  );
};
