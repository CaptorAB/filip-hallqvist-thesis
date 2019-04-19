import React from "react";
import {
  Pane,
  Small,
  TextInput,
  Text,
  Heading,
  TextInputField,
  withTheme
} from "evergreen-ui";
import { INSTRUMENT_NAMES } from "../constants";

import { Tabs, Tab } from "../Tabs/Tabs";

export const Parameters = withTheme(
  ({ theme, values, handleChange, setFieldValue }) => {
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
              <PortfolioParameters
                values={values}
                handleChange={handleChange}
              />
            </Tab>
            <Tab title="Simulation">
              <SimulationParameters
                values={values}
                handleChange={handleChange}
              />
            </Tab>
            <Tab title="Optimizer">
              <OptimizerParameters
                values={values}
                handleChange={handleChange}
              />
            </Tab>
            <Tab title="Reallocations">
              <ReallocationParameters
                values={values}
                handleChange={handleChange}
              />
            </Tab>
            <Tab title="Constraints">
              <ConstraintParameters
                values={values}
                handleChange={handleChange}
              />
            </Tab>
          </Tabs>
        </Pane>
      </>
    );
  }
);

const Parameter = props => (
  <TextInputField flexBasis="28%" marginX={8} {...props} />
);

export const ConstraintParameters = ({ values, handleChange }) => (
  <Pane display="flex" flexWrap="wrap" marginX={8}>
    <Pane>
      <Heading marginBottom={8}>Instrument allocation constraints</Heading>
      <Text>
        Define constraints on how large proportion of the portfolio might be
        invested in a particular instrument.
      </Text>
      <Pane is="table" marginX={-4} marginY={16}>
        <thead>
          <tr>
            <Text is="th">Instrument</Text>
            <Text is="th">Min</Text>
            <Text is="th">Max</Text>
          </tr>
        </thead>
        <tbody>
          {values.instrumentConstraints.map((allocationConstraint, index) => (
            <tr key={index}>
              <td>
                <Small>{INSTRUMENT_NAMES[index]}</Small>
              </td>
              <td>
                <TextInput
                  type="number"
                  name={`instrumentConstraints.${index}.0`}
                  value={values.instrumentConstraints[index][0]}
                  onChange={handleChange}
                  step="any"
                  min={0.0}
                  max={1.0}
                />
              </td>
              <td>
                <TextInput
                  type="number"
                  name={`instrumentConstraints.${index}.1`}
                  value={values.instrumentConstraints[index][1]}
                  onChange={handleChange}
                  step="any"
                  min={0.0}
                  max={1.0}
                />
              </td>
            </tr>
          ))}
        </tbody>
      </Pane>
    </Pane>
  </Pane>
);

export const ReallocationParameters = ({ values, handleChange }) => (
  <Pane display="flex" flexWrap="wrap" marginX={-8}>
    {values.transactionCosts.map((transactionCost, index) => (
      <Parameter
        key={index}
        name={`transactionCosts.${index}`}
        label={INSTRUMENT_NAMES[index]}
        placeholder="0.0"
        value={values.transactionCosts[index]}
        type="number"
        onChange={handleChange}
        step="any"
        min="0"
        max="1"
      />
    ))}
  </Pane>
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
        step="any"
        min="0"
        max="1"
      />
      <Parameter
        label="Initial funding ratio"
        hint="Initial funding ratio of the portfolio."
        placeholder="1.3"
        value={values.initialFundingRatio}
        name="initialFundingRatio"
        type="number"
        onChange={handleChange}
        step="any"
        min="0"
      />
      <Parameter
        label="Target funding ratio"
        hint="Penalize portfolios with funding ratio below this value."
        placeholder="1.3"
        value={values.targetFundingRatio}
        name="targetFundingRatio"
        type="number"
        onChange={handleChange}
        step="any"
        min="0"
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
        step="any"
        min="1"
        max="12"
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
        step="any"
        min="0"
      />
      <Parameter
        label="Elitism copies"
        hint="Keep a number of clones of the best individuals at iteration."
        placeholder="2"
        value={values.elitismCopies}
        name="elitismCopies"
        type="number"
        onChange={handleChange}
        step="any"
        min="0"
      />
      <Parameter
        label="Generations"
        hint="Iterate through this many generations."
        placeholder="100"
        value={values.generations}
        name="generations"
        type="number"
        onChange={handleChange}
        step="any"
        min="1"
      />
      <Parameter
        label="Mutation rate"
        hint="Probability of mutating a gene."
        placeholder="0.02"
        value={values.mutationRate}
        name="mutationRate"
        type="number"
        onChange={handleChange}
        step="any"
        min="0"
        max="1"
      />
      <Parameter
        label="Crossover rate"
        hint="Probability of crossing two individuals."
        placeholder="0.02"
        value={values.crossoverRate}
        name="crossoverRate"
        type="number"
        onChange={handleChange}
        step="any"
        min="0"
        max="1"
      />
    </Pane>
  );
};
