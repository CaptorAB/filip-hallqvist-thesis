import React from "react";
import { Pane, Heading, withTheme } from "evergreen-ui";

import { Tabs, Tab } from "../Tabs/Tabs";
import { Subscribe } from "unstated";
import { SimulatorContainer } from "../SimulatorContainer";
import { Plot } from "../Plot";
import { INSTRUMENT_NAMES } from "../constants";

export const Results = withTheme(({ theme }) => (
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
        <Heading size={200}>Results</Heading>
      </Pane>
    </Pane>
    <Pane display="flex" padding={16} flexDirection="column">
      <Tabs>
        <Tab title="Weights">
          <Pane display="flex" flexWrap="wrap" marginX={-8} marginY={-8}>
            <Subscribe to={[SimulatorContainer]}>
              {simulator =>
                simulator.state.weights.map((weight, index) => (
                  <Card
                    key={index}
                    title={INSTRUMENT_NAMES[index]}
                    value={weight.toFixed(2)}
                  />
                ))
              }
            </Subscribe>
          </Pane>
        </Tab>
        <Tab title="Metrics">
          <Pane display="flex" flexWrap="wrap" marginX={-8} marginY={-8}>
            <Subscribe to={[SimulatorContainer]}>
              {simulator => (
                <>
                  <Card
                    title="Fitness"
                    value={simulator.state.metrics.fitness.toFixed(2)}
                  />
                  <Card
                    title="Expected Return"
                    value={simulator.state.metrics.expectedReturn.toFixed(2)}
                  />
                  <Card
                    title="Expected Risk"
                    value={simulator.state.metrics.expectedRisk.toFixed(2)}
                  />
                </>
              )}
            </Subscribe>
          </Pane>
        </Tab>
        <Tab title="Scenarios">
          <Pane display="flex" flexWrap="wrap" marginX={-8} marginY={-8}>
            <Subscribe to={[SimulatorContainer]}>
              {simulator => (
                <>
                  <Card title="Expected Returns">
                    <Plot
                      data={[
                        {
                          x: simulator.state.finalReturns,
                          type: "histogram"
                        }
                      ]}
                    />
                  </Card>
                  <Card title="Best Scenario vs. Worst Scenario">
                    <Plot
                      data={[
                        {
                          y: simulator.state.bestPath,
                          type: "scatter",
                          name: "Best Scenario"
                        },
                        {
                          y: simulator.state.worstPath,
                          type: "scatter",
                          name: "Worst Scenario"
                        }
                      ]}
                    />
                  </Card>
                </>
              )}
            </Subscribe>
          </Pane>
        </Tab>
      </Tabs>
    </Pane>
  </>
));

const Card = ({ title, value, children, ...rest }) => (
  <Pane
    background="tint1"
    padding={16}
    flexBasis="15%"
    marginX={8}
    marginY={8}
    border="default"
    {...rest}
  >
    <Heading size={200} marginBottom={8}>
      {title}
    </Heading>
    {typeof value === "undefined" ? (
      children
    ) : (
      <Heading size={800}>{value}</Heading>
    )}
  </Pane>
);
