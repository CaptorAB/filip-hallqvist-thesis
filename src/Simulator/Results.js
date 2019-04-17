import React from "react";
import { Pane, Heading, withTheme } from "evergreen-ui";

import { Tabs, Tab } from "../Tabs/Tabs";
import { Subscribe } from "unstated";
import { SimulatorContainer } from "../SimulatorContainer";
import { PerformancePlot } from "./PerformancePlot";

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
        <Tab title="Metrics">
          <Pane display="flex" flexWrap="wrap" marginX={-8} marginY={-8}>
            <Subscribe to={[SimulatorContainer]}>
              {simulator => (
                <>
                  <Metric
                    title="Expected Return"
                    value={simulator.state.metrics.totalReturn.toFixed(2)}
                  />
                  <Metric
                    title="Expected Risk"
                    value={simulator.state.metrics.risk.toFixed(2)}
                  />
                </>
              )}
            </Subscribe>
          </Pane>
        </Tab>
        <Tab title="Performance">
          <Subscribe to={[SimulatorContainer]}>
            {simulator => (
              <PerformancePlot
                scenarios={simulator.state.scenarios}
                current={0}
                title="Domestic Equity"
              />
            )}
          </Subscribe>
        </Tab>
        <Tab title="Weights">Weights</Tab>
      </Tabs>
    </Pane>
  </>
));

const Metric = ({ title, value, children, ...rest }) => (
  <Pane
    background="tint1"
    padding={16}
    flexBasis="28%"
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
