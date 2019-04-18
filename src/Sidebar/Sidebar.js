import React from "react";
import {
  Pane,
  Heading,
  UnorderedList,
  ListItem,
  withTheme
} from "evergreen-ui";
import { Subscribe } from "unstated";
import { SimulatorContainer } from "../SimulatorContainer";

export const Sidebar = withTheme(({ theme, ...rest }) => (
  <Pane
    height="100%"
    overflowY="auto"
    width="200px"
    background="tint1"
    borderRight={`1px solid ${theme.colors.border.default}`}
    {...rest}
  >
    <Subscribe to={[SimulatorContainer]}>
      {simulator =>
        simulator.state.history.map((x, i) => (
          <Pane
            key={i}
            padding={16}
            borderBottom={`1px solid ${theme.colors.border.default}`}
          >
            <Heading size={200}>{x.title}</Heading>
            <Heading size={100}>{x.timestamp}</Heading>
            <UnorderedList marginBottom={0}>
              <Metric icon="symbol-triangle-up">
                {x.metrics.expectedReturn.toFixed(2)}
              </Metric>
              <Metric icon="symbol-triangle-down">
                {x.metrics.expectedRisk.toFixed(2)}
              </Metric>
            </UnorderedList>
          </Pane>
        ))
      }
    </Subscribe>
  </Pane>
));

const Metric = ({ children, ...rest }) => (
  <ListItem
    size={100}
    display="inline-block"
    marginTop={0}
    marginBottom={0}
    marginRight={32}
    paddingLeft={0}
    color="muted"
    fontSize={11}
    {...rest}
  >
    {children}
  </ListItem>
);
