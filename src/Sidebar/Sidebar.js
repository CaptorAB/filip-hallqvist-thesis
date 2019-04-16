import React from "react";
import { Pane, Heading, withTheme } from "evergreen-ui";

const items = [3, 2, 1];

export const Sidebar = withTheme(({ theme, ...rest }) => (
  <Pane
    minHeight="100%"
    width="200px"
    background="tint1"
    borderRight={`1px solid ${theme.colors.border.default}`}
    {...rest}
  >
    {items.map(x => (
      <Pane key={x} borderBottom={`1px solid ${theme.colors.border.default}`}>
        Backtest #{x}
      </Pane>
    ))}
  </Pane>
));
