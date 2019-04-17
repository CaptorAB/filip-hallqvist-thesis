import React from "react";
import { Pane, Heading, withTheme } from "evergreen-ui";

export const Header = withTheme(({ theme, ...rest }) => (
  <Pane
    display="flex"
    padding={16}
    background={theme.palette.blue.base}
    {...rest}
  >
    <Pane flex={1} alignItems="center" display="flex">
      <Heading color={theme.colors.border.default}>Capgen</Heading>
    </Pane>
    <Pane />
  </Pane>
));
