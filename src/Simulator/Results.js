import React from "react";
import { Pane, Heading, Button, withTheme } from "evergreen-ui";

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
  </>
));
