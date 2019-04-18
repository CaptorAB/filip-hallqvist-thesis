import React, { useState } from "react";
import { Tablist, Pane, Tab as EvergreenTab } from "evergreen-ui";

export const TabPane = ({ children, ...rest }) => (
  <Pane marginTop={16} {...rest}>
    {children}
  </Pane>
);

export const Tab = ({ children, ...rest }) => (
  <EvergreenTab {...rest}>{children}</EvergreenTab>
);

export const Tabs = ({ children }) => {
  const [selected, setSelected] = useState(0);
  return (
    <>
      <Tablist marginX={-4}>
        {React.Children.map(children, (x, i) => {
          const { title, children, ...rest } = x.props;
          return (
            <Tab
              key={i}
              is="a"
              isSelected={i === selected}
              onSelect={() => setSelected(i)}
              {...rest}
            >
              {title}
            </Tab>
          );
        })}
      </Tablist>
      {React.Children.map(
        children,
        (x, i) =>
          selected === i && <TabPane key={i}>{x.props.children}</TabPane>
      )}
    </>
  );
};
