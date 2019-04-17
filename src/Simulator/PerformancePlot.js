import React from "react";
import { Pane, Heading, withTheme } from "evergreen-ui";
import Plotly from "plotly.js-basic-dist";
import createPlotlyComponent from "react-plotly.js/factory";

const Plot = createPlotlyComponent(Plotly);

export const PerformancePlot = withTheme(
  ({ theme, scenarios, current, title }) => {
    const y = [];
    for (const scenario of scenarios) {
      y.push(scenario.priceChanges[current]);
    }

    return (
      <Pane>
        <Heading size={200} marginBottom={16}>
          {title}
        </Heading>
        <Pane height="160px">
          <Plot
            data={[
              {
                y,
                type: "scatter",
                mode: "lines",
                marker: { color: theme.palette.blue.base }
              }
            ]}
            layout={{
              margin: {
                l: 32,
                r: 32,
                b: 32,
                t: 32,
                pad: 0
              }
            }}
            useResizeHandler={true}
            style={{ width: "100%", height: "100%" }}
          />
        </Pane>
      </Pane>
    );
  }
);
