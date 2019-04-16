import React, { useState, useEffect } from "react";
import { Pane, Spinner } from "evergreen-ui";
import { Header } from "./Header/Header";
import { Sidebar } from "./Sidebar/Sidebar";
import { Simulator } from "./Simulator/Simulator";
import { libcapgen } from "./libcapgen";
import { Subscribe } from "unstated";
import { SimulatorContainer } from "./SimulatorContainer";

const HEADER_HEIGHT = "52px";

const App = () => {
  return (
    <Subscribe to={[SimulatorContainer]}>
      {simulator =>
        simulator.loading ? (
          <Spinner
            transform="translate(-50%, -50%)"
            position="absolute"
            left="50%"
            top="50%"
          />
        ) : (
          <Pane height="100vh" display="flex" flexDirection="column">
            <Header />
            <Pane display="flex" height={`calc(100% - ${HEADER_HEIGHT})`}>
              <Pane flex={0} height="100%">
                <Sidebar />
              </Pane>
              <Pane flex={1}>
                <Simulator />
              </Pane>
            </Pane>
          </Pane>
        )
      }
    </Subscribe>
  );
};

export default App;
