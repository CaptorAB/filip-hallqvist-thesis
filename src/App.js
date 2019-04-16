import React, { Component } from "react";
import { Pane, Heading } from "evergreen-ui";
import { Header } from "./Header/Header";
import { Sidebar } from "./Sidebar/Sidebar";
import { Simulator } from "./Simulator/Simulator";
import "./App.css";
import { Provider, Consumer } from "./Libcapgen";

const App = () => (
  <Provider>
    <Consumer>
      {({ loading, libcapgen }) =>
        loading ? (
          <Heading size={900}>Loading...</Heading>
        ) : (
          <Pane height="100vh" display="flex" flexDirection="column">
            <Header />
            <Pane display="flex" flex={1}>
              <Pane flex={0}>
                <Sidebar />
              </Pane>
              <Pane flex={1}>
                <Simulator />
              </Pane>
            </Pane>
          </Pane>
        )
      }
    </Consumer>
  </Provider>
);

export default App;
