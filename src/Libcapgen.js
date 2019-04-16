import React, { useState, useEffect } from "react";

const libcapgen = window.libcapgen();

const Context = React.createContext();

const INITIAL_STATE = {
  loading: true,
  libcapgen: null
};

export const Provider = ({ children }) => {
  const [state, setState] = useState(INITIAL_STATE);
  useEffect(() => {
    libcapgen.onRuntimeInitialized = () => {
      setState({
        loading: false,
        libcapgen
      });
    };
  }, []);
  return <Context.Provider value={state}>{children}</Context.Provider>;
};

export const Consumer = Context.Consumer;
