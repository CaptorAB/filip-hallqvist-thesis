import init from "../bin/libcapgen.js";
import wasm from "../bin/libcapgen.wasm";

const libcapgen = init({
  locateFile: path => (path.endsWith("libcapgen.wasm") ? wasm : path)
});

libcapgen.onRuntimeInitialized = () => {
  window.libcapgen = libcapgen;
};
