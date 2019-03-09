const webpack = require("webpack");
const path = require("path");

module.exports = {
  mode: "production",
  context: path.resolve(__dirname, "."),
  entry: "./app/index.js",
  output: {
    path: path.resolve(__dirname, "dist"),
    filename: "bundle.js"
  },
  node: {
    fs: "empty"
  },
  module: {
    rules: [
      {
        test: /bin\/libcapgen\.wasm$/,
        type: "javascript/auto",
        loader: "file-loader",
        options: {
          publicPath: "dist/"
        }
      }
    ]
  }
};
