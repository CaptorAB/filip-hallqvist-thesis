# Constrained Portfolio Optimization in Liability-Driven Investing

_(preliminary title)_

We formulate the classic portfolio optimization problem as a multi-stage stochastic programming problem, solve it using a genetic algorithm. The application allows the user to add allocation constraints, tune the parameters of the genetic algorithm, and edit the target funding ratio of the optimal portfolio.

The application is the result of a Master's Thesis, which (hopefully) will be published in June 2019.

## Getting Started

The application can be found [here](https://captorab.github.io/filip-hallqvist-thesis).

### Development

1. Clone this repository.
2. Make sure [Docker](https://docs.docker.com/install/) is installed
3. Pull an image for compiling WebAssembly: `docker pull trzeci/emscripten`
4. Install depencencies: `npm install`
5. Compile the WebAssembly module: `npm run docker`
6. Start the application: `npm start`
7. Fire up `localhost:3000` in your browser.

## Running the tests

There are tests written for both the genetic algorithm (C++) and the user interface (JavaScript).

Tests are run in CircleCI, but if you want to run the tests locally, you can do so by executing the `docker` script in `package.json`, and then executing the `npm run test:app` script.

```sh
npm run docker && npm run test:app
```

As previously mentioned, you need to have [Docker](https://docs.docker.com/install/) for it to work.

## Deployment

The application is currently hosted on GitHub Pages. Deployment is done easiest by executing the `deploy` npm script.

```sh
npm run deploy
```

**Remember to build the WebAssembly modules before deployment: `npm run docker`.**

## Built With

-   [Emscripten](https://emscripten.org/)
-   [Evergreen](https://evergreen.segment.com/)
-   [Catch2](https://github.com/catchorg/Catch2/)
-   [Random](https://github.com/effolkronium/random/)

## Versioning

We use [SemVer](http://semver.org/) for versioning.

## Authors

-   **Filip Hallqvist**

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Acknowledgments

Special thanks to **Holger Rootzén** for his constructive suggestions and insights during the research work. I am also particularly grateful for the assistance given by **Tor Nordqvist** and **Martin Karrin** at Captor for their invaluable feedback and help throughout the project.
