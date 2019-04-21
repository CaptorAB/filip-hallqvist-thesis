import React from "react";
import { Pane, Heading, Button, withTheme } from "evergreen-ui";
import { Subscribe } from "unstated";
import { Formik } from "formik";

import { Parameters } from "./Parameters";
import { Results } from "./Results";
import { SimulatorContainer } from "../SimulatorContainer";

export const Simulator = withTheme(({ theme, ...rest }) => (
  <Subscribe to={[SimulatorContainer]}>
    {simulator => (
      <Formik
        initialValues={{
          populationSize: 10,
          elitismCopies: 2,
          generations: 10,
          mutationRate: 0.02,
          crossoverRate: 0.02,
          steps: 4,
          riskAversion: 0.0,
          initialFundingRatio: 1.3,
          targetFundingRatio: 1.3,
          transactionCosts: [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
          ],
          instrumentConstraints: [
            [0.0, 1.0],
            [0.0, 1.0],
            [0.0, 1.0],
            [0.0, 1.0],
            [0.0, 1.0],
            [0.0, 1.0],
            [0.0, 1.0],
            [0.0, 1.0],
            [0.0, 1.0],
            [0.0, 1.0],
            [0.0, 1.0],
            [0.0, 1.0],
            [0.0, 1.0]
          ],
          marginConstraints: [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
          ]
        }}
        onSubmit={(values, { setSubmitting }) => {
          const instrumentConstraints = [
            ...values.instrumentConstraints.map(x => x[0]),
            ...values.instrumentConstraints.map(x => x[1])
          ];
          simulator.optimize({
            ...values,
            instrumentConstraints
          });
        }}
      >
        {({
          isSubmitting,
          handleSubmit,
          values,
          handleChange,
          setFieldValue
        }) => (
          <form onSubmit={handleSubmit}>
            <Pane display="flex" flexDirection="column" {...rest}>
              <Pane
                borderBottom={`1px solid ${theme.colors.border.default}`}
                display="flex"
                padding={16}
                alignItems="center"
              >
                <Pane flex={1}>
                  <Heading>Backtest</Heading>
                </Pane>
                <Pane>
                  <Button appearance="primary" iconBefore="play" type="submit">
                    Run
                  </Button>
                </Pane>
              </Pane>
              <Parameters
                values={values}
                handleChange={handleChange}
                setFieldValue={setFieldValue}
              />
              <Results />
            </Pane>
          </form>
        )}
      </Formik>
    )}
  </Subscribe>
));
