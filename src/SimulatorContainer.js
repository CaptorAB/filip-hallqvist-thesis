import { Container } from "unstated";
import { libcapgen } from "./libcapgen";
import { format } from "date-fns";

export class SimulatorContainer extends Container {
  state = {
    history: [],
    loading: false,
    metrics: {
      totalReturn: 0,
      risk: 0
    }
  };
  optimize = async (...args) => {
    const instance = await libcapgen();
    const result = instance.optimize(...args);
    const metrics = {
      totalReturn: result.totalReturn,
      risk: result.risk
    };
    this.setState({
      ...this.state,
      history: [
        {
          title: "Backtest",
          timestamp: format(new Date(), "yyyy-MM-dd HH:mm:ss"),
          metrics: { ...metrics }
        },
        ...this.state.history
      ],
      metrics: { ...metrics }
    });
  };
}
