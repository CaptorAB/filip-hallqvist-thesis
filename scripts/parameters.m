clear all;

global rho_min;  rho_min = -0.3;
global rho_step; rho_step = 0.05;
global rho_max;  rho_max = 0.3;
global sig_min;  sig_min = 1.0;
global sig_step; sig_step = 0.05;
global sig_max;  sig_max = 2.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Domestic Equity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prices = readmatrix('../data/omxs301d.csv');
closingPrices = prices(:, 2);
log_returns = log(closingPrices(2:size(closingPrices)) ./ closingPrices(1:size(closingPrices)-1));
dt = 1 / 365;

[best_rho, best_sig] = computeParams(log_returns, dt);
fprintf("Domestic Equity: rho = %.2f, sig = %.2f \n", best_rho, best_sig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute 
function [best_rho, best_sig] = computeParams(log_changes, dt)

  global rho_min;
  global rho_step;
  global rho_max;
  global sig_min;
  global sig_step;
  global sig_max;

  n_samples = 5000;
  x1 = log_changes;
  x2 = zeros(n_samples, 1);
  vol = sqrt(var(log_changes) / dt);
  best_p = 0.0;
  best_rho = 0.0;
  best_sig = 0.0;

  for rho = rho_min : rho_step : rho_max
    for sig = sig_min : sig_step : sig_max
      gam = computeGamma(rho, sig, vol, dt);

      for i = 1 : n_samples
        mu = (mean(log_changes) - gam) / dt;
        x2(i) = mu * dt - gam + vol * sqrt(dt) * sampleNln(rho, sig);
      end

      [h, p] = kstest2(x1, x2);
      if p > best_p
        best_p = p;
        best_rho = rho;
        best_sig = sig;
      end
    end
  end

end

% Sample nln
function u_norm = sampleNln(rho, sig)
  R = mvnrnd([0.0; 0.0], [sig, sig * rho; sig * rho, sig^2]);
  u = exp(0.5 * R(1)) * R(2);
  c1 = 1/2 * rho * sig * exp(1/8 * sig^2);
  c2 = exp(1/2 * sig^2) * (1 + rho^2 * sig^2 * (1 - 1/4 * exp(-1/4 * sig^2)));
  u_norm = (u - c1) / sqrt(c2);
end

% Gamma
function gamma = computeGamma(rho, sig, vol, dt)
  n_trials = 5000;
  gamma = zeros(n_trials, 1);
  c1 = 1/2 * rho * sig * exp(1/8 * sig^2);
  c2 = exp(1/2 * sig^2) * (1 + rho^2 * sig^2 * (1 - 1/4 * exp(-1/4 * sig^2)));
  for x = 1:100
      u_norm = (sampleNln(rho, sig) - c1) / sqrt(c2);
      gamma(x) = vol * sqrt(dt) * u_norm;
  end
  gamma = log(mean(exp(gamma)));
end