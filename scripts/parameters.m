clear all;

global rho_min;  rho_min = -0.3;
global rho_step; rho_step = 0.05;
global rho_max;  rho_max = 0.3;
global sig_min;  sig_min = 1.0;
global sig_step; sig_step = 0.05;
global sig_max;  sig_max = 2.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Credit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf("Credit...\n");
prices = readmatrix('../data/highyieldexcessreturn1d.csv');
closingPrices = prices(:, 2);
log_returns = log(closingPrices(2:size(closingPrices)) ./ closingPrices(1:size(closingPrices)-1));
dt = 1 / 12;

[mu, vol, gam, rho, sig, p] = computeGenericParams(log_returns, dt);
fprintf("mu = %.2f, vol = %.2f, gam = %.4f, rho = %.2f, sig = %.2f, p = %.2f \n\n", mu, vol, gam, rho, sig, p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alternative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf("Alternative...\n");
prices = readmatrix('../data/nhx1m.csv');
closingPrices = prices(:, 1);
log_returns = log(closingPrices(2:size(closingPrices)) ./ closingPrices(1:size(closingPrices)-1));
dt = 1 / 12;

[mu, vol, gam, rho, sig, p] = computeGenericParams(log_returns, dt);
fprintf("mu = %.2f, vol = %.2f, gam = %.4f, rho = %.2f, sig = %.2f, p = %.2f \n\n", mu, vol, gam, rho, sig, p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Real estate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf("Real estate...\n");
prices = readmatrix('../data/hox1m.csv');
closingPrices = prices(:, 1);
log_returns = log(closingPrices(2:size(closingPrices)) ./ closingPrices(1:size(closingPrices)-1));
dt = 1 / 12;

[mu, vol, gam, rho, sig, p] = computeGenericParams(log_returns, dt);
fprintf("mu = %.2f, vol = %.2f, gam = %.4f, rho = %.2f, sig = %.2f, p = %.2f \n\n", mu, vol, gam, rho, sig, p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global equity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf("Global equity...\n");
prices = readmatrix('../data/msciworld1d-filtered.csv');
closingPrices = prices(:, 2);
log_returns = log(closingPrices(2:size(closingPrices)) ./ closingPrices(1:size(closingPrices)-1));
dt = 1 / 365;

[mu, vol, gam, rho, sig, p] = computeGenericParams(log_returns, dt);
fprintf("mu = %.2f, vol = %.2f, gam = %.4f, rho = %.2f, sig = %.2f, p = %.2f \n\n", mu, vol, gam, rho, sig, p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Domestic Equity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf("Domestic Equity...\n");
prices = readmatrix('../data/omxs301d.csv');
closingPrices = prices(:, 2);
log_returns = log(closingPrices(2:size(closingPrices)) ./ closingPrices(1:size(closingPrices)-1));
dt = 1 / 365;

[mu, vol, gam, rho, sig, p] = computeGenericParams(log_returns, dt);
fprintf("mu = %.2f, vol = %.2f, gam = %.4f, rho = %.2f, sig = %.2f, p = %.2f \n\n", mu, vol, gam, rho, sig, p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [best_mu, best_vol, best_gam, best_rho, best_sig, best_p] = computeGenericParams(log_changes, dt)

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
  best_mu = 0.0;
  best_vol = 0.0;
  best_gam = 0.0;
  best_rho = 0.0;
  best_sig = 0.0;

  for rho = rho_min : rho_step : rho_max
    for sig = sig_min : sig_step : sig_max
      gam = computeGenericGamma(rho, sig, vol, dt);

      for i = 1 : n_samples
        mu = (mean(log_changes) - gam) / dt;
        x2(i) = mu * dt - gam + vol * sqrt(dt) * sampleNln(rho, sig);
      end

      [h, p] = kstest2(x1, x2);
      if p > best_p
        best_p = p;
        best_mu = mu;
        best_vol = vol;
        best_gam = gam;
        best_rho = rho;
        best_sig = sig;
      end
    end
  end

end

function nln_mean = nlnMean(rho, sig)
  nln_mean = 1/2 * rho * sig * exp(1/8 * sig^2);
end

function nln_var = nlnVar(rho, sig)
  nln_var = exp(1/2 * sig^2) * (1 + rho^2 * sig^2 * (1 - 1/4 * exp(-1/4 * sig^2)));
end

function u_stnd = standardize_nln(u, rho, sig)
  u_stnd = (u - nlnMean(rho, sig)) / sqrt(nlnVar(rho, sig));
end

function u_norm = sampleNln(rho, sig)
  R = mvnrnd([0.0; 0.0], [sig, sig * rho; sig * rho, sig^2]);
  u = exp(0.5 * R(1)) * R(2);
  u_norm = standardize_nln(u, rho, sig);
end

function gamma = computeForwardGamma(rho, sig, vol, dt, eigenvalues, eigenvectors)
  n_trials = 5000;
  gamma = zeros(n_trials, 1);
  j = 1; % 1-year forward rate
  for x = 1:n_trials
      epsilon = zeros(3, 1);
      for i = 1:3
        u = sampleNln(rho, sig);
        epsilon(i) = u * lambda * sqrt(dt) * eigenvectors(j, i);
      end
      gamma(x) = sum(epsilon);
  end
  gamma = log(mean(exp(gamma)));
end

function gamma = computeGenericGamma(rho, sig, vol, dt)
  n_trials = 1000;
  gamma = zeros(n_trials, 1);
  for x = 1:n_trials
      u = sampleNln(rho, sig);
      gamma(x) = vol * sqrt(dt) * u;
  end
  gamma = log(mean(exp(gamma)));
end