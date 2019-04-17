clear;

T = 100;

dailyParRates = dlmread("../data/swap1d.csv", ",", 1, 1);
dailyParRates = dailyParRates(:, 1:10) / 100; % Skip 12y, 15y, and 20y for now

dailyRows = size(dailyParRates, 1);

% Monthly par rates

monthlyParRates = dailyParRates(1:30:dailyRows, :);
cols = size(monthlyParRates, 2);
rows = size(monthlyParRates, 1);

% Discount factors

df = zeros(rows, cols);

for i = 1:rows
    s = 0;
    df(i, 1) = 1 / (1 + monthlyParRates(i, 1));
    
    for j = 2:cols
        s = s + df(i, j - 1);
        p = monthlyParRates(i, j);
        df(i, j) = (1 - p * s) / (1 + p);
    end
end

% Forward rates

f = zeros(rows, cols);
f_0 = 1;

for i = 1:rows
    f(i, 1) = - (log(df(i, 1)) - log(f_0));
    for j = 2:cols
        f(i, j) = - (log(df(i, j)) - log(df(i, j-1)));
    end
end

% Adjust for negative rate
f = f + 0.01;

% Forward rate changes

changes = real(log(f(2:rows, :) ./ f(1:rows-1, :)));

% Normalize
changes = normalize(changes); % TODO: Normalize assuming nln?

% PCA on changes: coeff pc vectors, score is new data, latent is variance
[coeff, score, latent] = pca(changes);

% Compute standard deviation
stds = sqrt(latent);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on;
for ir = 1:1
n_rands = 3;

rho = -0.1;
sigma = 1.7;
S_0 = 1.0;
T = 500;

u = zeros(T, n_rands);
for i = 1:T
    for j = 1:n_rands
        u(i, j) = nln(rho, sigma);
    end
end

mu_ir = mean(changes(:,ir));
sigma_ir = var(changes(:,ir));

u_norm = zeros(T, n_rands);

c1 = 1/2 * rho * sigma * exp(1/8 * sigma^2);
c2 = exp(1/2 * sigma^2) * (1 + rho^2 * sigma^2 * (1 - 1/4 * exp(-1/4 * sigma^2)));

for i = 1:n_rands
    u_norm(:, i) = (u(:, i) - c1) / sqrt(c2);
end

epsilon = zeros(T, 1);
F = zeros(T, 1);
S = zeros(T, 1);

for t = 1:T
    for i = 1:n_rands
        epsilon(t) = epsilon(t) + u_norm(t, i) * stds(ir) * coeff(ir, i) * sqrt(t);
    end
    
    gamma = zeros(10000, 1);
    for x = 1:10000
        for i = 1:n_rands
            gamma(x) = gamma(x) + u_norm(t, i) * stds(ir) * coeff(ir, i) * sqrt(t);
        end
    end
    gamma = log(mean(exp(gamma)));
    
    F(t) = S_0 * exp(mu_ir(ir) * sqrt(t));
    S(t) = F(t) * exp(-gamma + epsilon(t));
end

% generated_changes = diff(S)./S(1:size(S)-1);

plot(S);

% plot();
end % for k
