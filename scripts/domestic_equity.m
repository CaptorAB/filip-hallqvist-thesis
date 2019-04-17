clear;


prices = readmatrix('../data/omxs301d.csv');
closingPrices = prices(:, 2);
changes = diff(closingPrices)./closingPrices(1:size(closingPrices)-1);

rho = -0.02;
sigma = 1.1;
S_0 = 1.0;
T = 500;

u = zeros(T, 1);
for i = 1:T
    u(i) = nln(rho, sigma);
end

qqplot(changes);

mu_eq = mean(changes);
sigma_eq = var(changes);

c1 = 1/2 * rho * sigma * exp(1/8 * sigma^2);
c2 = exp(1/2 * sigma^2) * (1 + rho^2 * sigma^2 * (1 - 1/4 * exp(-1/4 * sigma^2)));

u_norm = (u - c1) / sqrt(c2);

epsilon = zeros(T, 1);
F = zeros(T, 1);
S = zeros(T, 1);

for t = 1:T
    epsilon(t) = sigma_eq * sqrt(t) * u_norm(t);
end

for t = 1:T
    epsilon(t) = sigma_eq * sqrt(t) * u_norm(t);
    
    gamma = zeros(100, 1);
    for x = 1:100
        u_normed = (nln(rho, sigma) - c1) / sqrt(c2);
        gamma(x) = sigma_eq * sqrt(t) * u_normed;
    end
    gamma = log(mean(exp(gamma)));
    
    F(t) = S_0 * exp(mu_eq * sqrt(t));
    S(t) = F(t) * exp(-gamma + epsilon(t));
end

generated_changes = diff(S)./S(1:size(S)-1);

% hold on;
figure
qqplot(generated_changes);
% plot(S);