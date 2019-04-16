function [u] = nln(rho, sig)

mu = [0 0];
sigma = [1 (rho * sig); (rho * sig) sig^2];
random = mvnrnd(mu, sigma);

epsilon = random(1);
eta = random(2);

u = exp(eta / 2) * epsilon;

end