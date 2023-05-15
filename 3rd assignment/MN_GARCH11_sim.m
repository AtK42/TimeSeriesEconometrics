function [sim_values, sim_sigma2] = MN_GARCH11_sim(len_path, mu, wghts, gamma0, gamma, Psi)
%--------------------------------------------------------------------------
% function to simulate from a mixed-normal GARCH(1,1) model
%--------------------------------------------------------------------------
% INPUT
% % len_path                integer, length of the path (T)
% % mu                      kx1 vector, mean vector
% % wghts                   kx1 vector, weights for the mixture
% % gamma0                  kx1 vector, long-term mean in GARCH model
% % gamma                   kx1 vector, GARCH-parameter
% % Psi                     kxk matrix, ARCH-parameter
%--------------------------------------------------------------------------
% OUTPUT
% % sim_values              len_path x 1 vector, simulated path
% % sim_sigma2              k x len_path vector, sigma^2 process
%--------------------------------------------------------------------------

k = length(mu);

% check if weights sum to one, otherwise adjust them
if sum(wghts) ~= 1
    wghts = wghts./sum(wghts);
end

% ensure that the rv has zero mean
mu_adj = mu;
mu_adj(end) = -sum(wghts(1:end-1)./wghts(end) .* mu(1:end-1));

% initialize variables
sim_values = zeros(len_path, 1);
sim_sigma2 = zeros(k, len_path);

% set initial values
sim_sigma2(:, 1) = gamma0;
sim_values(1) = wghts' * gamma0;

% simulation of the paths
for i_sim = 2:len_path
    sim_sigma2(:, i_sim) = gamma0 + gamma * sim_values(i_sim-1)^2 + Psi * sim_sigma2(:, i_sim-1);
    sim_values(i_sim) = wghts' * normrnd(mu_adj, sqrt(sim_sigma2(:, i_sim)));
end