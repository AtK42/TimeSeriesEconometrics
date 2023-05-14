% Here goes: Have a look at my mixed-normal GARCH construction, on page 478

% Your first task is to  * simulate *  this process. Use k=2 components, though
%   I might later ask you to try 3. You do NOT have to estimate the model, 
%   though it is also not difficult to set up the likelihood and estimate 
%   it. You just need to simulate it, and you can play around with the 
%   required parameters so that when you plot the resulting univariate time 
%   series, it should "look like actual daily stock returns".  Since you 
%   are simulating, you can simulate T=5000 observations -- or even more, 
%   all for free.

% The next part is to estimate a student-t-GARCH model based on your 
%   simulated data. Notice that the true DGP does NOT match the model you 
%   are estimating. That is okay -- that is reality. The student-t garch 
%   estimation is very easy to set up -- you just build on the codes I have
%   in my chapter. Obviously, the degrees of freedom of the student t gets 
%   jointly estimated with the garch parameters.  Surely, you can find 
%   toolboxes for matlab and python that also do it, ...)

% I think that is already a good start to the 3rd homework, and I will 
%   surely build on this, such as asking you to make similar VaR plots as 
%   I have on page 492. (So, you can do that now also, and thus I entice 
%   you to read some of chapter 11.)

%% 1) simulate a mixed-normal GARCH(1,1) construction
% set parameters of simulation
% % number of components of the mixed-normal model (corresponds to k)
n_comp = 3; % for k = 3, we need 14 free parameters
% % length of path (corresponds to T)
len_paths = 5e3;
% % number of paths
n_paths = 1e3;
% % number of lags for GARCH model, for this code fixed at (1,1)
%r=1; s=1;

% set the model parameters
% % mean
mu = [.22 -.31 .01]';
% % weights
wghts = [1, 3, 2]';
%wghts = 1/n_comp * ones(k, 1);
% % GARCH coefficients
gamma0 = .01 + zeros(3, 1);
gamma = [.1, .44, .26]';
Psi = [.4 .2 .3; .2 .8 .5; .3 .5 .9];

% checks
% % check no of components
if length(mu) ~= n_comp
    error(['mu is of the wrong dimension: should contain ', num2str(n_comp), ' elements but contains ', num2str(length(mu)), ' elements'])
%elseif length()
end
% % check if weights sum to one, otherwise adjust them
if sum(wghts) ~= 1
    wghts = wghts./sum(wghts);
end

% ensure that the rv has zero mean
mu_adj = mu;
mu_adj(end) = -sum(wghts(1:end-1)./wghts(end) .* mu(1:end-1));

% initialize variables
sim_values = zeros(len_paths, 1);
%sim_values = zeros(len_paths, n_paths);
sim_sigma2 = zeros(k, len_paths);

% set initial values
sim_sigma2(:, 1) = gamma0;
sim_values(1) = wghts' * gamma0;

% simulation of the paths
for i_sim = 2:len_paths
    sim_sigma2(:, i_sim) = gamma0 + gamma * sim_values(i_sim-1)^2 + Psi * sim_sigma2(:, i_sim-1);
    sim_values(i_sim) = wghts' * normrnd(mu_adj, sqrt(sim_sigma2(:, i_sim)));
end
% for i_path = 1:n_paths
%     for i_sim = 2:len_paths
%         sim_sigma = alpha_GARCH * sim_values(i_sim - 1, i_path) + beta_GARCH * sim_sigma();
%     end
%     sim_values(:, i_path) = 1;
% end

% Plot the simulated random variables
figure
plot(1:len_paths, sim_values);
title('simulated values');
xlabel('time'); ylabel('value')
xlim([0, len_paths])

%% 2) estimation student-t GARCH model

%% 3) VaR plots
