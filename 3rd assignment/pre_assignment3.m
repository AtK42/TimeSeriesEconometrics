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
len_path = 5e3;
% % number of paths
n_paths = 1e3;
% % number of lags for GARCH model, for this code fixed at (1,1)
%r=1; s=1;

rng(8, 'multFibonacci')
% set the model parameters
% % mean
% self chosen:
%mu = [.05, .075, .04]';
% from paper:
mu = [.164, -.153, -.865]';

% % weights
% self chosen:
%wghts = [15, 20, 4]';
%from paper:
wghts =  [.541, .433, .026]';

% % GARCH coefficients
% self chosen:
%gamma0 = [.01, .025, .03]';
% from paper:
gamma0 = [0, .012, .0332]';

% self chosen:
%gamma = [0.2, 0.85, 0.45]';
% from paper:
gamma = [.022, .197, 1.303]';

% self chosen:
%Psi = diag([.1,.3,.2]);
%Psi = toeplitz([.1 .3 .2]);
% from paper:
Psi = diag([.956, .835, .567])';

[sim_values, sim_sigma2] = MN_GARCH11_sim(len_path, mu, wghts, gamma0, gamma, Psi);

% plotting
subplot(2,1,1)
plot(sim_values, 'k-', 'LineWidth',.1);
title('simulated path')

subplot(2,1,2)
plot(sqrt(sim_sigma2(1,:)));
hold on
plot(sqrt(sim_sigma2(2,:)));
plot(sqrt(sim_sigma2(3,:)));
hold off
title('volatility process')

%% 2) estimation student-t GARCH model
est_model = estimate(garch('ARCHLags', 1, 'GARCHLags', 1, 'Distribution', 't'), sim_values);
%%
% initial parameter values
theta0 = [0.01, 0.1, 0.9, 8]; % [omega, alpha, beta, df]

% Define the log-likelihood function
logLikelihood = @(theta) -sum(log(tpdf(sim_values, theta(4)) .* sqrt(theta(2)/(theta(4)-2)) .* exp(-0.5*(sim_values.^2)./(theta(4)-2))));

% Optimize the log-likelihood function
options = optimset('MaxIter', 1000, 'MaxFunEvals', 1000);
theta_hat = fminsearch(logLikelihood, theta0, options);

% Extract the estimated parameters
omega_hat = theta_hat(1);
alpha_hat = theta_hat(2);
beta_hat = theta_hat(3);
df_hat = theta_hat(4);

% Display the estimated parameters
fprintf('Estimated GARCH(1,1) coefficients:\n');
fprintf('omega_hat = %.4f\n', omega_hat);
fprintf('alpha_hat = %.4f\n', alpha_hat);
fprintf('beta_hat = %.4f\n', beta_hat);
fprintf('df_hat = %.4f\n', df_hat);

%% 3) VaR plots
