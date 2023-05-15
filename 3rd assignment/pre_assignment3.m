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
mu = [.05, .075, .04]';
% % weights
wghts = [15, 20, 4]';
% % GARCH coefficients
gamma0 = [.01, .025, .03]';
gamma = [0.2, 0.85, 0.45]';
%Psi = diag([.1,.3,.2]);
Psi = toeplitz([.1 .3 .2]);

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
%% 3) VaR plots
