% The last part of this segment of the homework is to apply the parametric 
%   bootstrap to ALL 3 methods in the MA(1) case. Thus, you need to have 
%   these 3 methods programmed as perfect black-boxes. You use the 
%   parametric bootstrap to get 90% confidence intervals for the unknown 
%   MA parameter, and you need to *simulate the bootstrap method* in order 
%   to report the actual coverage obtained. I offhand presume B=400 
%   replications should be enough, but be smart and active --- recall the 
%   choice of B is ideally infinity, but is limited by computation time. 
%   Try to ascertain what value of B is adequate --- that means, small 
%   enough such that the results are stable.
 
% This is the usual bootstrap thing we did last semester (chapters 1 and 2 
%   of my statistics book) You use the SINGLE bootstrap, as opposed to 
%   double. And notice I said parametric, and not nonparametric. You are 
%   welcome of course to also do nonparametric, but it is harder, because 
%   you have to simulate MA processes with innovations taken from the 
%   filtered ones from your fitted model. I spare you this. So, parametric 
%   is easier, so just use that.
 
% You report the quality (meaning actual coverage) of nominal (single 
%   parametric bootstrap) 90% CIs for the MA(1) parameter, based on the 3 
%   methods: Durbin, my approximate method of estimation, and built into 
%   Matlab.

clear
clc
close all

% prep work
b_true = -.5; % freely assumed
T=100;
n_obs = T;% + p;
B = 400; % check if can be smaller or should be larger for reliable results
alpha = .1;

% initialize variables
vec_timeseries_Durbin = zeros(n_obs, B);
vec_timeseries_approxMLE = zeros(n_obs, B);
vec_timeseries_matlab = zeros(n_obs, B);
boot_MA1est_Durbin = zeros(B, 1);
boot_MA1est_approxMLE = zeros(B, 1);
boot_MA1est_matlab = zeros(B, 1);

% simulate MA(1) process and then estimate the parameter
% % for reproducibility
rng(42); seed = rng;

% % simulate MA(1) process with true parameter
timeseries = armasim(n_obs, 1, 0, b_true);

% % estimate MA parameter using the three methods
% %% Durbin
MA1est_Durbin = DurbinMA1959(timeseries, 1);
% %% approx MLE
MA1_temp = ma1(timeseries, 0); % 0 for approx MLE
MA1est_approxMLE = MA1_temp(1);
% %% matlab
MA1_temp = armax(timeseries, [0, 1]);
MA1est_matlab = MA1_temp.c(2:end);

for i = 1:B
    rng(seed.Seed + i);
    % simulation of MA process with ESTIMATED b
    vec_timeseries_Durbin(:,i) = armasim(n_obs, 1, 0, MA1est_Durbin);
    vec_timeseries_approxMLE(:,i) = armasim(n_obs, 1, 0, MA1est_approxMLE);
    vec_timeseries_matlab(:,i) = armasim(n_obs, 1, 0, MA1est_matlab);

    % estimation of bootstrapped timeseries
    boot_MA1est_Durbin(i) = DurbinMA1959(vec_timeseries_Durbin(:,i), 1);
    MA1_temp = ma1(vec_timeseries_approxMLE(:,i), 0);
    boot_MA1est_approxMLE(i) = MA1_temp(1);
    MA1_temp = armax(vec_timeseries_matlab(:,i), [0, 1]);
    boot_MA1est_matlab(i) = MA1_temp.c(2:end);
end

% calculate CI
ci_Durbin = quantile(boot_MA1est_Durbin, [alpha/2 1-alpha/2]);
ci_approxMLE = quantile(boot_MA1est_approxMLE, [alpha/2 1-alpha/2]);
ci_matlab = quantile(boot_MA1est_matlab, [alpha/2 1-alpha/2]);
