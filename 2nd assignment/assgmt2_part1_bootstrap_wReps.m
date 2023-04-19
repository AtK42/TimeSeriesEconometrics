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
b_true = -.6; % freely assumed
T=1000;
n_obs = T;
B = 400; % check if can be smaller or should be larger for reliable results
n_reps = 555;
alpha = .1;

% initialize variables
MA1est_Durbin = zeros(n_reps, 1);
MA1est_approxMLE = zeros(n_reps, 1);
MA1est_matlab = zeros(n_reps, 1);
vec_timeseries_Durbin = zeros(n_obs, B);
vec_timeseries_approxMLE = zeros(n_obs, B);
vec_timeseries_matlab = zeros(n_obs, B);
boot_MA1est_Durbin = zeros(B, n_reps);
boot_MA1est_approxMLE = zeros(B, n_reps);
boot_MA1est_matlab = zeros(B, n_reps);

% simulate MA(1) process and then estimate the parameter
% % for reproducibility
%rng(112); seed = rng;

parfor r = 1:n_reps
    % % simulate MA(1) process with true parameter
    timeseries = armasim(n_obs, 1, 0, b_true);
    
    % % estimate MA parameter using the three methods
    % %% Durbin
    MA1est_Durbin(r) = DurbinMA1959(timeseries, 1);
    % %% approx MLE
    MA1_temp = ma1(timeseries, 0); % 0 for approx MLE
    MA1est_approxMLE(r) = MA1_temp(1);
    % %% matlab
    MA1_temp = armax(timeseries, [0, 1]);
    MA1est_matlab(r) = MA1_temp.c(2:end);
    
    for i = 1:B
        %rng(seed.Seed + i);
        % delete variables
        %vec_timeseries_Durbin = []; vec_timeseries_approxMLE = []; vec_timeseries_matlab = [];
        % simulation of MA process with ESTIMATED b
        vec_timeseries_Durbin = armasim(n_obs, 1, 0, MA1est_Durbin(r));
        vec_timeseries_approxMLE = armasim(n_obs, 1, 0, MA1est_approxMLE(r));
        vec_timeseries_matlab = armasim(n_obs, 1, 0, MA1est_matlab(r));
    
        % estimation of bootstrapped timeseries
        boot_MA1est_Durbin(i,r) = DurbinMA1959(vec_timeseries_Durbin, 1);
        MA1_temp = ma1(vec_timeseries_approxMLE, 0);
        boot_MA1est_approxMLE(i,r) = MA1_temp(1);
        MA1_temp = armax(vec_timeseries_matlab, [0, 1]);
        boot_MA1est_matlab(i,r) = MA1_temp.c(2:end);
    end
    if mod(r, 2) == 0
        clc
        disp([num2str(r), ' out of ', num2str(n_reps), ' replications done (', num2str(round((r/n_reps)*100, 2)), '%).']);
    end
end
clc
disp('All done, enjoy!')

%% calculate CI
% get CI
ci_Durbin = quantile(boot_MA1est_Durbin, [alpha/2 1-alpha/2]);
ci_approxMLE = quantile(boot_MA1est_approxMLE, [alpha/2 1-alpha/2]);
ci_matlab = quantile(boot_MA1est_matlab, [alpha/2 1-alpha/2]);
% check which estimates are within the CI
cov_Durbin_vec = (b_true >= ci_Durbin(1,:) & (b_true <= ci_Durbin(2,:)));
cov_approxMLE_vec = (b_true >= ci_approxMLE(1,:) & (b_true <= ci_approxMLE(2,:)));
cov_matlab_vec = (b_true >= ci_matlab(1,:) & (b_true <= ci_matlab(2,:)));
% get coverage
cov_Durbin = mean(cov_Durbin_vec); disp(cov_Durbin);
cov_approxMLE = mean(cov_approxMLE_vec); disp(cov_approxMLE);
cov_matlab = mean(cov_matlab_vec); disp(cov_matlab);

%%
MA1est = [MA1est_Durbin, MA1est_approxMLE, MA1est_matlab];
boot_MA1est = zeros(B, n_reps, 3);
boot_MA1est(:,:,1) = boot_MA1est_Durbin;
boot_MA1est(:,:,2) = boot_MA1est_approxMLE;
boot_MA1est(:,:,3) = boot_MA1est_matlab;
%%
writematrix(MA1est, "MA1est.txt")
writematrix(boot_MA1est_Durbin, "boot_MA1est_Durbin.txt")
writematrix(boot_MA1est_approxMLE, "boot_MA1est_approxMLE.txt")
writematrix(boot_MA1est_matlab, "boot_MA1est_matlab.txt")
%%
struc_MA1est = struc(MA1est);
struc_boot_MA1est = struc(boot_MA1est);