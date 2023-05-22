% Here goes: Have a look at my mixed-normal GARCH construction, on page 478

% 1a):
% Your first task is to  * simulate *  this process. Use k=2 components, though
%   I might later ask you to try 3. You do NOT have to estimate the model, 
%   though it is also not difficult to set up the likelihood and estimate 
%   it. You just need to simulate it, and you can play around with the 
%   required parameters so that when you plot the resulting univariate time 
%   series, it should "look like actual daily stock returns".  Since you 
%   are simulating, you can simulate T=5000 observations -- or even more, 
%   all for free.

% 1b):  
% The next part is to estimate a student-t-GARCH model based on your
%   simulated data. Notice that the true DGP does NOT match the model you 
%   are estimating. That is okay -- that is reality. The student-t garch 
%   estimation is very easy to set up -- you just build on the codes I have
%   in my chapter. Obviously, the degrees of freedom of the student t gets 
%   jointly estimated with the garch parameters.  Surely, you can find 
%   toolboxes for matlab and python that also do it, ...)

% 1c:
% The one step ahead prediction of the t-GARCH model is a univariate 
%   Student t distribution with scale parameter coming from the GARCH 
%   recursion. Set this up (check for codes in my book), and do the 
%   following: Starting from t=250, and going up to t=999, you estimate the
%   t-GARCH model using data points 1 to t, and compute the predictive 
%   density for time t+1 based on data from time 1 to time t, and compute 
%   its 1% quantile (this latter quantile calculation of the Student t is 
%   built into Matlab, check tinv I think it is called). This is the Value 
%   at Risk (actually, when multiplied by the AUM, assets under management 
%   amount of money. We ignore this.) Compare it to the ACTUAL return 
%   (from your simulated data of course) at time t+1, and record a 1 if the
%   actual return is lower than the VaR, and zero otherwise. If the model 
%   is accurate, then approximately 1% of your 1000-250 recorded values 
%   will be a 1, and the rest zero. You report the percentage of 
%   "violations" (the number of ones). See Section 11.1 of the book.

%% 1a) simulate a mixed-normal GARCH(1,1) construction
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
grid()

subplot(2,1,2)
plot(sqrt(sim_sigma2(1,:)));
hold on
plot(sqrt(sim_sigma2(2,:)));
plot(sqrt(sim_sigma2(3,:)));
hold off
grid()
title('volatility process')

%% 1b) estimation student-t GARCH model
est_model = estimate(garch('ARCHLags', 1, 'GARCHLags', 1, 'Distribution', 't'), sim_values);
%
% THE FOLLOWING IS FROM CHAT GPT:
% % Set initial parameter values
% theta0 = [0.01, 0.1, 0.9, 10]; % [omega, alpha, beta, df]
% 
% % Define the log-likelihood function
% logLikelihood = @(theta) -sum(log(tpdf(sim_values, theta(4)) .* sqrt(theta(2)/(theta(4)-2)) .* exp(-0.5*(sim_values.^2)./(theta(4)-2))));
% 
% % Optimize the log-likelihood function
% options = optimset('MaxIter', 1000, 'MaxFunEvals', 1000);
% theta_hat = fminsearch(logLikelihood, theta0, options);
% 
% % Extract the estimated parameters
% omega_hat = theta_hat(1);
% alpha_hat = theta_hat(2);
% beta_hat = theta_hat(3);
% df_hat = theta_hat(4);
% 
% % Display the estimated parameters
% fprintf('Estimated GARCH(1,1) coefficients:\n');
% fprintf('omega_hat = %.4f\n', omega_hat);
% fprintf('alpha_hat = %.4f\n', alpha_hat);
% fprintf('beta_hat = %.4f\n', beta_hat);
% fprintf('df_hat = %.4f\n', df_hat);
%
% END CHAT GPT CODE
%% 1c) 99%-VaR calculations
alpha = .01;
t_start = 249;
t_end = 999;
delta = 1;
% from the book p.460, discussion following the formula for the APARCH 
%   model: "It turns out that, for daily financial asset returns data,
%   the likelihood is often relatively flat in parameter delta, with its 
%   maximum between one and two. We advocate just setting it equal to 
%   one."

% initialize variables
est_constant = zeros(t_end-t_start, 1);
est_GARCH = zeros(t_end-t_start, 1);
est_ARCH = zeros(t_end-t_start, 1);
est_dof = zeros(t_end-t_start, 1);
E = zeros(t_end-t_start, 1);
fore_sigma2 = zeros(t_end-t_start, 1);
violations = zeros(t_end-t_start, 1);

for t = t_start:t_end
    % estimate model and extract parameters
    est_model = estimate(garch('ARCHLags', 1, 'GARCHLags', 1, 'Distribution', 't'), sim_values(1:t), 'Display', 'off');
    est_constant(t) = est_model.Constant;
    est_GARCH(t) = cell2mat(est_model.GARCH);
    est_ARCH(t) = cell2mat(est_model.ARCH);
    est_dof(t) = est_model.Distribution.DoF;

    % set initial value of the sigma2 as the unconditional sigma2
    if t == t_start
        fore_sigma2(t_start) = 1/(1 - est_GARCH(t) - est_ARCH(t));
    end

    % recursive evaluation of the forecasted volatility term (book p. 467)
    gamma(t) = .5; % must be estimated somehow but not sure how
    E(t) = (abs(sim_values(t)) - gamma(t) * sim_values(t))^delta;
    fore_sigma2(t+1) = est_constant(t) + est_ARCH(t) * E(t) + est_GARCH(t) * fore_sigma2(t);

    % calculate the 1% quantile of the predictive density for t+1
    q = fore_sigma2(t+1) * tinv(alpha, est_dof(t)); % not sure whether this is correct
    violations(t+1) = (sim_values(t+1) < q);
end

mean(violations) % should be very close to 1%


















