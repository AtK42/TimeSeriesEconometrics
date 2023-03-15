%% 1: constant X matrix, null and al
% Via simulation, compare the theoretical distribution of the R^2 statistic, and the F statistic, 
%   to its empirically obtained one. Do so for a fixed X matrix of size 50 X 4, with a column of ones, 
%   and 3 regressors, each obtained as IID data from a normal distribution. The X matrix stays the same 
%   in all of your sim=1e5 simulations. What changes is, of course, the error term you add to the X*beta 
%   term. In order to assess the distribution of R2 and F under the null hypothesis, you need to take 
%   the beta vector to be (0, 0, 0, 0).
% What changes among your sim=1e5 runs is the residual vector --- in each of your simulations, the e vector 
%   is to be IID N(0, sigma^2), for which you please do 3 versions: sigma^2 = 1, 2, and 4. I suggest a nice 
%   graphical display, using, say, a kernel density of the 1e5 simulated values, overlaid with the true 
%   distributions of R2 and F under the null -- see page 9 of our book.
%
% Not done yet! Now via simulation only, show the distribution of R^2 and F when beta is (1,1,1,1). Here, 
%   the small sample distribution theory is more complicated, so we suffice ourselves with obtaining it 
%   via simulation. Make a kernel density using each of the 3 values of sigma^2 given above.
n_obs = 50;
n_regressors = 3;
n_sim = 1e5;
s2_vec = [1, 2, 4];
X_mat = [ones(n_obs, 1), normrnd(0, 1, n_obs, n_regressors)];
beta_vec_zero = zeros(1, 4)';
beta_vec_one = ones(1, 4)';

% set seed
rng('default')

% initialize variables
R2_zero1 = zeros(n_sim, 1); R2_one1 = zeros(n_sim, 1);
R2_zero2 = zeros(n_sim, 1); R2_one2 = zeros(n_sim, 1);
R2_zero3 = zeros(n_sim, 1); R2_one3 = zeros(n_sim, 1);
F_zero1 = zeros(n_sim, 1); F_one1 = zeros(n_sim, 1);
F_zero2 = zeros(n_sim, 1); F_one2 = zeros(n_sim, 1);
F_zero3 = zeros(n_sim, 1); F_one3 = zeros(n_sim, 1);

% simulate n_sim times the R2 and F statistic for different values of the
%   var of the e vector
for i = 1:n_sim
    % set seed
    rng(i);

    % % for beta vector (0,0,0,0)
    % get Y as prod X*beta (zero since beta is zero) plus some noise
    Y_mat_zero1 = X_mat * beta_vec_zero + normrnd(0, sqrt(s2_vec(1)), [n_obs, 1]);
    Y_mat_zero2 = X_mat * beta_vec_zero + normrnd(0, sqrt(s2_vec(2)), [n_obs, 1]);
    Y_mat_zero3 = X_mat * beta_vec_zero + normrnd(0, sqrt(s2_vec(3)), [n_obs, 1]);

    % calculate the beta hat for each of the 3 versions
    beta_hat_vec_zero1 = inv((X_mat' * X_mat)) * X_mat' * Y_mat_zero1;
    beta_hat_vec_zero2 = inv((X_mat' * X_mat)) * X_mat' * Y_mat_zero2;
    beta_hat_vec_zero3 = inv((X_mat' * X_mat)) * X_mat' * Y_mat_zero3;

    % calculate the Y hat for each of the 3 versions
    Y_hat_mat_zero1 = X_mat * beta_hat_vec_zero1;
    Y_hat_mat_zero2 = X_mat * beta_hat_vec_zero2;
    Y_hat_mat_zero3 = X_mat * beta_hat_vec_zero3;

    % calculate the realization of the (empirical) R2 r.v. 
    %   (R2 = ESS/TSS = 1 - RSS/TSS)
    R2_zero1(i) = 1 - sum((Y_mat_zero1 - Y_hat_mat_zero1).^2 / sum(Y_mat_zero1 - mean(Y_hat_mat_zero1)).^2);
    R2_zero2(i) = 1 - sum((Y_mat_zero2 - Y_hat_mat_zero2).^2 / sum(Y_mat_zero1 - mean(Y_hat_mat_zero1)).^2);
    R2_zero3(i) = 1 - sum((Y_mat_zero3 - Y_hat_mat_zero3).^2 / sum(Y_mat_zero1 - mean(Y_hat_mat_zero1)).^2);

    % calculate the (empirical) F statistic
    F_zero1(i) = (n_obs - n_regressors -1)/n_regressors * R2_zero1(i) /(1-R2_zero1(i));
    F_zero2(i) = (n_obs - n_regressors -1)/n_regressors * R2_zero2(i) /(1-R2_zero2(i));
    F_zero3(i) = (n_obs - n_regressors -1)/n_regressors * R2_zero3(i) /(1-R2_zero3(i));

    % % for beta vector (1,1,1,1)
    % get Y as prod X*beta (zero since beta is zero) plus some noise
    Y_mat_one1 = X_mat * beta_vec_one + normrnd(0, sqrt(s2_vec(1)), [n_obs, 1]);
    Y_mat_one2 = X_mat * beta_vec_one + normrnd(0, sqrt(s2_vec(2)), [n_obs, 1]);
    Y_mat_one3 = X_mat * beta_vec_one + normrnd(0, sqrt(s2_vec(3)), [n_obs, 1]);

    % calculate the beta hat for each of the 3 versions
    beta_hat_vec_one1 = inv((X_mat' * X_mat)) * X_mat' * Y_mat_one1;
    beta_hat_vec_one2 = inv((X_mat' * X_mat)) * X_mat' * Y_mat_one2;
    beta_hat_vec_one3 = inv((X_mat' * X_mat)) * X_mat' * Y_mat_one3;

    % calculate the Y hat for each of the 3 versions
    Y_hat_mat_one1 = X_mat * beta_hat_vec_one1;
    Y_hat_mat_one2 = X_mat * beta_hat_vec_one2;
    Y_hat_mat_one3 = X_mat * beta_hat_vec_one3;

    % calculate the realization of the (empirical) R2 r.v. 
    %   (R2 = ESS/TSS = 1 - RSS/TSS)
    R2_one1(i) = 1 - sum((Y_mat_one1 - Y_hat_mat_one1).^2);
    R2_one2(i) = 1 - sum((Y_mat_one2 - Y_hat_mat_one2).^2);
    R2_one3(i) = 1 - sum((Y_mat_one3 - Y_hat_mat_one3).^2);

    % calculate the (empirical) F statistic
    F_one1(i) = (n_obs - n_regressors -1)/n_regressors * R2_one1(i) /(1-R2_one1(i));
    F_one2(i) = (n_obs - n_regressors -1)/n_regressors * R2_one2(i) /(1-R2_one2(i));
    F_one3(i) = (n_obs - n_regressors -1)/n_regressors * R2_one3(i) /(1-R2_one3(i));
end

% [~,~,~,~,statistics] = regress(Y_mat, X_mat);
% R2(i) = statistics(1);
% F(i) = statistics(2);

%% plotting the simulated values against the theoretical ones for the first
%   part
x_vals_R2_plot = -3:.01:3;
x_vals_F_plot = -100:1:100;

theo_R2 = betapdf(x_vals_R2_plot, (n_regressors)/2, (n_obs-n_regressors-1)/2); % k includes beta0
theo_F = fpdf(x_vals_F_plot, n_regressors, n_obs-n_regressors-1);

emp_R2_zero1 = ksdensity(R2_zero1, x_vals_R2_plot); emp_R2_one1 = ksdensity(R2_one1, x_vals_R2_plot);
emp_R2_zero2 = ksdensity(R2_zero2, x_vals_R2_plot); emp_R2_one2 = ksdensity(R2_one2, x_vals_R2_plot);
emp_R2_zero3 = ksdensity(R2_zero3, x_vals_R2_plot); emp_R2_one3 = ksdensity(R2_one3, x_vals_R2_plot);

emp_F_zero1 = ksdensity(F_zero1, x_vals_F_plot); emp_F_one1 = ksdensity(F_one1, x_vals_F_plot);
emp_F_zero2 = ksdensity(F_zero2, x_vals_F_plot); emp_F_one2 = ksdensity(F_one2, x_vals_F_plot);
emp_F_zero3 = ksdensity(F_zero3, x_vals_F_plot); emp_F_one3 = ksdensity(F_one3, x_vals_F_plot);


figure
plot(x_vals_R2_plot, theo_R2, 'k--', ...
     x_vals_R2_plot, emp_R2_zero1, 'r-', ...
     x_vals_R2_plot, emp_R2_zero2, 'g-', ...
     x_vals_R2_plot, emp_R2_zero3, 'b-')
legend('theoretical R2', '\sigma^2 = 1', '\sigma^2 = 2', '\sigma^2 = 4', 'Location', 'northeast')
title('R^2, \beta = (0,0,0,0)')

figure
plot(x_vals_F_plot, theo_F, 'k--',...
     x_vals_F_plot, emp_F_zero1, 'r-', ...
     x_vals_F_plot, emp_F_zero2, 'g-', ...
     x_vals_F_plot, emp_F_zero3, 'b-')
legend('theoretical F', '\sigma^2 = 1', '\sigma^2 = 2', '\sigma^2 = 4', 'Location', 'northeast')
title('F, \beta = (0,0,0,0)')

figure
plot(x_vals_R2_plot, emp_R2_one1, 'r-', ...
     x_vals_R2_plot, emp_R2_one2, 'g-', ...
     x_vals_R2_plot, emp_R2_one3, 'b-')
legend('\sigma^2 = 1', '\sigma^2 = 2', '\sigma^2 = 4', 'Location', 'northeast')
title('R^2, \beta = (1,1,1,1)')

figure
plot(x_vals_F_plot, emp_F_one1, 'r-', ...
     x_vals_F_plot, emp_F_one2, 'g-', ...
     x_vals_F_plot, emp_F_one3, 'b--')
legend('\sigma^2 = 1', '\sigma^2 = 2', '\sigma^2 = 4', 'Location', 'northeast')
title('F, \beta = (1,1,1,1)')

%% 2:
% Fix the sample size at T=100. For the number of regressors k=2,3,...,10, simulate sim=1e4 replications 
%   of the R^2  statistic, based on a beta vector of all zeros, and X matrix obtained similar to above, 
%   i.e., IID normal. Now you use a different X matrix for each of the 1e4 replications, and of course 
%   also a different error vector, namely the same as above, IID N(0,sigma2), with sigma2 = 2. Now, for 
%   each k, you generate 1e4  R^2 values. You make a graphic containing 9 boxplots of the R2 values (obviously, 
%   one for each k). The goal is to see the behavior of  R^2 as the number of bogus (non-significant) 
%   regressors increases. Discuss what you find, and why it is occurring.
n_obs = 100; % same as T
k = 2:10;
n_regressors = 3;
n_sim = 1e4;

for i = 1:n_sim
    X_mat = [ones(n_obs, 1), randnorm(0, 1, n_obs, n_regressors)];
end
%% 3:
% Go and look at the regression model with AR(1) disturbances (and Gaussian innovations), equation lines 
%   5.1 and 5.2 in our book. For each of two sample sizes, namely T=20 and T=50, you are to simulate the 
%   process (say, 1e4 times), and estimate the parameters as described on page 225, end of section 5.1. 
%   Namely, OLS for the regression, followed by OLS for the AR(1) parameter a, followed by GLS for beta, 
%   followed by... until convergence, which I suggest you do until you achieve 5 points after the decimal. 
%   The true X matrix is a column of ones and 4 "random" regressors (you generate, as we did above in parts 
%   1 and 2, and not done new for each of the sim=1e4 times. What of course is different each time is the 
%   vector of U values -- you need to simulate an AR(1) time series model to get the error terms for the 
%   regression!
%
% So, once you can simulate this, do the simulation of the sim values, and report some kind of graphic or 
%   tables of the accuracy of beta and ar(1) parameter a. Now, you do this also for a grid of a-values, 
%   namely 0.0, 0.2, 0.4, 0.6, 0.8, and 0.99. Notice you do NOT need to compare to the exact MLE -- you do 
%   not do the exact MLE, but rather only the super-fast iterative method based on ols and gls.