% Now, continuing: That was all for MA(1). We now need to do the MA(q) 
%   model. We can use q=2 and just one parameter constellation, namely the 
%   choice from my book, page 300: b1 = −0.5 and b2 = −0.24. From the text 
%   in my book, it seems like an invertible model. Do a simulation for 
%   T=100 AND ALSO T=1000, comparing Durbin's method to my approximate 
%   method, and to the built in method of Matlab or python. Make a smart 
%   presentation, which could be a simple table with bias, variance, and 
%   MSE values. Or a graphic, up to you.

clear
clc
close all

q = 2;
b_true = [-0.5 -0.24];
T_vec = [100 1000];

% simulate MA(2) process and then estimate the parameters
% % simulation of MA process
rng(42);
vec_timeseries_100 = armasim(T_vec(1), 1, 0, b_true);
vec_timeseries_1000 = armasim(T_vec(2), 1, 0, b_true);

% % estimation of MA parameter
% %% Durbin
MA2est_Durbin_100 = DurbinMA1959(vec_timeseries_100, 2);
MA2est_Durbin_1000 = DurbinMA1959(vec_timeseries_1000, 2);

% %% approximate MLE
%[MA2_temp, ~, ~, ~, ~] = armareg(vec_timeseries_100, [], 0, 2, 0); % last 0 for approx MLE
%MA2est_approxMLE_100 = MA2_temp;

%[MA2_temp, ~, ~, ~, ~] = armareg(vec_timeseries_1000, [], 0, 2, 0); % last 0 for approx MLE
%MA2est_approxMLE_1000 = MA2_temp;

MA2temp = maq(vec_timeseries_100, 1, 2);
MA2est_approxMLE_100 = MA2temp(1:2);
MA2temp = maq(vec_timeseries_1000, 1, 2);
MA2est_approxMLE_1000 = MA2temp(1:2);

% %% bulilt-in MATLAB function
MA2_temp = armax(vec_timeseries_100, [0, 2]);
MA2est_matlab_100 = MA2_temp.c(2:end);

MA2_temp = armax(vec_timeseries_1000, [0, 2]);
MA2est_matlab_1000 = MA2_temp.c(2:end);

%% plotting
xlim_ = [(b_true(1) - .1), (b_true(2) + .1)];
ylim_ = [min([b_true, MA2est_Durbin_100', MA2est_approxMLE_100', MA2est_matlab_100])-.1, ...
        max([b_true, MA2est_Durbin_100', MA2est_approxMLE_100', MA2est_matlab_100])+.1];

horz = zeros(length(b_true), 1000);
horz(1,:) = linspace(b_true(1)-.02, b_true(1)+.02, 1000);
horz(2,:) = linspace(b_true(2)-.02, b_true(2)+.02, 1000);

% T = 100
figure('DefaultAxesFontSize', 14)
scatter(b_true, [mean(MA2est_Durbin_100(1),2)' mean(MA2est_Durbin_100(2),2)])
hold on
scatter(b_true, [mean(MA2est_approxMLE_100(1),2)' mean(MA2est_approxMLE_100(2),2)'], 'filled')
scatter(b_true, [mean(MA2est_matlab_100(1),2)' mean(MA2est_matlab_100(2),2)'], 'green', 'd')
plot(horz(1,:),ones(length(horz(1,:)),1) * b_true(1), 'black')
plot(horz(2,:),ones(length(horz(2,:)),1) * b_true(2), 'black')

xlim(xlim_); ylim(ylim_)
xlabel('true values'); ylabel('estimated values')
legend('Durbin', 'approx MLE', 'Matlab', 'true value', '', 'Location', 'south');
%legend('boxoff')
title(["True vs. estimated values of the MA(2)","parameter for T = 100"])

%%
% T = 1000
figure('DefaultAxesFontSize', 14)
scatter(b_true, [mean(MA2est_Durbin_1000(1),2)' mean(MA2est_Durbin_1000(2),2)])
hold on
scatter(b_true, [mean(MA2est_approxMLE_1000(1),2)' mean(MA2est_approxMLE_1000(2),2)'], 'filled')
scatter(b_true, [mean(MA2est_matlab_1000(1),2)' mean(MA2est_matlab_1000(2),2)'], 'green', 'd')
plot(horz(1,:),ones(length(horz(1,:)),1) * b_true(1), 'black')
plot(horz(2,:),ones(length(horz(2,:)),1) * b_true(2), 'black')

xlim(xlim_); ylim(ylim_)
xlabel('true values'); ylabel('estimated values')
legend('Durbin', 'approx MLE', 'Matlab', 'true value', '', 'Location', 'south');
%legend('boxoff')
title(["True vs. estimated values of the MA(2)","parameter for T = 1000"])

%%
% bias
biasb1_100_Durbin = mean(MA2est_Durbin_100(1),2)' - b_true(1);
biasb1_100_approxMLE = mean(MA2est_approxMLE_100(1),2)' - b_true(1);
biasb1_100_matlab = mean(MA2est_matlab_100(1),2)' - b_true(1);
biasb2_100_Durbin = mean(MA2est_Durbin_100(2),2)' - b_true(2);
biasb2_100_approxMLE = mean(MA2est_approxMLE_100(2),2)' - b_true(2);
biasb2_100_matlab = mean(MA2est_matlab_100(2),2)' - b_true(2);

biasb1_1000_Durbin = mean(MA2est_Durbin_1000(1),2)' - b_true(1);
biasb1_1000_approxMLE = mean(MA2est_approxMLE_1000(1),2)' - b_true(1);
biasb1_1000_matlab = mean(MA2est_matlab_1000(1),2)' - b_true(1);
biasb2_1000_Durbin = mean(MA2est_Durbin_1000(2),2)' - b_true(2);
biasb2_1000_approxMLE = mean(MA2est_approxMLE_1000(2),2)' - b_true(2);
biasb2_1000_matlab = mean(MA2est_matlab_1000(2),2)' - b_true(2);

% MSE
b1_true_100 = b_true(1) * ones(100, 1); b2_true_100 = b_true(2) * ones(100, 1);
b1_true_1000 = b_true(1) * ones(1000, 1); b2_true_1000 = b_true(2) * ones(1000, 1);

MSEb1_100_Durbin = mean( (MA2est_100_Durbin(1) - b1_true_100).^2, 1);
MSEb1_100_approxMLE = mean( (MA2est_100_approxMLE(1) - b1_true_100).^2, 1);
MSEb1_100_matlab = mean( (MA2est_100_matlab(1) - b1_true_100).^2, 1);
MSEb2_100_Durbin = mean( (MA2est_100_Durbin(2) - b2_true_100).^2, 1);
MSEb2_100_approxMLE = mean( (MA2est_100_approxMLE(2) - b2_true_100).^2, 1);
MSEb2_100_matlab = mean( (MA2est_100_matlab(2) - b2_true_100).^2, 1);

MSEb1_1000_Durbin = mean( (MA2est_1000_Durbin(1) - b1_true_1000).^2, 1);
MSEb1_1000_approxMLE = mean( (MA2est_1000_approxMLE(1) - b1_true_1000).^2, 1);
MSEb1_1000_matlab = mean( (MA2est_1000_matlab(1) - b1_true_1000).^2, 1);
MSEb2_1000_Durbin = mean( (MA2est_1000_Durbin(2) - b2_true_1000).^2, 1);
MSEb2_1000_approxMLE = mean( (MA2est_1000_approxMLE(2) - b2_true_1000).^2, 1);
MSEb2_1000_matlab = mean( (MA2est_1000_matlab(2) - b2_true_1000).^2, 1);