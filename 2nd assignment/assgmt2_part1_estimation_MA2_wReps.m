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
n_reps = 50;

vec_timeseries_100 = zeros(T_vec(1), n_reps);
vec_timeseries_1000 = zeros(T_vec(2), n_reps);
MA2est_Durbin_100 = zeros(n_reps, 2);
MA2est_Durbin_1000 = zeros(n_reps, 2);
MA2est_approxMLE_100 = zeros(n_reps, 2);
MA2est_approxMLE_1000 = zeros(n_reps, 2);
MA2est_matlab_100 = zeros(n_reps, 2);
MA2est_matlab_1000 = zeros(n_reps, 2);

% simulate MA(2) process and then estimate the parameters
%seed = rng(42);
T100 = T_vec(1);
T1000 = T_vec(2);
parfor r = 1:n_reps
    % % simulation of MA process
    vec_timeseries_100(:, r) = armasim(T100, 1, 0, b_true);
    vec_timeseries_1000(:, r) = armasim(T1000, 1, 0, b_true);

    % % estimation of MA parameter
    % %% Durbin
    MA2est_Durbin_100(r, :) = DurbinMA1959(vec_timeseries_100(:, r), 2)';
    MA2est_Durbin_1000(r, :) = DurbinMA1959(vec_timeseries_1000(:, r), 2)';
    
    % %% approximate MLE
    %[MA2_temp, ~, ~, ~, ~] = armareg(vec_timeseries_100, [], 0, 2, 0); % last 0 for approx MLE
    %MA2est_approxMLE_100 = MA2_temp;
    
    %[MA2_temp, ~, ~, ~, ~] = armareg(vec_timeseries_1000, [], 0, 2, 0); % last 0 for approx MLE
    %MA2est_approxMLE_1000 = MA2_temp;
    
    MA2temp = maq_approxMLE(vec_timeseries_100(:, r), 2);
    MA2est_approxMLE_100(r, :) = MA2temp(1:2)';
    MA2temp = maq_approxMLE(vec_timeseries_1000(:, r), 2);
    MA2est_approxMLE_1000(r, :) = MA2temp(1:2)';
    
    % %% bulilt-in MATLAB function
    MA2_temp = armax(vec_timeseries_100(:, r), [0, 2]);
    MA2est_matlab_100(r,:) = MA2_temp.c(2:end);
    MA2_temp = armax(vec_timeseries_1000(:, r), [0, 2]);
    MA2est_matlab_1000(r,:) = MA2_temp.c(2:end);

    if mod(r, 5) == 0
        clc
        disp([num2str(r), ' out of ', num2str(n_reps), ' replications done (', num2str(round((r/n_reps)*100, 2)), '%).']);
    end
end

%% plotting
xlim_ = [(b_true(1) - .1), (b_true(2) + .1)];
ylim_ = [min([b_true, mean(MA2est_Durbin_100), mean(MA2est_approxMLE_100), mean(MA2est_matlab_100)])-.1, ...
        max([b_true, mean(MA2est_Durbin_100), mean(MA2est_approxMLE_100), mean(MA2est_matlab_100)])+.1];

horz = zeros(length(b_true), 1000);
horz(1,:) = linspace(b_true(1)-.02, b_true(1)+.02, 1000);
horz(2,:) = linspace(b_true(2)-.02, b_true(2)+.02, 1000);

% T = 100
figure('DefaultAxesFontSize', 14)
scatter(b_true, [mean(MA2est_Durbin_100(:,2)) mean(MA2est_Durbin_100(:,1))])
hold on
scatter(b_true, [mean(MA2est_approxMLE_100(:,1)) mean(MA2est_approxMLE_100(:,2))], 'filled')
scatter(b_true, [mean(MA2est_matlab_100(:,1)) mean(MA2est_matlab_100(:,2))], 'green', 'd')
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
scatter(b_true, [mean(MA2est_Durbin_1000(:,2)) mean(MA2est_Durbin_1000(:,1))])
hold on
scatter(b_true, [mean(MA2est_approxMLE_1000(:,1)) mean(MA2est_approxMLE_1000(:,2))], 'filled')
scatter(b_true, [mean(MA2est_matlab_1000(:,1)) mean(MA2est_matlab_1000(:,2))], 'green', 'd')
plot(horz(1,:),ones(length(horz(1,:)),1) * b_true(1), 'black')
plot(horz(2,:),ones(length(horz(2,:)),1) * b_true(2), 'black')

xlim(xlim_); ylim(ylim_)
xlabel('true values'); ylabel('estimated values')
legend('Durbin', 'approx MLE', 'Matlab', 'true value', '', 'Location', 'south');
%legend('boxoff')
title(["True vs. estimated values of the MA(2)","parameter for T = 1000"])

%%
stat = {'Durbin'; 'approx MLE'; 'Matlab'};
% bias
biasb1_100_Durbin = mean(MA2est_Durbin_100(:,1)) - b_true(1);
biasb1_100_approxMLE = mean(MA2est_approxMLE_100(:,1)) - b_true(1);
biasb1_100_matlab = mean(MA2est_matlab_100(:,1)) - b_true(1);
biasb2_100_Durbin = mean(MA2est_Durbin_100(:,2)) - b_true(2);
biasb2_100_approxMLE = mean(MA2est_approxMLE_100(:,2)) - b_true(2);
biasb2_100_matlab = mean(MA2est_matlab_100(:,2)) - b_true(2);
bias_b1_100 = round([biasb1_100_Durbin; biasb1_100_approxMLE; biasb1_100_matlab], 4);
bias_b2_100 = round([biasb2_100_Durbin; biasb2_100_approxMLE; biasb2_100_matlab], 4);

biasb1_1000_Durbin = mean(MA2est_Durbin_1000(:,1)) - b_true(1);
biasb1_1000_approxMLE = mean(MA2est_approxMLE_1000(:,1)) - b_true(1);
biasb1_1000_matlab = mean(MA2est_matlab_1000(:,1)) - b_true(1);
biasb2_1000_Durbin = mean(MA2est_Durbin_1000(:,2)) - b_true(2);
biasb2_1000_approxMLE = mean(MA2est_approxMLE_1000(:,2)) - b_true(2);
biasb2_1000_matlab = mean(MA2est_matlab_1000(:,2)) - b_true(2);
bias_b1_1000 = round([biasb1_1000_Durbin; biasb1_1000_approxMLE; biasb1_1000_matlab], 4);
bias_b2_1000 = round([biasb2_1000_Durbin; biasb2_1000_approxMLE; biasb2_1000_matlab], 4);


% MSE
b1_true_100 = b_true(1) * ones(n_reps, 1); b2_true_100 = b_true(2) * ones(n_reps, 1);
b1_true_1000 = b_true(1) * ones(n_reps, 1); b2_true_1000 = b_true(2) * ones(n_reps, 1);

MSEb1_100_Durbin = mean( (MA2est_Durbin_100(:,1) - b1_true_100).^2, 1);
MSEb1_100_approxMLE = mean( (MA2est_approxMLE_100(:,1) - b1_true_100).^2, 1);
MSEb1_100_matlab = mean( (MA2est_matlab_100(:,1) - b1_true_100).^2, 1);
MSEb2_100_Durbin = mean( (MA2est_Durbin_100(:,2) - b2_true_100).^2, 1);
MSEb2_100_approxMLE = mean( (MA2est_approxMLE_100(:,2) - b2_true_100).^2, 1);
MSEb2_100_matlab = mean( (MA2est_matlab_100(:,2) - b2_true_100).^2, 1);
MSE_b1_100 = round([MSEb1_100_Durbin; MSEb1_100_approxMLE; MSEb1_100_matlab], 4);
MSE_b2_100 = round([MSEb2_100_Durbin; MSEb2_100_approxMLE; MSEb2_100_matlab], 4);

MSEb1_1000_Durbin = mean( (MA2est_Durbin_1000(:,1) - b1_true_1000).^2, 1);
MSEb1_1000_approxMLE = mean( (MA2est_approxMLE_1000(:,1) - b1_true_1000).^2, 1);
MSEb1_1000_matlab = mean( (MA2est_matlab_1000(:,1) - b1_true_1000).^2, 1);
MSEb2_1000_Durbin = mean( (MA2est_Durbin_1000(:,2) - b2_true_1000).^2, 1);
MSEb2_1000_approxMLE = mean( (MA2est_approxMLE_1000(:,2) - b2_true_1000).^2, 1);
MSEb2_1000_matlab = mean( (MA2est_matlab_1000(:,2) - b2_true_1000).^2, 1);
MSE_b1_1000 = round([MSEb1_1000_Durbin; MSEb1_1000_approxMLE; MSEb1_1000_matlab], 4);
MSE_b2_1000 = round([MSEb2_1000_Durbin; MSEb2_1000_approxMLE; MSEb2_1000_matlab], 4);

tab = table(stat, ...
            bias_b1_100, bias_b2_100, bias_b1_1000, bias_b2_1000, ...
            MSE_b1_100, MSE_b2_100, MSE_b1_1000, MSE_b2_1000);
disp(tab);
