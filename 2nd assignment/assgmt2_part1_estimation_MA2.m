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

% TODO: compare bias and MSE for parameter b

q = 2;
b_true = [-.5, -.24];
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
[MA2_temp, ~, ~, ~, ~] = armareg(vec_timeseries_100, [], 0, 2, 0); % last 0 for approx MLE
MA2est_approxMLE_100 = MA2_temp;

[MA2_temp, ~, ~, ~, ~] = armareg(vec_timeseries_1000, [], 0, 2, 0); % last 0 for approx MLE
MA2est_approxMLE_1000 = MA2_temp;

% %% bulilt-in MATLAB function
MA2_temp = armax(vec_timeseries_100, [0, 2]);
MA2est_matlab_100 = MA2_temp.c(2:end);

MA2_temp = armax(vec_timeseries_1000, [0, 2]);
MA2est_matlab_1000 = MA2_temp.c(2:end);



%% plotting

figure('DefaultAxesFontSize', 14)
scatter(-.9:.1:.9, mean(MA1est_Durbin,2)')
xticks(-.9:.1:.9)
xtickangle(90)
yticks(-.9:.1:.9)
hold on
scatter(-.9:.1:.9, mean(MA1est_approxMLE,2)', 'filled')
scatter(-.9:.1:.9, mean(MA1est_matlab,2)', 'green', 'd')


horz = zeros(length(b_true_vec), 1000);
for i = 1:length(b_true_vec)
    horz(i,:) = linspace(b_true_vec(i)-.05, b_true_vec(i)+.05, 1000);
end
plot(horz(1,:),ones(length(horz(1,:)),1) * b_true_vec(1), 'black')
ylim([-1 1]);xlim([-1 1]);xlabel('true values');ylabel('estimated values')
title(["True vs. estimated values of the MA(1)","parameter for different estimation methods"])
xticks(-.9:.1:.9)
hold on
for i = 2:length(b_true_vec)
    plot(horz(i,:), ones(length(horz(i,:)), 1) * b_true_vec(i), 'black')
end
%legend('Durbin', 'approx MLE', 'Matlab', 'true value', ...
%    'Location','northwest')
legend('Durbin', 'approx MLE', 'Matlab', '','','','','','','','','','','','','','','','','','','', ...
    'Location','northwest')


% plot bias
figure('DefaultAxesFontSize', 18)
plot(b_true_vec, mean(MA1est_Durbin,2)' - b_true_vec, 'r-', ...
     b_true_vec, mean(MA1est_approxMLE,2)' - b_true_vec, 'g:', ...
     b_true_vec, mean(MA1est_matlab,2)' - b_true_vec, 'b--', ...
     'linewidth',2)
legend('Durbin','approx MLE', 'Matlab', 'Location', 'north');
title('Bias');


% plot MSE
figure('DefaultAxesFontSize', 18)
true = kron(ones(n_reps, 1), b_true_vec);
plot(b_true_vec, mean((MA1est_Durbin' - true).^2,1), 'r-', ...
     b_true_vec, mean((MA1est_approxMLE' - true).^2,1), 'g:', ...
     b_true_vec, mean((MA1est_matlab' - true).^2,1), 'b--', ...
     'linewidth',2)
legend('Durbin','approx MLE', 'Matlab', 'Location', 'north');
title('MSE');


