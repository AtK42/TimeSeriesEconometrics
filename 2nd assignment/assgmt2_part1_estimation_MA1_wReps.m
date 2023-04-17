% Now for your simulations: For the MA(1) model, use T=100 observations, 
%   and a grid of true b values: -0.9, -0.8, ... 0, 0.1, up to 0.9. Use the
%   Durbin method in conjunction with simulation of the MA(1) process (it 
%   is easy to simulate an MA model, and I have codes in my book) to make 
%   nice performance *graphics* (not tables). Use your "smartz and common 
%   sense" to take it from here, as opposed to "uh, what exactly am I 
%   supposed to do?" Read carefully and THINK, and all should be clear.

% Next, please compare your findings above using a different method of 
%   MA(1) estimation, namely in my book, see my codes for the MLE and also 
%   the approximate MLE method. This is either in the chapter that does 
%   only MA(q) models, or in the chapter doing ARMA. You find it!

% You copy paste those codes and JUST USE the approximate method. Repeat: 
%   NOT the exact MLE method (because it takes too long to run, notably 
%   when done with the bootstrap.)  Do obviously the exact same simulation 
%   as above; make a graphic as before, but using my approx method, and 
%   compare performance. Got it?

% We are not done.
% Use a THIRD METHOD of parameter estimation in the MA(1) model, namely: 
%   In Matlab, google around, type help arma in Matlab, or recall it is 
%   written in my book (do a search for "matlab"). Matlab has a method 
%   built in, I think you need a particular toolbox to access it, but 
%   remember, with your uni license, you have access to all toolboxes. It 
%   is SUPER fast and works well. YOU do the research to figure out which 
%   function this is -- again, I think I list it in my book. I am not sure 
%   how Matlab's routine actually works --- you can attempt to read the 
%   help files and see if you can figure it out (and discuss it in your 
%   report). This is *optional* regarding figuring out how it works -- and 
%   I am not even sure the documentation is particularly clear anyway.

% So, as above: With this 3rd method, you make yet again some performance 
%   graphic, and OF COURSE you have a beautifully written, 
%   Shakespearean-English quality discussion of the comparison of the 3 
%   methods, in terms of bias and MSE of the estimators, along with 
%   computation time, and kick-ass graphics with lovely colors that make 
%   those gay-parade flags look boring.

% see listing 6.9 for code that simulates the exact and conditional MLE of
%   an MA(1) model and compares their bias and MSE for parameter b

% prep work
T=100;
%p = 20; % for AR estimates which are then used to estimate the MA parameters
n_obs = T;% + p;
b_true_vec = (-.9:.1:.9);
n_reps = 100;

% initialize variables
vec_timeseries = zeros(n_obs, length(b_true_vec));
MA1est_Durbin = zeros(length(b_true_vec), n_reps);
MA1est_approxMLE = zeros(length(b_true_vec), n_reps);
MA1est_approxMLE_alt = zeros(length(b_true_vec), n_reps);
MA1est_matlab = zeros(length(b_true_vec), n_reps);
time_Durbin = zeros(length(b_true_vec), 1);
time_approxMLE = zeros(length(b_true_vec), 1);
time_matlab = zeros(length(b_true_vec), 1);

total_start = tic;
% simulate MA(1) process and then estimate the parameter
for i = 1:length(b_true_vec)
    % simulation of MA process
    vec_timeseries(:,i) = armasim(n_obs, 1, b_true_vec(i), 0);
    
    Durbin_start = tic;
    for r = 1:n_reps
        % estimation of MA parameter
        % % Durbin
        MA1est_Durbin(i, r) = DurbinMA1959(vec_timeseries(:,i), 1);
    end
    time_Durbin(i) = toc(Durbin_start);

    approxMLE_start = tic;
    for r = 1:n_reps
        % % approximate MLE
        MA1_temp = ma1(vec_timeseries(:,i), 0); % 0 for approx MLE
        MA1est_approxMLE(i, r) = MA1_temp(1);
        MA1_temp = armaols(vec_timeseries(:,i),0,1);
        MA1est_approxMLE_alt(i, r) = MA1_temp(1);
    end
    time_approxMLE(i) = toc(approxMLE_start);

    matlab_start = tic;
    for r = 1:n_reps
        % % bulilt-in MATLAB function
        MA1_temp = armax(vec_timeseries(:,i), [0, 1]);
        MA1est_matlab(i, r) = MA1_temp.c(2:end);
    end
    time_matlab(i) = toc(matlab_start);
    disp(['done for b = ' num2str(b_true_vec(i))])
end
time_total = toc(total_start);


%% plotting

% Durbin_dev = (MA1est_Durbin) - b_true_vec';
% approxMLE_dev = MA1est_approxMLE - b_true_vec';
% matlab_dev = MA1est_matlab - b_true_vec';
% dev_cat = [Durbin_dev, approxMLE_dev, matlab_dev];


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


%%
% cats = categorical({'Durbin', 'MLE', 'Matlab'});
% tiledlayout(4,5)
% for i = 1:length(b_true_vec)
%     nexttile
%     scatter(cats, MA1est_cat(i,:), 25, 'filled')
%     ylim([-1 1])
%     yline(b_true_vec(i))
% end

% scatter(cats, MA1est_cat, 'filled')
% legend
