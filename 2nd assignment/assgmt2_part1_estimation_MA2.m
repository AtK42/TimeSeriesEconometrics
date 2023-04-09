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
MA2_temp = ma1(vec_timeseries_100, 0); % 0 for approx MLE
MA2est_approxMLE_100 = MA2_temp;

MA2_temp = ma1(vec_timeseries_1000, 0);
MA2est_approxMLE_1000 = MA2_temp(1);

% %% bulilt-in MATLAB function
MA2_temp = armax(vec_timeseries_100, [0, 2]);
MA2est_matlab_100 = MA2_temp.c(2:end);

MA2_temp = armax(vec_timeseries_1000, [0, 2]);
MA2est_matlab_1000 = MA2_temp.c(2:end);
