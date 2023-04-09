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

% initialize variables
vec_timeseries = zeros(n_obs, length(b_true_vec));
MA1est_Durbin = zeros(length(b_true_vec), 1);
MA1est_approxMLE = zeros(length(b_true_vec), 1);
MA1est_approxMLE_alt = zeros(length(b_true_vec), 1);
MA1est_matlab = zeros(length(b_true_vec), 1);

% simulate MA(1) process and then estimate the parameter
for i = 1:length(b_true_vec)
    % simulation of MA process
    vec_timeseries(:,i) = armasim(n_obs, 1, b_true_vec(i), 0);

    % estimation of MA parameter
    % % Durbin
    MA1est_Durbin(i) = DurbinMA1959(vec_timeseries(:,i), 1);

    % % approximate MLE
    MA1_temp = ma1(vec_timeseries(:,i), 0); % 0 for approx MLE
    MA1est_approxMLE(i) = MA1_temp(1);
    MA1_temp = armaols(vec_timeseries(:,i),0,1);
    MA1est_approxMLE_alt(i) = MA1_temp(1);

    % % bulilt-in MATLAB function
    MA1_temp = armax(vec_timeseries(:,i), [0, 1]);
    MA1est_matlab(i) = MA1_temp.c(2:end);

    disp(['done for b = ' num2str(b_true_vec(i))])
end





