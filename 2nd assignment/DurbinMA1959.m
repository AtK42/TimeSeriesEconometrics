%% computing the estimator of the MA(1) parameter
% First part is to please program his method of computing the estimator of 
%   the MA(1) parameter (I think it is just his equation 7), and the MA(q) 
%   parameters, this being, it seems, just equation 15. In both, as is 
%   clear from his article, the "a" terms are the AR estimated parameters, 
%   obtained from OLS: please see my chapter 4. So, easy peasy. For the 
%   MA(q) case, make a general program for arbitrary q, provided q is less 
%   than T, the number of observations in your time series -- your codes 
%   should check this! Below, we will use an MA(2) model, but please make 
%   your codes (i) specific for the q=1 case, and (ii) for the general q 
%   case.
 
% Notice that you may SKIP all the hard math and derivations in his paper! 
%   Thus, the assignment is very easy. I suggest you make ONE (perfect, 
%   beautiful, documented) computer program, called, say, DurbinMA1959 that
%   inputs a vector time series, and a value of q. If q=1, then you have an
%   IF THEN ELSE, and use his eq 7, which *obviously is a special case* of 
%   his eq 15, right? But the coding for q=1 is so much easier: So, a smart
%   computer science person would of course separate the two. And then, you
%   do the obvious -- you ensure that when you run the codes for the q>1 
%   case but in fact for q=1, you get exactly the same output. Be SURE you 
%   understand what I just said :-).
 
% Now for your simulations: For the MA(1) model, use T=100 observations, 
%   and a grid of true b values: -0.9, -0.8, ... 0, 0.1, up to 0.9. Use the
%   Durbin method in conjunction with simulation of the MA(1) process (it 
%   is easy to simulate an MA model, and I have codes in my book) to make 
%   nice performance *graphics* (not tables). Use your "smartz and common 
%   sense" to take it from here, as opposed to "uh, what exactly am I 
%   supposed to do?" Read carefully and THINK, and all should be clear.

%%% for MA(1) case
function MA_est = DurbinMA1959(vec_timeseries, q)
% function that estimates the coefficient of a MA(q) model
%--------------------------------------------------------------------------
% INPUT
% vec_timeseries                time series of a MA process
% q                             the order of the MA process
%--------------------------------------------------------------------------
% Output
% MA_est                        the estimated coefficient(s);
%                               real-valued for q = 1
%                               in R^q for q > 1
%--------------------------------------------------------------------------
    if nargin < 2
        warning('MA process and its assumed order is needed. q = 1 will be assumed')
        q = 1;
    end

    if q == 1
        a_vec_len = length(a_vec);
        num = a_vec(1:(a_vec_len-1))' * a_vec(2:a_vec_len);
        denom = sum(a_vec.^2);
        MA_est = - num / denom; 
    else
        if length(vec_timeseries) < q
            error('Please use a q smaller than the number of observations!')
        end
    end
end














































































































































































