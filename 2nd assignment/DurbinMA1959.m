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
%% estimation of the MA(q) parameters
function MA_est = DurbinMA1959(vec_timeseries, q)
% function that takes as input an MA(q) process and then estimates 
% the coefficient(s) of this MA(q) model
%--------------------------------------------------------------------------
% INPUT
% vec_timeseries                col vector, time series of an MA process
% q                             integer, the order of the MA process
%--------------------------------------------------------------------------
% OUTPUT
% MA_est                        the estimated coefficient(s);
%                               real-valued for q = 1
%                               in R^q for q > 1
%--------------------------------------------------------------------------
    if nargin < 2
        warning('MA process and its assumed order is needed. q = 1 will be assumed')
        q = 1;
    end

    % check that vec_timeseries is a col vec and get its length
    dim = size(vec_timeseries);
    if dim(1) == 1
        vec_timeseries = vec_timeseries';
    end
    %p = 20; % assuming an AR(p) process for the estimation of the MA params
    p = round(sqrt(length(vec_timeseries)));
    T = length(vec_timeseries);%-p;

    % estimate the AR parameters
    %Z_design_mat = zeros(T, p);
    %for i = 1:p
    %    Z_design_mat(:,i) = vec_timeseries(p+(1-i):T+p-i); % eq (6.31) in the book
    %end
    %a_vec = (Z_design_mat' * Z_design_mat) \ Z_design_mat' * vec_timeseries(p+1:end); % check eq (6.32) and take p = T, see also program listing 6.3
    %a_vec = [1; a_vec]; % set a_0 = 1 as suggested in Durbin's paper
    %a_vec_len = length(a_vec); % = k

    [~, a_vec] = yw(vec_timeseries, p);
    a_vec = [1; a_vec];

    % estimation for MA(1) case
    if q == 1
        %a_vec_len = length(a_vec);
        num = a_vec(1:(end-1))' * a_vec(2:end);
        denom = sum(a_vec.^2);
        MA_est = - num / denom; % eq 7 in Durbin(1959)
    
    % estimation for MA(q) case
    else
        % error checking
        if T < q
            error('Please use an MA order q smaller than the number of observations!')
        end
        LHS_mat = zeros(q, q);
        RHS_vec = zeros(q, 1);
        
        for i = 1:q
            for j = 1:q
                if i <= j
                    LHS_mat(i,j) = a_vec(1:(end-j))' * a_vec(j:end); % eq 15 in Durbin(1959)
                else
                    LHS_mat(i,j) = LHS_mat(j,i);
                end
            end
            RHS_vec(i) = a_vec(1:(end-i))' * a_vec(i:end);
        end
        MA_est = LHS_mat \ (-1 * RHS_vec);
    end
end

% code for computing (6.57) in the book (program listing 6.10, p. 306)
% T=length(y);
% p=round(sqrt(T));
% z=y(p+1:end);
% zl=length(z);
% Z=[];
% for i=1:p
%     Z=[Z y(p-i+1:p-i+zl)];
% end, a=inv(Z'*Z)*Z'*z;
% b=a(1)

% % function that computes the YW and least square estimator for an AR(p)
% %   model (see p. 292 in the book)
% % -------------------------------------------------------------------------
% function [ayw,aols]=yw(y,p)
% y=reshape(y,length(y),1); r=localsacf(y,p);
% v=[1; r(1:end-1)]; R=toeplitz(v,v); ayw=inv(R)*r;
% if nargout>1 % the o.l.s. estimator
%     z=y(p+1:end); zl=length(z); Z=[];
%     for i=1:p
%         Z=[Z y(p-i+1:p-i+zl)];
%     end
%     R=(Z'*Z) \ Z'; aols=R*z;
% end
% end
% 
% function acf=localsacf(x,imax) % computes the estimates of gamma_i
% T=length(x); a=zeros(imax,1);
% for i=1:imax
%     a(i)= sum(x(i+1:T) .* x(1:T-i) );
% end
% acf=a./sum(x.^2);
% end
% 
% 







































































































































































