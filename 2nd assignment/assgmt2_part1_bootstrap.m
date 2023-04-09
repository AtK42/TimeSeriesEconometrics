 
% The last part of this segment of the homework is to apply the parametric 
%   bootstrap to ALL 3 methods in the MA(1) case. Thus, you need to have 
%   these 3 methods programmed as perfect black-boxes. You use the 
%   parametric bootstrap to get 90% confidence intervals for the unknown 
%   MA parameter, and you need to *simulate the bootstrap method* in order 
%   to report the actual coverage obtained. I offhand presume B=400 
%   replications should be enough, but be smart and active --- recall the 
%   choice of B is ideally infinity, but is limited by computation time. 
%   Try to ascertain what value of B is adequate --- that means, small 
%   enough such that the results are stable.
 
% This is the usual bootstrap thing we did last semester (chapters 1 and 2 
%   of my statistics book) You use the SINGLE bootstrap, as opposed to 
%   double. And notice I said parametric, and not nonparametric. You are 
%   welcome of course to also do nonparametric, but it is harder, because 
%   you have to simulate MA processes with innovations taken from the 
%   filtered ones from your fitted model. I spare you this. So, parametric 
%   is easier, so just use that.
 
% You report the quality (meaning actual coverage) of nominal (single 
%   parametric bootstrap) 90% CIs for the MA(1) parameter, based on the 3 
%   methods: Durbin, my approximate method of estimation, and built into 
%   Matlab.
