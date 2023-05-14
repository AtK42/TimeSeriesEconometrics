% Here goes: Have a look at my mixed-normal GARCH construction, on page 478
% Your first task is to *simulate* this process. Use k=2 components, though
%   I might later ask you to try 3. You do NOT have to estimate the model, 
%   though it is also not difficult to set up the likelihood and estimate 
%   it. You just need to simulate it, and you can play around with the 
%   required parameters so that when you plot the resulting univariate time 
%   series, it should "look like actual daily stock returns".  Since you 
%   are simulating, you can simulate T=5000 observations -- or even more, 
%   all for free.

% The next part is to estimate a student-t-GARCH model based on your 
%   simulated data. Notice that the true DGP does NOT match the model you 
%   are estimating. That is okay -- that is reality. The student-t garch 
%   estimation is very easy to set up -- you just build on the codes I have
%   in my chapter. Obviously, the degrees of freedom of the student t gets 
%   jointly estimated with the garch parameters.  Surely, you can find 
%   toolboxes for matlab and python that also do it, ...)

% I think that is already a good start to the 3rd homework, and I will 
%   surely build on this, such as asking you to make similar VaR plots as 
%   I have on page 492. (So, you can do that now also, and thus I entice 
%   you to read some of chapter 11.)