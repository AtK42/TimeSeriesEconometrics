% function that computes the YW and least square estimator for an AR(p)
%   model (see p. 292 in the book)
% -------------------------------------------------------------------------
function [ayw,aols]=yw(y,p)
y=reshape(y,length(y),1); r=localsacf(y,p);
v=[1; r(1:end-1)]; R=toeplitz(v,v); ayw=inv(R)*r;
if nargout>1 % the o.l.s. estimator
    z=y(p+1:end); zl=length(z); Z=[];
    for i=1:p
        Z=[Z y(p-i+1:p-i+zl)];
    end
    R=(Z'*Z) \ Z'; aols=R*z;
end

function acf=localsacf(x,imax) % computes the estimates of gamma_i
T=length(x); a=zeros(imax,1);
for i=1:imax
    a(i)= sum(x(i+1:T) .* x(1:T-i) );
end
acf=a./sum(x.^2);



