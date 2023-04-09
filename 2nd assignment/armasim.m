function y = armasim(nobs,sig2,pv,qv)
% function that simulates nobs observations of an ARMA process
% -------------------------------------------------------------------------
% INPUT
% sig2                          integer, innovation variance
% pv                            vector, AR parameters
% qv                            vector, MA parameters
% -------------------------------------------------------------------------
% OUTPUT
% y                             ARMA process
% -------------------------------------------------------------------------
if nargin<4
    qv=[];
end

% some prep work
p=length(pv); q=length(qv);
pv=reshape(pv,1,p); qv=reshape(qv,1,q);
warmup=500;
e=sqrt(sig2)*randn(nobs+warmup,1);

% initialize variables
%init=0;
evec=zeros(q,1);
yvec=zeros(p,1);
y=zeros(nobs+warmup,1);

% simulate ARMA process
for i=1:nobs+warmup
    if p>0
        y(i) = y(i) + pv*yvec;
    end
    if q>0
        y(i) = y(i) + qv*evec;
    end
    y(i) = y(i) + e(i);

    % for AR part
    if p>1
        yvec(2:p)=yvec(1:p-1);
    end
    yvec(1)=y(i);

    % for MA part
    if q>1
        evec(2:q)=evec(1:q-1);
    end
    evec(1)=e(i);
end

% exclude values from warmup period
y=y(warmup+1:end);