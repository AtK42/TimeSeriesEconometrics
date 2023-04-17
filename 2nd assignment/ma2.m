function param=ma1(y)
% computes the MLE (conditional) of an MA(1) model

ylen=length(y);
y=reshape(y,ylen,1);
initvec=[0 0 std(y)]';
opt=optimset('Display','off','TolX',1e-4,'MaxIter',1000,'LargeScale','off'); % set display option to 'iter' for info on each iteration
[param,~,~]=fminunc(@condma1_,initvec,opt,y);

%b=param(1); littlesig=abs(param(2));
b1=param(1); b2=param(2); littlesig=abs(param(3));
if abs(b)>1
    b=1/b;
    littlesig=littlesig/abs(b);
end
param=[b littlesig]';


function [loglik,uvec]=condma1_(param,y)
ylen=length(y);
%uvec=zeros(ylen,1); pastu=0;
uvec=zeros(ylen,1); pastu=[0 0];
%b=param(1);
b1=param(1); b2 = param(2);

%sig=abs(param(2)); % this is NOT sigmaˆ2, but just (little) sigma.
sig = abs(param(3));

%if abs(b)>1
if abs(b1)>1
    %b=1/b; sig=sig/abs(b);
    b1=1/b1; b2=1/b2; sig=sig * abs(sqrt((1-b1^(-2)-b2^(-2)) / (1-b1^2-b2^2)));
end
for t=1:ylen
    %u=y(t)-b*pastu;
    u=y(t)-b1*pastu(1)-b2*pastu(2);
    %uvec(t)=u;
    uvec(t)=u;
    %pastu=u;
    pastu(2)=pastu(1);
    pastu(1)=u;
end
ll = - ylen * log(sig) - sum(uvec.^2)/(2*sig.^2);
loglik = -ll;
