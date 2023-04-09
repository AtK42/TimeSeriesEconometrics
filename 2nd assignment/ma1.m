function [param, stderr, resid]=ma1(y,exact)
% computes the MLE (exact or conditional) of an MA(1) model
% set exact to 1 to compute the exact MLE, else the conditional MLE is
%   computed
ylen=length(y);
y=reshape(y,ylen,1);
initvec=[0 std(y)]';
opt=optimset('Display','off','TolX',1e-4,'MaxIter',1000,'LargeScale','off'); % set display option to 'iter' for info on each iteration
if exact==1
    [param,loglik,exitflag]=fminunc(@exactma1_,initvec,opt,y);
else
    [param,loglik,exitflag]=fminunc(@condma1_,initvec,opt,y);
end

b=param(1); littlesig=abs(param(2));
if abs(b)>1
    b=1/b;
    littlesig=littlesig/abs(b);
end

param=[b littlesig]';
if nargout>1
    if exact==1
        H = -hessian(@exactma1_,param,y);
        stderr=sqrt(diag(inv(H)));
    else
        H = -hessian(@condma1_,param,y);
        stderr=sqrt(diag(inv(H)));
    end
end
if nargout>2
    if exact==1
        Sigma=ma1Sigma(b,ylen);
        SigInv=inv(Sigma);
        [V,D]=eig(0.5*(SigInv+SigInv'));
        W=sqrt(D); SigInvhalf = V*W*V';
        resid = SigInvhalf*y/littlesig;
    else
        [garb,uvec]=condma1_(param,y);
        resid=uvec/littlesig;
    end
end

function [loglik,uvec]=condma1_(param,y)
ylen=length(y);
uvec=zeros(ylen,1); pastu=0;
b=param(1);
sig=abs(param(2)); % this is NOT sigmaˆ2, but just (little) sigma.
if abs(b)>1
    b=1/b; sig=sig/abs(b);
end
for t=1:ylen
    u=y(t)-b*pastu;
    uvec(t)=u;
    pastu=u;
end
ll = - ylen * log(sig) - sum(uvec.^2)/(2*sig.^2);
loglik = -ll;

function loglik=exactma1_(param,y)
ylen=length(y);
b=param(1);
sig=abs(param(2));
if abs(b)>1
    b=1/b;
    sig=sig/abs(b);
end
Sigma=ma1Sigma(b,ylen); % varcov matrix, but not scaled by little sigma.
Vi=inv(Sigma);
detVi=det(Vi);
if detVi<=0
    loglik=abs(detVi+0.01)*1e4;
    return
end
ll = - ylen * log(sig) + 0.5*log(detVi) - y'*Vi*y/(2*sig^2);
loglik = -ll;

function Sigma=ma1Sigma(b,ylen) % construct the varcov in a primitive way:
Sigma=zeros(ylen,ylen);
v=(1+b^2);
for i=1:ylen
    Sigma(i,i)=v;
end
for i=1:(ylen-1)
    Sigma(i,i+1)=b;
    Sigma(i+1,i)=b;
end

