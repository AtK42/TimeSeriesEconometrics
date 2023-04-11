 function [MLE, stderr]=exactarp(y,p)
 n=length(y); y=reshape(y,n,1); initvec=[yw(y,p); std(y)];
 tol=1e-5; maxiter=200; show='none'; % 'iter','notify', or 'final'.
 options = optimset('Display',show,'TolX',tol,'Tolfun',tol, ...
 'MaxIter',maxiter,'LargeScale','off');
 [MLE,loglik,exitflag]=fminunc(@exactarp_,initvec,options,y,p);
 if nargout>1
 H = -hessian(@exactarp_,MLE,y,p);
 stderr=real(sqrt(diag(inv(H))));
 end

 function loglik=exactarp_(param,y,p)
 a=param(1:(end-1)); rr=roots([-a(end:-1:1); 1]); rootcheck=min(abs(rr));
 if rootcheck<=1, loglik=abs(0.999-rootcheck)*1e8; return, end
 sig=abs(param(end)); % this is not sigmaˆ2, but just (little) sigma.
 n=length(y); [Vi,detVi]=leeuwAR(a'); start=y(1:(p+1));
 lik0=0.5*log(detVi)-((n)/2)*log(2*pi*sig^2)-(0.5/sig^2)*start'*Vi*start;
 res=y((p+2):end); for i=1:p, res=res-a(i)*y(p+2-i:end-i); end
 ll=lik0-sum(res.^2)/2/sig^2; loglik = -ll;