k = [800 ; 400 ; 450; 290 ; 400];
R = randcorr(5);
S = diag([0.01;0.02;0.03;0.04;0.025]);
sigma = S*R*S;
mu = [0 ; 0 ; 0; 0; 0]%[0.01 ; -0.02 ; 0.03 ; 0.01 ; 0];

%Simulate 1000 returns for the 5 stocks
returns = [];
for i = 1:1000
    Z = mvnrnd(zeros(5,1),eye(5));
    D = diag(invGamRnd(k));
    returns(i, :) = mu + D*sqrtm(sigma)*Z';
end
figure;
plot(returns);
title('SMESTI distributed paths');
xlabel('t'); ylabel('value');
legend('1', '2', '3', '4', '5');


%Now do the second part of question 2 and estimate the CCC model
garchparams = [];
for i = 1:5
    garchparams(i,:) = babygarch(returns(:,i));
end

Q = corrcoef(returns);
R = diag(diag(Q).^(-1/2))*Q*diag(diag(Q).^(-1/2));

%R is the correlation matrix of the returns, and garchparams contains
%the elements of the 5 gaussian garch(1,1) models fit to the return series
garchparams
R



function value = invGamRnd(k)
value = [];
for i = 1:length(k)
v=k(i); Y = gamrnd(v/2, 1)/(v/2); value(i) = 1./Y;
end
end

 function [param,stderr,loglik,zvec] = babygarch(y)
 % normal-GARCH(1,1) with power=2. y is vector of log percentage returns
 initvec= [0 0.04 0.05 0.8];
 % mu c_0 c_1 d_1
 bound.lo= [-4 0 0 0 ];
 bound.hi= [ 4 0.5 1 1 ];
 bound.which=[ 1 1 1 1 ];
 opt=optimset('Display','None', 'Maxiter',500, 'TolFun',1e-6, ...
 'TolX',1e-6,'LargeScale','off');
 init=einschrk(initvec,bound);
 [pout,~,~,~,~,hess] = fminunc(@(param) like(param,y,bound),init,opt);
 [loglik,zvec]=like(pout,y,bound);
 V=pinv(hess)/length(y);
 [param,V]=einschrk(pout,bound,V); stderr=sqrt(diag(V));
 end

 function [loglik,zvec]=like(param,y,bound)
 param=einschrk(real(param),bound,999);
 meanterm=param(1); c0=param(2); c1=param(3); d1=param(4);
 e=y-meanterm; [zvec,sigvec]=ungarch(e,c0,c1,d1);
 K=sqrt(2*pi); ll = -0.5 * zvec.^2 - log(K) - log(sigvec);
 loglik = -mean(ll);
 end

 function [eout,sigvec]=ungarch(e,c0,c1,d1)
 sigvec=zeros(length(e),1); e2=e.^2;
 denom=1-c1-d1; if denom>0.001, sinit=c0/denom; else sinit=mean(e2); end
 einit=sinit;
 % do the recursion in sigvecˆdelta because it is faster
 sigvec(1)=c0+c1*einit+d1*sinit;
 for t=2:length(e), sigvec(t)=c0 + c1 *e2(t-1) + d1*sigvec(t-1); end
 sigvec=sigvec.^(1/2); eout=e./sigvec;
 end
