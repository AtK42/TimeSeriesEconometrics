function y=ma1sim(nobs,sigma,b)
u=sigma*randn(nobs,1);
y=zeros(nobs,1);
y(1)=u(1)+b*sigma*randn(1,1);
for i=2:nobs
    y(i) = u(i) + b*u(i-1);
end