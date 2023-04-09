function param = armaols(y,p,q)
% assumes zero mean stationary invertible ARMA(p,q)
% param = [AR terms, MA terms, sigma]
L=ceil(sqrt(length(y))); z=y(L+1:end);
Z=toeplitz(y(L:end-1),y(L:-1:1));
uhat=(eye(length(z)) - Z*inv(Z'*Z)*Z') * z; %#ok<*MINV>
sigmahat = std(uhat); yy=z-uhat; X=[]; m=max(p,q);
for i=1:p, X=[X z(m-i+1 : length(z)-i)]; end %#ok<*AGROW>
for i=1:q, X=[X uhat(m-i+1 : length(uhat)-i)]; end
yuse=yy(m+1:end); ARMAparam=inv(X'*X)*X'*yuse;
param=[ARMAparam ; sigmahat];