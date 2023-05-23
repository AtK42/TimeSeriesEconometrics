function [param, stderr, loglik, zvec] = babygarch(y)
    % normal-GARCH(1,1) with power=2. y is vector of log percentage returns
    initvec = [0.04 0.05 0.8 20]; % mu c_0 c_1 d_1 dof
    bound.lo = [0 0 0 0];
    bound.hi = [0.5 1 1 200];
    bound.which = [1 1 1 1];
    opt = optimset('Display', 'None', 'Maxiter', 500, 'TolFun', 1e-6, ...
        'TolX', 1e-6, 'LargeScale', 'off');
    init = einschrk(initvec, bound);
    [pout, ~, ~, ~, ~, hess] = fminunc(@(param) like(param, y, bound), init, opt);
    [loglik, zvec] = like(pout, y, bound);
    V = pinv(hess) / length(y);
    [param, V] = einschrk(pout, bound, V);
    stderr = sqrt(diag(V));
end

function [loglik, zvec] = like(param, y, bound)
    param = einschrk(real(param), bound, 999);
    %meanterm = param(1);
    c0 = param(1);
    c1 = param(2);
    d1 = param(3);
    dof = param(4);
    %e = y - meanterm;
    %[zvec, sigvec] = ungarch(e, c0, c1, d1);
    [zvec, sigvec] = ungarch(y, c0, c1, d1);
    %K = sqrt(2 * pi);
    ll = gammaln((dof+1)/2) - gammaln(dof / 2) - log(sqrt(dof * pi)) - log(sigvec) - ((dof + 1)/2) * log(1 + zvec.^2/dof);
    %ll = -0.5 * zvec.^2 - log(K) - log(sigvec);
    loglik = -mean(ll);
end

function [eout, sigvec] = ungarch(e, c0, c1, d1)
    sigvec = zeros(length(e), 1);
    e2 = e.^2;
    denom = 1 - c1 - d1;
    if denom > 0.001
        sinit = c0 / denom;
    else
        sinit = mean(e2);
    end
    einit = sinit;
    sigvec(1) = c0 + c1 * einit + d1 * sinit;
    for t = 2:length(e)
        sigvec(t) = c0 + c1 * e2(t-1) + d1 * sigvec(t-1);
    end
    sigvec = sigvec.^(1/2);
    eout = e ./ sigvec;
end

function [pout,Vout]=einschrk(pin,bound,Vin)
% [pout,Vout]=einschrk(pin,bound,Vin)
% if Vin specified, then pout is untransformed, otherwise pout is transformed
% M. Paolella, 1997

welche=bound.which;
if all(welche==0) % no bounds!
  pout=pin; 
  if nargin==3, Vout=Vin; end
  return
end

lo=bound.lo; hi=bound.hi;
if nargin < 3
  trans=sqrt((hi-pin) ./ (pin-lo));
  pout=(1-welche).* pin + welche .* trans;
  Vout=[];
else
  trans=(hi+lo.*pin.^2) ./ (1+pin.^2);
  pout=(1-welche).* pin + welche .* trans;
  % now adjust the standard errors
  trans=2*pin.*(lo-hi) ./ (1+pin.^2).^2;
  d=(1-welche) + welche .* trans; % either unity or delta method.
  J=diag(d);
  Vout = J * Vin * J;
end
end

