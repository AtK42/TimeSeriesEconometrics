function b_vec = maq_approxMLE(timeseries_vec,q)
n_obs = length(timeseries_vec);
timeseries_vec = reshape(timeseries_vec,n_obs,1);
initvec = zeros(q, 1);
opt=optimset('Display','off','TolX',1e-4,'MaxIter',1000,'LargeScale','off'); % set display option to 'iter' for info on each iteration

% minimize
b_vec = fminunc(@condmaq_,initvec,opt,timeseries_vec);

% change values back in case they were inverted before
for i = 1:q
    if abs(b_vec(i))>1
        b_vec(i)=1/b_vec(i);
    end
end


function loglik = condmaq_(param,timeseries_vec)
n_obs = length(timeseries_vec);
q = length(param);
uvec = zeros(n_obs+q, 1);

% check for b values outside the unity circle
for i = 1:q
    if abs(param(i))>1
        param(i)=1/param(i);
    end
end
for t=1:n_obs
    % invert order of u values to match order of timeseries
    pastu = flip(uvec(t:t+q-1));
    % calc latest error term
    uvec(t+q) = timeseries_vec(t) - param' * pastu;
end

% remove the initial values
uvec = uvec(q+1:end);

% calc loglik
loglik = sum(uvec.^2);
