%Part two of homework two in Topics in Timseries Analysis
T = 100

arpaic = [];
arpbic = [];
maqaic = [];
maqbic = [];
 for i = 1:10000
     simdata = armasim(T,2,[0.4 -0.5 -0.2], [-0.5 -0.24]);
     arpaic(i) = inf;
     arpbic(i) = inf;
     maqaic(i) = inf;
     maqbic(i) = inf;

     for j = 1:10
         mle_arp = exactarp(simdata,j);
         mle_maq = maq(simdata, 1, j);
        if aic(mle_arp(end)^2, j, T) < arpaic(i)
            arpaic(i) = j;
        end
        if bic(mle_arp(end)^2, j, T) < arpbic(i)
            arpbic(i) = j;
        end
        if aic(mle_maq(end)^2, j, T) < maqaic(i)
            maqaic(i) = j;
        end
        if bic(mle_maq(end)^2, j, T) < maqbic(i)
            maqbic(i) = j;
        end
     end
    i
 end

function [AIC] = aic(sigma_squared, z, T)
    AIC = log(sigma_squared) + (2*z)/T;
end
function [BIC] = bic(sigma_squared, z, T)
    BIC = log(sigma_squared) + z*log(T)/T;
end

