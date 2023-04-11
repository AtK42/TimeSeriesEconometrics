 function [param]=maq(y,exact,q)
 ylen=length(y); y=reshape(y,ylen,1); initvec=[zeros(1,q) std(y)]';
 opt=optimset('Display','off','TolX',1e-3,'MaxIter',100,'LargeScale','off');
 if exact==1, [param,loglik,exitflag]=fminunc(@exactmaq_,initvec,opt,y);
 end
 end
 function loglik=exactmaq_(param,y)
 ylen=length(y); b=param(1:length(param)-1); sig=abs(param(length(param)));
 q = length(param) - 1;
 Sigma=maqSigma(b,ylen, q); % varcov matrix, but not scaled by little sigma.
 Vi=inv(Sigma); detVi=det(Vi);
 if detVi<=0, loglik=abs(detVi+1)*1e4; return, end
 ll = - ylen * log(sig) + 0.5*log(detVi) - y'*Vi*y/(2*sig^2); loglik = -ll;
 end

 function Sigma=maqSigma(b,ylen, q); % construct the varcov in a primitive way:
 M=zeros(ylen,ylen); 
 for i = 1:ylen
    M(i,i) = 1;
    for j = i + 1:ylen
        if(j-i <= q)
            M(j,i) = b(j-i);
        end
    end
 end
 N1 = zeros(q, q);
 for i = 1:q
    for j = i:q
        N1(i,j) = b(q - (j-i));
    end
 end
 N = [N1 ; zeros(ylen-q, q)];
 Sigma = M*(M') + N*(N');

 end