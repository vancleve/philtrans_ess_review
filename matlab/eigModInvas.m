% Calculate leading eigenvalue (and derivative with respect to resident mutation rate)
% of product of jacobian matrices for modifier allele invasion
%
% Jeremy Van Cleve <vancleve@stanford.edu>
%
% 2009/08/13: creation date

function [ev dev d2ev] = eigModInvas(muRes, muMut, w11, w31, w12, w32, r, n)

% calculate allele frequencies and mean fitness in cycle
orbit = zeros(2*n,1);
wbar = zeros(2*n,1);
orbit(1) = AeqFreq(muRes, w11, w31, w12, w32, n);
wbarprod = 1;
for i=2:n+1
    orbit(i) = fMap(orbit(i-1), muRes, w11, w31);
    wbar(i-1) = w11*orbit(i-1) + w31*(1-orbit(i-1));
    wbarprod = wbarprod * wbar(i-1);
end
for i=n+2:2*n
    orbit(i) = fMap(orbit(i-1), muRes, w32, w12);
    wbar(i-1) = w32*orbit(i-1) + w12*(1-orbit(i-1));
    wbarprod = wbarprod * wbar(i-1);
end
wbar(2*n) = w11*orbit(2*n) + w31*(1-orbit(2*n));
wbarprod = wbarprod * wbar(2*n);

jprod = eye(2);
jprod2 = eye(2);
djprod = zeros(2);
d2jprod = zeros(2);
for i=1:n
    jprod = jac(orbit(i), muMut, w11, w31, r) * jprod;
    if (nargout > 2)
        d2jprod = jac(orbit(i), muRes, w11, w31, r) * d2jprod ...
                 + 2*djac(orbit(i), w11, w31, r) * djprod;
    end
    if (nargout > 1)
        djprod = jac(orbit(i), muRes, w11, w31, r) * djprod ...
                 + djac(orbit(i), w11, w31, r) * jprod2;
        jprod2 = jac(orbit(i), muRes, w11, w31, r) * jprod2;
    end
end
for i=n+1:2*n
    jprod = jac(orbit(i), muMut, w32, w12, r) * jprod;
    if (nargout > 2)
        d2jprod = jac(orbit(i), muRes, w32, w12, r) * d2jprod ...
                 + 2*djac(orbit(i), w32, w12, r) * djprod;
    end
    if (nargout > 1)
        djprod = jac(orbit(i), muRes, w32, w12, r) * djprod ...
                 + djac(orbit(i), w32, w12, r) * jprod2;
        jprod2 = jac(orbit(i), muRes, w32, w12, r) * jprod2;
    end
end

%if (any(isnan(jprod)))
%      disp(['muRes ' num2str(muRes) ' muMut ' num2str(muMut) ...
%            ' w11 ' num2str(w11) ' w31 ' num2str(w31) ' w12 ' num2str(w12) ' w32 ' num2str(w32) ...
%            ' r ' num2str(r) ' n ' num2str(n)]);
%end
% orbit'
% jprod
% eig(jprod)
ev = eig(jprod);
absev = abs(eig(jprod));
ev = ev(absev==max(absev));
if (length(ev) > 1)
    ev = ev(1);
end

if (nargout > 1)
    trjprod = trace(jprod);
    dertrjprod = trace(djprod);

    detjprod = det(jprod);
    if (abs(detjprod) < eps)
        detjprod = 0;
    end
    derdetjprod = -4 * n * detjprod / (1 - 2*muRes);
    dev = (ev * dertrjprod - derdetjprod) / (2 * ev - trjprod);
%     derdetjprod = -4 * n * (1 - 2*muRes)^(2*n-1) * (1 - r)^(2*n) * w11^n * w31^n * w12^n * w32^n / wbarprod^2
%     dev = (ev * dertrjprod - derdetjprod) / (2 * ev - trjprod)
    
    
    if (nargout > 2)
        der2trjprod = trace(d2jprod);
        der2detjprod = 8 * n * (2*n-1) * detjprod / (1 - 2*muRes)^2;
        d2ev = (2*dertrjprod*dev + der2trjprod*ev - 2*dev^2 - der2detjprod) / (2*ev - trjprod);
%         der2detjprod = 8 * n * (2*n-1) * (1 - 2*muRes)^(2*n-2) * (1 - r)^(2*n) * w11^n * w31^n * w12^n * w32^n / wbarprod^2
%         d2ev = (2*dertrjprod*dev + der2trjprod*ev - 2*dev^2 - der2detjprod) / (2*ev - trjprod)
    else
        d2ev = NaN;
    end
else
    dev = NaN;
end
ev = max(abs(eig(jprod)));

% calculate jacobian matrix for modifier allele invasion
function j = jac(x, mu, w1, w3, r)

j = [ (1-mu)*w1 - r*((1-mu)*w1 - mu*w3)*(1-x), mu*w3 + r*((1-mu)*w1 - mu*w3)*x;
      mu*w1 + r*((1-mu)*w3 - mu*w1)*(1-x), (1-mu)*w3 - r*((1-mu)*w3 - mu*w1)*x ] / (w1*x + w3*(1-x));

function dj = djac(x, w1, w3, r)

dj = [ -w1 + r*(1-x)*(w1+w3), w3 - r*x*(w1+w3);
       w1 - r*(1-x)*(w1+w3), -w3 + r*x*(w1+w3) ] / (w1*x + w3*(1-x));
  
% calculate one step iteration of allele frequency recursion
function f = fMap(x, mu, w1, w3)
    
f = ((1-mu)*w1*x + mu*w3*(1-x)) / (w1*x + w3*(1-x));