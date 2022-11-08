% Calculate one-locus equilibrium frequency for an allele in a fluctuating
% environment.
%
% Jeremy Van Cleve <vancleve@stanford.edu>
%
% 2008/08/13: creation date

function x = AeqFreq(mu, w11, w31, w12, w32, n)

fh = @(y) fNestedMapC(y, mu, w11, w31, w12, w32, n) - y;

initvals = linspace(0.1, 0.9, 5);
errorf = 1;

%disp(['mu ' num2str(mu) ' w11 ' num2str(w11) ' w31 ' num2str(w31) ' w12 ' num2str(w12) ' w32 ' num2str(w32) ...
%      ' n ' num2str(n)]);
%x = bisection(fh, 0, 1, 1e-6, 1e6);

if (sign(fh(0))*sign(fh(1)) < 0)
    x = fzero(fh, [0 1]);
else
    for i=1:5
        x = fzero(fh, initvals(i));
        if (x <= 1 && x >= 0)
            errorf = 0;
            break;
        end
    end
end

% if (errorf)
%     error(['no root for mu ' num2str(mu) ' w1 ' num2str(w1) ' w3 ' num2str(w3) ' n ' num2str(n)]);
% end

% calculate nested map
function f = fNestedMap(x, mu, w11, w31, w12, w32, n)

f = x;
for i=1:n
    f = fMap(f, mu, w11, w31);
end
for i=1:n
    f = fMap(f, mu, w32, w12);
end

% calculate one step iteration of allele frequency recursion
function f = fMap(x, mu, w1, w3)
    
f = ((1-mu)*w1*x + mu*w3*(1-x)) / (w1*x + w3*(1-x));
