% Calculate pairwise invasibility plot using leading eigen value of
% jacboian of mutation rate modifier
%
% Jeremy Van Cleve <vancleve@stanford.edu>
%
% 2009/10/05: creation date

function pip = eigPIP(muvals, w11, w31, w12, w32, r, n)

lenmuvals = length(muvals);
pip = zeros(lenmuvals, lenmuvals);

for i=1:lenmuvals
    for j=1:lenmuvals
        pip(i,j) = eigModInvas(muvals(j), muvals(i), w11, w31, w12, w32, r, n);
    end
end