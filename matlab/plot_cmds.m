muvals = [logspace(-3, log10(1/2), 101)];

s1 = 0.01;
s2 = 0.015;

% r = 0.0
plotandsave_pip(eigPIP(muvals, 1, 1-s1, 1, 1-s2, 0.0,50), 'pip_n50_s1_0-01_s2_0-015_r_0-0.pdf', muvals)
% r = 0.05
plotandsave_pip(eigPIP(muvals, 1, 1-s1, 1, 1-s2, 0.05, 50), 'pip_n50_s1_0-01_s2_0-015_r_0-05.pdf', muvals)
% r = 0.1
plotandsave_pip(eigPIP(muvals, 1, 1-s1, 1, 1-s2, 0.1, 50), 'pip_n50_s1_0-01_s2_0-015_r_0-1.pdf', muvals)
% r = 0.15
plotandsave_pip(eigPIP(muvals, 1, 1-s1, 1, 1-s2, 0.15, 50), 'pip_n50_s1_0-01_s2_0-015_r_0-15.pdf', muvals)