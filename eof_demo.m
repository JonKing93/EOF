% This is a demo for using EOF_Analysis.
% I would recommend running the demo section by section...
clearvars;
close all;

%% 1. Make Data
% Let's first create some data for data variables A-G. We'll make A and B
% strongly anticorrelated and with large magnitudes. C will be strongly 
% correlated with A, but will have a normal magnitude. 
% D and E will be moderately correlated.
% F and G will just be noise.

rng('default');   % Seed the random number generator so that the demo gives consistent results

% Make A, B, and C
bigMagnitude = 200;
A = bigMagnitude * (randn(1000,1));
B = bigMagnitude * (-0.9*A + 0.1*randn(1000,1));
C = 0.9*A + 0.1*randn(1000,1);

% Make C and D
D = randn(1000,1);
E = 0.3*D + 0.7*randn(1000,1);

% Make E and F
F = randn(1000,1);
G = randn(1000,1);

% Get the final data table and variable names
data = [A B C D E F G];
varNames = {'A','B','C','D','E','F','G'};

%% 2. Run the Analysis with the default settings.
s = EOF_Analysis( data, 'varNames', varNames, 'showProgress');

%% 3. Output plots

% Figure 1
% Before looking at any significance testing results, we need to make sure
% that the Rule N significance test converged to a stable significance
% threshold. Looking at Figure 1, we can see that this was the case. We can
% see that running Rule N with 100 or fewer iterations can cause some large
% variability in the significance threshold. However, adding more
% iterations decreases this variability, and tests with more than about 250
% iterations give very consistent results. We did 1000 iterations, so we
% can be confident in the robustness of the results of the Rule N test. 
% Also, if we work with similar datasets, we can expect >250 iterations to
% be sufficient.

% Figure 2
% This shows the actual significance level that would be tested for a
% different number of Monte Carlo iterations. We can see that EOF_Analysis
% is always conservative, and tests a slightly more strict significance
% level if the exact significance level for the number of iterations is not
% an integer. Here, we can see that anything with more than ~300 iterations
% will test at significance levels between 0.05 and 0.047. Ensuring that
% the true significance level of the test matches the desired significance
% level requires that p*MC is an integer. This is equivalent to requiring
% MC to be a multiple of 1/p. For p = 0.05, any multiple of 20 will work,
% and the default 1000 iterations test exactly at the p=0.05 level.

% Figure 3
% Here we see the results of the significance test. The first two EOF modes
% explain more variance than can be expected from AR(1) noise, so we find
% them to be significant. The remaining modes explain less than this
% criteria so they cannot be distinguished from AR(1) noise.

% Figure 4, 5
% Figures 4 and 5 show the two significant modes (also known as loadings).
% The first mode describes a pattern for which A and C are strongly,
% positively correlated, and B is strongly anti-correlated. The second mode
% describes a pattern for which D and E are strongly correlated.

% Figures 6, 7
% Figures 6 and 7 show the modes after rotating modes 1 and 2. Here, there
% there is very little change from the unrotated modes. For the purposes of
% this demo, our data series are constructed using random Gaussian noise
% and should be roughly orthogonal. For analyses with data series that are
% not already orthogonal, the rotated modes can be significantly different
% from their unrotated counterparts.
