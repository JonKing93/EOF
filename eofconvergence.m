function[ax] = eofconvergence(s)
%% Plots the convergence (or lack thereof) of the Rule N significance tests
%
% [ax] = plotRuleN(EOF)
% Plots the true confidence threshold and standardized significant
% eignvalue threshold for each iteration of the Monte Carlo Rule N
% significance test. Returns the axes handles for the plots.
%
%
% ----- Inputs -----
%
% EOF: A structure from an EOF_Analysis
%
%
% ----- Outputs -----
%
% ax: An array with the axes handles for the plots
%
%
% ----- Written By -----
% 
% Jonathan King, 2017, University of Arizona (jonking93@email.arizona.edu)

if ~isfield(s, 'MCsigExpVar') || ~isfield(s, 'MCtrue_p')
    warning('Insufficient data for the eofconvergence plot. Try running a convergence test...');
    return;
end

% Plot the true confidence threshold
figure();
plot(s.MCtrue_p)
xlabel('Number of Monte Carlo Iterations');
ylabel('True Significance Level of Significance Threshold.');
title('True Significance Level at Significance Thresholds for Monte Carlo Iterations ')
ax = gca;

% Plot the significant eigenvalues for each iteration
figure()
plot(zscore(s.MCsigExpVar));
title('Significance Threshold for Standaradized Random Explained Variances');
xlabel('Number of Monte Carlo Iterations');
ylabel('Standardized Explained Variance Significance Threshold');
ax = [ax; gca];

end