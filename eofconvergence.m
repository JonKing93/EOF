function[ax] = plotRuleN(EOF)
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

if ~isfield(EOF, 'iterTrueConf') || ~isfield(EOF, 'iterSigEigs')
    error('EOF does not contain convergence data. Rule N plots are not possible.');
end

% Plot the true confidence threshold
figure();
plot(EOF.iterTrueConf)
xlabel('Number of Monte Carlo Iterations');
ylabel('True Confidence Level of Eigenvalues');
title('True confidence level of iteratively selected significant eigenvalues.')
ax = gca;

% Plot the significant eigenvalues for each iteration
figure()
plot(zscore(EOF.iterSigEigs));
title('Standardized Significance Threshold for Random Eigenvalues');
xlabel('Number of Monte Carlo Iterations');
ylabel('Standardized Significance');
ax = [ax; gca];