function[ax] = eofconvergence(s)
%% Plots the evolution of the significance threshold for each iteration of Rule N.
%
% [ax] = eofconvergence(s)
% Plots the true confidence threshold and standardized significant
% eignvalue threshold for each iteration of the Monte Carlo Rule N
% significance test. Returns the axes handles for the plots.
%
%
% ----- Inputs -----
%
% s: An output structure from EOF_Analysis
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

if ~isfield(s, 'MCsigExpVar') || ~isfield(s, 'MCtrue_p') || ~isfield(s, 'MCnSig')
    warning('Insufficient data for the eofconvergence plot. Try running a convergence test...');
    return;
end

% Plot the true confidence threshold
figure();
plot(s.MCtrue_p)
xlabel('Number of Monte Carlo Iterations');
ylabel('True Significance Level of Significance Threshold.');
title('True Significance Level at Significance Thresholds for Monte Carlo Iterations ')
ax = [gca];

% Plot the significant eigenvalues for each iteration
figure()
plot(zscore(s.MCsigExpVar));
title('Significance Threshold for Standardized Random Explained Variances');
xlabel('Number of Monte Carlo Iterations');
ylabel('Standardized Explained Variance Significance Threshold');
if isfield(s, 'varNames')
    legend(s.varNames);
end
ax = [ax; gca]; 

% Plot the number of significant modes for each iteration
figure()
plot( s.MCnSig);
title('Number of significant leading modes for increasing Monte Carlo iterations');
xlabel('Number of Monte Carlo iterations');
ylabel('Number of significant leading modes');
ax = [ax; gca];

end