function[ax] = eofsignificance( s )
%% Plots the results of the Rule N significance test against EOF explained variances.
%
% [ax] = eofsignificance(s)
% Plots EOF explained variances against the Rule N explained variance
% significance threshold to test for EOF significance.
%
%
% ----- Inputs -----
%
% s: The output structure from EOF_Analysis.
%
%
% ----- Outputs -----
%
% ax: The axis handle for the plot.
%
%
% ----- Written By -----
% 
% Jonathan King, 2017, University of Arizona (jonking93@email.arizona.edu)

% Check that the required fields exist
if ~isfield(s, 'sigExpVar') || ~isfield(s, 'expVar') || ~isfield(s, 'p') || ~isfield(s,'nSig')
    warning('There is insufficient data for eofsignificance plot. Try running a significance test...');
    return;
end    

% Initialize a figure
figure();
hold on;

% Plot the significance threshold line
plot( s.sigExpVar, 'r-');

% Plot the actual explained variances
stem( s.expVar );
    
% Make everything look nice.
title(sprintf('Significance of EOF modes\r\nThe first %0.f modes are significant',s.nSig));
legend( sprintf('Significance Threshold, p < %0.2f', s.p), 'Explained Variance');
xlabel('EOF Mode (decreasing rank order)');
ylabel('Explained Variance');
xs = get(gca, 'XLim');
xs = [0 xs(2)+1];
ax = gca;
set(gca,'XLim',xs);
end