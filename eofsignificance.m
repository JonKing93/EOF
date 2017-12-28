function[f] = eofsignificance( s )
%% Makes a significance plot for the percent explained variance of data and
% Monte Carlo eigenvectors.
%
% [ax] = EOF_Sig_Plot(EOF)
% Makes a significance plot with all the eigenvectors. Returns the
% axis handle for the plot.
%
% [ax] = EOF_Sig_Plot(EOF, nEigs)
% Makes a significance plot with a user specified number of eigenvectors.
%
%
% ----- Inputs -----
%
% EOF: The output of the EOF_Analysis function
%
% nEigs: The number of eignevectors to show in the plot.
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
f = figure();
hold on;

% Plot the significance threshold line
plot( s.sigExpVar, 'r-');

% Plot the actual explained variances
stem( s.expVar );
    
% Make everything look nice.
title(sprintf('Significance of EOF modes\r\nThe first %0.f modes are significant',s.nSig));
legend( sprintf('Significance Threshold, p < %1.2f', s.p), 'Explained Variance');
xlabel('EOF Mode (decreasing rank order)');
ylabel('Explained Variance');
xs = get(gca, 'XLim');
xs = [0 xs(2)+1];
ax = gca;
set(gca,'XLim',xs);
end