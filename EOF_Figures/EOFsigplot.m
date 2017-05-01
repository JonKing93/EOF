function[ax] = EOFsigplot(EOF, varargin)
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

% Get the number of eigenvalues to plot
if isempty(varargin) 
    nEigs = size(EOF.randEigvals,2);
else
    [nEigs] = varargin{1};
end

figure(); clf; hold on;

% Plot the random eigenvalues at the tested threshold
plot( EOF.randEigvals( EOF.thresh, 1:nEigs), 'r-')

% Plot the data explained variances.
stem( EOF.expVar(1:nEigs));

% Make everything look nice
title('Significant Eigenvalues');
legend( sprintf('Minimum Significance Level, p = %1.2f', 1-EOF.conf), 'Eigenvalue Explained Variance');
xlabel('Eigenvalue Rank');
ylabel('Explained Variance');
xs = get(gca, 'XLim');
xs = [0 xs(2)+1];
ax = gca;
set(gca,'XLim',xs);