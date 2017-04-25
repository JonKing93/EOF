function[] = EOFsigplot(EOF, varargin)
%% Makes a significance plot for the percent explained variance of data and
% Monte Carlo eigenvectors.
%
% EOF_Sig_Plot(EOF)
%
% EOF_Sig_Plot(EOF, nEigs)
%
% ----- Inputs -----
%
% EOF: The output of the EOF_Analysis function
%
% *** Optional Inputs ***
%
% nEigs: The number of eignevectors to show in the plot. (Default = 10)


figure(); clf; hold on;

% Plot the random eigenvalues at the tested threshold
plot( EOF.randEigvals( EOF.thresh,:), 'r-')

% Plot the data explained variances.
stem( EOF.expVar);

title('Significant Eigenvalues');
legend( sprintf('Minimum Significance Level, p = %1.2f', 1-EOF.conf), 'Eigenvalue Explained Variance');
xlabel('Eigenvalue Rank');
ylabel('Explained Variance');
xs = get(gca, 'XLim');
xs = [0 xs(2)+1];
set(gca,'XLim',xs);