function[] = eofplot(s)
% A driver to create all the plots for an EOF_Analysis.

% Plot the results of the significance test
eofsignificance(s);

% Plots the loadings of the rotated and unrotated modes
eofloadings(s);

% Plot the Rule N convergence data
eofconvergence(s);
end

