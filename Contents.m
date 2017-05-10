% EOF
%
% This folder contains functions involved in EOF Analysis.
%
% The function EOF_Analysis runs an analysis in full, but the functions may
% also be used individually.
%
% Main Function:
%   EOF_Analysis  - Performs a full EOF Analysis of a data set. Uses a
%                   Rule N significance test and performs a Varimax rotation
%                   of significant modes.   
%
% Figure Functions:
%   EOFloadings  - Plots the loadings of data series on significant modes.
%   EOFsigplot   - A scree plot showing the modes that pass the significance test.
%   plotRuleN    - Plots the convergence (or lack thereof) of the Rule N Monte Carlo tests.
%
% Operational Functions:
%   simpleEOF    - Gets the EOF modes, and associated eigenvalues and explained variance for a data matrix
%   getSignals   - Gets the signals from a standardized dataset and its eigenvectors / modes.
%   scaleSignals - Scales signals to the standard deviation of the standardized data.
%   ruleN        - Runs a Rule N significance test on a data matrix and its eigenvalues / loadings.
%   scaleEigvecs - Scales the eigenvectors by the square root of the eigenvalues. Used for VARIMAX rotation.
%   varimaxRotation - Performs a VARIMAX rotation on a set of scaled eigenvectors and eigenvalues.
%
% Helper Functions:
%   hasNaN       - Tests if an N-dimensional array contains any NaN elements
%   posColSign   - Reverses the sign of a vector if the majority of elements are negative        
%   randNoiseSeries - Creates a random time series with specified noise properties 
%   quickSVD     - Runs an SVD decomposition on a matrix.
%
% References:
%   Deser, C., and M. L. Blackmon (1993), Surface climate variations over the
%   North Atlantic Oceanduring winter:  1900-1989, Journal of Climate,6(9),
%   1743-1753.
%
%   Principal Component Analysis in Metereology and Oceanography. Rudolph
%   Preisendorfer. Elsevier Science Publishers. New York, 1988.
%
% Written by: 
%   Jonathan King, University of Arizona (jonking93@email.arizona.edu)
%
% Acknowledgements:
%   This work based on assignments and material from the course "Spatiotemporal Data Analysis",
%   presented by Kevin Anchukaitis, University of Arizona, Fall 2016.
