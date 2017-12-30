% EOF
%
% This folder contains functions involved in EOF Analysis.
%
% The function EOF_Analysis runs an analysis in full, but the sub-functions
% may also be used individually.
%
% Main Function:
%   EOF_Analysis  - Performs a full EOF Analysis of a data set. Uses a
%                   Rule N significance test and performs a Varimax rotation
%                   of significant modes.   
%
% Figure Functions:
%   eofsignificance - Plots the results of the Rule N significance test,
%                     and shows significnant EOF modes.
%   eofloadings - Plots the loadings of data variables on significant and
%                 rotated EOF modes.
%   eofconvergence - Plots the change in significance threshold with
%                    increasing Monte Carlo iterations.
%
% Analysis Functions:
%   scaleSignals - Scale EOF signals (also known as PCs) for comparison
%                  with data in the EOF analysis matrix.
%
%   ruleN - Implements a "Rule N" Monte Carlo significance test in parallel or serial.
%   runRuleN - Generates a set of random explained variances for each
%              "Rule N" Monte Carlo iterations
%
%   eofSigThreshold - Determines the explained variance significance
%                     threshold based on "Rule N" explained variances.
%   ruleNConvergence - Calculates significance thresholds for successive
%                      Monte Carlo iterations to examine if Rule N
%                      converges to a stable threshold.
%
%   eofrotation - Implements the rotation of significant eof modes.
%
% Misc. Functions:
%   parseInputs - A general function to parse input flags and values.
%   progressbar - Implements a progress bar. (I DID NOT WRITE THIS. PLEASE SEE "progressbar_license.txt")
%   randNoiseSeries - Creates a random time series with specified noise properties 
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
%
%   "progressbar.m" by Steve Hoelzer. Please see "progressbar_license.txt"
%   for details.
