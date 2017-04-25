function[scaSignals] = scaleSignals(signals, eigVals)
%% Scales signals to the standard deviation of the standardized data.
%
% [sigScaled] = scaleSignals(signals, eigVals)
%
% ----- Inputs -----
%
% signals: a set of signals from an EOF analysis. Each column is one signal
%
% eigVals: A matrix containing the associated eigenvalues for each set of
%   signals Each column contains the eigenvalues for one set of signals. The
%   eigenvalues SHOULD NOT be normalized.
%
%
% ----- Outputs -----
%
% scaSignals: The set of scaled signals. These signals may be plotted
%   directly against the standardized dataset.

errCheck(signals, eigVals);

% Make the eigenvalues a row
if ~isrow(eigVals)
    eigVals = eigVals';
end

% Scale each signal by the corresponding eigenvalue
scaSignals = signals ./ sqrt(eigVals);
end

%%%%% Helper Functions %%%%%
function[nsignals] = errCheck(signals, eigVals)

% Ensure signals is a matrix
if ~ismatrix(signals)
    error('signals must be a matrix');
end

% Ensure eigVals is a vector
if ~isvector(eigVals)
    error('eigVals must be a vector');
end

% Ensure there are no NaNs
if hasNaN(signals)
    error('signals cannot contain NaN');
end
if hasNaN(eigVals)
    error('eigVals cannot contain NaN');
end

% Ensure eigVals are positive
if any(any( eigVals < 0))
    error('The eigenvalues must all be positive');
end

% Ensure the signals and eigenvalues align properly
[~, nsignals,] = size(signals);
[neig] = length(eigVals);
if nsignals ~= neig
    error('signals and eigVals do not have matching dimensions');
end
end
