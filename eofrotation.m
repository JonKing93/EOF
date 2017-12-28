function[rotModes, rotEigvals, rotExpVar, rotSignals, scaledRotSignals] = eofrotation(modes, eigvals, A, rotType)
%% Performs the EOF rotation and calculates associated values.

% Do the rotation
[rotModes, rotMatrix] = rotatefactors( modes, 'method', rotType );

% Get the rotated eigenvalues
rotEigvals = diag( diag(eigvals) * rotMatrix );

% Get the explained variance of each rotated mode
rotExpVar = (rotEigvals / sum(rotEigvals)) * sum(eigvals);

% Get the rotated signals
rotSignals = A * rotModes;

% Get the scaled, rotated signals
scaledRotSignals = scaleSignals( rotSignals, rotEigvals );

end