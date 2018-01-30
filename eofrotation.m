function[rotModes, rotEigvals, rotExpVar, rotSignals, scaledRotSignals] = eofrotation(modes, eigvals, A, rotType, expVar)
%% Performs the EOF rotation and calculates associated values.
%
% [rotModes, rotEigvals, rotExpVar, rotSignals, scaledRotSignals] = eofrotation(modes, eigvals, A, rotType)
%
%
% ----- Inputs -----
%
% modes: The modes from an eof analysis that should be rotated.
%
% eigvals: The eigenvalues corresponding to the rotating modes.
%
% A: The analysis matrix used to generate the rotating modes.
%
% rotType: A string flag for the type of rotation
%   'varimax': Varimax Rotation
%   'equamax': Equamax Rotation
%
%
% ----- Outputs -----
%
% rotModes: The rotated eof modes.
%
% rotEigvals: The eigenvalues for the rotated modes.
%
% rotExpVar: The variance explained by the rotated modes.
%
% rotSignals: The signals corresponding to the rotated modes.
%
% scaledRotSignals: The signals corresponding to the rotated modes scaled
%     to the analysis matrix.
%
%
% ----- Author -----
%
% Jonathan King, 2017, University of Arizona, jonking93@email.arizona.edu

% Do the rotation
[rotModes, rotMatrix] = rotatefactors( modes, 'method', rotType );

% Get the rotated eigenvalues
rotEigvals = diag( diag(eigvals) * rotMatrix );

% Get the explained variance of each rotated mode
rotExpVar = (rotEigvals / sum(rotEigvals)) * sum(expVar);

% Get the rotated signals
rotSignals = A * rotModes;

% Get the scaled, rotated signals
scaledRotSignals = scaleSignals( rotSignals, rotEigvals );

end