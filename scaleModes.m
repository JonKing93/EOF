function[scaledModes] = scaleModes(modes, eigVals)
%% Scales the EOF modes by the square root of the eigenvalues.
% A necessary step before EOF rotation.
%
% [scaledModes] = scaleModes(modes, eigVals)
%
%
% ----- Inputs -----
%
% modes: A set of EOF modes for Varimax rotation. Each column is a mode
%
% eigVals: A vector containing the eigenvalues associated with each mode.
%
%
% ----- Outputs -----
%
% scaledModes: The scaled EOF modes
%
%
% ----- Written By -----
% 
% Jonathan King, 2017, University of Arizona (jonking93@email.arizona.edu)

% Error check and make eigVals a row vector
[eigVals] = setup(modes, eigVals);

% Ensure compatibility with older versions
eigVals = repmat( eigVals, size(modes,1), 1);

% Scale the modes.
scaledModes = modes .* sqrt(eigVals);

end

%%%%% Helper Functions %%%%%
function[eigVals] = setup(modes, eigVals)

% Check eigVecs is a matrix
if ~ismatrix(modes)
    error('modes must be a matrix');
end

% Check eigVals is a vector
if ~isvector(eigVals)
    error('eigVals must be a vector');
end

% Ensure there are no NaNs
if hasNaN(modes)
    error('modes cannot contain NaN');
end
if hasNaN(eigVals)
    error('eigVals cannot contain NaN');
end

% Check that the dimensions of the eigvecs and eigvals align
[~, nModes] = size(modes);
[nEigvals] = size(eigVals);

if nEigvals ~= nModes
    error('The number of eigenvalues and modes do not align');
end

% Make eigVals a row vector
if ~isrow(eigVals)
    eigVals = eigVals';
end
end
    