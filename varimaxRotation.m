function[rotModes, varargout] = varimaxRotation(scaModes, varargin)
%% Performs a VARIMAX rotation on a set of scaled eigenvectors and eigenvalues.
%
% [rotModes, rotMatrix] = varimaxRotation(scaModes) returns a set of
%   eigenvectors rotated using the varimax criterion and the corresponding
%   rotation matrix
%
% [rotModes, rotEigvals, rotMatrix] = varimaxRotation(scaModes, eigVals)
%   also returns the rotated eigenvalues.
%
% [rotModes, rotEigvals, rotExpVar, rotMatrix] = varimaxRotation(scaModes, eigvals, expVar)
%   also returns the explained variance of each rotated mode.
%
%
% ----- Inputs -----
%
% scaModes: A set of modes scaled for varimax rotation by the square root
% of associated eigenvalues. Each column contains 1 mode.
%
% eigVals: The eigenvalues associated with each mode to be rotated.
%
% expVar: The explained variance of each mode to be rotated.
%
%
% ----- Outputs -----
%
% rotEigvecs: The rotated eigenvectors.
%
% rotEigvals: The rotated eigenvalues.
%
% rotExpVar: The explained variance of the rotated eigenvalues.
%
% rotMatrix: The matrix used to perform the eigenvector rotation.
%
%
% ----- Written By -----
% 
% Jonathan King, 2017, University of Arizona (jonking93@email.arizona.edu)

% Parse the inputs, do some error checking
[haveEigvals, haveExpVar, eigVals, expVar] = setup(scaModes, varargin{:});

% Perform the rotation
[rotModes, rotMatrix] = rotatefactors(scaModes);

% Calculate the rotated eigenvalues
if haveEigvals
    rotEigvals = diag(  diag(eigVals) * rotMatrix  );
    
    % Calculate rotated explained variance
    if haveExpVar
        rotExpVar = (rotEigvals / sum(rotEigvals)) * sum(expVar);
        varargout = {rotEigvals, rotExpVar, rotMatrix};
    else
        varargout = {rotEigvals, rotMatrix};
    end
else
    varargout = {rotMatrix};
end

end

%%%%% Helper Functions %%%%%
function[haveEigvals, haveExpVar, eigVals, expVar] = setup(scaModes, varargin)

% Set defaults
haveEigvals = false;
haveExpVar = false;
eigVals = [];
expVar=  [];

% Check for additional inputs
if length(varargin) == 1
    haveEigvals = true;
    eigVals = varargin{1};
elseif length(varargin) == 2
    haveEigvals = true;
    haveExpVar = true;    
    eigVals = varargin{1};
    expVar = varargin{2};
elseif length(varargin) > 2
    error('Too many inputs');
end

% Check that inputs are correctly formatted
if ~ismatrix(scaModes) || hasNaN(scaModes)
    error('scaModes must be a matrix and cannot contain NaN.');
end
if haveEigvals && (~isvector(eigVals) || hasNaN(eigVals) || any(eigVals<0)) 
    error('eigVals must be a vector of positive numbers and cannot contain NaN');
end
if haveExpVar && (~isvector(expVar) || hasNaN(expVar) || any(expVar<0))
    error('expVar must be a vector of positive numbers and cannot contain NaN');
end

% Ensure sizes align
nModes = size(scaModes, 2);
if haveEigvals
    nEigs = length(eigVals);
    if nEigs ~= nModes
        error('scaModes and eigVals must have the same number of modes/eigVals.');
    end
    
    if haveExpVar
        nExps = length(expVar);
        if nEigs ~= nExps
            errors('The number of explained variances must equal the number of eigenvalues');
        end
    end
end
end