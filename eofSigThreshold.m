function[varargout] = eofSigThreshold( randExpVar, p, expVar )
%% Gets the explained variance values required for significance at significance level p
%
% [sigExpVar, true_p] = getSigThreshold( randExpVar, p )
% Finds the significance threshold for explained variance values at a
% significance level p.
%
% [sigExpVar, true_p, nSig] = getSigThreshold( randExpVar, p, expVar )
% Also compares the significance threshold to a set of explained variances
% from an eof analysis and determines the number of significant leading
% modes.

% Quick error checking
if ~ismatrix(randExpVar) || any(isnan(randExpVar(:)))
    error('randExpVar must be a matrix and cannot contain NaN.');
elseif ~isfloat(p) || p<=0 || p>=1
    error('p must be a number on the open interval (0,1)');
elseif exist('expVar','var') && ( ~isvector(expVar) || numel(expVar)~=size(randExpVar,2) || any(isnan(expVar(:))) )
    error('expVar must be a vector with the same number of elements as columns in randExpVar. expVar may not contain NaN values.');
end

% Get the number of sets of random explained variance
nSets = size( randExpVar, 1);

% Sort the randomly generated explained variances
randExpVar = sort( randExpVar );

% Get the threshold index for significance <= p
threshold = ceil( nSets * (1-p) );

% Get the true significance of this threshold
true_p = 1 - threshold / nSets;

% Get the explained variances at this threshold
sigExpVar = randExpVar( threshold, : );

% If expVar was given, get the number of significant eof modes
if exist('expVar','var')
    nSig = find( expVar' < sigExpVar, 1, 'first') - 1;
    % Return the output
    varargout = {sigExpVar, true_p, nSig};
else
    varargout = {sigExpVar, true_p};
end
    
end