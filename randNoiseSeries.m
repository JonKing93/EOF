function[varargout] = randNoiseSeries(data, varargin)
%% Build a matrix of random time series with desired noise properties given
% an initial time series vector or matrix, and scales to the standard deviation of
% the original data.
%
% [randSeries] = randNoiseSeries(ts, ...)
% Constructs an artificial noisy time series based on the properties of a time series ts.
%
% [randSeries] = randNoiseSeries(dataMatrix, ...)
% If Data is an m x n matrix, returns an m x n matrix with n random series. 
% Each of the n columns in randSeries uses the properties of the analogous
% column in "matrix".
%
% [..., noiseArgs] = randNoiseSeries(data, 'red')
% Uses red, AR(1) noise properties and returns the desired properties for
% future calls to randNoiseSeries.
%
% [..., noiseArgs] = randNoiseSeries(data, 'white')
% Uses white, Gaussian noise. Returns the desired properties for future
% calls to randNoiseSeries.
%
% [...] = randNoiseSeries(..., 'noScaling')
% Does not scale the random series to data standard deviation. All random
% series will have mean = 0, an variance = 1.
%
% [randSeries] = randNoiseSeries( noiseArgs )
% Uses output from a past call to generate a new set of noisy time series.
% This method is recommended for Monte Carlo processes to avoid the
% creation of unecessarily large matrices.
%
%
% ----- Inputs -----
%
% ts: A time series. Must be a vector.
%
% dataMatrix: A 2D matrix. Each column is a time series.
%
% noiseArgs: A structure with possible fields...
%   'scale': The scaling matrix
%   'ar1': AR1 properties
%   'size': output size
%   'noise': The noise type
%
%
% ----- Output -----
%
% randSeries: The matrix of randomly generated, noisy time series. Each column
%   of the matrix is a separate time series.
%
% noiseArgs: A structure with possible fields...
%   'scale': The scaling matrix
%   'ar1': AR1 properties
%   'size': output size
%   'noise': The noise type
%
%
% ----- Written By -----
%
% Jonathan King, 2017, University of Arizona, jonking93@email.arizona.edu

% Parse inputs, determine noise properties, error check
[sData, ar1, scale, noise] = setup( data, varargin{:} );

% White noise
if isnan(ar1)
    randSeries = randn( sData );

% Red noise
else    
    % Preallocate
    randSeries = NaN( sData );
    
    % Initialize the first row
    randSeries(1,:) = randn(1,sData(2));
    
    % Calculate autocorrelation through the matrix, add white noise
    for k = 1:sData(1)-1
        randSeries(k+1,:) = ( ar1 .* randSeries(k,:) ) + randn(1, sData(2));
    end
end

% Standardize
randSeries = zscore(randSeries);

% Scaling
if ~isnan(scale)
    randSeries = randSeries .* scale;
end

% Create noise args
if nargout == 1
    varargout = {randSeries};
elseif nargout == 2
    noiseArgs = struct('size',sData,'noise',noise,'ar1',ar1,'scale',scale);
    varargout = {randSeries, noiseArgs};
else
    error('Too many outputs');
end

    
end

%%%%% Helper Function %%%%%
function[sData, ar1, scale, noise] = setup(data, varargin)

    % Set defaults
    scale = NaN;
    ar1 = NaN;

    % If noiseArgs was provided, get the properties
    if nargin == 1
        if isstruct( data )
            if ~isfield(data,'noise') || ~isfield(data,'size') || ~isfield(data,'ar1') || ~isfield(data, 'scale')
                error('noiseArgs structure must contain ''noise'', ''size'', ''ar1'', and ''scale'' fields');
            end
            noise = data.noise;
            sData = data.size;
            ar1 = data.ar1;
            scale = data.scale;
            if ~ismember(noise, {'white','red'})
                error('Unreognized noise type');
            elseif (~isscalar(scale) && ~isequal(size(scale),sData)) || (isscalar(scale) && ~isnan(scale) && ~isequal(size(scale),sData))
                error('The scaling field in noiseArgs must be the same size as the ''size'' field dimensions');
            elseif numel(sData)~=2
                error('The output size in noiseArgs must be 2D');
            elseif (~isscalar(ar1) && (~isvector(ar1) || numel(ar1)~=sData(2))) || (isscalar(ar1) && ~isnan(ar1) && numel(ar1)~=sData(2))
                error('ar1 in noiseArgs must be a vector containing the same number of elements as columns in the output series');
            end
            if ~isrow(ar1)
                ar1 = ar1';
            end
        else
            error('noiseArgs must be a structure');
        end
        
    % Otherwise determine data noise properties from data
    elseif nargin < 4
        if ~ismember(varargin{1}, {'white','red'})
            error('Unrecognized noise type');
        end
        noise = varargin{1};
        if ~ismatrix(data)
            error('data must be a matrix');
        end

        % Get the size of the output series
        sData = size(data);

        % Get red noise properties
        if strcmpi(noise, 'red')
            ar1 = diag( corr( data(1:end-1,:), data(2:end,:) ) )';
        end

        % Get scaling properties
        if nargin == 2
            scale = repmat( std(data), [sData(1), 1]);
        elseif ~strcmpi( varargin{2}, 'noScaling')
            error('Unrecognized input');
        end
    else
        error('Unrecognized Inputs');
    end
end