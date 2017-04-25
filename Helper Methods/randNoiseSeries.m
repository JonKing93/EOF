function[randSeries] = randNoiseSeries(noiseType, Data, varargin)
%% Build a matrix of random time series with desired noise properties given
% an initial time series or matrix, and scales to the standard deviation of
% the original data.
%
% [randSeries] = randNoiseSeries(ts, nSeries, noiseType)
% If Data is a time series, constructs a matrix with nSeries artificial
% time series with the desired noise properties.
%
% [randSeries] = randNoiseSeries(matrix, noiseType)
% If Data is an m x n matrix, returns a matrix with n random series. Each
% of the n series contains the desired noise properties of the analogous
% series in Data.
%
% [...] = randNoiseSeries(..., 'noScaling')
% Does not scale the random series to the data standard deviation. All
% random series will have mean = 0, and variance = 1.
%
% ----- Inputs -----
%
% noiseType: The type of noise to introduce into the random time series
%   'white': white Gaussian noise
%   'red': lag-1 autocorrelated red noise
%
% Data: A time series or matrix of time series on which to base the random
%       series. If Data is a matrix, each column contains a time series.
%
% nSeries: If data is a time series, the number of random time series to generate
%
%
% ----- Output -----
%
% randSeries: The matrix of randomly generated, noisy time series. Each column
%   of the matrix is a separate time series.

[nSeries, scaling] = parseInputs(Data, varargin{:});

% If data is a row vector, make into a column vector
if isrow(Data)
    Data = Data';
end

% Get the length of the series
[lts] = size(Data, 1);

% White noise
if strcmpi(noiseType, 'white')
   randSeries = randn(lts, nSeries);
    
% Red noise
elseif strcmpi(noiseType, 'red')
    
    % Preallocate
    randSeries = NaN(lts, nSeries);
    
    % Get the lag-1 autocorrelations
    ar1 = diag(  corr( Data(1:end-1,:), Data(2:end,:) )  )';
    
    % Initialize the first row
    randSeries(1,:) = randn(1,nSeries);
    
    % Calculate autocorrelation through the matrix, add white noise
    for k = 1:lts-1
        randSeries(k+1,:) = ( ar1 .* randSeries(k,:) );
        randSeries(k+1,:) = randSeries(k+1,:) + randn(1, nSeries);
    end
    
% Othr noise type
else
    error('Unrecognized noiseType')
end

% Standardize the random series so that scaling is correct.
randSeries = zscore(randSeries);

% Scale to data standard deviation if desired
if scaling
    randSeries = randSeries .* std(Data);
end    
    
end

% ----- Helper Functions -----
function[nSeries, scaling] = parseInputs(Data, varargin)
inArgs = varargin;

% Set defaults
scaling = true;

% Check dataType
if isvector(Data)
    dataType = 'vector';
    
    % Ensure the number of series is given
    if ~isempty(inArgs) && isscalar(inArgs{1}) && (inArgs{1} > 0)
        nSeries = inArgs{1};
    else
        error('When Data is a vector, a positive integer must specify the number of random series to generate');
    end
elseif ismatrix(Data)
    dataType = 'matrix';
    nSeries = size(Data,2);    
else
    error('Data must be a vector or a matrix');
end

if ~isempty(inArgs)
    % Get the inputs
    for k = 1:length(inArgs)
        arg = inArgs(k);

        if strcmpi(dataType,'vector') && k==1
            % Do nothing, we already assigned nSeries
        elseif strcmpi(arg, 'noScaling')
            scaling = false;
        else
            error('Unrecognized Input');
        end
    end
end
end