function[eigVals, modes, expVar, Datax0, C] = simpleEOF(Data, matrix, varargin)
%% Gets the EOF modes and loadings of a data matrix.
% 
% [eigVals, modes, expVar, Datax0, C] = simpleEOF(Data, matrix)
%
% [...] = simpleEOF(Data, matrix, 'econ')
% Uses an economy sized svd decomposition rather than the full svd.
%
% [...] = simpleEOF(Data, matrix, 'svds', nEigs)
% Uses the svds decomposition and determines the first nEigs eigenvalues.
%
% [...] = simpleEOF(Data, matrix, 'svd')
% Performs the normal svd.
%
%
% ----- Inputs -----
% 
% Data: A 2D data matrix. Each column corresponds to a particular time
%   series. Data cannot contain NaN entries.
%
% matrix: The desired analysis matrix.
%       'cov': Covariance matrix -- Minimizes variance along EOFs
%       'corr': Correlation matrix -- Minimizes relative variance along
%               EOFs. Useful for data series with significantly
%               different magnitudes.
%
%
% ----- Outputs -----
%
% eigVals: The eigenvalues of the analysis matrix. These determine the
%       significance of each mode.
%
% modes: The EOF modes. These are the eigenvectors of the analysis 
%       matrix. Each column is one mode.
%
% expVar: The amount of data variance explained by each mode. Equivalent to
%       the normalized eigenvalues.
%
% Datax0: The standardized or detrended data matrix
%
% C: The analysis matrix for the PCA. This may be the standardized data
%       matrix, its covariance matrix, or its correlation matrix, as appropriate.

% Get the inputs
[svdArgs] = parseInputs(varargin{:});

% Error check
errCheck(Data, matrix);

% Standardize Data
if strcmp(matrix, 'corr')
    Datax0 = zscore(Data);
else    % matrix = 'cov'
    Datax0 = detrend(Data, 'constant');
end

% Get the analysis matrix for decomposition.
C = cov(Datax0);

% Run SVD(S)
[eigVals, modes] = quickSVD(C, svdArgs{:});

% Get the expalined variance = normalized eigenvalues
expVar = eigVals / sum(var(Datax0));

end


%%%%% Helper Functions %%%%%
function[svdArgs] = parseInputs(varargin)
inArgs = varargin;

% Set the default
svdArgs = {'svd'};

if ~isempty( inArgs)
    isSvdsArg = false;
    for k = 1:length(inArgs)
        arg = inArgs{k};
        
        if isSvdsArg
            if isscalar(arg)
                svdArgs = {'svds',arg};
            else
                error('The nEigs parameter must be a positive integer');
            end
        elseif strcmpi(arg, 'econ')
            svdArgs = {'svd', 'econ'};            
        elseif strcmpi(arg, 'svd')
            % Do nothing
        elseif strcmpi(arg, 'svds')
            if length(inArgs) >= k+1
                isSvdsArg = true;
            else
                error('The svds flag must be followed by the nEigs parameter');
            end
        else
            error('Unrecognized Input');
        end
    end
end
end

function[] = errCheck(Data, matrix)
%% Ensure data matrix is 2D    
if ~ismatrix(Data)
    error('Data must be a 2D matrix');
end

% Ensure data does not contain NaNs
if hasNaN(Data)
    error('Data cannot contain NaNs');
end

% Matrix is recognized
if ~any( strcmpi(matrix, {'corr','cov'}) )
    error('Unrecognized matrix');
end
end