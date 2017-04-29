function[nSig, randEigSort, thresh, trueConf, varargout] = ...
    ruleN(Data, matrix, normEigvals, MC, noiseType, pval, varargin)
%% Runs a Rule N significance test on a data matrix and its eigenvalues.
%
% [nSig, randEigSort, thresh, trueConf, iterTrueConf, iterSigEigs] = ...
%    ruleN(Data, matrix, eigVals, MC, noiseType, pval)
% Runs a Rule N significance test on a dataset and saves Monte Carlo
% convergence data.
%
% [...] = ruleN(..., showProgress)
% Choose whether to display the current Monte Carlo iteration against the total number of
% simulations.
%
% [nSig, randEigSort, thresh, trueConf] = ruleN(..., convergeTest)
% Choose whether to include or block the recording of the Monte Carlo 
% iteration convergence. Blocking may speed runtime for large analyses, but
% causes a loss of information.
%
% [...] = ruleN(..., 'svds', 'econ')
% Uses the economy sized svds decomposition during Rule N.
%
% [...] = ruleN(..., 'svd', nModes')
% Uses the svds decomposition to get eigenvalues for the first nModes
% modes.
%
% 
% ----- Inputs -----
%
% Data: A 2D data matrix. Each column corresponds to a series of
%   observations.
%
% matrix: The desired analysis matrix.
%       'cov': Covariance matrix -- Minimizes variance along EOFs
%       'corr': Correlation matrix -- Minimizes relative variance along
%               EOFs. Often useful for data series with significantly
%               different magnitudes.
%       'none': Perform svd directly on data matrix.
%
% normEigvals: The normalized eigenvalues of the analysis matrix. Generally
%       equivalent to the explained variance of EOF modes.
%
% MC: The number of Monte Carlo iterations to perform
%
% noiseType: 
%   'white':    white noise
%   'red':      lag-1 autocorrelated red noise with Gaussian white noise.
%
% pval: The significance level desired for the test to pass. Must be on the
%       interval (0 1).
%
% showProgress: A flag for displaying the current Monte Carlo iteration number
%       'showProgress' -- Displays the current Monte Carlo iteration number
%       'noProgress' (Default) -- Does not display the current Monte Carlo number
%
% convergeTest: A flag to save or ignore Monte Carlo convergence information
%       'testConverge' (Default) -- Saves the Monte Carlo convergence data
%       'noConvergeTest' -- Does not save Monte Carlo convergence data
%
%            
% ----- Outputs -----
%
% lastSigNum: The number of eigenvalues that pass rule N
%
% randEigSort: The matrix of random, normalized, sorted eigenvalues
% 
% thresh: The index of the threshold row in randEigSort of which data
%       eigenvalues must exceed to remain significant.
%
% trueConf: The true confidence level of this threshold
%
% iterSigEigs: The set of random eigenvalues that data eigenvalues must
%       exceed to remain significant after each additional Monte Carlo iteration.
%
% iterTrueConf: The true significance level of the threshold eigenvalues
%       after each additional Monte Carlo iteration.
%
%
% ----- Written By -----
% 
% Jonathan King, 2017, University of Arizona (jonking93@email.arizona.edu)

% Inputs and error checking
[showProgress, testConverge, svdArgs] = parseInputs(varargin{:});
errCheck(Data, normEigvals, MC, pval)

% Preallocate output
[~, n] = size(Data);
randEigvals = NaN(MC,n);
if testConverge
    iterSigEigs = NaN(MC, n);
    iterTrueConf = NaN(MC, 1);
else
    iterSigEigs = [];
    iterTrueConf = [];
end

% Run Rule N...
for k = 1:MC
    
    % Display progress if desired
    if showProgress
        fprintf('Running Monte Carlo simulation: %i / %i\r\n', k, MC);
    end
    
    % Create a random matrix with the desired noise properties, scaled to data standard deviation.
    g = randNoiseSeries(noiseType, Data);
    
    % Run an EOF analysis on the random matrix
    [randEig, ~] = simpleEOF(g, matrix, svdArgs);
    
    % Normalize the eigenvalues
    randEig = randEig ./ sum(randEig);
    
    % Store the random eigenvalues
    randEigvals(k,:) = randEig;
    
    % If testing Monte Carlo convergence...
    if testConverge
        % Sort the current set of random eigenvalues
        randEigvals = sort(randEigvals);
        
        % Calculate the current confidence level threshold
        thresh = ceil(k * (1-pval));
        iterTrueConf(k) = thresh / k;
        
        % Get the set of values on the confidence interval
        iterSigEigs(k,:) = randEigvals(thresh,:);
    end
    
end

% Sort the eigenvalues when there is no convergence test
if ~testConverge
    randEigSort = sort(randEigvals);
else
    randEigSort = randEigvals;
    varargout = cell(2,1);
    varargout{1} = iterSigEigs;
    varargout{2} = iterTrueConf;
end


% Calculate the confidence level threshold and its true confidence level
thresh = ceil( MC * (1-pval) );
trueConf = thresh / MC;

% Find the significant values
for k = 1:n
    if normEigvals(k) <= randEigSort(thresh, k)
        nSig = k-1;
        break;
    end
end

end

%%%%% Helper Functions %%%%%
function[showProgress, testConvergence, svdArgs] = parseInputs(varargin)
inArgs = varargin;

% Set defaults
showProgress = false;
testConvergence = true;
svdArgs = {'svd'};

% Get input values
if ~isempty(inArgs)
    isSvdsArg = false;
    
    % Get each input
    for k = 1:length(inArgs)
        arg = inArgs{k};
        
        
        if isSvdsArg
            if isscalar(arg) || strcmpi(arg,'econ')
                svdArgs = {'svds', arg};
            else
                error('The svds flag must be followed by nEigs or the ''econ'' flag');
            end
        elseif strcmpi(arg, 'showProgress')
            showProgress = true;
        elseif strcmpi(arg, 'blockProgress')
            % Do nothing
        elseif strcmpi(arg, 'noConvergeTest')
            testConvergence = false;
        elseif strcmpi(arg, 'testConverge')
            % Do nothing
        elseif strcmpi(arg, 'svd')
            % Do nothing
        elseif strcmpi(arg, 'svds')
            if length(inArgs) >= k+1
                isSvdsArg = true;
            else
                error('The svds flag must be followed by nEigs or the ''econ'' flag');
            end
        else
            error('Unrecognized Input');
        end
    end
end
end

function[] = errCheck(Data, eigVals, MC, pval)
if hasNaN(Data)
    error('Data cannot contain NaN');
elseif hasNaN(eigVals)
    error('eigVals cannot contain NaN');
elseif MC < 1
    error('The Monte Carlo number must be a positive integer');
elseif pval<=0 || pval>=1
    error('The p values must be on the interval (0,1)');
end
end
    