function[randExpVar] = ruleN(Data, matrix, MC, noise, varargin)
%% Runs a Rule N significance test on a data matrix and its eigenvalues.
%
% [randExpVar] = ruleN(Data, matrix, MC, noise)
% Runs a Rule N significance test using MC iterations. Tests for
% significance <= p against white Gaussian, or red AR(1) noise. Records
% Monte Carlo convergence data.
%
% ruleN(..., 'parallel')
% To speed up runtime, runs the Monte Carlo process in parallel using 
% MATLAB's default settings. If the Distributed Computing Toolbox is not 
% licensed, reverts to serial computing. 
%
% ruleN(..., 'showProgress')
% Displays a progress bar showing the progress through the Monte Carlo 
% iterations. This option is disabled when computing in parallel.
%
% ruleN(..., 'estimateRuntime')
% Estimates total runtime based on the number of iterations and available
% workers in the computing pool.
%
% ruleN(..., 'pcaArgs', {'param1',val1,'param2',val2...})
% Uses alternative arguments for the pca function. See the MATLAB
% documentation on "pca" for more details.
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
%       'raw': Perform svd directly on data matrix.
%
% MC: The number of Monte Carlo iterations to perform
%
% noise: A flag for the noise type
%   'white':    white Gaussian noise (mean = 0, variance = 1)
%   'red':      AR(1) red noise with added Gaussian white noise.
%
%            
% ----- Outputs -----
%
% randExpVar: The set of random explained variances.
%
%
% ----- Written By -----
% 
% Jonathan King, 2017, University of Arizona (jonking93@email.arizona.edu)

% Parse inputs and error checking
[parallel, showProgress, estimateRuntime, pcaArgs] = setup( Data, matrix, MC, noise, varargin{:} );

% Preallocate Eigenvalue array
[~, nCols] = size(Data);
randExpVar = NaN(MC,nCols);

% Initialize the parallel pool if computing in parallel
nWorkers = 1;
if parallel
    pool = gcp;
    nWorkers = pool.NumWorkers;
end

% Get the noise Arguments
[~, noiseArgs] = randNoiseSeries( Data, noise );

% Initialize the progress bar if displaying
if showProgress
    progressbar(0);
end

% Run the first iteration. Estimate runtime if desired.
if estimateRuntime
    startTime = tic;
end
randExpVar(1,:) = ruleNStep(noiseArgs, matrix, pcaArgs);
if estimateRuntime
    time = toc(startTime);
    time = time*MC/nWorkers;
    h = time / 360;
    m = mod(time,360)/60;
    s = mod(time,60);
    fprintf('Estimated runtime: %0.f hour(s), %0.f minute(s), %0.f seconds \r\n',h,m,s);
end

% Run Rule N...
if parallel          % In parallel
    parfor k = 2:MC
        randExpVar(k,:) = ruleNStep(noiseArgs, matrix, pcaArgs);
    end
    
else                 % In serial
    for k = 2:MC
        randExpVar(k,:) = ruleNStep(noiseArgs, matrix, pcaArgs);
        
        % Update the progress bar if displaying
        if showProgress
            progressbar(k/MC);
        end
    end
end

end

%%%%% Helper Functions %%%%%
function[parallel, showProgress, estimateRuntime, pcaArgs] = setup( Data, matrix, MC, noise, varargin )
% This does some initial error checking and setup for Rule N

% Parse the inputs
[parallel, showProgress, estimateRuntime, pcaArgs] = parseInputs( varargin,...
    {'parallel','showProgress','estimateRuntime','pcaArgs'},... % Flags
    {false,false,false,{}}, {'b','b','b',{}} );      % Defaults and switches

% Error checking
if ~ismatrix(Data)
    error('Data must be a 2D matrix');
elseif ~any( strcmpi( matrix, {'corr','cov','raw'} ) )
    error('Unrecognized matrix');
elseif ~any( strcmpi( noise, {'red','white'}))
    error('Unrecognized noise type');
elseif ~isfloat(MC) || MC < 1 || mod(MC,1)~=0
    error('MC must be a positive integer.');
elseif any( strcmpi( 'Centered', pcaArgs ) )
    error('The ''Centered'' option for pcaArgs is forbidden. RuleN implements this separately.');
end

% If analyzing raw matrix, remove centering from pca
if strcmpi(matrix, 'raw')
    pcaArgs = [pcaArgs, 'Centered', false];
end

% If running in parallel, check for the Distributed Computing Toolbox
if parallel
    if ~license('test', 'Distrib_Computing_Toolbox')
        warning( sprintf('The Distributed Computing Toolbox is not licensed on this computer.\r\nCannot run in parallel. Reverting to default serial mode...')); %#ok<SPWRN>
        parallel = false;
    end
end

% Disable 'showProgress' for parallel computing
if parallel && showProgress
    warning('Cannot display progress for parallel computations.');
    showProgress = false;
end

end

function[randExpVar] = ruleNStep(noiseArgs, matrix, pcaArgs)
% This runs a single iteration of the Monte Carlo Rule N process

% Create a random matrix with the desired noise properties, scaled to data standard deviation.
g = randNoiseSeries(noiseArgs);
    
% Get the analysis matrix
if strcmp(matrix, 'corr')
    g = zscore(g);
elseif strcmp(matrix, 'cov')
    g = detrend(g, 'constant');
end
    
% Run a pca (eof analysis) on the random matrix and get explained variance
[~,~,~,~, randExpVar] = pca(g, pcaArgs{:});

end   