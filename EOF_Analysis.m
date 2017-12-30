function[s] = EOF_Analysis(Data, varargin)
%% Performs an EOF Analysis (also known as PC Analysis) with data normalization, significance testing, rotation, and plotting.
%
% [s] = EOF_Analysis(Data)
% Normalizes data and performs an EOF analysis. Applies a 1000-iteration
% Monte Carlo "Rule N" procedure to EOF modes (also known as loadings) to 
% test for significance (p<0.05) against an AR(1) red-noise background. 
% Records convergence information for the Monte Carlo procedure. Rotates 
% significant modes according to the varimax criterion. Returns a structure
% 's' with analysis output and metadata. Plots analysis data.
%
% EOF_Analysis(..., 'noplot')
% Supresses the output plots.
%
% EOF_Analysis(..., 'varNames', {'Name1','Name2',...})
% Includes variable names in plots.
% 
% ANALYSIS OPTIONS:
% Need something other than the default analysis? Try these options...
% 
% EOF_Analysis(..., 'matrix', matrix)
% Specifies whether to perform the analysis on normalized, centered, or raw
% data. By default, EOF_Analysis normalizes data. Normalizing data causes
% each mode to minimize RELATIVE covariance (i.e. correlation). It is most
% useful for examining shared correlations or for data vectors with very
% different magnitudes. Centering data causes each mode to minimize TOTAL
% covariance between data vectors. Analyzing raw data minimizes total
% variance in the original dataset.
%
% EOF_Analysis(..., 'MC', MC)
% Specifies the number of Monte Carlo iterations to use during significance
% testing. A greater number of iterations increases statistical robustness
% but also raises runtime.
%
% EOF_Analysis(..., 'noiseType', noise)
% Specifies whether to perform the significance test against a white or a
% AR(1) red-noise process. By default, EOF_Analysis uses a red significance
% test. Use a white significance test (Gaussian noise, variance=1, mean=0)
% for processes with minimal autocorrelation. Use a red test for processes
% with high autocorrelation or long memory. A red test is generally more
% stringent and applicable to many geophysical / geological processes.
%
% EOF_Analysis(..., 'p', p)
% Specify the significance level cutoff to use for significance testing. By
% default, EOF_Analysis uses p<0.05 to assign significance.
%
% EOF_Analysis(..., 'noSigTest')
% Blocks the Monte Carlo Rule-N significance test.
%
% EOF(..., 'pcaArgs', {'param1',val1,'param2',val2...})
% Implement alternative arguments for the pca function (e.g. alternative
% algorithms, treatment of NaN values, analysis of only a few modes, weighted
% data vectors, etc.). See the MATLAB documentation on "pca" for more details. 
% EOF_Analysis handles data normalization and centering separately from the
% "pca" function, so the 'Centering' argument has been disabled.
%
% EOF(..., 'nRotate', nRot)
% Specifies the number of EOF modes to rotate. This will rotate the leading
% modes.
%
% EOF(..., 'rotType', rotType)
% Specify the type of rotation to use (Varimax or Equamax). By default,
% EOF_Analysis uses Varimax rotation.
%
% RUNTIME OPTIONS:
% Huge dataset? Monte Carlo processes taking forever? Try some of these
% settings...
%
% [s] = EOF_Analysis(..., 'estimateRuntime')
% Displays an estimate of total runtime for the Monte Carlo process.
%
% [s] = EOF_Analysis(..., 'showProgress')
% Displays a progress bar with percent completion and estimated remaining
% time for Rule-N. Not available for parallel computations.
%
% [s] = EOF_Analysis(..., 'parallel')
% If possible, runs the Monte Carlo significance tests in parallel. This
% can significantly speed up runtime for calculations longer than several
% minutes. If the "Distributed Computing Toolbox" is not licensed, notifies
% the user and proceeds in serial.
%
% [s] = EOF_Analysis(..., 'noConvergeTest')
% A flag to block the test for Monte Carlo convergence. This may modestly
% improve runtime for very large significance tests.
%
%
% ----- Inputs -----
%
% Data: A 2D data matrix. Each column corresponds to a particular data
%       series. Data may only contain numeric entries.
%
% matrix: A flag for the desired analysis matrix.
%       'corr' (Default): Correlation matrix. Minimize RELATIVE variance on
%                         EOF modes. Best for data series with significantly
%                         different magnitudes or where correlation is of interest.
%       'cov': Covariance matrix. Minimizes TOTAL variance along EOF modes.
%              Recommended for analyses of aggregate changes to a dataset.
%       'raw': Perform analysis directly on data matrix.
%  
% MC: The number of Monte Carlo iterations used in the Rule N significance
%       test. Must be a positive integer.
%
% noise: A flag for the noise used in the Rule N significance test
%       'red' (Default): AR(1) noise with added white noise
%       'white': white Gaussian noise. (mean = 0, variance = 1)
%
% p: The significance level that the significance test should pass. Must be
%    a positive number on the interval (0, 1).
%
% nRotate: The number of modes to rotate. Must be between 0 and the number
%          of data variables
%
% rotType: The type of rotation to use
%       'varimax' (Default): Varimax Rotation
%       'equamax': Equamax Rotation
%
%
% ----- Outputs -----
%
% s: A structure containing some or all of the following fields
%
%   varNames: The names of the data variables.
%
%   A: The analysis matrix. This is the standardized, detrended or raw data matrix.
%
%   modes: (Also known as loadings). Each column contains the coefficients
%          for one EOF signal. Each mode is an eigenvector of the analysis 
%          matrix. They record how strongly each data vector is associated
%          with each EOF signal.
%
%   signals: The signal for each EOF mode. Each column contains a signal.Signals are the imprint of each
%       mode on the original data series, also known as scores or EOF 
%       time series. Each column is one signal.
%
%   eigVals: A vector with the eigenvalues of the analysis matrix sorted in
%            descending order. Each eigenvalue corresponds to the strength
%            of the associated EOF mode.
%
%   expVar: The variance explained by each mode. Equivalent to the
%           normalized eigenvalues.
%
%   scaledSignals: The signals scaled to the analysis data matrix. Allows
%                  for direct comparison of signals with the original
%                  standardized/centered/raw data series.
%
%   randEigvals: The set of random, normalized eigenvalues generated during
%       the Rule N significance test. Each row contains the eigenvalues
%       at a particular confidence interval.
%
%   p: The significance level used for the Rule N test.
%
%   sigExpVar: The explained variance threshold for significance. The
%       explained variance of EOF modes must be greater than or equal to 
%       this threshold to be significant.
%
%   true_p: The actual significance level of the Rule N test. (Depending on
%   the number of Monte Carlo iterations, the significance level of the
%   Rule N test may be slightly higher than the user-specified value.
%
%   nSig: The number of leading eof modes that pass the Rule N significance test.
%
%   MCsigExpVar: The explained variance significance threshold at each
%       successive Monte Carlo iteration.
%
%   MCtrue_p: The true significance level being tested at each successive
%       Monte Carlo iteration.
%
%   rotModes: The rotated eof modes.
%
%   rotEigvals: The eigenvalues for the rotated modes.
%
%   rotExpVar: The variance explained by the rotated modes.
%
%   rotSignals: The signals corresponding to the rotated modes.
%
%   scaledRotSignals: The signals corresponding to the rotated modes scaled
%       to the analysis matrix.
%
%   metadata: Information to support replication of the analysis.
%
%
% ----- Written By -----
% 
% Jonathan King, 2017, University of Arizona (jonking93@email.arizona.edu)

% Parse inputs, error checking, setup for the analysis
[matrix, pcaArgs, sigTest, ruleNArgs, p, convergeTest, nRotate, rotType, plotting, varNames] = setup( Data, varargin{:} );% Setup the analysis, do some error checking

% Declare the intial structure
s = struct();
if ~isempty(varNames)
    s.varNames = varNames;
end

% Get the analysis matrix
if strcmp(matrix, 'corr')
    s.A = zscore(Data);
elseif strcmp(matrix, 'cov')
    s.A = detrend(Data, 'constant');
else
    s.A = Data;
    pcaArgs = [pcaArgs, 'Centered', false];
end

% Run the initial EOF on the Data
[s.modes, s.signals, s.eigvals, ~, s.expVar] = pca(s.A, pcaArgs{:});

% Scale the signals to the standardized data
s.scaledSignals = scaleSignals( s.signals, s.eigvals );

% If testing significance...
if sigTest
    % Do the Rule N generating process
    s.randExpVar = ruleN( Data, matrix, ruleNArgs{:});
    
    % Record the significance level used
    s.p = p;
    
    % Test for significance at the desired significance level
    [s.sigExpVar, s.true_p, s.nSig] = eofSigThreshold(s.randExpVar, p, s.expVar);
    
    % Record Monte Carlo convergence data.
    if convergeTest
        [s.MCsigExpVar, s.MCtrue_p] = ruleNConvergence( s.randExpVar, s.p );
    end
end

% Determine how many modes to rotate if the user did not override the rotation number
if isnan(nRotate) && isfield(s, 'nSig')
    nRotate = s.nSig;
end
    
% Rotation must be performed on 2+ modes or there is no change...
if ~isnan(nRotate) && nRotate > 1
    % Perform the rotation
    [s.rotModes, s.rotEigvals, s.rotExpVar, s.rotSignals] = ...
        eofrotation( s.modes(:,1:nRotate), s.eigvals(1:nRotate), s.A, rotType );
end
    
% Add metadata for the analysis
s.metadata = [{'matrix';'MC';'noiseType';'p';'pcaArgs';'Rotated Modes';'Rotation Type'},...
    {matrix; ruleNArgs{1}; ruleNArgs{2}; s.p; pcaArgs; nRotate; rotType}];

% Plot the output if desired
if plotting
    eofconvergence(s);
    eofsignificance(s);
    eofloadings(s);
end

end

%%%%% Helper Functions %%%%%
function[matType, pcaArgs, sigTest, ruleNArgs, p,convergeTest, nRotate, rotType, plotting, varNames] = setup( Data, varargin )

% Parse the inputs
[plotting, matType, MC,  noise,      p,  sigTest,  pcaArgs,  nRotate,  rotType,  estimateRuntime,  showProgress,  parallel,   convergeTest,     varNames] = parseInputs(varargin,...
{'noplot','matrix','MC','noiseType','p','noSigTest','pcaArgs','nRotate','rotType','estimateRuntime','showProgress','parallel','noConvergeTest','varNames'},...   % The string flags
{ true,  'corr',  1000,  'red',   0.05,  true,       {},       NaN,   'varimax',     false,           false,       false,       true,            {}},...       % The default values
{'b',{'corr','cov','raw'},{},{'red','white'},{},'b',   {},      {},{'varimax','equamax'},'b',            'b',         'b',       'b',           {}} );   % Switches

% Ensure Data is a matrix
if ~ismatrix(Data)
    error('Data must be a 2D matrix.');
end

% Disable the 'Centered' option for pca. (Centering is implemented
% separately by EOF_Analysis
if any( strcmpi( 'Centered', pcaArgs ) )
    error('The ''Centered'' option for pcaArgs is forbidden. EOF_Analysis implements this separately.');
end

% Ensure that the number of rotated modes is valid if specified by the user
if ~isnan(nRotate)
    nCols = size(Data,2);
    if ~isfloat(nRotate) || nRotate<0 || nRotate>nCols || mod(nRotate,1)~=0
        error('The number of rotated modes must be an integer on the closed interval from 0 to the number of data vectors.');
    end
end

% Set the rule N argument list
ruleNArgs = {MC, noise, 'pcaArgs', pcaArgs};
if estimateRuntime
    ruleNArgs = [ruleNArgs, 'estimateRuntime'];
end
if showProgress
    ruleNArgs = [ruleNArgs, 'showProgress'];
end
if parallel
    ruleNArgs = [ruleNArgs, 'parallel'];
end

end