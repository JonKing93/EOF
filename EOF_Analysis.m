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
% ANALYSIS OPTIONS:
% Need something other than the default analysis? Try these options...
% 
% EOF_Analysis(..., 'matrix', matrix)
% Specifies whether to perform the analysis on normalized, centered, or raw
% data. By default, EOF_Analysis normalizes data. Normalizing data causes
% each mode to minimize RELATIVE covariance. It is most useful for
% examining shared correlations or for data vectors with very different
% magnitudes. Centering data causes each mode to minimize TOTAL variance
% between centered data vectors. Analyzing raw data minimizes total
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
% test. Use a white significance test (random noise, variance=1, mean=0)
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
% EOF(..., 'equamax')
% Uses an equamax rotation instead of the default varimax rotation.
%
% RUNTIME OPTIONS:
% Huge dataset? Monte Carlo processes taking forever? Try some of these
% settings...
%
% [s] = EOF_Analysis(..., 'estimateRuntime')
% Displays an estimate of total runtime for the Monte Carlo process.
%
% [s] = EOF_Analysis(..., 'showProgress')
% Displays the percent completion of the Monte Carlo process.
%
% [s] = EOF_Analysis(..., 'parallel')
% If possible, runs the Monte Carlo significance tests in parallel. This
% can significantly speed up runtime for calculations longer than several
% minutes. If the "Distributed Computing Toolbox" is not licensed, notifies
% the user and proceeds in serial.
%
% [s] = EOF_Analysis(..., 'noConvergeTest')
% A flag to block the test for Monte Carlo convergence. This may modestly
% improve runtime, but all convergence data will be lost.
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
%       'white': white Gaussian noise. (mean = 0, variance = 1)
%       'red': lag-1 autocorrelated noise with added white noise
%
% p: The significance level that the significance test should pass. Must be
%    a positive number on the interval (0, 1).
%
%
% ----- Outputs -----
%
% s: A structure containing the following fields
%
%   A: The analysis matrix. This is the standardized, detrended or raw data matrix.
%
%   modes: (Also known as loadings). Each column contains the coefficients
%          for one EOF signal. Each mode is an eigenvector of the analysis 
%          matrix. They record how strongly each data vector is associated
%          with each EOF signal.
%
%   eigVals: A vector with the eigenvalues of the analysis matrix sorted in
%            descending order. Each eigenvalue corresponds to the strength
%            of the associated EOF mode.
%
%   expVar: The variance explained by each mode. Equivalent to the
%           normalized eigenvalues.
%
%   signals: The signal for each EOF mode. Each column contains a signal.Signals are the imprint of each
%       mode on the original data series, also known as scores or EOF 
%       time series. Each column is one signal.
%
%   scaledSignals: The signals scaled to the analysis data matrix. Allows
%                  for direct comparison of signals with the
%                  standardized/centered/raw data series.
%
%   nSig: The number of modes that pass the Rule N significance test. 
%
%   randEigvals: The set of random, normalized eigenvalues generated during
%       the Rule N significance test. Each row contains the eigenvalues
%       at a particular confidence interval.
%
%   thresh: The index of the threshold row in randEigvals that the data 
%       eigenvalues must exceed in order to pass the significance test.
%
%   conf: The confidence level for the significance tests. 
%
%   trueConf: The true confidence level of the threshold row.
%
%   iterTrueConf: The true confidence level of the threshold row after
%       each iteration of the Monte Carlo simulations.
%
%   iterSigEigs: The set of eigenvalues that the data values must exceed
%       for significance after each successive Monte Carlo iteration.
%
%   scaModes: The scaled modes used for VARIMAX rotation. Modes are scaled 
%       by the square root of the loadings.
%
%   rotModes: The VARIMAX rotated modes.
%
%   rotEigvals: The eigenvalues for the rotated modes.
%
%   rotExpVar: The variance explained by the rotated loadings.
%
%   rotSignals: The signals corresponding to the rotated modes.
%
%   scaRotSignals: The scaled signal for each rotated mode.
%
%   rotMatrix: The rotation matrix used to rotate the significant modes.
%
%   metadata: Information concerning the settings used for the analysis.
%       Contains: matrix, MC, noisetype, pval, and any additional flags.
%
%
% ----- Written By -----
% 
% Jonathan King, 2017, University of Arizona (jonking93@email.arizona.edu)

% Parse Inputs, all error checking will occur in called functions
% Need to write

% Declare the intial structure
s = struct();

% Get the analysis matrix
if strcmp(matType, 'corr')
    s.A = zscore(Data);
elseif strcmp(matType, 'cov')
    s.A = detrend(Data, 'constant');
end

% Run the initial EOF on the Data
[s.modes, s.signals, s.eigvals, ~, s.expVar] = pca(s.A, pcaArgs{:});

% Scale the signals to the standardized data
s.scaledSignals = scaleSignals( s.signals, s.eigvals );

% If testing significance...
if sigTest
    % ... Run Rule N
%     ruleN( Data, matrix, s.expVar, MC, noisetype, pval, showProgress, guessRuntime)
    % Need to write.
end

% Determine how many modes to rotate if the user did not override the
% rotation number
if isnan(nRotate)     % User did not override rotation number
    nRotate = s.nSig;
end
    
% Rotation must be performed on 2+ modes or there is no change...
if nRotate > 1
    % Perform the rotation
    [s.rotModes, s.rotEigvals, s.rotExpVar, s.rotSignals] = ...
        eofrotation( s.modes(:,1:nRotate), s.eigvals(1:nRotate), s.A, rotType );
end
    
% Add metadata for the analysis
s.metadata = [{'matrix';'MC';'noiseType';'p';'pcaArgs';'Rotated Modes';'Rotation Type'},...
    {matType; MC; noise; p; pcaArgs; nRotate; rotType}];

% Plot the output if desired
if plotting
    eofplot(s);
end

end

%%%%% Helper Functions %%%%%
function[MC, noiseType, pval, showProgress, svdArgs, blockMC, convergeTest, guessRuntime] = parseInputs(varargin)
inArgs = varargin;

% Set defaults
MC = NaN;
noiseType = NaN;
pval = NaN;
showProgress = 'blockProgress';
guessRuntime = 'noRuntime';
svdArgs = {'svd'};
blockMC = false;
convergeTest = true;

% Get input values
if ~isempty(inArgs)
    
    % Set flag switches to false
    isSvdsArg = false;
    isNoiseType = false;
    isPval = false;
    
    % Get each input
    for k = 1:length(inArgs)
        arg = inArgs{k};
        
        % Set MC if using rule N
        if (k==1) && ~strcmpi(arg, 'noSigTest')
            if length(inArgs) < 3
                error('Insufficient parameters given for Rule N test');
            else
                MC = arg;
                isNoiseType = true;
            end    
            
        % Get noisetype
        elseif isNoiseType
            noiseType = arg;
            isPval = true;
            isNoiseType = false;
            
        % Get the p value
        elseif isPval
            pval = arg;
            isPval = false;    
            
        % Get svds Args
        elseif isSvdsArg
            if isscalar(arg)
                svdArgs = {'svds', arg};
            else
                error('The svds flag must be followed by nModes');
            end
            
        % Decide whether to do the economy sized decomposition
        elseif strcmpi(arg, 'econ')
            svdArgs = {'svd', 'econ'};
            
        % Decide whether to show the MC iteration
        elseif strcmpi(arg, 'showProgress') 
            showProgress = 'showProgress';
            
        % Decide whether to guess the runtime
        elseif strcmpi(arg, 'estimateRuntime')
            guessRuntime = 'estimateRuntime';
            
        % Decide whether to perform the ruleN significance test
        elseif strcmpi(arg, 'noSigTest')
            blockMC = true;
            
        % Decide whether to save data on Monte Carlo convergence
        elseif strcmpi(arg, 'noConvergeTest')
            convergeTest = false;
            
        % Note that the next input is the svd Arg
        elseif strcmpi(arg, 'svds')
            if length(inArgs) >= k+1
                isSvdsArg = true;
            else
                error('The svds flag must be followed by nEigs or the ''econ'' flag');
            end
            
        % Anything else
        else
            error('Unrecognized Input');
        end
    end
    
% Must have at least some inputs
else
    error('Insufficient inputs');
end
end             
