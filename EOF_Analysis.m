function[s] = EOF_Analysis(Data, matrix, varargin)
%% Performs a full EOF Analysis of a data set.
% 
% [s] = EOF_Analysis(Data, matrix, MC, noiseType, pval)
% Finds the EOF modes, explained variances, signals, and signals scaled to
% the standardized dataset. Performs a Rule N significance test on the modes
% and rotates the significant modes according to VARIMAX criterion. Returns
% all calculated values in a structure, s.
%
% [s] = EOF_Analysis(Data, matrix, 'noSigTest')
% A flag to block the Rule N significance testing. The returned structure
% will not contain any fields requiring the significance test.
%
% [s] = EOF_Analysis(..., 'showProgress')
% Displays the current Monte Carlo iteration onscreen.
%
% [s] = EOF_Analysis(..., 'econ')
% Performs the economy sized svd decomposition rather than the full svd.
%
% [s] = EOF_Analysis(..., 'svds', nModes)
% Uses the svds decomposition and determines the first nModes modes.
%
% [s] = EOF_Analysis(..., 'noConvergeTest')
% A flag to block the test for Monte Carlo convergence. The returned
% structure will contain neither the iterSigLevel nor iterTrueSig fields.
%
%
% ----- Inputs -----
%
% Data: A 2D data matrix. Each column corresponds to a particular data
%       series. Data may only contain numeric entries.
%
% matrix: The desired analysis matrix.
%       'cov': Covariance matrix -- Minimizes variance along EOFs
%       'corr': Correlation matrix -- Minimizes relative variance along
%               EOFs. Often useful for data series with significantly
%               different magnitudes.
%       'none': Perform svd directly on data matrix.
%  
% MC: The number of Monte Carlo iterations used in the Rule N significance
%       test
%
% noiseType: The noise used in the Rule N significance test
%       'white': white Gaussian noise
%       'red': lag-1 autocorrelated noise with added white noise
%
% pval: The significance level that the significance test should pass.      
%
%
% ----- Outputs -----
%
% s: A structure containing the following fields
%
%   Datax0: The standardized or detrended data matrix
%
%   C: The analysis matrix. The covariance or correlation matrix of Datax0,
%       or the original dataset.
%
%   eigVals: A vector with the eigenvalues of the analysis matrix. Each 
%       eigenvalue corresponds to an EOF mode. The larger the eigenvalue
%       magnitude, the more data variance explained by the mode.
%
%   modes: The EOF modes. These are the eigenvectors of the analysis 
%       matrix, also known as loadings. Each column is one mode.
%
%   expVar: The data variance explained by each mode. Equivalent to the
%       normalized eigenvalues.
%
%   signals: The signal for each EOF mode. Signals are the imprint of each
%       mode on the original data series, also known as scores or EOF 
%       time series. Each column is one signal.
%
%   scaSignals: The signals scaled to the standardized data matrix. Allows
%       for quick comparison of signals with the standardized data series.
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
[MC, noiseType, pval, showProgress, svdArgs, blockMC, convergeTest] = parseInputs(varargin{:});

% Declare the intial structure
s = struct();

% Run the initial EOF on the Data
[s.eigVals, s.modes, s.expVar, s.Datax0, s.C] = simpleEOF(Data, matrix, svdArgs{:});

% Get the signals
s.signals = getSignals(s.Datax0, s.modes);

% Scale the signals to the standardized data
s.scaledSignals = scaleSignals(s.signals, s.eigVals);

% If testing significance...
if ~blockMC

    % Run Rule N...
    if convergeTest                 % ... with convergence testing
        [s.nSig, s.randEigvals, s.thresh, s.trueConf, s.iterSigEigs, s.iterTrueConf] = ...
                ruleN(Data, matrix, s.expVar, MC, noiseType, pval, svdArgs, showProgress);
    
    else                            % ... without convergence testing
        [s.nSig, s.randEigvals, s.thresh, s.trueConf] = ...
            ruleN(Data, matrix, s.expVar, MC, noiseType, pval, svdArgs, showProgress, 'noConvergeTest');
    end
    
    % Perform a Varimax rotation
    if s.nSig < 2  % Less than 2 significant modes, rotation is not possible
        % Do nothing
        
    else 
        % Scale the eigenvectors
        s.scaModes = scaleModes( s.modes(:,1:s.nSig), s.eigVals(1:s.nSig));

        % Rotate the scaled eigenvectors
        [s.rotModes, s.rotEigvals, s.rotExpVar, s.rotMatrix] = ...
            varimaxRotation( s.scaModes(:,1:s.nSig), s.eigVals(1:s.nSig), s.expVar(1:s.nSig) );

        % Get the rotated signals
        s.rotSignals = getSignals(s.Datax0, s.rotModes);

        % Get the scaled rotated signals
        s.scaRotSignals = scaleSignals( s.rotSignals, s.rotEigvals);
    end
end

% Add metadata
s.conf = 1-pval;
s.metadata = [{'matrix';'MC';'noiseType';'svdArgs';'sigTest';'convergeTest'},...
    {matrix; MC; noiseType; svdArgs; ~blockMC; convergeTest}];
end

%%%%% Helper Functions %%%%%
function[MC, noiseType, pval, showProgress, svdArgs, blockMC, convergeTest] = parseInputs(varargin)
inArgs = varargin;

% Set defaults
MC = NaN;
noiseType = NaN;
pval = NaN;
showProgress = 'blockProgress';
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
