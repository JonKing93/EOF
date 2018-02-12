function[varargout] = mcthreshold( mcVals, p, varargin )
%% Gets significance or confidence thresholds for a Monte Carlo process
%
% [mcThresh, true_p] = mcthreshold( mcVals, p )
% Finds the threshold for the p significance level (e.g. p<0.05) for a 
% Monte Carlo process and calculates the true significance level tested.
%
% [mcThresh, true_c] = mcthreshold( mcVals, c, 'conf' )
% Finds the threshold for the c confidence interval (e.g. 95% confidence)
% for a Monte carlo process and determines the true confidence interval
% tested.
%
% [upThresh, lowThresh, (true_p / true_c)] = mcthreshold( ..., '2tail' )
% Performs a centered, 2-tailed significance or confidence test. Upper and
% lower tails are stored in the last dimension of the output, respectively.
% 
% [...] = mcthreshold( ..., 'dim', d )
% Specifies the Monte carlo dimension. By default, mcthreshold assumes the
% first dimension holds successive monte carlo iterations.
%
% [...] = mcthreshold( ..., 'converge')
% Calculates the threshold values and significance levels after each
% successive Monte Carlo iteration to test the process for convergence to a
% stable value.


% Parse inputs. Error Check.
[conf, conf2, isconf, d, converge] = setup(mcVals, p, varargin{:});

% Get the initial size / dimensionality
sInitial = size(mcVals);

% If the Monte Caro dimension is not first, permute it to first
valDims = 1: max( ndims(mcVals), d);
if d~=1
    valDims(d) = 1;
    valDims(1) = d;
    
    mcVals = permute( mcVals, valDims );
end

% Get the number of Monte Carlo iterations
sVals = size(mcVals);
nMC = sVals(1);

% Determine the number of data subsets used in the analysis
kstart = nMC;
if converge
    kstart = 1;
end

% Preallocate 
nElts = nMC-kstart+1;
nConf = length(conf);

trueConf = NaN(nElts,nConf);
upperThresh = NaN( [nConf, nElts sVals(2:end)] );
if ~isnan(conf2)
    lowerThresh = upperThresh;
end

% For each data subset
row = 1;
for k = kstart : nMC

    % Sort the random MC variables
    mcVals(1:k,:) = sort( mcVals(1:k,:), 1 );

    % Get the threshold index for significance. Round conservatively.
    upperIndex = ceil( k .* conf );

    % Get the true confidence interval for the test
    trueConf(row,:) = upperIndex ./ k;

    % Get the threshold values
    upperThresh(:,row,:) = mcVals(upperIndex,:);


    % If this is a two tailed test...
    if ~isnan(conf2)

        % Get the lower threshold. Round conservatively.
        lowIndex = max( floor( k * conf2 ), 1);

        % Get the threshold values
        lowerThresh(:,row,:) = mcVals(lowIndex,:);

        % Get the 2-tailed true confidence
        trueConf(row,:) = (upperIndex - floor(k.*conf2)) ./ k;
    end
    
    % Increment the row
    row = row+1;
end

% Permute to original form with P as the final dimension
upperThresh = permOutput( upperThresh, valDims, sInitial );
if ~isnan(conf2)
    lowerThresh = permOutput(lowerThresh, valDims, sInitial);
end

% If the user provided significance levels, convert from c to p
if ~isconf
    trueConf = 1-trueConf;
end

% Set the output
if isnan(conf2)
    varargout = {upperThresh, trueConf};
else
    varargout = {upperThresh, lowerThresh, trueConf};
end

end

function[thresh] = permOutput( thresh, valDims, sInitial )
    % Restore output dimensionality
    thresh = permute(thresh, [1 1+valDims]);
    % Add p to final dimension
    if sInitial(end) == 1          % Column vector, add p to d2
        thresh = permute( thresh, [2 1] );
    else
        thresh = permute(thresh, [2:max(valDims+1),1]);
    end
end

function[conf, conf2, isconf, d, converge] = setup(mcVals, conf, varargin)
        
    % Parse the inputs
    [isconf, tail2, d, converge] = parseInputs( varargin, {'conf','2tail','dim','converge'}, {false, false, 1, false}, {'b','b',{},'b'} );

    % Check that p/c is on the correct interval.
    if ~isvector(conf) || any(conf<=0) || any(conf>=1)
        error('''p'' or ''c'' must be a vector of values on the open interval (0,1)');
    end

    % Convert to confidence value if this is a p-value
    if ~isconf
        conf = 1-conf;
    end

    % If this is a 2 tailed test, get the centered tails
    conf2 = NaN;
    if tail2
        conf2 = (1-conf)/2;
        conf = conf + conf2;
    end

    % Warn about NaN sorting
    if any( isnan( mcVals(:) ) )
        warning('mcVals contains NaN elements. This may affect sorting.');
    end
    
    % Check that d is an allowed dimension
    if ~isscalar(d) || d<=0 || mod(d,1)~=0
        error('d must be a positive, integer scalar.');
    end
end
