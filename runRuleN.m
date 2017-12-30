function randExpVar = runRuleN(k, MC, noise, Data, matrix, pcaArgs, estimateRuntime, showProgress, nWorkers)
%% Generates random explained variances within a for or parfor loop in ruleN.
%
% randExpVar = runRuleN(k, MC, noise, Data, matrix, pcaArgs, estimateRuntime, showProgress, nWorkers)
%
%
% ----- Author -----
%
% Jonathan King, 2017, University of Arizona, jonking93@email.arizona.edu
     
% Begin recording runtime if an estimate is desired
if estimateRuntime && k==1
    startTime = tic;
end

% Create a random matrix with the desired noise properties, scaled to data standard deviation.
g = randNoiseSeries(noise, Data);
    
% Get the analysis matrix
if strcmp(matrix, 'corr')
    g = zscore(g);
elseif strcmp(matrix, 'cov')
    g = detrend(g, 'constant');
end
    
% Run a pca (eof analysis) on the random matrix and get explained variance
[~,~,~,~, randExpVar] = pca(g, pcaArgs{:});
    
% Report runtime estimate
if estimateRuntime && k==1
    time = toc(startTime);
    time = time*MC/nWorkers;
    h = time / 360;
    m = mod(time,360)/60;
    s = mod(time,60);
    fprintf('Estimated runtime: %0.f hour(s), %0.f minute(s), %0.f seconds \r\n',h,m,s);
end
    
% Update the progress bar if displaying
if showProgress
    progressbar(k/MC);
end

end   