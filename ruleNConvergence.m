function[MCsigExpVar, MCtrue_p] = ruleNConvergence( randExpVar, p )
%% Returns Monte Carlo convergence data for a significance level.

% Preallocate the output
MCsigExpVar = NaN( size( randExpVar) );
MCtrue_p = NaN( size(randExpVar,1), 1);

% For each MC iteration...
for k = 1:MC
    % Get the associated significance threshold and true significance level
    [MCsigExpVar(k,:), MCtrue_p(k)] = eofSigThreshold( randExpVar(1:k,:), p );
end

end