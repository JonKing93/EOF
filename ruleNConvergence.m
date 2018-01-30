function[MCsigExpVar, MCtrue_p] = ruleNConvergence( randExpVar, p )
%% Returns Monte Carlo convergence data for a significance level.
%
% [MCsigExpVar, MCtrue_p] = ruleNConvergence( randExpVar, p )
%
%
% ----- Inputs -----
%
% randExpVar: the random explained variances from a Rule N process
%
% p: The desired significance level for the test. p must be on the open
%   interval (0, 1).
%
%
% ----- Outputs -----
%
% MCsigExpVar: The explained variance significance threshold at each
%     successive Monte Carlo iteration.
%
% MCtrue_p: The true significance level being tested at each successive
%     Monte Carlo iteration.
%
%
% ----- Author -----
%
% Jonathan King, 2017, University of Arizona, jonking93@email.arizona.edu

% Preallocate the output
MCsigExpVar = NaN( size( randExpVar) );
MCtrue_p = NaN( size(randExpVar,1), 1);

% For each MC iteration...
for k = 1: size(randExpVar,1)
    % Get the associated significance threshold and true significance level
    [MCsigExpVar(k,:), MCtrue_p(k)] = eofSigThreshold( randExpVar(1:k,:), p );
end

end