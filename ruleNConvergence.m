function[MCsigExpVar, MCtrue_p, MCnSig] = ruleNConvergence( randExpVar, p, expVar )
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
% expVar: The explained variance of data EOF modes.
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
% MCnSig: The number of significant modes at each successive Monte Carlo
%   iteration.
%
% ----- Author -----
%
% Jonathan King, 2017, University of Arizona, jonking93@email.arizona.edu

% Get the significance threshold at each MC iteration
[MCsigExpVar, MCtrue_p] = mcthreshold( randExpVar, p, 'converge' );

% Preallocate
nMC = size(MCsigExpVar,1);
MCnSig = NaN( nMC, 1);

% Find the number of significant leading modes for each iteration
for k = 1:nMC
    MCnSig(k) = find( expVar <= MCsigExpVar(k,:), 1, 'first') - 1;
end

end