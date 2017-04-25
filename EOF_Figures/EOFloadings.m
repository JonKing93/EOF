function[] = EOFloadings(EOF, names, rotType)
%% Plots the loadings for the elements on an EOF
%
% -- In --
%
% EOF: the structure from an EOF_Analysis
%
% names: The names of each of the elements
%
% rotType: a flag to use rotated vs unrotated modes
%   'rotated': Use the rotated modes
%   'unrotated': Use the unrotated modes

if EOF.nSig > 0
    
    % Choose between rotated or unrotated modes
    if strcmpi( rotType, 'rot') || strcmpi(rotType, 'rotated')
        loadings = EOF.rotModes;
        rotString = 'Rotated';
        
    elseif strcmpi(rotType,'unrot') || strcmpi(rotType, 'unrotated') || isempty(rotType)
        loadings = EOF.modes;
        rotString = 'Unrotated';
        
    else
        error('Unrecognized rotType');
    end
    
    % For each significant loading
    for k = 1:EOF.nSig
        figure(); clf; hold on;
        
        % Sort the loadings in ascending order
        [~, iL] = sort(loadings(:,k));
        
        % Make a bar plot of the loading
        bar( loadings(iL, k) );
        
        % Plot the names on the loadings
        set(gca, 'XTick', 1:size(loadings,1), 'XTickLabel', names(iL));
        title( sprintf('%s Mode %i', rotString, k) );
        
    end
end