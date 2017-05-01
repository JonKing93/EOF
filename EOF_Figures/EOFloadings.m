function[ax] = EOFloadings(EOF, names, rotType, varargin)
%% Plots the loadings of data elements for EOFs
%
% [ax] = EOFloadings(EOF, names, rotType)
% Plots the loadings for data elements on significant EOFs. Returns the
% plot axes.
%
% [ax] = EOFloadings(EOF, names, rotType, 'all')
% Plots the loadings for data elements on all EOFs
%
% [ax] = EOFloadings(EOF, names, rotType, iEOF)
% Plots the loadings for data elements on user-selected EOFs
%
%
% ----- Inputs -----
%
% EOF: the structure from an EOF_Analysis
%
% names: A string array of the names of each of the Data elements in the 
%       order they are given in EOF_Analysis.
%
% rotType: a flag to use rotated vs unrotated modes
%   'rot' OR 'rotated': Use the rotated modes
%   'unrot' OR 'unrotated: Use the unrotated modes
%
% iEOF: An array containing the indices of the desired EOFS
%       (e.g. iEOF = [1 2 5] will plot the first, second, and fifth EOF loadings)
%
%
% ----- Outputs -----
%
% ax: The axes handle for the plot
%
%
% ----- Written By -----
% 
% Jonathan King, 2017, University of Arizona (jonking93@email.arizona.edu)

% Parse Inputs, error check
[rotString, loadings, iEOF] = parseInputs(EOF, rotType, varargin{:});
errCheck(loadings, names, iEOF);

% Get the EOF loadings to plot if not specified
if strcmpi(iEOF, 'all')
    iEOF = 1:size(loadings, 2);
elseif isnan(iEOF)
    iEOF = 1:EOF.nSig;
end

if ~isempty(iEOF)
   
    % For each significant loading
    for k = iEOF
        figure(); clf; hold on;
        
        % Sort the loadings in ascending order
        [~, iL] = sort(loadings(:,k));
        
        % Make a bar plot of the loading
        bar( loadings(iL, k) );
        
        % Plot the names on the loadings
        set(gca, 'XTick', 1:size(loadings,1), 'XTickLabel', names(iL));
        title( sprintf('%s Mode %i', rotString, k) );
        
    end
else
    
end
end

% ----- Helper Methods -----
function[rotString, loadings, iEOF] = parseInputs(EOF, rotType, varargin)

% string for rotated vs unrotated modes
if strcmpi( rotType, 'rot') || strcmpi(rotType, 'rotated')
    loadings = EOF.rotModes;
    rotString = 'Rotated';

elseif strcmpi(rotType,'unrot') || strcmpi(rotType, 'unrotated') || isempty(rotType)
    loadings = EOF.modes;
    rotString = 'Unrotated';

else
    error('Unrecognized rotType');
end

% Optional inputs
inArgs = varargin;

% Set defaults
iEOF = NaN;

if ~isempty(inArgs)
    
    if length(inArgs) == 1
        if strcmpi( inArgs{1}, 'all')
            iEOF = 'all';
        elseif isvector( inArgs{1} )
            iEOF = inArgs{1};
        else
            error('Unrecognized input');
        end        
    else
        error('Too many input arguments');
    end
end
end
   
function[] = errCheck(loadings, names, iEOF)
[nEls, nLoads] = size(loadings); 

if length(names) ~= nEls
    error('Number of names given does not match the number of data columns');
elseif ~ischar(iEOF) && ~isstring(iEOF) && ~any(isnan(iEOF))
    if any(iEOF > nLoads)
        error('EOF index exceeds the number of EOFs');
    elseif any(iEOF<0 | mod(iEOF,1)~=0)
        error('iEOF must be a vector of positive integers');
    end
end
end
