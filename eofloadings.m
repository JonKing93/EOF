function[ax] = eofloadings(s)
%% Plots the loadings of significant EOF modes and rotated modes
%
% [ax] = eofloadings(s)
% Plots the loadings for all significant and rotated modes.
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
% ax: An array of the axes handles.
%
%
% ----- Written By -----
% 
% Jonathan King, 2017, University of Arizona (jonking93@email.arizona.edu)

% Check for required fields
if ~isfield(s,'nSig') || ~isfield(s, 'modes')
    warning('Insufficient data to plot significant loadings. Try running a significance test.');
    return;
elseif s.nSig == 0
    warning('None of the EOF modes are significant. Not plotting loadings...');
    return;
end
if isfield(s, 'rotModes')
    rot = size(s.rotModes,2);
else
    rot = 0;
end

% Create variable names if required
if isfield(s, 'varNames')
    varNames = s.varNames;
else
    nVars = size( s.modes,1);
    varNames = cell( nVars, 1);
    for k = 1:nVars
        varNames{k} = ['Var',num2str(k)];
    end
end

% For each plot
ax = [];
nPlots = s.nSig + rot;
for k = 1:nPlots
    % Initialize a figure
    figure();
    hold on;

    % Get the current set of loadings and some information for the title
    if k <= s.nSig
        loading = s.modes(:,k);
        tString = 'Unrotated';
        tNum = k;
    else
        loading = s.rotModes(:,k-s.nSig);
        tString = 'Rotated';
        tNum = k-s.nSig;
    end

    % Sort the loadings in ascending order
    [~, iLoad] = sort( loading );

    % Make a bar plot of the loading
    bar( loading( iLoad ) );

    % Add variable names
    set( gca, 'XTick', 1:length(loading), 'XTickLabel', varNames(iLoad) );

    % Add a title
    title( sprintf('%s Mode %0.f', tString, tNum) );

    % Add the axis handle to the output array
    ax = [ax; gca];
end