function[ax] = eofloadings(s)
%% Plots the loadings of significant EOF modes and rotated modes.
%
% [ax] = eofloadings(s)
% Plots the loadings for all significant and rotated modes.
%
% ----- Inputs -----
%
% s: The output structure from an EOF_Analysis
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
    fprintf( sprintf('Insufficient data to plot significant loadings. Try running a significance test.\r\n') );
    return;
elseif s.nSig == 0
    fprintf( sprintf('None of the EOF modes are significant. Not plotting loadings...\r\n') );
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