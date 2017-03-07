%% load fit fe structures across repeats and view a bunch of edges
% Brent McPherson
% 20161220
%

% will need to loop over subjects

% using HCP 105115 prob lmax10 as example
load data/hcp_105115_prob_lmax10_rep02.mat

% load LiFE fit fg
load /N/dc2/projects/lifebid/code/ccaiafa/Caiafa_Pestilli_paper2015/Results/ETC_Dec2015/Single_TC/105115/fe_structure_105115_STC_run01_SD_PROB_lmax10_connNUM02.mat

% fix fiber group, because that is necessary
wbFG = feGet(fe, 'fg acpc');

% for currently loaded subject
for ii = 1:length(pconn)
    
    % check if there are enough fibers to plot, otherwise skip
    if length(pconn{ii}.end.nzfibs) < 50
        continue
    end
    
    % create title
    roi1 = strrep(pconn{ii}.roi1, '_label', '');
    roi1 = strrep(roi1, '_', '.');
    roi2 = strrep(pconn{ii}.roi2, '_label', '');
    roi2 = strrep(roi2, '_', '.');
    fgName = ['Index: ' num2str(ii) '; ' roi1 '-to-' roi2];
    
    % pull indices from pconn for all fibers
    % add if statement to find all / non-zero fibers
    fgex{ii} = fgCreate('name', fgName, 'colorRgb', [0 1 0], 'fibers', {wbFG.fibers{pconn{ii}.end.nzfibs}}');
    
    % clean outlier fibers
    fgcx{ii} = mbaComputeFibersOutliers(fgex{ii}, 3, 3, 100, 'mean');
    
    % plot command
    bsc_quickEdgeCheck(fgex{ii}, fgcx{ii}, fgName);
    
    pause;
    
    % save a simple png of plot
    %figName = ['edge_' num2str(ii) '_of_' roi1 '_' roi2];
    %print(['figs/quickEdges/' figName '.png'], '-dpng', gcf);
    %close all;

end



