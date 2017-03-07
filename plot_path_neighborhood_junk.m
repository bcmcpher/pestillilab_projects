function plot_path_neighborhood()
%% plot path neighborhood of an edge
% modified from Cesar's code by Brent
%

%% load data

% using HCP 105115 prob lmax10 as example
load data/hcp_105115_prob_lmax10_rep02.mat

% load LiFE fit fg
load /N/dc2/projects/lifebid/code/ccaiafa/Caiafa_Pestilli_paper2015/Results/ETC_Dec2015/Single_TC/105115/fe_structure_105115_STC_run01_SD_PROB_lmax10_connNUM02.mat

% fix fiber group, because that is necessary
wbFG = feGet(fe, 'fg acpc');

anatomy  = niftiRead(fe.path.anatomy);

%% build plotting paramaters

% switch to clean tract before render - false/true
doCleaning = true; 

% pick connection
conn = 20; 

% create title
roi1 = strrep(pconn{conn}.roi1, '_label', '');
roi1 = strrep(roi1, '_', '.');
roi2 = strrep(pconn{conn}.roi2, '_label', '');
roi2 = strrep(roi2, '_', '.');
fgName = ['Index: ' num2str(conn) '; ' roi1 '-to-' roi2];

% generate a series of costum colors for the tracts
colors     = {[.1 .25 .65],[.75 .25 .1]};
viewCoords = [90,0];
slice      = [-1 0 0];
proportion_to_show = .01;
threshold_length = 10;
     
%% make the figures

for conn = conn_set
    
    %% STC TENSOR
 
    %run01
    fe_structureFile = deblank(ls(char(fullfile(dataRootPathSTC,'105115',strcat('fe_structure_*',subject,'*_STC_','*run01','_tensor__','*',conn,'.mat'))))); 
    TractsFile = deblank(ls(char(fullfile(dataRootPathSTC,strcat('fe_structure_*',subject,'*_STC_','*run01_500000','_tensor__','*',conn,'_TRACTS.mat')))));
    
    load(fe_structureFile);
    load(TractsFile); 

    ind_tracts1 = find(classification.index==conn); 
    
    % We clean all major tracs
    switch doCleaning
        case {true}
            for iTract = conn
                fprintf('\n Cleaning %s ...',classification.names{iTract})
                [fg_tract, keep] = mbaComputeFibersOutliers(fascicles(iTract),3,3);
            end
            
            saveCleanTag = '_CLEAN';
            ind_tracts1 = ind_tracts1(keep);
        otherwise
            saveCleanTag = '_';
    end
    
    ind_nnz = find(fe.life.fit.weights);
    ind1 = ind_nnz(ind_tracts1);
    ind_tracts2 = feGet(fe,'Path Neighborhood',ind1);
    fg          = feGet(fe,'fibers img');
    fg_pathn    = fgExtract(fg,ind_tracts2,'keep');
    clear fg
    roicoords = feGet(fe,'coords from fibers',ind1); 
    fg_pathn  = feClipFibersToVolume(fg_pathn,roicoords,3);
    fg_pathn.fibers  = mbaFiberSplitLoops(fg_pathn.fibers);
    xform     = feGet(fe,'img2acpcxform');
    fg_pathn  = dtiXformFiberCoords(fg_pathn,xform,'acpc');

    % Prepare the plot of tract and neighborhood
    fg{1}    = fg_tract;
    fg_pathn = rmfield(fg_pathn,'coordspace');
    fg_pnplot = fg_pathn;
    
    c = 1;
    for ii = 1:length(fg_pathn.fibers)
       if length(fg_pathn.fibers{ii}) > threshold_length
           fibers{c} = fg_pathn.fibers{ii};
        c = c + 1;
       end
    end
    fg_pathn.fibers = fibers; clear fibers
    fibs_indx = randsample(1:length(fg_pathn.fibers),round(length(fg_pathn.fibers)*proportion_to_show));
    fg_pnplot.fibers = fg_pnplot.fibers(fibs_indx);
    fg{2}    = fg_pnplot;
 
    fig_name      = sprintf('anat_tracts_%s_%s',tract_names,strcat('fe_structure_',subject,'_STC_','run01','_tensor__',conn_set{1}));    
    [fig_h, ~, ~] = plotFascicles(fg, colors, slice, anatomy, viewCoords, fig_name, [1 2]);
    feSavefig(fig_h,'verbose','yes','figName',[fig_name, 'leftSAG'],'figDir',dataOutputPath,'figType','jpg');
    close all; drawnow
 
    %% STC PROB/STREAM lmax2...lmax12
    for alg = alg_set
        for lparam = lparam_set
          
            %run01
            fe_structureFile = deblank(ls(char(fullfile(dataRootPathSTC,'105115',strcat('fe_structure_*',subject,'*_STC_','*run01','*',alg,'*',lparam,'*',conn,'.mat'))))); 
            TractsFile = deblank(ls(char(fullfile(dataRootPathSTC,strcat('fe_structure_*',subject,'*_STC_','*run01','*',alg,'*',lparam,'*',conn,'_TRACTS.mat')))));
            load(fe_structureFile);
            load(TractsFile); 
            ind_tracts1 = find(classification.index==conn); 
    
             % We clean all major tracs
             switch doCleaning
                 case {true}
                     for iTract = conn
                         fprintf('\n Cleaning %s ...',classification.names{iTract})
                         [fg_tract, keep] = mbaComputeFibersOutliers(fascicles(iTract),3,3);
                     end
                     
                     saveCleanTag = '_CLEAN';
                     ind_tracts1 = ind_tracts1(keep);
                 otherwise
                     saveCleanTag = '_';
             end
             
             ind_nnz = find(fe.life.fit.weights);
             ind1    = ind_nnz(ind_tracts1);
             ind_tracts2 = feGet(fe,'Path Neighborhood',ind1);
             fg          = feGet(fe,'fibers img');
             fg_pathn    = fgExtract(fg,ind_tracts2,'keep');
             clear fg
             roicoords = feGet(fe,'coords from fibers',ind1);
             fg_pathn = feClipFibersToVolume(fg_pathn,roicoords,3);
             fg_pathn.fibers = mbaFiberSplitLoops(fg_pathn.fibers);
             xform     = feGet(fe,'img2acpcxform');
             fg_pathn  = dtiXformFiberCoords(fg_pathn,xform,'acpc');
             
             % Prepare the plot of tract and neighborhood
             fg{1}    = fg_tract;
             fg_pathn = rmfield(fg_pathn,'coordspace');
             fg_pnplot = fg_pathn;
             
             c = 1;
             for ii = 1:length(fg_pathn.fibers)
                 if length(fg_pathn.fibers{ii}) > threshold_length
                     fibers{c} = fg_pathn.fibers{ii};
                     c = c + 1;
                 end
             end
             fg_pathn.fibers = fibers; clear fibers
             fibs_indx = randsample(1:length(fg_pathn.fibers),round(length(fg_pathn.fibers)*proportion_to_show));
             fg_pnplot.fibers = fg_pnplot.fibers(fibs_indx);
             fg{2}    = fg_pnplot;
             
             % plot ARC and PN
             fig_name      = char(strcat('anat_tracts_',tract_names,'_fe_structure_', ...
                 subject,'_STC_run01_', ...
                 alg ,'_', ...
                 lparam,'_', ...
                 conn) );
             [fig_h, ~, ~] = plotFascicles(fg, colors, slice, anatomy, viewCoords, fig_name, [1 2]);
             feSavefig(fig_h,'verbose','yes','figName',[fig_name, 'leftSAG'],'figDir',dataOutputPath,'figType','jpg');
             close all; drawnow
             
             % plot ARC
             fig_name      = char(strcat('anat_tracts_','ARC','_fe_structure_', ...
                 subject,'_STC_run01_', ...
                 alg ,'_', ...
                 lparam,'_', ...
                 conn) );
             [fig_h, ~, ~] = plotFascicles(fg, colors, slice, anatomy, viewCoords, fig_name, [1]);
             feSavefig(fig_h,'verbose','yes','figName',[fig_name, 'leftSAG'],'figDir',dataOutputPath,'figType','jpg');
             close all; drawnow
             
             % plot PN
             fig_name      = char(strcat('anat_tracts_','PN','_fe_structure_', ...
                 subject,'_STC_run01_', ...
                 alg ,'_', ...
                 lparam,'_', ...
                 conn) );
             [fig_h, ~, ~] = plotFascicles(fg, colors, slice, anatomy, viewCoords, fig_name, [2]);
             feSavefig(fig_h,'verbose','yes','figName',[fig_name, 'leftSAG'],'figDir',dataOutputPath,'figType','jpg');
             close all; drawnow             
             
        end
    end
end

end %% Main function

% Local functions to plot the tracts
function [fig_h, light_h, brain_h] = plotFascicles(fascicles, color, slice, anatomy, viewCoords, fig_name,tracts_to_clean)
fig_h = figure('name',fig_name,'color','k');
set(gca,'visible','on','ylim',[-108 69],'xlim',[-75 75],'zlim',[-45 78])
brain_h = mbaDisplayBrainSlice(anatomy, slice,gca);
hold on
for iFas  = 1:length(tracts_to_clean)
    hold on
    [~, light_h] = mbaDisplayConnectome(fascicles{ tracts_to_clean(iFas) }.fibers,fig_h,color{ tracts_to_clean(iFas) },'single');
    delete(light_h)
end
view(viewCoords(1),viewCoords(2))
light_h = camlight('right');
lighting phong;

%set(fig_h,'Units','normalized', 'Position',[0.5 .2 .4 .8]);
drawnow

end

% Local functions to plot the tracts
function [fig_h, light_h] = plotFasciclesNoAnat(fascicles, color, viewCoords, fig_name,tracts_to_clean)
fig_h = figure('name',fig_name,'color','k');
hold on
set(gca,'visible','off','ylim',[-108 69],'xlim',[-75 75],'zlim',[-45 78])
for iFas  = 1:length(tracts_to_clean)
    [~, light_h] = mbaDisplayConnectome(fascicles{ tracts_to_clean(iFas) }.fibers,fig_h,color{ tracts_to_clean(iFas) },'single');
    delete(light_h)
end
view(viewCoords(1),viewCoords(2))
light_h = camlight('right');
lighting phong;
%set(fig_h,'Units','normalized', 'Position',[0.5 .2 .4 .8]);
drawnow

end

function colors = getColors

% % Prepare the colors for plotting, Left HM warm, Right HM cool
% numColRes = 12;
% allColors = 1:1:numColRes;
% colormaps = {'spri ng','summer','autumn','winter','bone'};
% for iMap = 1:length(colormaps)
% %figure(iMap)
% for iFas  = 1:length(allColors)
% color{iMap, iFas} = getSmoothColor(allColors(iFas),numColRes,colormaps{iMap});
% %plot(iFas,1,'o','markerfacecolor',color{iMap, iFas},'markersize',16); hold on
% %text(iFas-.1,1,sprintf('%i',iFas),'color','w')
% end
% end
% colors = {color{1,10}, color{1,10}, color{1,5}, color{1,5}, ...
%     color{2,5}, color{2,5},   color{3,1}, color{3,1}, ...
%     color{5,8}, color{5,10},  color{3,6}, color{3,6}, ...
%     color{4,6},  color{4,6},  color{4,9},  color{4,9}, ...
%     color{2,3},  color{2,3},  color{4,3},  color{4,3}};
%
colors = {[233,150,122]./255, [233,150,122]./255, ... % Salmon
    [255,215,0]./255,   [255,215,0]./255, ... % Gold
    [64,224,208]./255,  [64,224,208]./255, ...% Turquise
    [255,99,71]./255,   [255,99,71]./255,  ...% Tomato
    [220 220 220]./255, [220 220 220]./255, ...% Gainsboro
    [220,20,60]./255,   [220,20,60]./255,   ...
    [221,160,221]./255, [221,160,221]./255, ...
    [199,21,133]./255,  [199,21,133]./255, ...
    [230,230,250]./255, [230,230,250]./255, ...
    [100,149,237]./255, [100,149,237]./255};

%figure
%for iFas  = 1:length(colors)
%plot(iFas,1,'o','markerfacecolor',colors{iFas},'markersize',16); hold on
%text(iFas-.1,1,sprintf('%i',iFas),'color','w')
%end
end

function figDir = feSavefig(h,varargin)
% Saves a figure for publishing purpose.
%
%  function figDir = savefig(h,varargin)
%
% INPUTS: Must be pairs, e.g., 'name',value
%
%   figName  - the name of the figure file.
%              Default: sprintf('feFig_%s',get(h,'Name'))
%   figDir   - Name of the subfolder where to save this figure.
%              Default: current dir ('.')
%   figType  - fig, eps, jpg, png
%              Default: 'png' fastest smallest file, low resolution.
%   verbose  - display on screen the print command being invoked.
%
% NOTES:
%   This function invokes a series of print commands:
%   print(h, '-cmyk', '-painters','-depsc2','-tiff','-r500', '-noui', 'figureName')
%
% EXAMPLE:
%   feSavefig(figureHandle,'verbose','yes','figName','myfig','figDir','/path/to/fig/folder/');
%
%
% Copyright (2016), Franco Pestilli, Indiana University, pestillifranco@gmail.com.

% set up default variables:
figName           = sprintf('feFig_%s',get(h,'Name')); % the name of the figure file
figDir            = '.'; % the subfolder where to save this figure
figType           = 'png';
verbose           = 'yes'; % 'yes', 'no', display on screen what it is going on

if ~isempty(varargin)
    if mod(length(varargin),2), error('varargin must be pairs'); end
    for ii=1:2:(length(varargin)-1)
        eval(sprintf('%s = ''%s'';',varargin{ii}, varargin{ii+1}));
    end
end

% make the figure dir if it does not exist:
if ~isdir(figDir), mkdir(figDir);end

% Create a print command that will save the figure
switch figType
    case {'png'}
        printCommand = ...
            sprintf('print(%s, ''-painters'',''-dpng'', ''-noui'', ''%s'')', ...
            num2str(h.Number),fullfile(figDir,figName));
        
    case {'jpg'}
        printCommand = ...
            sprintf('print(%s, ''-djpeg95'',''-r500'', ''-noui'', ''%s'')', ...
            num2str(h.Number),fullfile(figDir,figName));
        
    case {'eps'}
        printCommand = ...
            sprintf('print(%s, ''-cmyk'', ''-painters'',''-depsc2'',''-tiff'',''-r500'' , ''-noui'', ''%s'')', ...
            num2str(h.Number),fullfile(figDir,figName));
        
    case {'fig'}
        printCommand = ...
            sprintf('saveas(%s, ''%s'', ''fig'')', num2str(h.Number),figName);
        
    otherwise
        error('[%s] Cannot save figure with type set to: %s', mfilename,figType)
end

if strcmpi(verbose,'yes')
    fprintf('[%s] Saving eps figure %s/%s\nUsing command: \n%s...\n', ...
        mfilename,figDir,figName,printCommand);
end

% Do the printing here:
eval(printCommand);

% Delete output if it was nto requested
if (nargout < 1), clear figDir;end

end
