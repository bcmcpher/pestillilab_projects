%% Brent McPherson
% 20161030
% plot anatomy of edges
%

%% load data

% using HCP 105115 prob lmax10 as example
load data/hcp_105115_prob_lmax10_rep02.mat

% load LiFE fit fg
load /N/dc2/projects/lifebid/code/ccaiafa/Caiafa_Pestilli_paper2015/Results/ETC_Dec2015/Single_TC/105115/fe_structure_105115_STC_run01_SD_PROB_lmax10_connNUM02.mat

% fix fiber group, because that is necessary
wbFG = feGet(fe, 'fg acpc');

% point to anatomy to display
anat = niftiRead(fe.path.anatomy);

% load an FA map - doesn't exist yet...
famp = niftiRead('/N/dc2/projects/lifebid/2t1/HCP/105115/fibers_new/dwi_data_b2000_aligned_trilin_fa.nii.gz');

% or just use dt6
%dt6p = '/N/dc2/projects/lifebid/HCP/Dan/105115/dt6_b2000trilin/dtiInit/dt6.mat';

%% make a plot of fiber groups

% because why wouldn't I need to parallelize this
parpool;

% find stongest connections in all connections...
emdm = emat(:,:,14);

% plot labeled matrix because I give up
pltLabel = strrep(roiNames, '_label.nii.gz', '');
pltLabel = strrep(pltLabel, '_', '.');
figure('Position', [350 150 1250 1100]);
colormap('hot');
imagesc(log(emdm));
axis('square'); axis('equal'); axis('tight');
%y = colorbar('FontSize', 29, 'Ticks', [0, 25, 50], 'TickLabels', {0, 25, 50});
caxis([-2 3]);
%ylabel(y, 'LiFE Measure');
%set(gca, 'XTickLabel', pltLabel, 'YTickLabel', pltLabel, 'XTick', [1:68], 'YTick', [1:68], 'XTickLabelRotation', 90);
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1], 'LineWidth', 2);
line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1], 'LineWidth', 2);
line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1], 'LineWidth', 2);

% subset to lower diagonal so ranks aren't duplicated
emdl = tril(emdm);

% sort and find x / y indices in matrix for strongest connections
[ vals, mInd ] = sort(emdl(:), 1, 'descend');
[x, y] = ind2sub(size(emdl), mInd);

% look, I'm not insane
% show the first 10 x, y coordinates for highest values
[x(1:10) y(1:10)]
% look up the highest value in the matrix
emdm(x(2), y(2))
% check that the highest value is the same in sort as matrix
vals(2)
% all checks out...

% create x / y index key
indx = nchoosek(1:68, 2);
tmp = 1:length(indx);
tmp = tmp';
indx = [indx tmp];
clear tmp

% go from x / y in matrix to index
xi = 20;
yi = 40;
indx((indx(:,1) == xi & indx(:,2) == yi), 3)
indx((indx(:,1) == yi & indx(:,2) == xi), 3)
% this can't be right...
% it turns out the stongest connection in the network is one really highly
% weighted fiber. well so much for that.

% go from index to x / y in matrix
val = 2171; % 2199
indx(indx(:,3) == val, 1:2)

% % run virtual lesion again to make sure I'm not insane
% [ ewVL, ewoVL ] = feComputeVirtualLesion(fe, pconn{2171}.end.fibers);
% vlOut = feComputeEvidence(ewoVL, ewVL);
% vlOut.em.mean

% or just look at 1
%conn = 1;

% quick plot a bunch of fibers to quickly inspect anatomy
for ii = 1:length(pconn)
    
    % check if there are enough fibers to plot, otherwise skip
    if length(pconn{ii}.end.fibers) < 500
        continue
    end
    
    % create title
    roi1 = strrep(pconn{ii}.roi1, '_label', '');
    roi1 = strrep(roi1, '_', '.');
    roi2 = strrep(pconn{ii}.roi2, '_label', '');
    roi2 = strrep(roi2, '_', '.');
    fgName = ['Index: ' num2str(ii) '; ' roi1 '-to-' roi2];
    
    % pull indices from pconn for non-zero fibers
    fgex{ii} = fgCreate('name', fgName, 'colorRgb', [0 1 0], 'fibers', {wbFG.fibers{pconn{ii)}.end.fibers}}');
    fgcx{ii} = mbaComputeFibersOutliers(fgex{ii}, 4, 4, 100, 'mean');
    bsc_quickEdgeCheck(fgex{ii}, fgcx{ii}, fgName);
    pause;
end
% create good renders of indices (1000 min): 3, 1720, 1725
% create good renders of indices (750  min): 1730, 2201
% create good renders of indices (500  min): 509, 746, 1962, 2144, 
% 2069 - rh.insula to rh.superiorparietal has 2 bundles as 1 edge 

% connections I want rendered pretty
conn = [3 1720 1725 1730 2201 509 746 1962 2144];

% create tract profile
for ii = 1:length(conn)
    
    % fix ROI names for plotting / labeling
    roi1 = strrep(pconn{ii}.roi1, '_label', '');
    roi1 = strrep(roi1, '_', '.');
    roi2 = strrep(pconn{ii}.roi2, '_label', '');
    roi2 = strrep(roi2, '_', '.');
    
    % pull indices from pconn for non-zero fibers
    fgName = [roi1 '-to-' roi2];
    fgex{ii} = fgCreate('name', fgName, 'colorRgb', [0 1 0], 'fibers', {wbFG.fibers{pconn{conn(ii)}.end.fibers}}');
    
    % clean fibers with mba
    fgcx{ii} = mbaComputeFibersOutliers(fgex{ii}, 4, 4, 100, 'mean');
    %fgcx{ii} = fgex{ii};
    
    % find stats on FA map w/ vistasoft
    fgst{ii} = dtiComputeDiffusionPropertiesAlongFG(fgcx{ii}, famp, [], [], 100);
    
    % y axis labels
    uplw = minmax(fgst{ii});
    mdpt = ((max(fgst{ii})-min(fgst{ii})) / 2) + min(fgst{ii});
    ybnd = [uplw(1) mdpt uplw(2)];
    
    % plot the FA tract profile
    figure('Position', [950 605 1350 645]); hold on;
    %title([ 'Tract Profile between ' roi1 ' and ' roi2 ]);
    plot(fgst{ii}, 'LineWidth', 4);
    set(gca, 'XTick', [0:50:100], 'YTick', ybnd, 'FontSize', 29);
    %xlabel([ 'Fiber Nodes (' roi1 ' to ' roi2 ')']);
    xlabel('Nodes Along Fiber', 'FontSize', 29);
    ylabel('FA Value', 'FontSize', 29);
    hold off;
    %print(gcf, '-cmyk', '-painters', '-depsc2', '-tiff', '-r500', '-noui', [ 'sfn2016/anatomy/edge-stats-' sprintf('%02d', ii)]);
    
    close all;
end

% create rendereing of tract
for ii = 1:length(conn)
    
    % render FA tract profile
    core = dtiComputeSuperFiberRepresentation(fgcx{ii}, [], 100);
    coords  = core.fibers{1};
    radius  = 3;    
    color   = fgst{ii}; %dtiComputeDiffusionPropertiesAlongFG(fgst{ii}, famp);
    subdivs = 20;
    cmap = 'jet';
    crange = [0.20 0.50];
    
    % not the best plot to show these
    %lightH = AFQ_RenderTractProfile(coords, radius, color, subdivs, cmap, crange);
    
    % how do I plot the tract profile w/ fibers? necessary?
    %tmp = AFQ_CreateTractProfile();
    %AFQ_RenderFibers(fgcx{ii}, 'tractprofile', tmp);
    
    % create figure and axes label point
    ax1 = axes;
    
    % load a brain slice and display fibers
    fh1 = mbaDisplayBrainSlice(anat, [-1 0 0], ax1);
    fh1 = mbaDisplayConnectome(fgcx{ii}.fibers, fh1, [0 1 0], 'single');
    
    pause;
    
    % this doesn't work displaying things from the right...
    
    % render superfiber w/ mba
    %[ fh1, lh1] = mbaDisplayConnectome(core.fibers, fh1, [0 0 1], 'single', 'jet', 1.0, 1);
    
    % manipulate views
%     view(-90,0); % left view
%     lh1 = camlight('left');
%     zoom(2);
%     ylim([-80, 75]); % crop front / back
%     zlim([ -40, 80]); % crop top / bottom
    
    %feSavefig(fhNum,'verbose','yes', ...
    %    'figName',fig.names{iview}, ...
    %    'figDir',fullFigureOutDir, ...
    %    'figType','jpg');
    
%     % clean out plotting space
%     delete(fh1);
%     delete(ax1);
%     delete(lh1);
%     delete(lightH);
    
end


%% because all of these have to be manipulated manually and differently...

ii = 8;

ax = axes;
%fh = mbaDisplayBrainSlice(anat, [0 10 0], ax); % [saggital coronal axial]
fh = mbaDisplayBrainSlice(anat, [-1 0 0], ax);
[ fh, lh ] = mbaDisplayConnectome(fgcx{ii}.fibers, fh, [0 1 0], 'single');

% the good view rotations
view(0, 15);
view(-45, -15); % good coronal slice


% manipulate views
view(-90,0); % left view
lh = camlight('left');
zoom(2);
ylim([-80, 75]); % crop front / back
zlim([ -40, 80]); % crop top / bottom

%% coronal slice for poster;
ii = 1;
ax = axes;
fh = mbaDisplayBrainSlice(anat, [0 10 0], ax);
[ fh, lh ] = mbaDisplayConnectome(fgcx{ii}.fibers, fh, [0 1 0], 'single');
view(-15, -15);
zoom(2);
% still need to manually pan, but not that bad all things considered

% tract profile - build axes
uplw = minmax(fgst{ii});
mdpt = ((max(fgst{ii})-min(fgst{ii})) / 2) + min(fgst{ii});
ybnd = [round(uplw(1),2) round(mdpt,2) round(uplw(2),2)];
% plot the FA tract profile
figure('Position', [950 605 1350 645]); hold on;
plot(fgst{ii}, 'LineWidth', 4, 'Color', [0 1 0]);
set(gca, 'XTick', [0:50:100], 'YTick', ybnd, 'FontSize', 29);
xlabel('Nodes Along Fiber', 'FontSize', 29);
ylabel('FA Value', 'FontSize', 29);
hold off;
%print(gcf, '-cmyk', '-painters', '-depsc2', '-tiff', '-r500', '-noui', [ 'sfn2016/anatomy/edge-stats-' sprintf('%02d', ii)]);

%% another coronal view
ii = 2;
ax = axes;
fh = mbaDisplayBrainSlice(anat, [0 15 0], ax);
[ fh, lh ] = mbaDisplayConnectome(fgcx{ii}.fibers, fh, [0 1 0], 'single');
view(-15, -15);
zoom(2);

% tract profile - build axes
uplw = minmax(fgst{ii});
mdpt = ((max(fgst{ii})-min(fgst{ii})) / 2) + min(fgst{ii});
ybnd = [round(uplw(1),2) round(mdpt,2) round(uplw(2),2)];
% plot the FA tract profile
figure('Position', [950 605 1350 645]); hold on;
plot(fgst{ii}, 'LineWidth', 4, 'Color', [0 1 0]);
set(gca, 'XTick', [0:50:100], 'YTick', ybnd, 'FontSize', 29);
xlabel('Nodes Along Fiber', 'FontSize', 29);
ylabel('FA Value', 'FontSize', 29);
hold off;

%% sagital slice, motor/premotor
ii = 6;
ax = axes;
fh = mbaDisplayBrainSlice(anat, [-1 0 0], ax);
[ fh, lh ] = mbaDisplayConnectome(fgcx{ii}.fibers, fh, [1 0 0], 'single');
view(-90,0);
zoom(2);
lh = camlight(75, 0);

% tract profile - build axes
uplw = minmax(fgst{ii});
mdpt = ((max(fgst{ii})-min(fgst{ii})) / 2) + min(fgst{ii});
ybnd = [round(uplw(1),2) round(mdpt,2) round(uplw(2),2)];
% plot the FA tract profile
figure('Position', [950 605 1350 645]); hold on;
plot(fgst{ii}, 'LineWidth', 4, 'Color', [1 0 0]);
set(gca, 'XTick', [0:50:100], 'YTick', ybnd, 'FontSize', 29);
xlabel('Nodes Along Fiber', 'FontSize', 29);
ylabel('FA Value', 'FontSize', 29);
hold off;

%% last sagital slice
ii = 8;
ax = axes;
fh = mbaDisplayBrainSlice(anat, [-1 0 0], ax);
[ fh, lh ] = mbaDisplayConnectome(fgcx{ii}.fibers, fh, [0 0 1], 'single');
view(90, 18);

% tract profile - build axes
uplw = minmax(fgst{ii});
mdpt = ((max(fgst{ii})-min(fgst{ii})) / 2) + min(fgst{ii});
ybnd = [round(uplw(1),2) round(mdpt,2) round(uplw(2),2)];
% plot the FA tract profile
figure('Position', [950 605 1350 645]); hold on;
plot(fgst{ii}, 'LineWidth', 4, 'Color', [0 0 1]);
set(gca, 'XTick', [0:50:100], 'YTick', ybnd, 'FontSize', 29);
xlabel('Nodes Along Fiber', 'FontSize', 29);
ylabel('FA Value', 'FontSize', 29);
hold off;



