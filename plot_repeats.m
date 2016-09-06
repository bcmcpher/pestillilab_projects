%% load subject data

% build data structures
[ avg, mat, fns, pval, fdat ] = feMergeRepeats('stn', 'FP', 'prob', '10');

% load roi labels
load('roilabels.mat');

% Frobenius norm - Cesar
d = norm(A(:)-B(:));

% what matrix index to use for plots
mt = 1;

%% figure 1

% figs of interest: 1 5 13 14
% figs of importance: 1 2 3 4 5 6 7 8 13 14

% figure labels
plab = {'Fiber Count', 'Fiber Density', 'Fiber Length', 'Fiber Density x Length', ...
        'Weighted Fiber Count', 'Weighted Fiber Density', 'Weighted Fiber Length', 'Weighted Fiber Density x Length', ...
        'Sum of Weights', 'Weights / Count', 'Weights / Density', 'Weights / Length', ...
        'Strength of Evidence', 'Earth Movers Distance', 'Jeffery''s Divergence', 'Kullback-Leibler'};

% axis labels
pyax = {'Number of Fibers', 'Density of Fibers', 'Length of Fibers', 'Density x Length', ...
        'NZ Number of Fibers', 'NZ Density of Fibers', 'NZ Length of Fibers', 'NZ Density x Length', ...
        'Strength of Evidence', 'Strength of Evidence', 'Strength of Evidence', 'Strength of Evidence', ...
        'Strength of Evidence', 'Strength of Evidence', 'Strength of Evidence', 'Strength of Evidence'};

% loop to create and save figures
for ii = 1:size(avg, 2)
    fh = figure;
    colormap('hot');
    imagesc(avg{ii}.emat.mean);
    axis('square'); axis('equal'); axis('tight');
    title(plab{ii});
    xlabel('FS DK Regions');
    ylabel('FS DK Regions');
    y = colorbar;
    ylabel(y, pyax{ii});
    line([34.5 34.5], [0 68.5], 'Color', [0 0 1]);
    line([0 68.5], [34.5 34.5], 'Color', [0 0 1]);
    line([68.5 0], [68.5 0], 'Color', [0 0 1]);
    
    % save plots
    fname = ['figs/fig' sprintf('%02d', ii) '.png'];
    print(fh, fname, '-dpng');
    close all;
end

% load a deterministic avg for 6 bar plots... 
% ...


fh = figure;
colormap('hot');
imagesc(avg{1}.emat.mean);
axis('square'); axis('equal'); axis('tight');
title('Mean Fiber Count');
xlabel('FS DK Regions');
ylabel('FS DK Regions');
y = colorbar;
ylabel(y, 'Number of Fibers');
line([34.5 34.5], [0 68.5], 'Color', [0 0 1]);
line([0 68.5], [34.5 34.5], 'Color', [0 0 1]);
line([68.5 0], [68.5 0], 'Color', [0 0 1]);

fh = figure;
colormap('hot');
imagesc(avg{5}.emat.mean);
axis('square'); axis('equal'); axis('tight');
title('Mean Non-Zero Fiber Count');
xlabel('FS DK Regions');
ylabel('FS DK Regions');
y = colorbar;
ylabel(y, 'Number of Fibers');
line([34.5 34.5], [0 68.5], 'Color', [0 0 1]);
line([0 68.5], [34.5 34.5], 'Color', [0 0 1]);
line([68.5 0], [68.5 0], 'Color', [0 0 1]);

fh = figure;
colormap('hot');
imagesc(avg{13}.emat.mean);
axis('square'); axis('equal'); axis('tight');
title('Virtual Lesion - SOE');
xlabel('FS DK Regions');
ylabel('FS DK Regions');
y = colorbar;
ylabel(y, 'Strength of Evidence');
line([34.5 34.5], [0 68.5], 'Color', [0 0 1]);
line([0 68.5], [34.5 34.5], 'Color', [0 0 1]);
line([68.5 0], [68.5 0], 'Color', [0 0 1]);

fh = figure;
colormap('hot');
imagesc(avg{14}.emat.mean);
axis('square'); axis('equal'); axis('tight');
title('Virtual Lesion - EMD');
xlabel('FS DK Regions');
ylabel('FS DK Regions');
y = colorbar;
ylabel(y, 'Strength of Evidence');
line([34.5 34.5], [0 68.5], 'Color', [0 0 1]);
line([0 68.5], [34.5 34.5], 'Color', [0 0 1]);
line([68.5 0], [68.5 0], 'Color', [0 0 1]);


%% plot subject data

% plot the repeats individual for a matrix
fh = figure;
for kk = 1:size(mat{mt}, 3)
    subplot(3, 4, kk);
    colormap('hot');
    imagesc(mat{mt}(:,:,kk));
    axis('square'); axis('equal'); axis('tight');
    xlabel('FS DK Regions');
    ylabel('FS DK Regions');
    y = colorbar;
    ylabel(y, 'Strength of Connection');
end

clear fh kk y

%% plot all average matrices

% label matrix
plab = {'Fiber Count', 'Fiber Density', 'Fiber Length', 'Fiber Density x Length', ...
        'Weighted Fiber Count', 'Weighted Fiber Density', 'Weighted Fiber Length', 'Weighted Fiber Density x Length', ...
        'Sum of Weights', 'Weights / Count', 'Weights / Density', 'Weights / Length', ...
        'Strength of Evidence', 'Earth Movers Distance', 'Jeffery''s Divergence', 'Kullback-Leibler'};

% plot figure
fh = figure('Position', [325 25 1550 1175]);
for kk = 1:size(avg, 2)
    subplot(4, 4, kk);
    colormap('hot');
    imagesc(avg{kk}.emat.mean);
    axis('square'); axis('equal'); axis('tight');
    title([plab{kk} ' Mean']);
    xlabel('FS DK Regions');
    ylabel('FS DK Regions');
    y = colorbar;
    ylabel(y, 'Strength of Connection');
end

fh = figure('Position', [325 25 1550 1175]);
for kk = 1:size(avg, 2)
    subplot(4, 4, kk);
    colormap('winter');
    imagesc(avg{kk}.emat.std);
    axis('square'); axis('equal'); axis('tight');
    title([plab{kk}, ' STD']);
    xlabel('FS DK Regions');
    ylabel('FS DK Regions');
    y = colorbar;
    ylabel(y, 'Strength of Connection');
end

%% inspect average figure

fh = figure('Position', [650 525 1000 475]);
ax1 = subplot(1, 2, 1);
colormap(ax1, 'hot');
imagesc(avg{mt}.emat.mean);
axis('square'); axis('equal'); axis('tight');
title('Connectivity');
xlabel('FS DK Regions');
ylabel('FS DK Regions');
y = colorbar;
ylabel(y, 'Strength of Connection');

ax2 = subplot(1, 2, 2);
colormap(ax2, 'winter');
imagesc(avg{mt}.emat.std);
axis('square'); axis('equal'); axis('tight');
title('Standard Deviation');
xlabel('FS DK Regions');
ylabel('FS DK Regions');
y = colorbar;
ylabel(y, 'Strength of Connection');

clear y ax1 ax2

%% the dream (impossible) combined plot

% create data
ldat = log(tril(avg{mt}.emat.mean));
udat = triu(avg{mt}.emat.std);

fh = figure('Position', [1000 550 825 625]);
ax1 = axes;
colormap(ax1, 'hot');
im1 = imagesc(ldat);
%alpha(im1, 0.75);
axis('square'); axis('equal'); axis('tight');
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', [], 'NextPlot', 'add');
title('Average Structural Connectivity and Standard Deviation');
xlabel('FreeSurfer Desikan-Killiany Regions');
y1 = colorbar(ax1, 'Location', 'westoutside');
ylabel(y1, 'Strength of Connection');
ya = colorbar(ax1, 'Location', 'eastoutside', 'visible', 'off');
caxis(ax1, [0 max(max(ldat))]);

ax2 = axes;
axis(ax2, 'off');
colormap(ax2, 'winter');
im2 = imagesc(udat);
%alpha(im2, 0.25);
axis('square'); axis('equal'); axis('tight');
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', [], 'NextPlot', 'add');
title('Average Structural Connectivity and Standard Deviation');
xlabel('FreeSurfer Desikan-Killiany Regions');
y2 = colorbar(ax2, 'eastoutside');
ylabel(y2, 'Standard Deviation of Connection');
yb = colorbar(ax2, 'Location', 'westoutside', 'visible', 'off');
caxis(ax2, [0 max(max(udat))]);
linkaxes([ax1, ax2]);

% even bother?
line([34.5 34.5], [0 68.5], 'Color', [0 0 1]);
line([0 68.5], [34.5 34.5], 'Color', [0 0 1]);
line([68.5 0], [68.5 0], 'Color', [0 0 1]);

%% practical combined plot

% combine data
pmat = tril(avg{mt}.emat.mean) + triu(avg{mt}.emat.std);

% build combined colormap (?)
% not really happening...

% plot figure
fh = figure('Position', [850 400 850 600]);
colormap('hot');
imagesc(log(pmat));
axis('square'); axis('equal'); axis('tight');
title('Average Structural Connectivity and Standard Deviation');
xlabel('FreeSurfer Desikan-Killiany Regions');
ylabel('FreeSurfer Desikan-Killiany Regions');
y1 = colorbar();
ylabel(y1, 'Log Strength of Connection');
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
line([34.5 34.5], [0 68.5], 'Color', [1 1 1], 'LineWidth', 1.5);
line([0 68.5], [34.5 34.5], 'Color', [1 1 1], 'LineWidth', 1.5);
line([68.5 0], [68.5 0], 'Color', [1 1 1], 'LineWidth', 1.5);

%% bar plots


br = avg{mt}.fns.nrm.mean.str;
ebr = avg{mt}.fns.nrm.std.str;

fh = figure('Position', [530 700 1720 550]);
hold on;
bar(br);
title('Average of Repeats - Connection Strength');
ylabel('Strength');
set(gca, 'XTick', 1:length(roiLabel), 'XTickLabel', roiLabel, 'XTickLabelRotation', 45);
errorbar(br, ebr, '.');

%% spider plots

% merge data from objects
dat = [ fns{mt}.nrm.smwrld'; fns{mt}.nrm.mcoef'; fns{mt}.nrm.trans'; ...
        fns{mt}.nrm.glEff'; fns{mt}.nrm.stat'; fns{mt}.nrm.mod' ];

% axis labels
lglb = {'Small Worldness', 'Mean Clustering Coefficient', 'Transitivity', ... 
        'Global Efficiency', 'Community Structure', 'Modularity'};

% line labels (legend)
legLab = {'rep01', 'rep02', 'rep03', 'rep04', 'rep05', 'rep06', 'rep07', 'rep08', 'rep09', 'rep10'};
    
% make plot
titl = 'Repeated Global Network Measures';
spider(dat, titl, [], lglb, legLab);

%% bar plots of global network measures

dat = [ fns{mt}.nrm.smwrld'; fns{mt}.nrm.mcoef'; fns{mt}.nrm.trans'; ...
        fns{mt}.nrm.glEff'; fns{mt}.nrm.stat'; fns{mt}.nrm.mod' ];
    
br = mean(dat, 2);
ebr = std(dat, 0, 2);

bar(br); hold on;
title('Average of Repeats - Global Measures');
set(gca, 'XTick', 1:length(lglb), 'XTickLabel', lglb, 'XTickLabelRotation', 45);
errorbar(br, ebr, '.');






%% online example of 2 colormaps

%create data
[X,Y,Z] = peaks(30);
Zprime = del2(Z);
contourmin = min(Z(:));
contourmax = max(Z(:));
pcolormin = min(Zprime(:));
pcolormax = max(Zprime(:));

%create figure and store handle
hF = figure;
%create axes for pcolor and store handle
hAxesP = axes;
%set colormap for pcolor axes
colormap(hAxesP,gray);

%plot pcolor for gradient
pcolorPlot = pcolor(X,Y,Zprime);
set(pcolorPlot,'FaceColor','interp','EdgeColor','interp');
%create color bar and set range for color
cbP = colorbar(hAxesP,'Location','westoutside');
caxis(hAxesP,[pcolormin pcolormax]);

%create axes for the countourm axes
hAxesCM = axes;
%set visibility for axes to 'off' so it appears transparent
axis(hAxesCM,'off')
%set colormap for contourm axes
colormap(hAxesCM,cool);

%plot contourm
contourmPlot = contourm(X,Y,Z,20);
%create color bar and set range for color
cbCM = colorbar(hAxesCM,'Location','eastoutside');
caxis(hAxesCM,[contourmin contourmax]);

%link the two overlaying axes so they match at all times to remain accurate
linkaxes([hAxesP,hAxesCM]);





% create summary plot
plab = {'Fiber Count', 'Fiber Density', 'Fiber Length', 'Fiber Density x Length', ...
        'Weighted Fiber Count', 'Weighted Fiber Density', 'Weighted Fiber Length', 'Weighted Fiber Density x Length', ...
        'Sum of Weights', 'Weights / Count', 'Weights / Density', 'Weights / Length', ...
        'Strength of Evidence', 'Earth Movers Distance', 'Jeffery''s Divergence', 'Kullback-Leibler'};

fh = figure;
for kk = 1:size(emat, 3)
    subplot(4, 4, kk);
    colormap('hot');
    imagesc(log(emat(:,:,kk)));
    title(plab{kk});
    xlabel('FS DK Regions');
    ylabel('FS DK Regions');
    axis('square'); axis('equal'); axis('tight');
    set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
    y = colorbar;
    ylabel(y, 'Strength of Connection');
end
