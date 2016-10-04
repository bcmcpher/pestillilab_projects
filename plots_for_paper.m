%% plots for paper
% Brent McPherson
% 20160613
%

%% Figure 1

% load data
f01_detr = feMergeRepeats('hcp', '105115', 'detr', '2');
f01_prob = feMergeRepeats('hcp', '105115', 'prob', '10');
indx = 2;

% raw data

figure('Position', [700 430 770 630]);
colormap('hot');
imagesc(f01_detr{indx}.emat.mean);
axis('square'); axis('equal'); axis('tight');
title('Deterministic Network');
xlabel('FS DK Regions');
ylabel('FS DK Regions');
y = colorbar;
ylabel(y, 'Number of Streamlines');
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1]);
line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1]);
line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1]);

figure('Position', [700 430 770 630]);
colormap('hot');
imagesc(f01_prob{indx}.emat.mean);
axis('square'); axis('equal'); axis('tight');
title('Probabilistic Network');
xlabel('FS DK Regions');
ylabel('FS DK Regions');
y = colorbar;
ylabel(y, 'Number of Streamlines');
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1]);
line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1]);
line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1]);

% log data

figure('Position', [700 430 770 630]);
colormap('hot');
imagesc(log(f01_detr{indx}.emat.mean));
axis('square'); axis('equal'); axis('tight');
title('Deterministic Network');
xlabel('FS DK Regions');
ylabel('FS DK Regions');
y = colorbar;
ylabel(y, 'Log Number of Streamlines');
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1]);
line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1]);
line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1]);

figure('Position', [700 430 770 630]);
colormap('hot');
imagesc(log(f01_prob{indx}.emat.mean));
axis('square'); axis('equal'); axis('tight');
title('Probabilistic Network');
xlabel('FS DK Regions');
ylabel('FS DK Regions');
y = colorbar;
ylabel(y, 'Log Number of Streamlines');
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1]);
line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1]);
line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1]);

% testing template

indx = 1;
figure('Position', [700 430 770 630]);
colormap('hot');
imagesc(f01_prob{indx}.emat.mean > 0);
axis('square'); axis('equal'); axis('tight');
%title('Probabilistic Network');
%xlabel('FS DK Regions');
%ylabel('FS DK Regions');
%y = colorbar;
%ylabel(y, 'Number of Streamlines');
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1]);
line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1]);
line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1]);

% multiplot w/ std

mod = 1;
figure;

f1 = subplot(2, 2, 1);
colormap(f1, 'hot');
imagesc(f01_detr{mod}.emat.mean);
axis('square'); axis('equal'); axis('tight');
title('Deterministic Network');
xlabel('FS DK Regions');
ylabel('FS DK Regions');
y = colorbar;
ylabel(y, 'Number of Streamlines');
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1]);
line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1]);
line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1]);

f2 = subplot(2, 2, 2);
colormap(f2, 'winter');
imagesc(f01_detr{mod}.emat.std);
axis('square'); axis('equal'); axis('tight');
title('Deterministic Standard Deviation');
xlabel('FS DK Regions');
ylabel('FS DK Regions');
y = colorbar;
ylabel(y, 'Number of Streamlines');
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
line([34.5 34.5], [0.5 68.5], 'Color', [1 0 0]);
line([0.5 68.5], [34.5 34.5], 'Color', [1 0 0]);
line([68.5 0.5], [68.5 0.5], 'Color', [1 0 0]);

f3 = subplot(2, 2, 3);
colormap(f3, 'hot');
imagesc(f01_prob{mod}.emat.mean);
axis('square'); axis('equal'); axis('tight');
title('Probabilistic Network');
xlabel('FS DK Regions');
ylabel('FS DK Regions');
y = colorbar;
ylabel(y, 'Number of Streamlines');
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1]);
line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1]);
line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1]);

f4 = subplot(2, 2, 4);
colormap(f4, 'winter');
imagesc(f01_prob{mod}.emat.std);
axis('square'); axis('equal'); axis('tight');
title('Probabilistic Standard Deviation');
xlabel('FS DK Regions');
ylabel('FS DK Regions');
y = colorbar;
ylabel(y, 'Number of Streamlines');
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
line([34.5 34.5], [0.5 68.5], 'Color', [1 0 0]);
line([0.5 68.5], [34.5 34.5], 'Color', [1 0 0]);
line([68.5 0.5], [68.5 0.5], 'Color', [1 0 0]);

%% Figure 2

% recreate LiFE paper figures

%% Figure 3

mInd = [1, 5, 13, 14];
mLabel = {'Fiber Count', 'NonZero FC', 'SOE', 'EMD'};

mInd = [1, 2, 3, 4, 5, 6, 7, 8, 13, 14];
mLabel = {'Fiber Count', 'Fiber Density', 'Fiber Length', 'Fiber Len*Den', ...
    'NonZero FC', 'NonZero FD', 'NonZero FL', 'NonZero FDFL', 'SOE', 'EMD'};

% pull network wide statistics from output
glbVars = {'smwrld', 'mcoef', 'trans', 'glEff', 'stat', 'mod'};
glbLabs = {'Small Worldness', 'Mean Clustering Coefficient', 'Transitivity', ... 
           'Global Efficiency', 'Community Structure', 'Modularity'};

% path stems to pull variables from
dmstem = 'f01_detr{mInd(ii)}.fns.nrm.mean.';
dsstem = 'f01_detr{mInd(ii)}.fns.nrm.std.';
pmstem = 'f01_prob{mInd(ii)}.fns.nrm.mean.';
psstem = 'f01_prob{mInd(ii)}.fns.nrm.std.';

% create data to plot
for ii = 1:length(mInd)
    for jj = 1:length(glbVars)
        mdDat(jj, ii) = eval([dmstem glbVars{jj}]); 
        sdDat(jj, ii) = eval([dsstem glbVars{jj}]); 
        mpDat(jj, ii) = eval([pmstem glbVars{jj}]); 
        spDat(jj, ii) = eval([psstem glbVars{jj}]); 
    end
end

% plot global network data
for ii = 1:length(glbVars)
    figure;
    bar(mdDat(ii,:));
    hold on;
    errorbar(mdDat(ii,:), sdDat(ii,:), '.');
    set(gca, 'XTick', 1:length(mInd), 'XTickLabel', mLabel, 'XTickLabelRotation', 45);
    hold off;
end

for ii = 1:length(glbVars)
    figure;
    bar(mpDat(ii,:));
    hold on;
    errorbar(mpDat(ii,:), spDat(ii,:), '.');
    set(gca, 'XTick', 1:length(mInd), 'XTickLabel', mLabel, 'XTickLabelRotation', 45);
    hold off;
end

% all in 1
figure;
for ii = 1:length(glbVars)
    
    % deterministic
    subplot(2, 6, ii);
    bar(mdDat(ii,:));
    title(glbLabs{ii});
    hold on;
    errorbar(mdDat(ii,:), sdDat(ii,:), '.');
    set(gca, 'XTick', 1:length(mInd), 'XTickLabel', mLabel, 'XTickLabelRotation', 45);
    hold off;
    
    % probabilistic
    subplot(2, 6, ii+length(glbVars));
    bar(mdDat(ii,:));
    title(glbLabs{ii});
    hold on;
    errorbar(mdDat(ii,:), sdDat(ii,:), '.');
    set(gca, 'XTick', 1:length(mInd), 'XTickLabel', mLabel, 'XTickLabelRotation', 45);
    hold off;
    
end



figure;
bar(mdDat);
hold on;
title('Smallworldness');
errorbar(mdDat, sdDat, '.');
set(gca, 'XTick', 1:length(mInd), 'XTickLabel', mLabel, 'XTickLabelRotation', 45);
hold off;

figure;
bar(mpDat);
hold on;
title('Smallworldness');
errorbar(mpDat, spDat, '.');
set(gca, 'XTick', 1:length(mInd), 'XTickLabel', mLabel, 'XTickLabelRotation', 45);
hold off;

clear ii mInd mLabel mdDat sdDat mpDat spDat

% Indices of summary data = 1, 2, 3, 4, 5, 6, 7, 8, 13, 14
% 1, 5, 13, 14

% plot multiple metrics before LiFE (top row) and after LiFE removes
% non-zero fibers (bottom row)
figure('Position', [140 270 1800 960]);
iter = 1;
for ii = 1:2
    for jj = 1:4
        subplot(2, 4, iter);
        colormap('hot');
        imagesc(f01_detr{iter}.emat.mean);
        axis('square'); axis('equal'); axis('tight');
        %title('Deterministic Network');
        %xlabel('FS DK Regions');
        %ylabel('FS DK Regions');
        y = colorbar;
        %ylabel(y, 'Number of Streamlines');
        set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
        line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1]);
        line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1]);
        line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1]);
        iter = iter + 1;
    end
end

% probabilistic of above
figure;
iter = 1;
for ii = 1:2
    for jj = 1:4
        subplot(2, 4, iter);
        colormap('hot');
        imagesc(f01_prob{iter}.emat.mean);
        axis('square'); axis('equal'); axis('tight');
        title('Probabilistic Network');
        xlabel('FS DK Regions');
        ylabel('FS DK Regions');
        y = colorbar;
        ylabel(y, 'Number of Streamlines');
        set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
        line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1]);
        line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1]);
        line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1]);
        iter = iter + 1;
    end
end

labs = {'Number of Streamlines', 'Streamline Density', 'Streamline Length', 'Streamline Den*Len', ...
    'Number of Streamlines', 'Streamline Density', 'Streamline Length', 'Streamline Den*Len'};

saveName = {'count', 'density', 'length', 'denlen', 'nzcnt', 'nzden', 'nzlen', 'nzdenlen'};
for ii = 1:8
    figure('Position', [700 430 770 630]);
    colormap('hot');
    imagesc(detr02{ii}.emat.mean);
    axis('square'); axis('equal'); axis('tight');
    %y = colorbar;
    %ylabel(y, labs{ii});
    set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
    line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1]);
    line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1]);
    line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1]);
    saveas(gcf, ['~/Desktop/figs_fall2016/unlabeled_detr_02_' saveName{ii} '.png']);
    close(gcf);
end


%% show average and individual plots

% load multiple data sets
[ detr02, idetr02, sdetr02 ] = feMergeRepeats('hcp', '105115', 'detr', '2');
[ detr10, idetr10, sdetr10 ] = feMergeRepeats('hcp', '105115', 'detr', '10');
[ prob02, iprob02, sprob02 ] = feMergeRepeats('hcp', '105115', 'prob', '2');
[ prob10, iprob10, sprob10 ] = feMergeRepeats('hcp', '105115', 'prob', '10');

%% Smith scatter plot

% compare detr / prob

% up front values
edge_index = triu(true(size(zeros(68,68))));
indx = 13;

% create x-axes
xa = detr02{1}.emat.mean(edge_index);
xb = detr02{3}.emat.mean(edge_index);
xc = detr02{5}.emat.mean(edge_index);
xd = detr02{7}.emat.mean(edge_index);

% pull unique node values
% USE VARIANCE INSTEAD OF STANDARD DEVIATION (FANO FACTOR)
mu1 = detr02{indx}.emat.mean(edge_index);
sg1 = detr02{indx}.emat.std(edge_index);

mu2 = prob02{indx}.emat.mean(edge_index);
sg2 = prob02{indx}.emat.std(edge_index);

% find means greater than 1
gt1 = mu1 > 1;
gt2 = mu2 > 1;

% create average
cove1 = sg1 ./ mu1;
cove2 = sg2 ./ mu2;

% weighted mean coefficient of variation
wmcove1 = sum(2.^(log(mu1(gt1))) .* cove1(gt1)) / sum(2.^(log(mu1(gt1)))); 
wmcove2 = sum(2.^(log(mu2(gt2))) .* cove2(gt2)) / sum(2.^(log(mu2(gt2)))); 

% basic plot
figure('Position', [450 650 1000 450]); hold on;
scatter(log(xc(gt1)), cove1(gt1), 10, [ 1 0 0 ], 'filled');
scatter(log(xc(gt2)), cove2(gt2), 10, [ 0 0 1 ], 'filled');
title('Consistency of SOE Relative to Streamline Density');
xlabel('Log Fiber Count');
ylabel('CoV - Coefficient of Variation');
rln1 = refline([ 0 wmcove1 ]);
rln1.Color = [ 1 0 0 ];
rln2 = refline([ 0 wmcove2 ]);
rln2.Color = [ 0 0 1 ];
legend('Deterministic', 'Probabilistic');

% compare count / nz-count / emd

% up front values
edge_index = triu(true(size(zeros(68,68))));

% create x-axes
xa = detr02{1}.emat.mean(edge_index);
xb = detr02{5}.emat.mean(edge_index);

xc = detr02{3}.emat.mean(edge_index);
xd = detr02{7}.emat.mean(edge_index);

% pull unique node values
mu1 = detr02{1}.emat.mean(edge_index);
sg1 = detr02{1}.emat.std(edge_index);

mu2 = prob02{5}.emat.mean(edge_index);
sg2 = prob02{5}.emat.std(edge_index);

mu3 = prob02{14}.emat.mean(edge_index);
sg3 = prob02{14}.emat.std(edge_index);

% find means greater than 1
gt1 = mu1 > 1;
gt2 = mu2 > 1;
gt3 = mu3 > 1;

% create average
cove1 = sg1 ./ mu1;
cove2 = sg2 ./ mu2;
cove3 = sg3 ./ mu3;

% weighted mean coefficient of variation
wmcove1 = sum(2.^(log(mu1(gt1))) .* cove1(gt1)) / sum(2.^(log(mu1(gt1)))); 
wmcove2 = sum(2.^(log(mu2(gt2))) .* cove2(gt2)) / sum(2.^(log(mu2(gt2)))); 
wmcove3 = sum(2.^(log(mu3(gt3))) .* cove3(gt3)) / sum(2.^(log(mu3(gt3)))); 

% basic plot
figure('Position', [450 650 1000 450]); hold on;
scatter(log(xb(gt1)), cove1(gt1), 10, [ 1 0 0 ], 'filled');
scatter(log(xb(gt2)), cove2(gt2), 10, [ 0 1 0 ], 'filled');
scatter(log(xb(gt3)), cove3(gt3), 10, [ 0 0 1 ], 'filled');
title('Comparison of CoV Estimates Between Methods Relative to Streamline Density');
xlabel('Log Non-Zero Fiber Count');
ylabel('CoV - Coefficient of Variation');
rln1 = refline([ 0 wmcove1 ]);
rln1.Color = [ 1 0 0 ];
rln2 = refline([ 0 wmcove2 ]);
rln2.Color = [ 0 1 0 ];
rln3 = refline([ 0 wmcove3 ]);
rln3.Color = [ 0 0 1 ];
legend('Fiber Count', 'Non-Zero Fiber Count', 'EMD');

% 3d scatter plot w/ length

% up front values
edge_index = triu(true(size(zeros(68,68))));
indx = 13;

% create x-axes
xa = detr02{5}.emat.mean(edge_index);
xb = detr02{7}.emat.mean(edge_index);

% pull unique node values
mu1 = detr02{indx}.emat.mean(edge_index);
sg1 = detr02{indx}.emat.std(edge_index);

mu2 = prob02{indx}.emat.mean(edge_index);
sg2 = prob02{indx}.emat.std(edge_index);

% find means greater than 1 (mu# > 1)
% this only plots edges w/ greater than 0 log-scale fibers
gt1 = log(xa) > 0;
gt2 = log(xa) > 0;

% create average
cove1 = sg1 ./ mu1;
cove2 = sg2 ./ mu2;

% weighted mean coefficient of variation
wmcove1 = sum(2.^(log(mu1(gt1))) .* cove1(gt1)) / sum(2.^(log(mu1(gt1)))); 
wmcove2 = sum(2.^(log(mu2(gt2))) .* cove2(gt2)) / sum(2.^(log(mu2(gt2)))); 

% basic plot
figure; hold on;
scatter3(log(xa(gt1)), xb(gt1), cove1(gt1), 10, [ 1 0 0 ], 'filled');
scatter3(log(xa(gt2)), xb(gt2), cove2(gt2), 10, [ 0 0 1 ], 'filled');
title('Number of Non-Zero Fibers X Non-Zero Fiber Length X CoV');
xlabel('Log Non-Zero Fiber Count');
ylabel('Non-Zero Fiber Length');
zlabel('CoV - Coefficient of Variation')
legend('Deterministic', 'Probabilistic');

%% ICC of matrix

% create data for count, non-zero count, and EMD
for ii = 1:10
    tmpc = idetr02{1}(:,:,ii);
    tmpn = idetr02{1}(:,:,ii);
    tmpl = idetr02{1}(:,:,ii);
    count(:,ii) = reshape(tmpc(edge_index), 1, [])';
    nzcnt(:,ii) = reshape(tmpn(edge_index), 1, [])';
    lfemd(:,ii) = reshape(tmpl(edge_index), 1, [])';
end

clear tmpc tmpn tmpl

% full stats on ICC
[ icc_r, icc_lb, icc_ub, icc_f, icc_df1, icc_df2, icc_p ] = ICC(count, '1-1', 0.01, 0.5);
clear icc_r icc_lb icc_ub icc_f icc_df1 icc_df2 icc_p

% run ICC estimates
icc_count = ICC(count, '1-1', 0.001, 1);
icc_nzcnt = ICC(nzcnt, '1-1', 0.001, 1);
icc_lfemd = ICC(lfemd, '1-1', 0.001, 1);

%% matrix plots

% plot average fiber count and 3 example individual data sets
fh = figure('Position', [320 350 1100 800]);

% big matrix
subplot(3, 4, [1, 2, 3, 5, 6, 7, 9, 10, 11]);
colormap('hot');
imagesc(log(detr02{1}.emat.mean));
axis('square'); axis('equal'); axis('tight');
title('Average Deterministic Network');
%xlabel('FS DK Regions');
%ylabel('FS DK Regions');
y = colorbar;
ylabel(y, 'Log Number of Streamlines');
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1]);
line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1]);
line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1]);

subplot(3, 4, 4);
colormap('hot');
imagesc(log(edetr02{1}(:,:,1)));
axis('square'); axis('equal'); axis('tight');
title('Example Deterministic Networks');
%xlabel('FS DK Regions');
%ylabel('FS DK Regions');
%y = colorbar;
%ylabel(y, 'Number of Streamlines');
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1]);
line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1]);
line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1]);

subplot(3, 4, 8);
colormap('hot');
imagesc(log(edetr02{1}(:,:,3)));
axis('square'); axis('equal'); axis('tight');
%title('Example Deterministic Network');
%xlabel('FS DK Regions');
%ylabel('FS DK Regions');
%y = colorbar;
%ylabel(y, 'Number of Streamlines');
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1]);
line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1]);
line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1]);

subplot(3, 4, 12);
colormap('hot');
imagesc(log(edetr02{1}(:,:,5)));
axis('square'); axis('equal'); axis('tight');
%title('Example Deterministic Network');
%xlabel('FS DK Regions');
%ylabel('FS DK Regions');
%y = colorbar;
%ylabel(y, 'Number of Streamlines');
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1]);
line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1]);
line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1]);

% look at all 10 repeats- no great difference
fh = figure;
for ii = 1:10
    subplot(3, 4, ii)
    colormap('hot');
    imagesc(log(iprob10{1}(:,:,ii)));
    axis('square'); axis('equal'); axis('tight');
    title('Average Deterministic Network');
    set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
    line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1]);
    line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1]);
    line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1]);
end

%% plot statistics

% across x axis: value of global metric w/ error bars for single subject
% across y axis: jittered subject index
% point shapes : tracking method
% point colors : data set

hcp = {'105115' '110411' '111312' '113619'};
stn = {'FP' 'HT' 'KK' 'MP'};

% create plotting data
for ii = 1:4
    
    % load a subject
    [ hcpVal, ~, hcpRaw ] = feMergeRepeats('hcp', hcp{ii}, 'prob', '10');
    [ stnVal, ~, stnRaw ] = feMergeRepeats('stn', stn{ii}, 'prob', '10');
    
    % grab value
    hcpDat(ii, 1) = hcpVal{1}.fns.nrm.mean.smwrld;
    stnDat(ii, 1) = stnVal{1}.fns.nrm.mean.smwrld;
    
    % pull intermediary values to plot standard error of mean (SEM)
    ht1 = hcpVal{1}.fns.nrm.std.smwrld;
    ht2 = sqrt(size(hcpRaw{1}.nrm.smwrld, 1));
    
    st1 = stnVal{1}.fns.nrm.std.smwrld;
    st2 = sqrt(size(stnRaw{1}.nrm.smwrld, 1));
    
    % create SEM
    ht3 = ht1 / ht2;
    st3 = st1 / st2;
    
    % put error values in matrix
    hcpDat(ii, 2) = hcpDat(ii, 1) - ht3;
    hcpDat(ii, 3) = hcpDat(ii, 1) + ht3;

    stnDat(ii, 2) = stnDat(ii, 1) - st3;
    stnDat(ii, 3) = stnDat(ii, 1) + st3;
    
    % clean up
    clear ht1 ht2 ht3 hcpVal hcpRaw st1 st2 st3 stnVal stnRaw
    
end

%% loop across parameters

% load data - edit function to add different parameters
dat = Gen_data_for_plots();

% copy of parameter set to parse data object
lparam_set = {'2', '4', '6', '8', '10', '12'}; 
alg_set = {'detr', 'prob'};
glbVars = {'smwrld', 'mcoef', 'trans', 'glEff', 'stat', 'mod'};
nodVars = {'str', 'deg', 'btw', 'eigv', 'pcoef', 'modz'};
connMod = {2, 6, 14};

% labels for parameters
glbLabel = {'Small Worldness', 'Mean Clustering Coefficient', 'Transitivity', ... 
           'Global Efficiency', 'Community Structure', 'Modularity'};
nodLabel = {'Strength', 'Degree', 'Betweenness', 'Eigenvector Centrality', ...
            'Partial Coefficiency', 'Modularity'};
modLabel = {'Fiber Density', 'NonZero Density', 'EMD'};

% build color map 
%tcol = parula(8);
c1 = colormap(parula(64));
c2 = colormap(autumn(64));
tcol = [ c2([49 40 28 4],:); c1([4 15 27 31],:) ];
close all;

% jitter for subject and parameter, combined in loop
sjit = linspace(-0.30, 0.30, length(dat));
ljit = linspace(-0.025, 0.025, length(lparam_set));

% global parameter value to plot

for imod = 1:length(glbLabel)
    
    modIndex = imod;
    
    figure('Position', [948 420 560 825]); hold on;
    ylim([ 0.50 3.50 ]);
    % for every subject loaded in a data structure
    for ii = 1:length(dat)
        % for every model
        for jj = 1:length(connMod)
            % for every tracking type
            for kk = 1:length(alg_set)
                alg = alg_set{kk};
                switch alg
                    case 'detr'
                        shape = 'o';
                    case 'prob'
                        shape = 's';
                end
                % for every lmax
                for ll = 1:length(lparam_set)
                    param = lparam_set{ll};
                    % combine subject and lmax jitter for more obvious offsets
                    jitter = sjit(ii) + ljit(ll);
                    dlab = [alg '_' param];
                    pt = eval(['dat{' num2str(ii) '}.' dlab '.glob{' num2str(jj) '}(' num2str(modIndex) ',:)' ]);
                    plot(pt(1), jj+jitter, shape, 'Color', tcol(ii,:), 'MarkerFaceColor', tcol(ii,:), 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerSize', 10);
                    plot([pt(2) pt(3)], [jj jj]+jitter,  'LineWidth', 2, 'Color', [0 0 0]);
                end
            end
            % add tenosr points here
        end
    end
    % plot labels
    xlabel(glbLabel{modIndex});
    set(gca, 'YTick', [1 2 3], 'YTickLabel', {});
    if imod == 1
        ylabel('Connectome Models');
        set(gca, 'YTick', [1 2 3], 'YTickLabel', modLabel);
    end
    hold off;
    
end



% combine data sets
pdat = [hcpDat; stnDat];


% catch values for plotting
jity = [0.85; 0.95; 1.05; 1.15]; % y value to be jittered for subject plotting
jity = [jity; jity+1];

%val = prob10{1}.fns.nrm.mean.smwrld; % value on y-axis to be plotted
%err = prob10{1}.fns.nrm.std.smwrld / sqrt(size(sprob10{1}.nrm.smwrld, 1)); % distance from y value for error bar

figure;
xlim([0 4]);
ylim([0 2.5]);
hold on;

% plot average point
scatter(pdat(:,1), jity, 100, [ 0 0 1 ], 'o', 'LineWidth', 2);

% plot horizontal error bar
%plot([val-err val+err], [jity jity],  'LineWidth', 2, 'Color', [0 0 1]);
for ii = 1:8
    plot([pdat(ii,2) pdat(ii,3)], [jity(ii) jity(ii)],  'LineWidth', 2, 'Color', [0 0 1]);
end

% labeling - need to scale figure too
title('Small-Worldness Estimate');
xlabel('Small-Worldness');
ylabel('Subject Index');
        