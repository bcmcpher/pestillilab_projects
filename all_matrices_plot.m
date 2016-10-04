%% plots for paper
% Brent McPherson
% 20160926
%

%% load and configure data for plotting

% load averaged / individual / statistics for 1 subject across repeats for
% 2x2 tracking/lmax parameter combinations
[ detr02, idetr02, sdetr02 ] = feMergeRepeats('hcp', '105115', 'detr', '2');
[ detr10, idetr10, sdetr10 ] = feMergeRepeats('hcp', '105115', 'detr', '10');
[ prob02, iprob02, sprob02 ] = feMergeRepeats('hcp', '105115', 'prob', '2');
[ prob10, iprob10, sprob10 ] = feMergeRepeats('hcp', '105115', 'prob', '10');

%% plot labeled matrix

% isolate data set for plotting here
% log transform if desired
mat = log(idetr02{1}(:, :, 3));

% plot connectivity matrix w/ labels
figure('Position', [700 430 770 630]);
colormap('hot');
imagesc(mat);
axis('square'); axis('equal'); axis('tight');
title('Probabilistic Network');
xlabel('FS DK Regions');
ylabel('FS DK Regions');
y = colorbar;
%ylabel(y, 'Number of Streamlines');
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1]);
line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1]);
line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1]);

%% plot plain matrix

% isolate data set for plotting here
% log transform if desired
mat = log(idetr02{1}(:, :, 3));

% plot connectivity matrix w/o labels
figure('Position', [700 430 770 630]);
colormap('hot');
imagesc(mat);
axis('square'); axis('equal'); axis('tight');
title('Deterministic Network');
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1]);
line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1]);
line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1]);

%% save a bunch of unlabeled average matrices

% global labels
mInd = [1, 2, 3, 4, 5, 6, 7, 8, 13, 14];
mLabel = {'Fiber Count', 'Fiber Density', 'Fiber Length', 'Fiber Len*Den', ...
    'NonZero FC', 'NonZero FD', 'NonZero FL', 'NonZero FDFL', 'SOE', 'EMD'};
mFile = {'fcnt', 'fden', 'flen', 'flnd', 'nzfc', 'nzfd', 'nzfl', 'nzfx', 'lsoe', 'lemd'};

for ii = 1:length(mInd)
    
    % grab data
    mat = log(detr02{mInd(ii)}.emat.mean);
    
    % generic plot
    figure('Position', [700 430 770 630]);
    colormap('hot');
    imagesc(mat);
    axis('square'); axis('equal'); axis('tight');
    set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
    line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1]);
    line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1]);
    line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1]);
    
    % reduced margins in saved figures
    set(gca,'position',[0 0 1 1],'units','normalized');
    
    % save plot    
    fnam = ['detr02_' mFile{ii} '_avg_log_matrix.png'];
    fout = fullfile('/N/dc2/projects/lifebid/HCP/Brent/cogs610/matlab/figs', fnam);
    saveas(gcf, fout, 'png');
    close all;
    
end

%% save a bunch of unlabeled individual matrices

% individual matrices saved

% detr02: 3, 6, 10
individ = [3, 6, 10];

% prob10: 1, 6, 10
individ = [1, 6, 10];

% global labels
mInd = [1, 2, 3, 4, 5, 6, 7, 8, 13, 14];
mLabel = {'Fiber Count', 'Fiber Density', 'Fiber Length', 'Fiber Len*Den', ...
    'NonZero FC', 'NonZero FD', 'NonZero FL', 'NonZero FDFL', 'SOE', 'EMD'};
mFile = {'fcnt', 'fden', 'flen', 'flnd', 'nzfc', 'nzfd', 'nzfl', 'nzfx', 'lsoe', 'lemd'};

for ii = 1:length(mInd)
    for jj = 1:length(individ)
        
        % grab data
        mat = log(iprob10{mInd(ii)}(:, :, individ(jj)));
        
        % generic plot
        figure('Position', [700 430 770 630]);
        colormap('hot');
        imagesc(mat);
        axis('square'); axis('equal'); axis('tight');
        set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
        line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1]);
        line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1]);
        line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1]);
        
        % reduce margins in saved figure
        set(gca,'position',[0 0 1 1],'units','normalized');
        
        % save plot
        idout = sprintf('%02d', individ(jj));
        fnam = ['prob10_' mFile{ii} '_rep' idout '_log_matrix.png'];
        fout = fullfile('/N/dc2/projects/lifebid/HCP/Brent/cogs610/matlab/figs', fnam);
        saveas(gcf, fout, 'png');
        close all;
        
    end
end

