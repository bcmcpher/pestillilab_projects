%% SFN poster matrices
% Brent McPherson
% 20161007
%

% will go through all subjects
% 2 HCP, 2 STN avg lmax12, lmax8, lmax4, lmax2 prob / detr + tens

% other 2 HCP, 2 STN all 10 individual repeats for lmax10 prob

% load average counts
[ detr02, idetr02, sdetr02 ] = feMergeRepeats('hcp', '105115', 'detr', '2');
[ detr10, idetr10, sdetr10 ] = feMergeRepeats('hcp', '105115', 'detr', '10');
[ prob02, iprob02, sprob02 ] = feMergeRepeats('hcp', '105115', 'prob', '2');
[ prob10, iprob10, sprob10 ] = feMergeRepeats('hcp', '105115', 'prob', '10');
[ tens02, itens02, stens02 ] = feMergeRepeats('hcp', '105115', 'tens', '2');


%% development of plots

% isolate data set for plotting here
% log transform if desired
mat = log(tens02{1}.emat.mean);
out = 'tens10';

% plot connectivity matrix w/ labels
figure('Position', [700 430 770 630]);
colormap('hot');
imagesc(mat);
axis('square'); axis('equal'); axis('tight');
%title('Probabilistic Network');
%xlabel('FS DK Regions');
%ylabel('FS DK Regions');
y = colorbar('FontSize', 29, 'Ticks', [0, 4, 8], 'TickLabels', {0, 4, 8});
caxis([0 8]);
ylabel(y, 'Log Number of Streamlines');
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1], 'LineWidth', 2);
line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1], 'LineWidth', 2);
line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1], 'LineWidth', 2);
print(gcf, '-cmyk', '-painters', '-depsc2', '-tiff', '-r500', '-noui', [ 'sfn2016/labeled_' out]);
%saveas(gcf, [ 'sfn2016/labeled_' out '.eps'], 'epsc');
%saveas(gcf, [ 'sfn2016/labeled_' out '.png'], 'png');

% plot connectivity matrix w/o labels
figure('Position', [700 430 770 630]);
colormap('hot');
imagesc(mat);
axis('square'); axis('equal'); axis('tight');
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
caxis([0 8]);
line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1], 'LineWidth', 2);
line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1], 'LineWidth', 2);
line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1], 'LineWidth', 2);
%print(gcf, '-cmyk', '-painters', '-depsc2', '-tiff', '-r500', '-noui', [ 'sfn2016/test_' out ]); 
print(gcf, '-cmyk', '-painters', '-dpng', '-r500', '-noui', [ 'sfn2016/test_' out ]); 

%% load up 2 of 3T data sets across multiple parameters and save down figures

dset = {'hcp', 'hcp', 'stn', 'stn'};
subj = {'110411', '111312', 'FP', 'MP'};
%dset = {'stn', 'stn'};
%subj = {'FP', 'MP'};
lmax = {'2', '4', '8', '12'};
trck = {'prob', 'detr'};
cscl = {[0 15], [0 15], [0 1], [0 1]};

for ii = 1:length(subj)
    for jj = 1:length(trck)
        for kk = 1:length(lmax)
            
            % load subject
            tmp = feMergeRepeats(dset{ii}, subj{ii}, trck{jj}, lmax{kk});
            
            % pull average count matrix for parameters
            mat = tmp{14}.emat.mean;

            % plot the thing
            figure('Position', [700 430 770 630]);
            colormap('hot');
            imagesc(mat);
            axis('square'); axis('equal'); axis('tight');
            set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
            caxis(cscl{ii});
            line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1], 'LineWidth', 2);
            line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1], 'LineWidth', 2);
            line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1], 'LineWidth', 2);
            
            % save the thing
            fname = [ 'sfn2016/' dset{ii} '_' subj{ii} '_' trck{jj} '_lmax' lmax{kk}, '_emd_average' ];
            print(gcf,  '-cmyk', '-painters', '-dpng', '-r500', '-noui', fname);
            close all;
            clear tmp mat fname;
        end
    end
    
    % repeate load of subject for tensor
    tmp = feMergeRepeats(dset{ii}, subj{ii}, 'tens', '2');
    
    % pull average count matrix for parameters
    mat = tmp{14}.emat.mean;
    
    % plot the thing
    figure('Position', [700 430 770 630]);
    colormap('hot');
    imagesc(mat);
    axis('square'); axis('equal'); axis('tight');
    set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
    caxis(cscl{ii});
    line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1], 'LineWidth', 2);
    line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1], 'LineWidth', 2);
    line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1], 'LineWidth', 2);
            
    % save the thing
    fname = [ 'sfn2016/' dset{ii} '_' subj{ii} '_tens_lmax2_emd_average' ];
    print(gcf, '-cmyk', '-painters', '-dpng', '-r500', '-noui', fname);
    close all;
    clear tmp mat fname;
    
end

% in imagemagic to crop background out of images
% mogrify -trim +repage *.png

clear ii jj kk

%% save 10 repeats (prob10) individually for 4 subjects (2 HCP, 2 STN)

%dset = {'hcp', 'hcp', 'hcp', 'stn', 'stn'};
%subj = {'105115', '110411', '111312', 'FP', 'MP'};

%dset = {'stn', 'stn'};
%subj = {'FP', 'MP'};

for ii = 1:length(subj)
    for jj = 1:10
        
        % load subject
        [~, tmp] = feMergeRepeats(dset{ii}, subj{ii}, 'prob', '10');
        
        % pull average count matrix for parameters
        %mat = log(tmp{14}(:,:,jj));
        mat = tmp{14}(:,:,jj);
        
        % plot the thing
        figure('Position', [700 430 770 630]);
        colormap('hot');
        imagesc(mat);
        axis('square'); axis('equal'); axis('tight');
        set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
        caxis(cscl{ii});
        line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1], 'LineWidth', 2);
        line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1], 'LineWidth', 2);
        line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1], 'LineWidth', 2);
        
        % save the thing
        fname = [ 'sfn2016/' dset{ii} '_' subj{ii} '_prob_lmax10_rep' sprintf('%02d', jj) '_emd' ];
        print(gcf, '-cmyk', '-painters', '-dpng', '-r500', '-noui', fname);
        close all;
        clear tmp mat fname;
    end
end

clear ii jj

%% stanford and hcp are on very different emd scales

[ adetr02, aidetr02, asdetr02 ] = feMergeRepeats('hcp', '105115', 'detr', '2');
[ bdetr02, bidetr02, bsdetr02 ] = feMergeRepeats('stn', 'FP', 'detr', '2');

%% make 4 of each individual one to layer for main plot

mods = {2, 6, 14};
modl = {'dens', 'nzdn', 'lemd'};
mscl = {[ 0 1200 ], [ 0 2 ], [ 0 1000000 ]};

% load subject
[~, dat] = feMergeRepeats('hcp', '105115', 'detr', '2');
        
for ii = 1:length(mods)
    for jj = 1:4
        
        % pull average count matrix for parameters
        mat = dat{mods{ii}}(:,:,jj);
        
        % plot the thing
        figure('Position', [700 430 770 630]);
        colormap('hot');
        imagesc(mat);
        axis('square'); axis('equal'); axis('tight');
        set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
        %caxis(mscl{ii});
        line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1], 'LineWidth', 2);
        line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1], 'LineWidth', 2);
        line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1], 'LineWidth', 2);
        
        % save the thing
        fname = [ 'sfn2016/' 'hcp_105115_detr_lmax02_rep' sprintf('%02d', jj) '_' modl{ii}];
        print(gcf, '-cmyk', '-painters', '-dpng', '-r500', '-noui', fname);
        close all;
        clear mat fname;
    end
end

clear ii jj


