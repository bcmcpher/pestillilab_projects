function [ fh ] = average_plots()
%% plots for Paolo
%

% parameters to loop over
lparam_set = {'8'}; 
alg_set = {'detr', 'prob'};
connMod = {1, 2};
connRng = {[0 8], [-9 -2]};
connLab = {'count', 'density'};

% subjects to loop over
stn = {'FP','HT','KK','MP'};
hc3 = {'105115','110411','111312','113619'};
hc7 = {'108323','109123','131217','910241'};
subjs = [stn hc3 hc7];

fcrng = [0 8];
dnrng = [-9 -2];

% for every subject
for ii = 1:length(subjs)
    subj = subjs{ii};
    
    % for every tracking type
    for jj = 1:length(alg_set)
        alg = alg_set{jj};
        
        % for every parameter
        for kk = 1:length(lparam_set)
            param = lparam_set{kk};
            
            % pick a subject, load individual data, and save plots
            switch subj
                case stn
                    tmp = feMergeRepeats('stn', subj, alg, param);
                    plotMatrix(tmp{1}.emat.mean, fcrng, [ 'paolo/stn_' subj '_' alg '_lmax' param '_count' ]);
                    close all;
                    plotMatrix(tmp{2}.emat.mean, dnrng, [ 'paolo/stn_' subj '_' alg '_lmax' param '_density' ]);
                    close all;
                case hc3
                    tmp = feMergeRepeats('hcp', subj, alg, param);
                    plotMatrix(tmp{1}.emat.mean, fcrng, [ 'paolo/hcp_' subj '_' alg '_lmax' param '_count' ]);
                    close all;
                    plotMatrix(tmp{2}.emat.mean, dnrng, [ 'paolo/hcp_' subj '_' alg '_lmax' param '_density' ]);
                    close all;
                case hc7
                    tmp = feMergeRepeats('7T', subj, alg, param);
                    plotMatrix(tmp{1}.emat.mean, fcrng, [ 'paolo/hc7_' subj '_' alg '_lmax' param '_count' ]);
                    close all;
                    plotMatrix(tmp{2}.emat.mean, dnrng, [ 'paolo/hc7_' subj '_' alg '_lmax' param '_density' ]);
                    close all;
            end
        end
    end
    
    % grab tensor plots / data, too
    switch subj
        case stn
            tmp = feMergeRepeats('stn', subj, 'tens', '2');
            plotMatrix(tmp{1}.emat.mean, fcrng, [ 'paolo/stn_' subj '_tens_count' ]);
            close all;
            plotMatrix(tmp{2}.emat.mean, dnrng, [ 'paolo/stn_' subj '_tens_density' ]);
            close all;
        case hc3
            tmp = feMergeRepeats('hcp', subj, 'tens', '2');
            plotMatrix(tmp{1}.emat.mean, fcrng, [ 'paolo/hcp_' subj '_tens_count' ]);
            close all;
            plotMatrix(tmp{2}.emat.mean, dnrng, [ 'paolo/hcp_' subj '_tens_density' ]);
            close all;
        case hc7
            tmp = feMergeRepeats('7T', subj, 'tens', '2');
            plotMatrix(tmp{1}.emat.mean, fcrng, [ 'paolo/hc7_' subj '_tens_count' ]);
            close all;
            plotMatrix(tmp{2}.emat.mean, dnrng, [ 'paolo/hc7_' subj '_tens_density' ]);
            close all;
    end
end

clear ii jj kk

%% all in 1 plots

%for every model
for hh = 1:length(connMod)
    
    % for every tracking type
    for ii = 1:length(alg_set)
        alg = alg_set{ii};
        
        % for every parameter
        for jj = 1:length(lparam_set)
            param = lparam_set{jj};
            
            fh = figure(); ax = gca;
            for kk = 1:length(subjs)
                subj = subjs{kk};
                
                % pick a subject, load individual data, and save plots
                switch subj
                    case stn
                        tmp = feMergeRepeats('stn', subj, alg, param);
                        subplot(3, 4, kk);
                        plotMatrixOnly(tmp{connMod{hh}}.emat.mean, connRng{hh});
                    case hc3
                        tmp = feMergeRepeats('hcp', subj, alg, param);
                        subplot(3, 4, kk);
                        plotMatrixOnly(tmp{connMod{hh}}.emat.mean, connRng{hh});
                    case hc7
                        tmp = feMergeRepeats('7T', subj, alg, param);
                        subplot(3, 4, kk);
                        plotMatrixOnly(tmp{connMod{hh}}.emat.mean, connRng{hh});
                end
            end
            
            % add global colorbar
            %h = colorbar;
            
            %             hp4 = get(subplot(2,2,4),'Position');
            %             colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.1  hp4(2)+hp4(3)*2.1])
            %
            %             set(h, 'Position', [.8314 .11 .0581 .8150])
            %             for i=1:3
            %                 pos=get(ax(i), 'Position');
            %                 set(ax(i), 'Position', [pos(1) pos(2) 0.85*pos(3) pos(4)]);
            %             end
            %
            % save file
            print(gcf, '-cmyk', '-painters', '-dpng', '-r500', '-noui', ['paolo/All_' alg '_lmax' param '_' connLab{hh}]);
            close all;
        end
        
        %         % grab tensor plots / data, too
        %         switch subj
        %             case stn
        %                 tmp = feMergeRepeats('stn', subj, 'tens', '2');
        %
        %             case hc3
        %                 tmp = feMergeRepeats('hcp', subj, 'tens', '2');
        %
        %             case hc7
        %                 tmp = feMergeRepeats('7T', subj, 'tens', '2');
        %
        %          end
    end
end

end

%% internal plot function

function [fh] = plotMatrix(inMat, crng, fout)

% create color scale
mp = (crng(2) + crng(1)) / 2;
apts = [crng(1) mp crng(2)];
albs = {crng(1), mp, crng(2)};

% do the plot thing
fh = figure;
colormap('hot');
imagesc(log(inMat));
axis('square'); axis('equal'); axis('tight');
caxis(crng);
colorbar('Ticks', apts, 'TickLabels', albs);
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1], 'LineWidth', 2);
line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1], 'LineWidth', 2);
line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1], 'LineWidth', 2);

% save file
print(gcf, '-cmyk', '-painters', '-dpng', '-r500', '-noui', fout);

close all;

end

function [fh] = plotMatrixOnly(inMat, crng)

% create color scale
mp = (crng(2) + crng(1)) / 2;
apts = [crng(1) mp crng(2)];
albs = {crng(1), mp, crng(2)};

% do the plot thing
colormap('hot');
imagesc(log(inMat));
axis('square'); axis('equal'); axis('tight');
caxis(crng);
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1], 'LineWidth', 2);
line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1], 'LineWidth', 2);
line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1], 'LineWidth', 2);

end