%% plot global / node based network repeated stats across all lmax/tracking types
% Brent McPherson
% 20160905
%

% load data - edit function to add different parameters
dat = Gen_data_for_plots();

% copy of parameter set to parse data object - from top of Gen_data_for_plots()
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

%% build color map - just hard code numbers once I have 7T values

% % this is the order and subsets Cesar pulled on the subjects
% HCP_subject_set = {'111312','105115','113619','110411'};
% STN_subject_set = {'KK_96dirs_b2000_1p5iso','FP_96dirs_b2000_1p5iso','HT_96dirs_b2000_1p5iso','MP_96dirs_b2000_1p5iso'};
% HCP7T_subject_set = {'108323','131217','109123','910241'};

% % this is my order
%sbj_set = {'105115', '110411', '111312', '113619'};
%sbj_set = {'FP', 'HT', 'KK', 'MP'};
%sbj_set = {'108323','109123','131217','910241'};

% switch color_type
%     case 'cold'
%         c = [c1([1 3 6 9],:) ];
%     case 'medium'
%         c = [c1([12 16 19 23],:) ];
%     case 'hot'
%         c = [c2([32 25 13 5],:)];
% end
% from 2 color maps

c1 = colormap(parula(64));
c2 = colormap(autumn(64));

% in the end, hardcode tcol w/ all colors in the order I use the subjects
%tcol = [ c1([3 9 1 6],:); c1([16 19 12 23],:); c2([32 13 25 5],:); ];
tcol = [ c1([3 9 1 6],:); c1([34 44 54 64],:); c2([32 13 25 5],:); ];
clear c1 c2;
close all;

% % hard coded colors to subjects
% tcol = [ ...
% 0.2123 0.2138 0.6270; % 105115
% 0.0117 0.3875 0.8820; % 110411
% 0.2081 0.1663 0.5292; % 111312
% 0.1707 0.2919 0.7792; % 113619
% 0.0779 0.5040 0.8384; % FP
% 0.0641 0.5570 0.8240; % HT
% 0.0329 0.4430 0.8720; % KK
% 0.0239 0.6287 0.8038; % MP
% 1.0000 0.4921      0; % 108323
% 1.0000 0.1905      0; % 109123
% 1.0000 0.3810      0; % 131217
% 1.0000 0.0635      0; % 910241
% ];

%% global network statistics

% jitter for subject and parameter, combined in loop
% will need to modify once 7T subjects added
sjit = linspace(-0.35, 0.35, length(dat));
ljit = linspace(-0.025, 0.025, length(lparam_set));

% make a plot of each glbVars for every subject/parameter combo
for imod = 1:length(glbLabel)
    % create empty figure
    figure('Position', [950 -170 1000 1100]); hold on;
    ylim([ 0.50 3.50 ]); % adjust ylim based on how many connMods are used
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
                    try
                        pt = eval(['dat{' num2str(ii) '}.' dlab '.glob{' num2str(jj) '}(' num2str(imod) ',:)' ]);
                        plot(pt(1), jj+jitter, shape, 'Color', tcol(ii,:), 'MarkerFaceColor', tcol(ii,:), 'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 18);
                        plot([pt(2) pt(3)], [jj jj]+jitter,  'LineWidth', 2, 'Color', [0 0 0]);
                    catch
                        %display([dat{ii}.subj '_' dlab ' not plotted']);
                    end
                end
            end
            % plot tensor points
            pt = eval(['dat{' num2str(ii) '}.tens.glob{' num2str(jj) '}(' num2str(imod) ',:)' ]);
            plot(pt(1), jj+sjit(ii), 'diamond', 'Color', tcol(ii,:), 'MarkerFaceColor', tcol(ii,:), 'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 18);
            plot([pt(2) pt(3)], [sjit(ii) sjit(ii)]+jj,  'LineWidth', 2, 'Color', [0 0 0]);
        end
    end
    % plot labels
    title(glbLabel{imod});
    
    % calculate x axis labels
    xtcks = get(gca, 'XTick');
    mxtck = max(xtcks);
    xpnts = [0 mxtck/2 mxtck];
    
    % adjust x / y axis
    set(gca, 'YTick', [1 2 3], 'YTickLabel', modLabel, 'YTickLabelRotation', 90, ...
        'XTick', xpnts, 'XTickLabel', xpnts, 'XLim', [0 mxtck]);
    
    % set fontsize for minimum visibility on poster
    set(gca, 'FontAngle', 'oblique', 'FontSize', 29.37);
    
    hold off;
    clear xtcks mxtck xpnts
end

clear ii jj kk ll imod alg shape param dlab jitter pt 

%% node based statistics

% can't really tell what anything is, poor discriminability between param / subj...

% for 1st node statistic
nodMod = 1;

% to better stagger subjects, make a really wide plot - will have to manually add more when 7T are done
% doesn't really improve resolution, still hard to see
sjit(1,:) = 1:8:544;
sjit(2,:) = 2:8:544;
sjit(3,:) = 3:8:544;
sjit(4,:) = 4:8:544;
sjit(5,:) = 5:8:544;
sjit(6,:) = 6:8:544;
sjit(7,:) = 7:8:544;
sjit(8,:) = 8:8:544;

% nodes along x-axis all the values as y
figure; hold on;
% for every subject loaded in a data structure
for ii = 1:length(dat)
    % for every model
    for jj = 1:length(connMod) % IS THIS OVERLAYING 3 MODELS ON THE SAME PLOT?
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
            for ll = 1:length(lparam_set) % perhaps vary shape by lmax? - more discernable difference that way?
                param = lparam_set{ll};
                % combine subject and lmax jitter for more obvious offsets
                dlab = [alg '_' param];
                pt = eval(['dat{' num2str(ii) '}.' dlab '.node{' num2str(nodMod) ',' num2str(jj) '}' ]);
                plot(sjit(ii,:), pt(:,1), shape, 'Color', tcol(ii,:), 'MarkerFaceColor', tcol(ii,:), 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerSize', 5);
                plot([sjit(ii,:); sjit(ii,:)], [pt(:,2)'; pt(:,3)'],'LineWidth', 1, 'Color', [0 0 0]);  
            end
        end
        % plot tensor points
        %pt = eval(['dat{' num2str(ii) '}.tens.glob{' num2str(jj) '}(' num2str(imod) ',:)' ]);
        %plot(pt(1), jj+sjit(ii), 'diamond', 'Color', tcol(ii,:), 'MarkerFaceColor', tcol(ii,:), 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerSize', 10);
        %plot([pt(2) pt(3)], [sjit(ii) sjit(ii)]+jj,  'LineWidth', 2, 'Color', [0 0 0]);
    end
end
% plot labels
ylabel(['Node ' nodLabel{nodMod}]);
load('roilabels.mat');
xlabel('Network Nodes');
set(gca, 'XTick', 4.5:8:545, 'XTickLabel', roiLabel, 'XTickLabelRotation', 45);
hold off;

