%% load and compare pre/post LiFE measures between data sets to see scaling of SNR

% left precentral-postcentral (motor strip): pconn{517}
% right precentral-postcentral (motor strip): pconn{1962}

% left parstriagularis-parsorbitalis: pconn{200}
% right parstriagularis-parsorbitalis: pconn{1815}

% load all the subjects
subjects = {'hcp_105115', 'hcp_110411', 'hcp_111312', 'hcp_113619', ...
            'stn_FP', 'stn_MP', 'stn_HT', 'stn_KK', ...
            '7T_108323', '7T_109123', '7T_131217', '7T_910241'};

reps = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10'};
        
% multiple outputs for pre/post LiFE
pre_dens = zeros(length(subjects), 4, 10);
post_dens = zeros(length(subjects), 4, 10);

% I didn't save down EMD vectors...
life_emd = zeros(length(subjects), 4, 10);

% new checks
chck01 = zeros(length(subjects), 4, 10);
chck02 = zeros(length(subjects), 4, 10);
chck03 = zeros(length(subjects), 1, 10);
chck04 = zeros(length(subjects), 1, 10);

% create index into EMD matrix b/c I was dumb before
pairs = nchoosek(1:68, 2);

% for every subject
for subj = 1:length(subjects)
    
    % for every repeat
    for jj = 1:length(reps)
        
        % load data
        load(strcat('data/', subjects{subj}, '_prob_lmax8_rep', reps{jj}, '.mat'), 'pconn');
        load(strcat('data/', subjects{subj}, '_prob_lmax8_rep', reps{jj}, '.mat'), 'emat');
        load(strcat('data/', subjects{subj}, '_prob_lmax8_rep', reps{jj}, '.mat'), 'out');
                
        % left
        pre_dens(subj, 1, jj) = (2 * size(pconn{517}.end.fibers, 1)) / size(unique(pconn{517}.roisvx, 'rows'), 1);
        pre_dens(subj, 2, jj) = (2 * size(pconn{200}.end.fibers, 1)) / size(pconn{517}.roisvx, 1);
        
        % right
        pre_dens(subj, 3, jj) = (2 * size(pconn{1962}.end.fibers, 1)) / size(unique(pconn{1962}.roisvx, 'rows'), 1);
        pre_dens(subj, 4, jj) = (2 * size(pconn{1815}.end.fibers, 1)) / size(unique(pconn{1962}.roisvx, 'rows'), 1);
        
        % left
        post_dens(subj, 1, jj) = (2 * size(pconn{517}.end.nzfibs, 1)) / size(unique(pconn{517}.roisvx, 'rows'), 1);
        post_dens(subj, 2, jj) = (2 * size(pconn{200}.end.nzfibs, 1)) / size(unique(pconn{1962}.roisvx, 'rows'), 1);
        
        % right
        post_dens(subj, 3, jj) = (2 * size(pconn{1962}.end.nzfibs, 1)) / size(unique(pconn{517}.roisvx, 'rows'), 1);
        post_dens(subj, 4, jj) = (2 * size(pconn{1815}.end.nzfibs, 1)) / size(unique(pconn{1962}.roisvx, 'rows'), 1);
        
        % grab life EMD output
        life_emd(subj, 1, jj) = emat(pairs(517, 1), pairs(517, 2), 14);
        life_emd(subj, 2, jj) = emat(pairs(200, 1), pairs(200, 2), 14);
        life_emd(subj, 3, jj) = emat(pairs(1962, 1), pairs(1962, 2), 14);
        life_emd(subj, 4, jj) = emat(pairs(1815, 1), pairs(1815, 2), 14);
        
        % checks
        
        % voxel size
        chck01(subj, 1, jj) = size(unique(pconn{517}.roisvx, 'rows'), 1);
        chck01(subj, 2, jj) = size(unique(pconn{200}.roisvx, 'rows'), 1);
        chck01(subj, 3, jj) = size(unique(pconn{1962}.roisvx, 'rows'), 1);
        chck01(subj, 4, jj) = size(unique(pconn{1815}.roisvx, 'rows'), 1);
        
        % streamline count
        chck02(subj, 1, jj) = size(pconn{517}.end.fibers, 1);
        chck02(subj, 2, jj) = size(pconn{200}.end.fibers, 1);
        chck02(subj, 3, jj) = size(pconn{1962}.end.fibers, 1);
        chck02(subj, 4, jj) = size(pconn{1815}.end.fibers, 1);
        
        % all ROI voxel size
        for kk = 1:length(out)
            chck03(subj, jj) = chck03(subj) + size(unique(out{kk}.vxroi, 'rows'), 1);
        end
        
        % all streamline count
        for ll = 1:length(pconn)
            chck04(subj, jj) = chck04(subj) + size(pconn{ll}.end.fibers, 1);
        end
 
    end
end

clear subj jj kk ll pconn emat out

%% plot the data

% get the colors
c1 = colormap(parula(64));
c2 = colormap(autumn(64));
tcol = [ c1([3 9 1 6],:); c1([34 44 54 64],:); c2([32 13 25 5],:); ];
clear c1 c2;
close all;

% build x-coordinates in a dumb way
xcoord1 = [1, 1, 1, 1, 3, 3, 3, 3, 5, 5, 5, 5];
xcoord2 = [2, 2, 2, 2, 4, 4, 4, 4, 6, 6, 6, 6];
xcoord3 = [];

xlabel = {'Pre-LiFE Density', 'Post-LiFE Density'};

% build a jitter
% - left / right
% - between observations
lrjit = 0.15;
objit = linspace(-0.05, 0.05, 10);

prex = 1 + objit;
pstx = 2 + objit;

% plots

% % pre/post-density, motorstrip
% figure; hold on;
% title('Motorstrip - Density pre/post LiFE');
% xlim([0.50 2.50]); %ylim(log10([0.001 .18]));
% for ii = 1:8
%     
%     for jj = 1:10
%         
%         % prelife density
%         plot(prex(jj) - lrjit, log10(pre_dens(ii, 1, jj)), 'o', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
%         plot(prex(jj) + lrjit, log10(pre_dens(ii, 3, jj)), 's', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
%         
%         % postlife density
%         plot(pstx(jj) - lrjit, log10(post_dens(ii, 1, jj)), 'o', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
%         plot(pstx(jj) + lrjit, log10(post_dens(ii, 3, jj)), 's', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
%         
%     end
% end
% set(gca, 'XTick', [1 2], 'XTickLabel', xlabel); clear ii jj

% dens / life, motorstrip
figure; hold on;
title('Motorstrip - Density compared to LiFE');
xlim([0.50 2.50]); %ylim([-0.50 5]);
for ii = 1:8
    
    for jj = 1:10
        
        % prelife density
        plot(prex(jj) - lrjit, log10(pre_dens(ii, 1, jj)), 'o', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
        plot(prex(jj) + lrjit, log10(pre_dens(ii, 3, jj)), 's', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
        
        % life EMD
        plot(pstx(jj) - lrjit, log10(life_emd(ii, 1, jj)), 'o', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
        plot(pstx(jj) + lrjit, log10(life_emd(ii, 3, jj)), 's', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
        
    end
end
set(gca, 'XTick', [1 2], 'XTickLabel', {'Pre-LiFE Density', 'LiFE EMD'}); clear ii jj

% pre/post-density, small
figure; hold on;
title('ParsTri-ParsOp - Density pre/post LiFE');
xlim([0.50 2.50]); ylim([0 .08]);
for ii = 1:8
    
    for jj = 1:10
        
        % prelife density
        plot(prex(jj) - lrjit, pre_dens(ii, 2, jj), 'o', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
        plot(prex(jj) + lrjit, pre_dens(ii, 4, jj), 's', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
        
        % postlife density
        plot(pstx(jj) - lrjit, post_dens(ii, 2, jj), 'o', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
        plot(pstx(jj) + lrjit, post_dens(ii, 4, jj), 's', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
        
    end
end
set(gca, 'XTick', [1 2], 'XTickLabel', xlabel); clear ii jj

% dens / life, motorstrip
figure; hold on;
title('ParsTri-ParsOp - Density compared to LiFE');
xlim([0.50 2.50]); ylim([-0.50 5]);
for ii = 1:8
    
    for jj = 1:10
        
        % prelife density
        plot(prex(jj) - lrjit, pre_dens(ii, 2, jj), 'o', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
        plot(prex(jj) + lrjit, pre_dens(ii, 4, jj), 's', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
        
        % life EMD
        plot(pstx(jj) - lrjit, life_emd(ii, 2, jj), 'o', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
        plot(pstx(jj) + lrjit, life_emd(ii, 4, jj), 's', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
        
    end
end
set(gca, 'XTick', [1 2], 'XTickLabel', {'Pre-LiFE Density', 'LiFE EMD'}); clear ii jj

%% original plots 

% pre-density, motorstrip
figure; hold on;
xlim([0.25 6.5]); ylim([0 .04]);
for ii = 1:12
    plot(xcoord1(ii), pre_dens(ii, 1), 'o', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
    plot(xcoord2(ii), pre_dens(ii, 3), 'o', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
end
set(gca, 'XTickLabel', xlabel); clear ii


% post-density, motorstrip
figure; hold on;
xlim([0.25 6.5]); ylim([0 .04]);
for ii = 1:12
    plot(xcoord1(ii), post_dens(ii, 1), 'o', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
    plot(xcoord2(ii), post_dens(ii, 3), 'o', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
end
set(gca, 'XTickLabel', xlabel); clear ii

% EMD, motorstrip
figure; hold on;
xlim([0.25 6.5]); ylim([-0.5 5.5]);
for ii = 1:12
    plot(xcoord1(ii), life_emd(ii, 1), 'o', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
    plot(xcoord2(ii), life_emd(ii, 3), 'o', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
end
set(gca, 'XTickLabel', xlabel); clear ii

% pre-density, small
figure; hold on;
xlim([0.25 6.5]); ylim([0 .07]);
for ii = 1:12
    plot(xcoord1(ii), pre_dens(ii, 2), 'o', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
    plot(xcoord2(ii), pre_dens(ii, 4), 'o', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
end
set(gca, 'XTickLabel', xlabel); clear ii

% post-density, small
figure; hold on;
xlim([0.25 6.5]); ylim([0 .012]);
for ii = 1:12
    plot(xcoord1(ii), post_dens(ii, 2), 'o', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
    plot(xcoord2(ii), post_dens(ii, 4), 'o', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
end
set(gca, 'XTickLabel', xlabel); clear ii

% EMD, small
figure; hold on;
xlim([0.25 6.5]); ylim([-.5 5.5]);
for ii = 1:12
    plot(xcoord1(ii), life_emd(ii, 2), 'o', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
    plot(xcoord2(ii), life_emd(ii, 4), 'o', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
end
set(gca, 'XTickLabel', xlabel); clear ii



%% plot development
figure; hold on;
xlim([0.25 6.5]);
ylim([0 .04]);

% for every subject - switch between pre_dens / post_dens / life_emd
for ii = 1:12
    plot(xcoord1(ii), pre_dens(ii, 2), 'o', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
    plot(xcoord2(ii), pre_dens(ii, 4), 'o', 'color', 'black', 'MarkerFaceColor', tcol(ii, :));
end

% set plot axis labels
set(gca, 'XTickLabel', xlabel);

clear ii












