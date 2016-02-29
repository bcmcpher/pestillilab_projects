%% make plots w/ BCT of EM / Strength of Evidence matrices
% Brent McPherson
% 20151019
%

%load ../../mrtrix/matrix_20160119.mat

%% BCT Toolbox Guesses

% set proportional / absolute thresholding
em_tha = threshold_absolute(em, 7);
em_thp = threshold_proportional(em, 0.75);

%em_tha = log(em);
%em_tha(isinf(em_tha)) = 0;
%em_tha(em_tha < 0) = 0;

% weighted conversion
em_bin = weight_conversion(em_tha, 'binarize');
em_nrm = weight_conversion(em_tha, 'normalize');
em_len = weight_conversion(em_tha, 'lengths');
em_fix = weight_conversion(em_tha, 'autofix');

% clustering / strength / density of em_tha network
em_deg = degrees_und(em_tha);
em_str = strengths_und(em_tha);
em_den = density_und(em_tha);

% topological measure - for binary undirected network
em_top = gtom(em_bin, 1);
em_top = gtom(em_bin, 2);

% neighborhood overlap
em_nei = edge_nei_overlap_bu(em_bin);

% matching index
em_match = matching_ind_und(em_bin);

% rentian scaling - NEED SPATIAL COORDIANTES OF MATRIX NODES
%[ em_rnt_node, em_rnt_edge  ] = rentian_scaling(em, XYZ, n);
% XYZ is a 3xlength(labs) matrix of spatial coordinates of nodes
% n is partitions to compute to estimate the scaling; more is better

% clustering coefficient - must used normalized (0-1) data
em_ccoef = clustering_coef_wu(em_nrm);

% clustering_coef_wu_sign has alternative algorithms for calculating outs
[ em_Cpos, em_Cneg, em_Tpos, em_Tneg ] = clustering_coef_wu_sign(em_tha, 1); % Onnela et al.
[ em_Cpos, em_Cneg, em_Tpos, em_Tneg ] = clustering_coef_wu_sign(em_tha, 2); % Zhang & Horvath
[ em_Cpos, em_Cneg, em_Tpos, em_Tneg ] = clustering_coef_wu_sign(em_tha, 3); % Constantini & Perugini

% rich club / transitivity coefficient
em_rc = rich_club_wu(em_tha);
em_tr = transitivity_wu(em_tha);

% global / local binary efficiency
em_bgl = efficiency_bin(em_bin, 0);
em_blc = efficiency_bin(em_bin, 1);

% connected components
[ em_cmp, em_cmp_size ] = get_components(em_bin);

% global / local weighted efficiency
em_glb = efficiency_wei(em_tha, 0);
em_loc = efficiency_wei(em_tha, 1);

% community structure
em_lov = community_louvain(em_tha);
em_mod = modularity_und(em_tha);

% S score - not sure it matters...
[ em_sscore, em_ssize ] = score_wu(em_tha, 0.5);

% distance and path length
[ em_dis, em_ded ] = distance_wei(em_len);

% node betweenness / centrality metrics
em_btw = betweenness_wei(em_tha);
[ em_edge_btw, em_edge_cvc ] = edge_betweenness_wei(em_tha);

% self-referential connectedness
em_eig = eigenvector_centrality_und(em_tha);

% network summaries using community affiliation vector
em_pcoef = participation_coef(em_tha, em_lov);
em_modz = module_degree_zscore(em_tha, em_lov);

% reorder EM matrix - 3 methods
[ em_reord, em_reorI, em_cost ] = reorder_matrix(em_tha, 'line', 1);

[ em_ind, em_ord ] = reorder_mod(em_tha, em_lov);

[ em_MATr, em_MATi, em_MATc ] = reorderMAT(em_tha, 1000, 'line');

% backbone of  connectivity
[ em_Tree, em_Ptree ] = backbone_wu(em_tha, 2);

%% BCT Primary Plot

% 1. Threshold more strictly and binarize
% 2. Find and recreate a matrix (as best I can) from Olaf
% 3. Figure out how to evaluate all of this

% pull raw EMD to binary matrix
% normalize relabeled binary matrix

% make copies for sanity
tmp1 = em;
tmp2 = em_bin;

% make a subset vector and counter
tmp3 = tmp1(tmp2 == 1);
indx = 1;
out = zeros(size(tmp1, 1), size(tmp2, 2));

% because subsetting by a logical matrix would be too easy...
% symmetric, so it doesn't matter row / column?

% for each row
for yy = 1:size(tmp2, 1)
    % for each column
    for zz = 1:size(tmp2, 2)
        % assign the value and iterate if binary  matrix is 1
        if em_bin(yy, zz) == 1;
            out(yy, zz) = tmp3(indx);
            indx = indx + 1;
        % else assign 0
        else 
            out(yy, zz) = 0;
        end
    end
end

% normalize remaining weights
out = threshold_absolute(em, 8);

% find specific labels and use just those
% left / right labels
% EMD colorbar label
% very detailed identification of some connections

h = figure('name','Thresheld EMD Connectome Matrix', 'color','w', 'Position', [618 -221 1190 907]);
colormap('hot');
hax = axes;
imagesc(out);
hold on;
%title('Thresheld EMD Connectome Matrix');
%set(gca, 'YTick', 1:1:length(labs), 'YTickLabel', labs,'tickdir','out');
set(gca, 'XTick', [], 'YTick', []);
set(gcf, 'color', 'w');
%xlabel('FreeSurfer Define Anatomical Regions');
%ylabel('FreeSurfer Define Anatomical Regions');
axis('square'); axis('equal'); axis('tight');
y = colorbar;
ylabel(y, 'Connection Evidence (emd)');
SP = size(out, 1) / 2;
line([SP SP], get(hax, 'YLim'), 'Color', [1 1 1]);
line(get(hax, 'XLim'), [SP SP], 'Color', [1 1 1]);
line([SP*2 0], [SP*2 0], 'Color', [1 1 1]);
xlab = title('Left Hemisphere   Right Hemisphere', 'FontWeight', 'normal');
ylab = ylabel('Right Hemisphere   Left Hemisphere');
set(hax, 'XLim', [ 0 length(labs) ], 'YLim', [ 0 length(labs) ], ...
    'FontName', 'Helvetica', 'FontAngle', 'italic', 'FontSize', 24, 'FontWeight', 'normal', ...
    'XTickLabel', [], 'YTickLabel', []);

% save both when writing out
feSavefig(h, 'figType', 'eps', 'figName', 'em_thr8_20151026.eps','figDir','/N/dc2/projects/lifebid/HCP/Brent');
feSavefig(h, 'figType', 'png', 'figName', 'em_thr8_20151026.png','figDir','/N/dc2/projects/lifebid/HCP/Brent');

%% redundant figures

% raw EMD matrix
figure('name','EMD Connectome Matrix')
colormap('hot');
imagesc(em); axis('square'); axis('equal'); axis('tight');
set(gca, 'YTick', 1:1:length(labs), 'YTickLabel', labs);
colorbar

% absolute threshold EMD
figure('name','Abs. Thresheld EMD Connectome Matrix')
colormap('hot');
imagesc(em_tha); axis('square'); axis('equal'); axis('tight');
set(gca, 'YTick', 1:1:length(labs), 'YTickLabel', labs);
colorbar

% proportional threshold EMD
figure('name','Proportional Thresheld EMD Connectome Matrix')
colormap('hot');
imagesc(em_thp); axis('square'); axis('equal'); axis('tight');
set(gca, 'YTick', 1:1:length(labs), 'YTickLabel', labs);
colorbar

% binary EMD
figure('name','Binary EMD Connectome Matrix')
colormap('hot');
imagesc(em_bin); axis('square'); axis('equal'); axis('tight');
set(gca, 'YTick', 1:1:length(labs), 'YTickLabel', labs);
colorbar

% normalized EMD
figure('name','Normalized EMD Connectome Matrix')
colormap('hot');
imagesc(em_nrm); axis('square'); axis('equal'); axis('tight');
set(gca, 'YTick', 1:1:length(labs), 'YTickLabel', labs);
colorbar

% length of EMD
figure('name','Length EMD Connectome Matrix')
colormap('hot');
imagesc(em_len); axis('square'); axis('equal'); axis('tight');
set(gca, 'YTick', 1:1:length(labs), 'YTickLabel', labs);
colorbar

% fixed (?) EMD
figure('name','Fixed EMD Connectome Matrix')
colormap('hot');
imagesc(em_fix); axis('square'); axis('equal'); axis('tight');
set(gca, 'YTick', 1:1:length(labs), 'YTickLabel', labs);
colorbar

% S-Score of matrix
figure('name','S-Score EMD Connectome Matrix')
colormap('hot');
imagesc(em_sscore); axis('square'); axis('equal'); axis('tight');
set(gca, 'YTick', 1:1:length(labs), 'YTickLabel', labs);
colorbar

% Weighted Path based on Distance
figure('name','Weighted Path - Distance - EMD Connectome Matrix')
colormap('hot');
imagesc(em_dis); axis('square'); axis('equal'); axis('tight');
set(gca, 'YTick', 1:1:length(labs), 'YTickLabel', labs);
colorbar

% Weighted Path based on number of edges
figure('name','Weighted Path - Edges - EMD Connectome Matrix')
colormap('hot');
imagesc(em_ded); axis('square'); axis('equal'); axis('tight');
set(gca, 'YTick', 1:1:length(labs), 'YTickLabel', labs);
colorbar

% Edge betweenness in EMD matrix
figure('name','Edge Betweenness EMD Connectome Matrix')
colormap('hot');
imagesc(em_edge_btw); axis('square'); axis('equal'); axis('tight');
set(gca, 'YTick', 1:1:length(labs), 'YTickLabel', labs);
colorbar

% reordered EMD
figure('name','Reordered EMD Connectome Matrix')
colormap('hot');
imagesc(em_reord); axis('square'); axis('equal'); axis('tight');
set(gca, 'YTick', 1:1:length(labs), 'YTickLabel', labs(em_reorI));
colorbar

% reoder EM matrix w/ community affiliation
figure('name','Reordered Community Weighted EMD Connectome Matrix')
colormap('hot');
imagesc(em_ord_log); axis('square'); axis('equal'); axis('tight');
%set(gca, 'YTick', 1:1:length(labs), 'YTickLabel', labs(em_ind));
colorbar

% reorderMAT of EM
figure('name','ReorderedMAT EMD Connectome Matrix')
colormap('hot');
imagesc(em_MATr); axis('square'); axis('equal'); axis('tight');
set(gca, 'YTick', 1:1:length(labs), 'YTickLabel', labs(em_MATi));
colorbar

% backbone tree of EM 
figure('name','Backbone Tree EMD Connectome Matrix')
colormap('hot');
imagesc(em_Tree); axis('square'); axis('equal'); axis('tight');
set(gca, 'YTick', 1:1:length(labs), 'YTickLabel', labs);
colorbar

% backbone tree + strongest connections up to average
figure('name','Backbone Plus EMD Connectome Matrix')
colormap('hot');
imagesc(em_Ptree); axis('square'); axis('equal'); axis('tight');
set(gca, 'YTick', 1:1:length(labs), 'YTickLabel', labs);
colorbar
