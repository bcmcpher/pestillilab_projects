function [ emat, nmat, fns, pval, fdat ] = feMatrixFromTensor_v01(fe, anatPath, outPath, nclust, cachePath) 
%% development of tensor based matrix construction
% Brent McPherson
% 20160716
% 
% fe: path to fit FE structure with evaluated life model
% anatPath: path to subject directory in FreeSurfer $SUBJECTS_DIR - need to
%           run mri_annotation2label & mri_label2vol for all ROIs
% outPath: path with file name to save output .mat object
% nclust: number of cores to use with parallel toolbox
% cachePath: where temporary space for parallel will be (important in Karst/BR2, jobs fail if it uses /tmp
% 
% tic; [ z1, z2, z3, z4, z5 ] = feMatrixFromTensor_v01('/N/dc2/projects/lifebid/code/ccaiafa/Caiafa_Pestilli_paper2015/Results/ETC_Dec2015/Single_TC/105115/fe_structure_105115_STC_run01_SD_PROB_lmax10.mat', '/N/dc2/projects/lifebid/2t1/anatomy/105115/', '~/feMatrixFromTensor_test/test_output.mat', 16, '~/feMatrixFromTensor_test'); toc;
%

%% pull all ROI NIfTI files of FreeSurfer labels

% path to FreeSurfer output
%roiNames = dir(fullfile(anatPath, 'label', '*_label.nii.gz'));

% % create cell array of file names
% roiNames = {roiNames(:).name}';
% 
% % subset left / right indices of 34 DK labels to use
% roiIndex = [[88:90, 92, 94:123] [88:90, 92, 94:123] + 123];
% roiNames = roiNames(roiIndex);
% 
% % sort DK labels by lobes
% roiOrder = [ 29, 28, 3, 19, 20, 21, 13, 15, 25, 17, 6, 27, 2, 10, ... % frontal
%              30, 8, 32, 23, 26, 24, 11, ... % parietal
%              9, 16, 31, 1, 7, 34, 5, 33, 18, ... % temporal
%              12, 14, 4, 22 ]; % occipital
%              
% % roiOrder = [ 29, 28, 3, 19, 20, 21, 13, 15, 25, 17, 6, ... % frontal
% %              30, 8, 32, 23, 26, ... % parietal
% %              9, 16, 31, 1, 7, 34, 5, 33, 18, ... % temporal
% %              12, 14, 4, 22, ... % occipital
% %              27, 2, 24, 11, 10]; % cingulate (frontal, frontal, parietal, parietal, insula)
% 
% % create sorted left / right order for matrix
% roiOrder = [roiOrder [roiOrder+34]]';
% 
% % order ROIs
% roiNames = roiNames(roiOrder);
% 
% % create ROI labels
% roiLabel = strrep(roiNames, '_label.nii.gz', '');
% roiLabel = strrep(roiLabel, '_', '.');

% need to be numeric indices to pull from single file
roiNames = {...
    'lh_superiorfrontal_label.nii.gz';          % 1028
    'lh_rostralmiddlefrontal_label.nii.gz';     % 1027
    'lh_caudalmiddlefrontal_label.nii.gz';      % 1003
    'lh_parsopercularis_label.nii.gz'; 	    % 1018
    'lh_parsorbitalis_label.nii.gz'; 	    % 1019
    'lh_parstriangularis_label.nii.gz'; 	    % 1020
    'lh_lateralorbitofrontal_label.nii.gz';     % 1012
    'lh_medialorbitofrontal_label.nii.gz'; 	    % 1014
    'lh_precentral_label.nii.gz'; 		    % 1024
    'lh_paracentral_label.nii.gz'; 		    % 1017
    'lh_frontalpole_label.nii.gz'; 		    % 1032
    'lh_rostralanteriorcingulate_label.nii.gz'; % 1026
    'lh_caudalanteriorcingulate_label.nii.gz';  % 1002
    'lh_insula_label.nii.gz'; 		    % 1035
    'lh_superiorparietal_label.nii.gz'; 	    % 1029
    'lh_inferiorparietal_label.nii.gz'; 	    % 1008
    'lh_supramarginal_label.nii.gz'; 	    % 1031
    'lh_postcentral_label.nii.gz'; 		    % 1022
    'lh_precuneus_label.nii.gz'; 		    % 1025
    'lh_posteriorcingulate_label.nii.gz'; 	    % 1023
    'lh_isthmuscingulate_label.nii.gz'; 	    % 1010
    'lh_inferiortemporal_label.nii.gz'; 	    % 1009
    'lh_middletemporal_label.nii.gz'; 	    % 1015
    'lh_superiortemporal_label.nii.gz'; 	    % 1030
    'lh_bankssts_label.nii.gz'; 		    % 1001
    'lh_fusiform_label.nii.gz'; 		    % 1007
    'lh_transversetemporal_label.nii.gz'; 	    % 1034
    'lh_entorhinal_label.nii.gz'; 		    % 1006
    'lh_temporalpole_label.nii.gz'; 	    % 1033
    'lh_parahippocampal_label.nii.gz'; 	    % 1016
    'lh_lateraloccipital_label.nii.gz'; 	    % 1011
    'lh_lingual_label.nii.gz'; 		    % 1013
    'lh_cuneus_label.nii.gz'; 		    % 1005
    'lh_pericalcarine_label.nii.gz';	    % 1021
    'rh_superiorfrontal_label.nii.gz'; 	    % 2028
    'rh_rostralmiddlefrontal_label.nii.gz';	    % 2027
    'rh_caudalmiddlefrontal_label.nii.gz';	    % 2003
    'rh_parsopercularis_label.nii.gz';	    % 2018
    'rh_parsorbitalis_label.nii.gz';	    % 2019
    'rh_parstriangularis_label.nii.gz';	    % 2020
    'rh_lateralorbitofrontal_label.nii.gz';	    % 2012
    'rh_medialorbitofrontal_label.nii.gz';	    % 2014
    'rh_precentral_label.nii.gz';		    % 2024
    'rh_paracentral_label.nii.gz';		    % 2017
    'rh_frontalpole_label.nii.gz';		    % 2032
    'rh_rostralanteriorcingulate_label.nii.gz'; % 2026
    'rh_caudalanteriorcingulate_label.nii.gz';  % 2002
    'rh_insula_label.nii.gz';		    % 2035
    'rh_superiorparietal_label.nii.gz';	    % 2029
    'rh_inferiorparietal_label.nii.gz';	    % 2008
    'rh_supramarginal_label.nii.gz';	    % 2031
    'rh_postcentral_label.nii.gz';		    % 2022
    'rh_precuneus_label.nii.gz';		    % 2025
    'rh_posteriorcingulate_label.nii.gz';	    % 2023
    'rh_isthmuscingulate_label.nii.gz';	    % 2010
    'rh_inferiortemporal_label.nii.gz';	    % 2009
    'rh_middletemporal_label.nii.gz';	    % 2015
    'rh_superiortemporal_label.nii.gz';	    % 2030
    'rh_bankssts_label.nii.gz';		    % 2001
    'rh_fusiform_label.nii.gz';		    % 2007
    'rh_transversetemporal_label.nii.gz';	    % 2034
    'rh_entorhinal_label.nii.gz';		    % 2006
    'rh_temporalpole_label.nii.gz';		    % 2033
    'rh_parahippocampal_label.nii.gz';	    % 2016
    'rh_lateraloccipital_label.nii.gz';	    % 2011
    'rh_lingual_label.nii.gz';		    % 2013
    'rh_cuneus_label.nii.gz';		    % 2005
    'rh_pericalcarine_label.nii.gz'};	    % 2021

% create plain name label
roiLabels = strrep(roiNames, '_label.nii.gz', '');
roiLabels = strrep(roiLabels, '_', '.');

%% start processing

display('Loading fe Structure...');

% load data
load(fe);

% start parallel pool for fefgGet
display(['Opening Parallel Pool with ', num2str(nclust), ' Cores...']);

% create parallel cluster object
c = parcluster;

% set number of cores from arguments
c.NumWorkers = nclust;

% set temporary cache directory
t = tempname(cachePath);

% make cache dir
OK = mkdir(t);

% check and set cachedir location
if OK
    % set local storage for parpool
    c.JobStorageLocation = t;
end

% start parpool
parpool(c, nclust, 'IdleTimeout', 720);

% convert all fiber nodes to voxel coord indices in tensor
display('Finding Nodes...');
fibNodeIndices = fefgGet(fe.fg, 'nodes2voxels', fe.roi.coords); % ~2 hours
fibLength = fefgGet(fe.fg, 'length'); % instantly

%% for loop to find fibers in each ROI

display('Finding ROI fibers...');

% for every roi 
for ii = 1:length(roiNames)
    
    % put ROI name in output
    out{ii}.name = roiNames{ii};
    
    % load ROI
    roi = dtiImportRoiFromNifti(char(fullfile(anatPath, 'label', roiNames{ii})));
    
    % grab size of ROI fibers in output object
    out{ii}.vxsize = size(roi.coords, 1);
    
    % transform roi coordinates into image space for tensor
    vx.coords = mrAnatXformCoords(fe.life.xform.acpc2img, roi.coords);
    
    % round xformed coords like dtiExportRoiToNifti
    vx.coords = ceil(vx.coords);
    % CAUSES OVERLAP IN ROIS AND SUBSEQUENTLY FIBERS
    % USE APARC+ASEG IMAGE DIRECTLY - next improvement
    % OR COMPROMISE FIBER ASSIGNMENT - fix this post hoc
    
    % keep coordinates to speed up fiber check
    out{ii}.vxroi = vx.coords;
    
    % find the transformed ROI coordinates in the fe object and return the
    % indices. These are the fe.roi.coords' indices to find fibers in tensor
    % the coordinate logical (c)index
    vx.cindex = ismember(fe.roi.coords, vx.coords, 'rows');
    
    % convert logical of fe.roi.coords to indices of fe.roi.coords for
    % subsetting fibers in sptensor
    vx.index = find(vx.cindex);
    
    % create sub-index of sptensor for fibers to keep
    [ inds1, ~ ] = find(fe.life.M.Phi(:, vx.index, :));
    
    % This failed one job - why?
    % pull unique intersecting fiber indices into output object
    if size(inds1) > 0
        out{ii}.int.fibers = unique(inds1(:, 3));
    else
        out{ii}.int.fibers = [];
    end
    
    % pull fiber node indices for subset of fibers intersecting ROI
    vx.nodes = fibNodeIndices(out{ii}.int.fibers);
    
    % for every fiber intersecting the ROIs, find the 2 end points
    % able to back-track end points by +/- indices into endpoint cell-array
    
    for jj = 1:length(vx.nodes)
        vx.ep(jj, 1) = vx.nodes{jj}(1);
        vx.ep(jj, 2) = vx.nodes{jj}(length(vx.nodes{jj}));
    end
    
    % redo the logic above w/ endpoint indices
    
    % fiber endpoints of indices found within ROI
    % must do separately to maintain logical order
    vx.ep1indx = ismember(vx.ep(:,1), vx.index, 'rows');
    vx.ep2indx = ismember(vx.ep(:,2), vx.index, 'rows');
    
    % convert logical of endpoint indices to index value of int.fibers
    vx.lep1 = find(vx.ep1indx);
    vx.lep2 = find(vx.ep2indx);
    
    % find unique end point indices
    vx.uep = unique([vx.lep1; vx.lep2]);
    
    % combine unique fiber endpoint indices to id fiber indices
    % index into out{ii}.int.fibers to pull right fiber indices
    out{ii}.end.fibers = out{ii}.int.fibers(vx.uep);
    
    % id non-endpoint intersecting fibers
    out{ii}.nid.fibers = setdiff(out{ii}.int.fibers, out{ii}.end.fibers);
    
    % additional info (fiber length / fiber weights)
    
    % catch output number of nodes in each fiber
    out{ii}.int.lengths = fibLength(out{ii}.int.fibers);
    out{ii}.end.lengths = fibLength(out{ii}.end.fibers);
    out{ii}.nid.lengths = fibLength(out{ii}.nid.fibers);
    
    % weights of fibers
    out{ii}.int.weights = fe.life.fit.weights(out{ii}.int.fibers);
    out{ii}.end.weights = fe.life.fit.weights(out{ii}.end.fibers);
    out{ii}.nid.weights = fe.life.fit.weights(out{ii}.nid.fibers);
    
    % weighted fibers
    
    % weighted fiber intersections
    nzint = out{ii}.int.weights > 0;
    out{ii}.int.nzfibs = out{ii}.int.fibers(nzint);
    out{ii}.int.nzleng = out{ii}.int.lengths(nzint);
    out{ii}.int.nzwght = out{ii}.int.weights(nzint);
    
    % weighted endpoint intersections
    nzend = out{ii}.end.weights > 0;
    out{ii}.end.nzfibs = out{ii}.end.fibers(nzend);
    out{ii}.end.nzleng = out{ii}.end.lengths(nzend);
    out{ii}.end.nzwght = out{ii}.end.weights(nzend);
    
    % weighted non-endpoint intersections
    nzend = out{ii}.nid.weights > 0;
    out{ii}.nid.nzfibs = out{ii}.nid.fibers(nzend);
    out{ii}.nid.nzleng = out{ii}.nid.lengths(nzend);
    out{ii}.nid.nzwght = out{ii}.nid.weights(nzend);
    
    clear roi vx inds1 jj nzint nzend 
end

% clean up workspace
clear ii roi vx inds1 jj nzint nzend

% redo sanity check to see if the correction below is still needed
% % basic (and lazy) sanity checks
% for ii = 1:68
%     % intersect of endpoint in raw should be 100% - all endpoint fibers
%     % should intersect, but very few do.
%     tst1(ii) = length(intersect(out{ii}.int.fibers, out{ii}.end.fibers));
%     tst2(ii) = length(out{ii}.int.fibers);
%     tst3(ii) = length(out{ii}.end.fibers);
%     tst4(ii) = length(out{ii}.nid.fibers);
% end
% 
% tst1 == tst3 % are all endpoints intersecting as well?
% tst2 > tst1  % are there fewer endpoints than intersections?
% (tst2 - tst3) == tst4 % the difference between intersect and end is the same as nid

%% for every unique combinations' intersection, find fiber indices

display('Finding each unique pairs intersections...');

% create combinations
pairs = nchoosek(1:length(roiNames), 2);

% for every unique pair, find intersecting fibers
for ii = 1:length(pairs)
    
    % pull output regions into tmp objects for clarity of code
    reg1 = out{pairs(ii, 1)};
    reg2 = out{pairs(ii, 2)};
    
    % grab roi names
    pconn{ii}.roi1 = strrep(reg1.name, '.nii.gz', '');
    pconn{ii}.roi2 = strrep(reg2.name, '.nii.gz', '');
    
    % do the rois intersect
    pconn{ii}.intersect = size(intersect(reg1.vxroi, reg2.vxroi, 'rows'), 1) > 0;
    
    % grab roi coords
    % this is a lot of info - necessary?
    pconn{ii}.roi1vx = size(reg1.vxroi, 1);
    pconn{ii}.roi2vx = size(reg2.vxroi, 1);
    pconn{ii}.roisvx = [reg1.vxroi; reg2.vxroi];
    
    % find fiber intersections
    [ pconn{ii}.int.fibers, ind, ~ ] = intersect(reg1.int.fibers, reg2.int.fibers);
    pconn{ii}.int.lengths = reg1.int.lengths(ind);
    pconn{ii}.int.weights = reg1.int.weights(ind);
    
    % find non-zero fiber intersections
    [ pconn{ii}.int.nzfibs, ind, ~ ] = intersect(reg1.int.nzfibs, reg2.int.nzfibs);
    pconn{ii}.int.nzleng = reg1.int.nzleng(ind);
    pconn{ii}.int.nzwght = reg1.int.nzwght(ind);
    
    % find endpoint intersections
    [ pconn{ii}.end.fibers, ind, ~ ] = intersect(reg1.end.fibers, reg2.end.fibers);
    pconn{ii}.end.lengths = reg1.end.lengths(ind);
    pconn{ii}.end.weights = reg1.end.weights(ind);
    
    % find non-zero endpoint intersections
    [ pconn{ii}.end.nzfibs, ind, ~ ] = intersect(reg1.end.nzfibs, reg2.end.nzfibs);
    pconn{ii}.end.nzleng = reg1.end.nzleng(ind);
    pconn{ii}.end.nzwght = reg1.end.nzwght(ind);
    
    % find endpoint to non-endpoint intersections
    [ pconn{ii}.nid1.fibers, ind, ~ ] = intersect(reg1.end.fibers, reg2.nid.fibers);
    pconn{ii}.nid1.lengths = reg1.end.lengths(ind);
    pconn{ii}.nid1.weights = reg1.end.weights(ind);
    
    % find non-zero endpoint to non-endpoint intersections
    [ pconn{ii}.nid1.nzfibs, ind, ~ ] = intersect(reg1.end.nzfibs, reg2.nid.nzfibs);
    pconn{ii}.nid1.nzleng = reg1.end.lengths(ind);
    pconn{ii}.nid1.nzwght = reg1.end.weights(ind);
    
    % find endpoint to non-endpoint intersections
    [ pconn{ii}.nid2.fibers, ind, ~ ] = intersect(reg2.end.fibers, reg1.nid.fibers);
    pconn{ii}.nid2.lengths = reg2.end.lengths(ind);
    pconn{ii}.nid2.weights = reg2.end.weights(ind);
    
    % find non-zero endpoint to non-endpoint intersections
    [ pconn{ii}.nid2.nzfibs, ind, ~ ] = intersect(reg2.end.nzfibs, reg1.nid.nzfibs);
    pconn{ii}.nid2.nzleng = reg2.end.lengths(ind);
    pconn{ii}.nid2.nzwght = reg2.end.weights(ind);
    
end

clear reg1 reg2

display('Removing duplicate fibers...');

% remove fibers identified multiple times within each pair
comb = nchoosek(1:size(pairs, 1), 2);

for ii = 1:size(comb, 1)
    [ ~, pt1, pt2 ] = setxor(pconn{comb(ii,1)}.int.fibers, pconn{comb(ii,2)}.int.fibers);
    pconn{comb(ii,1)}.int.fibers = pconn{comb(ii,1)}.int.fibers(pt1);
    pconn{comb(ii,1)}.int.lengths = pconn{comb(ii,1)}.int.lengths(pt1);
    pconn{comb(ii,1)}.int.weights = pconn{comb(ii,1)}.int.weights(pt1);
    pconn{comb(ii,2)}.int.fibers = pconn{comb(ii,2)}.int.fibers(pt2);
    pconn{comb(ii,2)}.int.lengths = pconn{comb(ii,2)}.int.lengths(pt2);
    pconn{comb(ii,2)}.int.weights = pconn{comb(ii,2)}.int.weights(pt2);
    
    [ ~, pt1, pt2 ] = setxor(pconn{comb(ii,1)}.int.nzfibs, pconn{comb(ii,2)}.int.nzfibs);
    pconn{comb(ii,1)}.int.nzfibs = pconn{comb(ii,1)}.int.nzfibs(pt1);
    pconn{comb(ii,1)}.int.nzleng = pconn{comb(ii,1)}.int.nzleng(pt1);
    pconn{comb(ii,1)}.int.nzwght = pconn{comb(ii,1)}.int.nzwght(pt1);
    pconn{comb(ii,2)}.int.nzfibs = pconn{comb(ii,2)}.int.nzfibs(pt2);
    pconn{comb(ii,2)}.int.nzleng = pconn{comb(ii,2)}.int.nzleng(pt2);
    pconn{comb(ii,2)}.int.nzwght = pconn{comb(ii,2)}.int.nzwght(pt2);
    
    % find indices of each non overlapping end fiber and keep just the unique indices
    [ ~, pt1, pt2 ] = setxor(pconn{comb(ii,1)}.end.fibers, pconn{comb(ii,2)}.end.fibers);
    pconn{comb(ii,1)}.end.fibers = pconn{comb(ii,1)}.end.fibers(pt1);
    pconn{comb(ii,1)}.end.lengths = pconn{comb(ii,1)}.end.lengths(pt1);
    pconn{comb(ii,1)}.end.weights = pconn{comb(ii,1)}.end.weights(pt1);
    pconn{comb(ii,2)}.end.fibers = pconn{comb(ii,2)}.end.fibers(pt2);
    pconn{comb(ii,2)}.end.lengths = pconn{comb(ii,2)}.end.lengths(pt2);
    pconn{comb(ii,2)}.end.weights = pconn{comb(ii,2)}.end.weights(pt2);
    
    [ ~, pt1, pt2 ] = setxor(pconn{comb(ii,1)}.end.nzfibs, pconn{comb(ii,2)}.end.nzfibs);
    pconn{comb(ii,1)}.end.nzfibs = pconn{comb(ii,1)}.end.nzfibs(pt1);
    pconn{comb(ii,1)}.end.nzleng = pconn{comb(ii,1)}.end.nzleng(pt1);
    pconn{comb(ii,1)}.end.nzwght = pconn{comb(ii,1)}.end.nzwght(pt1);
    pconn{comb(ii,2)}.end.nzfibs = pconn{comb(ii,2)}.end.nzfibs(pt2);
    pconn{comb(ii,2)}.end.nzleng = pconn{comb(ii,2)}.end.nzleng(pt2);
    pconn{comb(ii,2)}.end.nzwght = pconn{comb(ii,2)}.end.nzwght(pt2);
    
    % find indices of each non overlapping nid1 fiber and keep just the unique indices
    [ ~, pt1, pt2 ] = setxor(pconn{comb(ii,1)}.nid1.fibers, pconn{comb(ii,2)}.nid1.fibers);
    pconn{comb(ii,1)}.nid1.fibers = pconn{comb(ii,1)}.nid1.fibers(pt1);
    pconn{comb(ii,1)}.nid1.lengths = pconn{comb(ii,1)}.nid1.lengths(pt1);
    pconn{comb(ii,1)}.nid1.weights = pconn{comb(ii,1)}.nid1.weights(pt1);
    pconn{comb(ii,2)}.nid1.fibers = pconn{comb(ii,2)}.nid1.fibers(pt2);
    pconn{comb(ii,2)}.nid1.lengths = pconn{comb(ii,2)}.nid1.lengths(pt2);
    pconn{comb(ii,2)}.nid1.weights = pconn{comb(ii,2)}.nid1.weights(pt2);
    
    [ ~, pt1, pt2 ] = setxor(pconn{comb(ii,1)}.nid1.nzfibs, pconn{comb(ii,2)}.nid1.nzfibs);
    pconn{comb(ii,1)}.nid1.nzfibs = pconn{comb(ii,1)}.nid1.nzfibs(pt1);
    pconn{comb(ii,1)}.nid1.nzleng = pconn{comb(ii,1)}.nid1.nzleng(pt1);
    pconn{comb(ii,1)}.nid1.nzwght = pconn{comb(ii,1)}.nid1.nzwght(pt1);
    pconn{comb(ii,2)}.nid1.nzfibs = pconn{comb(ii,2)}.nid1.nzfibs(pt2);
    pconn{comb(ii,2)}.nid1.nzleng = pconn{comb(ii,2)}.nid1.nzleng(pt2);
    pconn{comb(ii,2)}.nid1.nzwght = pconn{comb(ii,2)}.nid1.nzwght(pt2);
    
    % find indices of each non overlapping nid2 fiber and keep just the unique indices
    [ ~, pt1, pt2 ] = setxor(pconn{comb(ii,1)}.nid2.fibers, pconn{comb(ii,2)}.nid2.fibers);
    pconn{comb(ii,1)}.nid2.fibers = pconn{comb(ii,1)}.nid2.fibers(pt1);
    pconn{comb(ii,1)}.nid2.lengths = pconn{comb(ii,1)}.nid2.lengths(pt1);
    pconn{comb(ii,1)}.nid2.weights = pconn{comb(ii,1)}.nid2.weights(pt1);
    pconn{comb(ii,2)}.nid2.fibers = pconn{comb(ii,2)}.nid2.fibers(pt2);
    pconn{comb(ii,2)}.nid2.lengths = pconn{comb(ii,2)}.nid2.lengths(pt2);
    pconn{comb(ii,2)}.nid2.weights = pconn{comb(ii,2)}.nid2.weights(pt2);
    
    [ ~, pt1, pt2 ] = setxor(pconn{comb(ii,1)}.nid2.nzfibs, pconn{comb(ii,2)}.nid2.nzfibs);
    pconn{comb(ii,1)}.nid2.nzfibs = pconn{comb(ii,1)}.nid2.nzfibs(pt1);
    pconn{comb(ii,1)}.nid2.nzleng = pconn{comb(ii,1)}.nid2.nzleng(pt1);
    pconn{comb(ii,1)}.nid2.nzwght = pconn{comb(ii,1)}.nid2.nzwght(pt1);
    pconn{comb(ii,2)}.nid2.nzfibs = pconn{comb(ii,2)}.nid2.nzfibs(pt2);
    pconn{comb(ii,2)}.nid2.nzleng = pconn{comb(ii,2)}.nid2.nzleng(pt2);
    pconn{comb(ii,2)}.nid2.nzwght = pconn{comb(ii,2)}.nid2.nzwght(pt2);
end

display('Create output matrices...');

% emat is endpoints (fixed), imat is intersection (old) 
% [e|i]mat will be 16:
% count/dens/leng/denLeng
% nzLiFE of count/dens/leng/denLeng
% wghtSum/wghtWght/wghtDens/wghtLeng
% SOE/EMD/VL1/VL2
% 
% nmat is 1 asym-matrix of endpoints to non-endpoint intersection
%

% run virtual - catch relevant output; ~ 3 hours
% run virtual - catch relevant output; ~ 3 hours
display('Running Virtual Lesion...');

for ii = 1:length(pconn)
    try
        [ ewVL, ewoVL ] = feComputeVirtualLesion(fe, pconn{ii}.end.nzfibs);
        efevl{ii} = feComputeEvidenceFix(ewoVL, ewVL);
        
    catch  
        efevl{ii}.s.mean = 0;
        efevl{ii}.em.mean = 0;
        efevl{ii}.j.mean = 0;
        efevl{ii}.kl.mean = 0;
    end

end

% create empty output matrix
emat = zeros(size(roiNames, 1), size(roiNames, 1), 16);
nmat = zeros(size(roiNames, 1), size(roiNames, 1), 1);

% create matrix outputs
for ii = 1:length(pairs)
    
    % 01. fiber count
    efc = length(pconn{ii}.end.fibers);
    emat(pairs(ii, 1), pairs(ii, 2), 1) = efc;
    emat(pairs(ii, 2), pairs(ii, 1), 1) = efc;
    
    % 02. fiber density
    efc = 2*efc / size(pconn{ii}.roisvx, 1);
    emat(pairs(ii, 1), pairs(ii, 2), 2) = efc;
    emat(pairs(ii, 2), pairs(ii, 1), 2) = efc;
    
    % 03. fiber length
    efl = (1/efc) * sum(pconn{ii}.end.lengths);
    emat(pairs(ii, 1), pairs(ii, 2), 3) = efl;
    emat(pairs(ii, 2), pairs(ii, 1), 3) = efl;
    
    % 04. fiber density * length
    efdl = (2 / size(pconn{ii}.roisvx, 1)) * sum(1 / pconn{ii}.end.lengths);
    emat(pairs(ii, 1), pairs(ii, 2), 4) = efdl;
    emat(pairs(ii, 2), pairs(ii, 1), 4) = efdl;
    
    % 05. nz fiber count
    wefc = length(pconn{ii}.end.nzfibs);
    emat(pairs(ii, 1), pairs(ii, 2), 5) = wefc;
    emat(pairs(ii, 2), pairs(ii, 1), 5) = wefc;
    
    % 06. nz fiber density
    wefd = 2*wefc / size(pconn{ii}.roisvx, 1);
    emat(pairs(ii, 1), pairs(ii, 2), 6) = wefd;
    emat(pairs(ii, 2), pairs(ii, 1), 6) = wefd;
    
    % 07. nz fiber length
    wefl = (1/wefc) * sum(pconn{ii}.end.nzleng);
    emat(pairs(ii, 1), pairs(ii, 2), 7) = wefl;
    emat(pairs(ii, 2), pairs(ii, 1), 7) = wefl;
    
    % 08. nz fiber density * length
    wefdl = (2 / size(pconn{ii}.roisvx, 1)) * sum(1 / pconn{ii}.end.nzleng);
    emat(pairs(ii, 1), pairs(ii, 2), 8) = wefdl;
    emat(pairs(ii, 2), pairs(ii, 1), 8) = wefdl;
    
    % 09. weights sum    
    wesum = sum(pconn{ii}.end.nzwght);
    emat(pairs(ii, 1), pairs(ii, 2), 9) = wesum;
    emat(pairs(ii, 2), pairs(ii, 1), 9) = wesum;
    
    % 10. weighted count
    wewgh = (1/wefc) * sum(pconn{ii}.end.nzwght);
    emat(pairs(ii, 1), pairs(ii, 2), 10) = wewgh;
    emat(pairs(ii, 2), pairs(ii, 1), 10) = wewgh;
    
    % 11. weighted density
    weden = (1/wefd) * sum(pconn{ii}.end.nzwght);
    emat(pairs(ii, 1), pairs(ii, 2), 11) = weden;
    emat(pairs(ii, 2), pairs(ii, 1), 11) = weden;
    
    % 12. weighted length
    wwefl = (1/wefl) * sum(pconn{ii}.end.nzwght);
    emat(pairs(ii, 1), pairs(ii, 2), 12) = wwefl;
    emat(pairs(ii, 2), pairs(ii, 1), 12) = wwefl;
    
    % 13. SOE
    emat(pairs(ii, 1), pairs(ii, 2), 13) = efevl{ii}.s.mean;
    emat(pairs(ii, 2), pairs(ii, 1), 13) = efevl{ii}.s.mean;
    
    % 14. EMD
    emat(pairs(ii, 1), pairs(ii, 2), 14) = efevl{ii}.em.mean;
    emat(pairs(ii, 2), pairs(ii, 1), 14) = efevl{ii}.em.mean;
    
    % 15. Jeffery's Divergence    
    emat(pairs(ii, 1), pairs(ii, 2), 15) = efevl{ii}.j.mean;
    emat(pairs(ii, 2), pairs(ii, 1), 15) = efevl{ii}.j.mean;
    
    % 16. Kullback-Leibler Divergence
    emat(pairs(ii, 1), pairs(ii, 2), 16) = efevl{ii}.kl.mean;
    emat(pairs(ii, 2), pairs(ii, 1), 16) = efevl{ii}.kl.mean;
    
    % create matrix output for endpoint to non-endpoint connectivity
    nmat(pairs(ii, 1), pairs(ii, 2), 1) = length(pconn{ii}.nid1.fibers);
    nmat(pairs(ii, 2), pairs(ii, 1), 1) = length(pconn{ii}.nid2.fibers);
    
    % clear temporary objects
    clear efc efd efl efdl wefc wefd wefl wefdl wesum wewgh weden wwefl

end

% fix zeros / nan
for ii = 1:size(emat, 3)

    elen = emat(:,:,ii);
    elen(isnan(elen)) = 0;
    emat(:,:,ii) = elen;
    
    clear elen
end

clear elen

esoe = emat(:,:,13);
esoe(esoe < 0) = 0;
emat(:,:,13) = esoe;

ejd = emat(:,:,15);
ejd(ejd < 0) = 0;
emat(:,:,15) = ejd;

ekl = emat(:,:,16);
ekl(ekl < 0) = 0;
emat(:,:,16) = ekl;

clear esoe ejd ekl

% % create summary plot
% plab = {'Fiber Count', 'Fiber Density', 'Fiber Length', 'Fiber Density x Length', ...
%         'Weighted Fiber Count', 'Weighted Fiber Density', 'Weighted Fiber Length', 'Weighted Fiber Density x Length', ...
%         'Sum of Weights', 'Weights / Count', 'Weights / Density', 'Weights / Length', ...
%         'Strength of Evidence', 'Earth Movers Distance', 'Jeffery''s Divergence', 'Kullback-Leibler'};
% 
% fh = figure;
% for kk = 1:size(emat, 3)
%     subplot(4, 4, kk);
%     colormap('hot');
%     imagesc(log(emat(:,:,kk)));
%     title(plab{kk});
%     xlabel('FS DK Regions');
%     ylabel('FS DK Regions');
%     axis('square'); axis('equal'); axis('tight');
%     set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
%     y = colorbar;
%     ylabel(y, 'Strength of Connection');
% end
% 
% figure;
% imagesc(nmat);
% colormap('hot');
% xlabel('FS DK Regions');
% ylabel('FS DK Regions');
% title('Endpoint to Non-Endpoint Intersections');
% axis('square'); axis('equal'); axis('tight');
% set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
% y = colorbar;
% ylabel(y, 'Strength of Connection');

% calculate network statistics
display('Calculating Network Statistics...');

for ii = 1:size(emat, 3)
    [ fns{ii}, pval{ii}, fdat{ii} ] = feFiberNetworkStats(emat(:,:,ii), 100, 0.75, 18);
end

display('Saving data...');

save(outPath, 'emat', 'out', 'pconn', 'nmat', 'fns', 'pval', 'fdat', 'roiNames', 'roiLabels');

end

