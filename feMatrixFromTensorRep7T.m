function [ emat, nmat, fns, pval, fdat, imat, time ] = feMatrixFromTensorRep7T(subj, dmdl, lmax, rep, nclust) 
%% development of tensor based matrix construction
% Brent McPherson
% 20160710
% 
% [ z1, z2, z3, z4, z5, z6, z7 ] = feMatrixFromTensorRep7T(2, 'detr', '2', 1, 16);
tic;

% subject index - 1:10
% HCP Data Order: '102311', '108323', '109123', '131217', '158035', '200614', '352738', '573249', '770352', '910241'
% Current HCP Data Order: '108323', '109123', '131217', '910241'
file.name.indx = subj;

% run 'stn' or 'hcp' data
%dataGroup = dgrp;

% run 'prob', 'detr', or 'tens' connectomes
dataModel = dmdl;

% define lmax: 2, 4, 6, 8
lmx = ['lmax' num2str(lmax)];

%% create paths to data

display('Building File Structure...')

% define subject ID(s)
%file.hcp.name.dir = {'105115', '110411', '111312', '113619'};
%file.hc7.name.dir = {'102311', '108323', '109123', '131217', '158035', '200614', '352738', '573249', '770352', '910241'};
file.hc7.name.dir = {'108323', '109123', '131217', '910241'};
%file.hcp.name.params = {'SD_PROB', 'SD_STREAM', 'wm_tensor_'};
file.hc7.name.params = {'SD_PROB', 'SD_STREAM', 'tensor'};
file.name.lmax = lmx;

% output directory for saved network results
file.out.dir = '/N/dc2/projects/lifebid/HCP/Brent/cogs610/7t_reps_data';

% build subject names
%file.hcp.name.subj = file.hcp.name.dir;
file.hc7.name.subj = file.hc7.name.dir{file.name.indx};

% define paths to folders
file.dir.path = '/N/dc2/projects/lifebid/code/ccaiafa/Caiafa_Pestilli_paper2015/Results/ETC_Dec2015/Single_TC/';
file.dir.h7path = '/N/dc2/projects/lifebid/HCP7';
%file.stn.dir.anat = fullfile('/N/dc2/projects/lifebid/2t1/anatomy/', file.stn.name.subj);
%file.hcp.dir.anat = fullfile('/N/dc2/projects/lifebid/2t1/anatomy/', file.hcp.name.subj);
file.hc7.dir.anat = fullfile(file.dir.h7path, file.hc7.name.subj, 'anatomy');

% anat file names / repeats
file.name.anat = 'T1w_acpc_dc_restore_1.05.nii.gz';
file.name.atlas = 'aparc+aseg.nii.gz';
file.name.iter = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10'};

% % hcp file names
% file.hcp.name.prob.sngl = strcat('fe_structure_', file.hcp.name.subj(file.name.indx), '_STC_run01_', file.hcp.name.params{1},'_', file.name.lmax, '.mat');
% file.hcp.name.detr.sngl = strcat('fe_structure_', file.hcp.name.subj(file.name.indx), '_STC_run01_', file.hcp.name.params{2},'_', file.name.lmax, '.mat');
% file.hcp.name.tens.sngl = strcat('fe_structure_', file.hcp.name.subj(file.name.indx), '_STC_run01_', file.hcp.name.params{3}, '.mat');

% hcp 7t reps
file.hc7.name.prob.reps = strcat('fe_structure_', file.hc7.name.subj, '_STC_run01_', file.hc7.name.params{1},'_', file.name.lmax, '_connNUM', file.name.iter, '.mat');
file.hc7.name.detr.reps = strcat('fe_structure_', file.hc7.name.subj, '_STC_run01_', file.hc7.name.params{2},'_', file.name.lmax, '_connNUM', file.name.iter, '.mat');
file.hc7.name.tens.reps = strcat('fe_structure_', file.hc7.name.subj, '_STC_run01_', file.hc7.name.params{3},'_', '_connNUM', file.name.iter, '.mat');

% fullfile path to files / inputs / outputs
file.hc7.path.anat = char(fullfile(file.hc7.dir.anat, file.name.anat));
file.hc7.path.prob = char(fullfile(file.dir.path, strcat('7t_', file.hc7.name.subj), file.hc7.name.prob.reps{rep}));
file.hc7.path.detr = char(fullfile(file.dir.path, strcat('7t_', file.hc7.name.subj), file.hc7.name.detr.reps{rep}));
file.hc7.path.tens = char(fullfile(file.dir.path, strcat('7t_', file.hc7.name.subj), file.hc7.name.tens.reps{rep}));
%file.dir.anat = fullfile(file.hc7.dir.anat, file.name.anat);
file.hc7.path.anat = fullfile(file.hc7.dir.anat, file.name.anat);
file.hc7.path.atlas = fullfile(file.hc7.dir.anat, file.name.atlas);

%% pull all ROI NIfTI files of FreeSurfer labels

%display('Loading ROIs...');

% pull all freesurfer label files
% switch dataGroup
%     case 'stn'
%         %roiNames = dir(fullfile(file.stn.dir.anat{file.name.indx}, 'label', '*_label.nii.gz'));
%         file.dir.anat = file.stn.dir.anat{file.name.indx};
%     case 'hcp'
%         %roiNames = dir(fullfile(file.hcp.dir.anat{file.name.indx}, 'label', '*_label.nii.gz'));
%         
%     otherwise
%         error('No valid data group selected');
% end

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
    'lh-superiorfrontal.nii.gz';          % 1028
    'lh-rostralmiddlefrontal.nii.gz';     % 1027
    'lh-caudalmiddlefrontal.nii.gz';      % 1003
    'lh-parsopercularis.nii.gz'; 	    % 1018
    'lh-parsorbitalis.nii.gz'; 	    % 1019
    'lh-parstriangularis.nii.gz'; 	    % 1020
    'lh-lateralorbitofrontal.nii.gz';     % 1012
    'lh-medialorbitofrontal.nii.gz'; 	    % 1014
    'lh-precentral.nii.gz'; 		    % 1024
    'lh-paracentral.nii.gz'; 		    % 1017
    'lh-frontalpole.nii.gz'; 		    % 1032
    'lh-rostralanteriorcingulate.nii.gz'; % 1026
    'lh-caudalanteriorcingulate.nii.gz';  % 1002
    'lh-insula.nii.gz'; 		    % 1035
    'lh-superiorparietal.nii.gz'; 	    % 1029
    'lh-inferiorparietal.nii.gz'; 	    % 1008
    'lh-supramarginal.nii.gz'; 	    % 1031
    'lh-postcentral.nii.gz'; 		    % 1022
    'lh-precuneus.nii.gz'; 		    % 1025
    'lh-posteriorcingulate.nii.gz'; 	    % 1023
    'lh-isthmuscingulate.nii.gz'; 	    % 1010
    'lh-inferiortemporal.nii.gz'; 	    % 1009
    'lh-middletemporal.nii.gz'; 	    % 1015
    'lh-superiortemporal.nii.gz'; 	    % 1030
    'lh-bankssts.nii.gz'; 		    % 1001
    'lh-fusiform.nii.gz'; 		    % 1007
    'lh-transversetemporal.nii.gz'; 	    % 1034
    'lh-entorhinal.nii.gz'; 		    % 1006
    'lh-temporalpole.nii.gz'; 	    % 1033
    'lh-parahippocampal.nii.gz'; 	    % 1016
    'lh-lateraloccipital.nii.gz'; 	    % 1011
    'lh-lingual.nii.gz'; 		    % 1013
    'lh-cuneus.nii.gz'; 		    % 1005
    'lh-pericalcarine.nii.gz';	    % 1021
    'rh-superiorfrontal.nii.gz'; 	    % 2028
    'rh-rostralmiddlefrontal.nii.gz';	    % 2027
    'rh-caudalmiddlefrontal.nii.gz';	    % 2003
    'rh-parsopercularis.nii.gz';	    % 2018
    'rh-parsorbitalis.nii.gz';	    % 2019
    'rh-parstriangularis.nii.gz';	    % 2020
    'rh-lateralorbitofrontal.nii.gz';	    % 2012
    'rh-medialorbitofrontal.nii.gz';	    % 2014
    'rh-precentral.nii.gz';		    % 2024
    'rh-paracentral.nii.gz';		    % 2017
    'rh-frontalpole.nii.gz';		    % 2032
    'rh-rostralanteriorcingulate.nii.gz'; % 2026
    'rh-caudalanteriorcingulate.nii.gz';  % 2002
    'rh-insula.nii.gz';		    % 2035
    'rh-superiorparietal.nii.gz';	    % 2029
    'rh-inferiorparietal.nii.gz';	    % 2008
    'rh-supramarginal.nii.gz';	    % 2031
    'rh-postcentral.nii.gz';		    % 2022
    'rh-precuneus.nii.gz';		    % 2025
    'rh-posteriorcingulate.nii.gz';	    % 2023
    'rh-isthmuscingulate.nii.gz';	    % 2010
    'rh-inferiortemporal.nii.gz';	    % 2009
    'rh-middletemporal.nii.gz';	    % 2015
    'rh-superiortemporal.nii.gz';	    % 2030
    'rh-bankssts.nii.gz';		    % 2001
    'rh-fusiform.nii.gz';		    % 2007
    'rh-transversetemporal.nii.gz';	    % 2034
    'rh-entorhinal.nii.gz';		    % 2006
    'rh-temporalpole.nii.gz';		    % 2033
    'rh-parahippocampal.nii.gz';	    % 2016
    'rh-lateraloccipital.nii.gz';	    % 2011
    'rh-lingual.nii.gz';		    % 2013
    'rh-cuneus.nii.gz';		    % 2005
    'rh-pericalcarine.nii.gz'};	    % 2021

% create plain name label
roiLabel = strrep(roiNames, '.nii.gz', '');
roiLabel = strrep(roiLabel, '-', '.');

% label values in aparc+aseg file that correspond to ROIs
% roiIndices = {...
%     1028;
%     1027;
%     1003;
%     1018;
%     1019;
%     1020;
%     1012;
%     1014;
%     1024;
%     1017;
%     1032;
%     1026;
%     1002;
%     1035;
%     1029;
%     1008;
%     1031;
%     1022;
%     1025;
%     1023;
%     1010;
%     1009;
%     1015;
%     1030;
%     1001;
%     1007;
%     1034;
%     1006;
%     1033;
%     1016;
%     1011;
%     1013;
%     1005;
%     1021;
%     2028;
%     2027;
%     2003;
%     2018;
%     2019;
%     2020;
%     2012;
%     2014;
%     2024;
%     2017;
%     2032;
%     2026;
%     2002;
%     2035;
%     2029;
%     2008;
%     2031;
%     2022;
%     2025;
%     2023;
%     2010;
%     2009;
%     2015;
%     2030;
%     2001;
%     2007;
%     2034;
%     2006;
%     2033;
%     2016;
%     2011;
%     2013;
%     2005;
%     2021};

%% start processing

display('Loading fe Structure...');

% load data
% load saved fe structure
switch dataModel
    case 'prob'
        load(file.hc7.path.prob);
    case 'detr'
        load(file.hc7.path.detr);
    case 'tens'
        load(file.hc7.path.tens);
    otherwise
        error('No valid connectome model selected');
end

% start parallel pool for fefgGet
display(['Opening Parallel Pool with ', num2str(nclust), ' Cores']);

% create parallel cluster object
c = parcluster;

% set number of cores from arguments
c.NumWorkers = nclust;

% set temporary cache directory
t = tempname('/N/dc2/projects/lifebid/HCP/Brent/cogs610/matlab/parCache');

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

% % load fe anatomy
% feAnat = readFileNifti(fe.path.anatomy);
% 
% % load aparc+aseg.nii.gz
% aparc = readFileNifti(file.hc7.path.atlas);
% 
% % find transform from aparc space to fe anatomy space
% 
% % apply transform to image for alignment to ACPC from fe anatomy
% acpc_aparc = niftiApplyXform(aparc, fe.life.xform.acpc2img);
% % not working... can't permute / resample image
% 
% % go ROI by ROI - like I currently am
% [ x, y, z ] = ind2sub(size(aparc.data), find(aparc.data == roiIndices{1}));
% 
% % merge into ROI
% roi = [x, y, z];
% 
% % convert coordinates to acpc space
% acpc_roi = mrAnatXformCoords(fe.life.xform.acpc2img, roi);
% acpc_roi = ceil(acpc_roi);
% 
% % keep only uniqe voxel indices
% acpc_roi = unique(acpc_roi, 'rows');
% 
% % check if transformed ROI is in fe.roi
% sum(ismember(fe.roi.coords, acpc_roi, 'rows'))
% % not working... indices way out of range

% iterate over ROI values to find indices of intersecting fibers

% for every roi 
for ii = 1:length(roiNames)
    
    % put ROI name in output
    out{ii}.name = roiNames{ii};
    
    % load ROI
    roi = dtiImportRoiFromNifti(char(fullfile('/N/dc2/projects/lifebid/HCP/Brent/7t_rerun', file.hc7.name.subj, 'label', roiNames{ii})));
    
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
    
    % pull unique intersecting fiber indices into output object
    if size(inds1, 2) > 2
        out{ii}.int.fibers = unique(inds1(:, 3));
    else
        out{ii}.int.fibers = [];
    end
    
    % pull fiber node indices for subset of fibers intersecting ROI
    vx.nodes = fibNodeIndices(out{ii}.int.fibers);
    
    % for every fiber intersecting the ROIs, find the 2 end points
    % able to back-track end points by +/- indices into endpoint cell-array
    
    % there are fibers
    if ~isempty(vx.nodes)
        for jj = 1:length(vx.nodes)
            vx.ep(jj, 1) = vx.nodes{jj}(1);
            vx.ep(jj, 2) = vx.nodes{jj}(length(vx.nodes{jj}));
        end
    else % there are not fibers, set to zero
        vx.ep(1, 1) = 0;
        vx.ep(1, 2) = 0;
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
parfor ii = 1:length(pairs)
    
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

% CHECK IF THIS IS NECESSARY BEFORE RUNNING - MAY BE RESOLVED BY USING SINGLE ROI FILE
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

% create empty output matrix
emat = zeros(size(roiNames, 1), size(roiNames, 1), 16);
imat = zeros(size(roiNames, 1), size(roiNames, 1), 16);
nmat = zeros(size(roiNames, 1), size(roiNames, 1), 1);

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
        [ iwVL, iwoVL ] = feComputeVirtualLesion(fe, pconn{ii}.int.nzfibs);
        ifevl{ii} = feComputeEvidenceFix(iwoVL, iwVL);
        
    catch 
        
        ifevl{ii}.s.mean = 0;
        ifevl{ii}.em.mean = 0;
        ifevl{ii}.j.mean = 0;
        ifevl{ii}.kl.mean = 0;
        
    end
    
    try
        [ ewVL, ewoVL ] = feComputeVirtualLesion(fe, pconn{ii}.end.nzfibs);
        efevl{ii} = feComputeEvidenceFix(ewoVL, ewVL);
        
    catch  
        efevl{ii}.s.mean = 0;
        efevl{ii}.em.mean = 0;
        efevl{ii}.j.mean = 0;
        efevl{ii}.kl.mean = 0;
    end
    
    %clear iwVL iwoVL
    %clear ewVL ewoVL
    
end

% create empty output matrix
emat = zeros(size(roiNames, 1), size(roiNames, 1), 16);
imat = zeros(size(roiNames, 1), size(roiNames, 1), 16);
nmat = zeros(size(roiNames, 1), size(roiNames, 1), 1);

% create matrix outputs
for ii = 1:length(pairs)
    
    % 01. fiber count
    ifc = length(pconn{ii}.int.fibers);
    imat(pairs(ii, 1), pairs(ii, 2), 1) = ifc;
    imat(pairs(ii, 2), pairs(ii, 1), 1) = ifc;

    efc = length(pconn{ii}.end.fibers);
    emat(pairs(ii, 1), pairs(ii, 2), 1) = efc;
    emat(pairs(ii, 2), pairs(ii, 1), 1) = efc;
    
    % 02. fiber density
    ifd = 2*ifc / size(pconn{ii}.roisvx, 1);
    imat(pairs(ii, 1), pairs(ii, 2), 2) = ifd;
    imat(pairs(ii, 2), pairs(ii, 1), 2) = ifd;

    efc = 2*efc / size(pconn{ii}.roisvx, 1);
    emat(pairs(ii, 1), pairs(ii, 2), 2) = efc;
    emat(pairs(ii, 2), pairs(ii, 1), 2) = efc;
    
    % 03. fiber length
    ifl = (1/ifc) * sum(pconn{ii}.int.lengths);
    imat(pairs(ii, 1), pairs(ii, 2), 3) = ifl;
    imat(pairs(ii, 2), pairs(ii, 1), 3) = ifl;

    efl = (1/efc) * sum(pconn{ii}.end.lengths);
    emat(pairs(ii, 1), pairs(ii, 2), 3) = efl;
    emat(pairs(ii, 2), pairs(ii, 1), 3) = efl;
    
    % 04. fiber density * length
    ifdl = (2 / size(pconn{ii}.roisvx, 1)) * sum(1 / pconn{ii}.int.lengths);
    imat(pairs(ii, 1), pairs(ii, 2), 4) = ifdl;
    imat(pairs(ii, 2), pairs(ii, 1), 4) = ifdl;

    efdl = (2 / size(pconn{ii}.roisvx, 1)) * sum(1 / pconn{ii}.end.lengths);
    emat(pairs(ii, 1), pairs(ii, 2), 4) = efdl;
    emat(pairs(ii, 2), pairs(ii, 1), 4) = efdl;
    
    % 05. nz fiber count
    wifc = length(pconn{ii}.int.nzfibs);
    imat(pairs(ii, 1), pairs(ii, 2), 5) = wifc;
    imat(pairs(ii, 2), pairs(ii, 1), 5) = wifc;

    wefc = length(pconn{ii}.end.nzfibs);
    emat(pairs(ii, 1), pairs(ii, 2), 5) = wefc;
    emat(pairs(ii, 2), pairs(ii, 1), 5) = wefc;
    
    % 06. nz fiber density
    wifd = 2*wifc / size(pconn{ii}.roisvx, 1);
    imat(pairs(ii, 1), pairs(ii, 2), 6) = wifd;
    imat(pairs(ii, 2), pairs(ii, 1), 6) = wifd;

    wefd = 2*wefc / size(pconn{ii}.roisvx, 1);
    emat(pairs(ii, 1), pairs(ii, 2), 6) = wefd;
    emat(pairs(ii, 2), pairs(ii, 1), 6) = wefd;
    
    % 07. nz fiber length
    wifl = (1/wifc) * sum(pconn{ii}.int.nzleng);
    imat(pairs(ii, 1), pairs(ii, 2), 7) = wifl;
    imat(pairs(ii, 2), pairs(ii, 1), 7) = wifl;

    wefl = (1/wefc) * sum(pconn{ii}.end.nzleng);
    emat(pairs(ii, 1), pairs(ii, 2), 7) = wefl;
    emat(pairs(ii, 2), pairs(ii, 1), 7) = wefl;
    
    % 08. nz fiber density * length
    wifdl = (2 / size(pconn{ii}.roisvx, 1)) * sum(1 / pconn{ii}.int.nzleng);
    imat(pairs(ii, 1), pairs(ii, 2), 8) = wifdl;
    imat(pairs(ii, 2), pairs(ii, 1), 8) = wifdl;

    wefdl = (2 / size(pconn{ii}.roisvx, 1)) * sum(1 / pconn{ii}.end.nzleng);
    emat(pairs(ii, 1), pairs(ii, 2), 8) = wefdl;
    emat(pairs(ii, 2), pairs(ii, 1), 8) = wefdl;
    
    % 09. weights sum
    wisum = sum(pconn{ii}.int.nzwght);
    imat(pairs(ii, 1), pairs(ii, 2), 9) = wisum;
    imat(pairs(ii, 2), pairs(ii, 1), 9) = wisum;
    
    wesum = sum(pconn{ii}.end.nzwght);
    emat(pairs(ii, 1), pairs(ii, 2), 9) = wesum;
    emat(pairs(ii, 2), pairs(ii, 1), 9) = wesum;
    
    % 10. weighted count
    wiwgh = (1/wifc) * sum(pconn{ii}.int.nzwght);
    imat(pairs(ii, 1), pairs(ii, 2), 10) = wiwgh;
    imat(pairs(ii, 2), pairs(ii, 1), 10) = wiwgh;
    
    wewgh = (1/wefc) * sum(pconn{ii}.end.nzwght);
    emat(pairs(ii, 1), pairs(ii, 2), 10) = wewgh;
    emat(pairs(ii, 2), pairs(ii, 1), 10) = wewgh;
    
    % 11. weighted density
    widen = (1/wifd) * sum(pconn{ii}.int.nzwght);
    imat(pairs(ii, 1), pairs(ii, 2), 11) = widen;
    imat(pairs(ii, 2), pairs(ii, 1), 11) = widen;
    
    weden = (1/wefd) * sum(pconn{ii}.end.nzwght);
    emat(pairs(ii, 1), pairs(ii, 2), 11) = weden;
    emat(pairs(ii, 2), pairs(ii, 1), 11) = weden;
    
    % 12. weighted length
    wwifl = (1/wifl) * sum(pconn{ii}.int.nzwght);
    imat(pairs(ii, 1), pairs(ii, 2), 12) = wwifl;
    imat(pairs(ii, 2), pairs(ii, 1), 12) = wwifl;
    
    wwefl = (1/wefl) * sum(pconn{ii}.end.nzwght);
    emat(pairs(ii, 1), pairs(ii, 2), 12) = wwefl;
    emat(pairs(ii, 2), pairs(ii, 1), 12) = wwefl;
    
    % 13. SOE
    imat(pairs(ii, 1), pairs(ii, 2), 13) = ifevl{ii}.s.mean;
    imat(pairs(ii, 2), pairs(ii, 1), 13) = ifevl{ii}.s.mean;
    
    emat(pairs(ii, 1), pairs(ii, 2), 13) = efevl{ii}.s.mean;
    emat(pairs(ii, 2), pairs(ii, 1), 13) = efevl{ii}.s.mean;
    
    % 14. EMD
    imat(pairs(ii, 1), pairs(ii, 2), 14) = ifevl{ii}.em.mean;
    imat(pairs(ii, 2), pairs(ii, 1), 14) = ifevl{ii}.em.mean;
    
    emat(pairs(ii, 1), pairs(ii, 2), 14) = efevl{ii}.em.mean;
    emat(pairs(ii, 2), pairs(ii, 1), 14) = efevl{ii}.em.mean;
    
    % 15. Jeffery's Divergence
    imat(pairs(ii, 1), pairs(ii, 2), 15) = ifevl{ii}.j.mean;
    imat(pairs(ii, 2), pairs(ii, 1), 15) = ifevl{ii}.j.mean;
    
    emat(pairs(ii, 1), pairs(ii, 2), 15) = efevl{ii}.j.mean;
    emat(pairs(ii, 2), pairs(ii, 1), 15) = efevl{ii}.j.mean;
    
    % 16. Kullback-Leibler Divergence
    imat(pairs(ii, 1), pairs(ii, 2), 16) = ifevl{ii}.kl.mean;
    imat(pairs(ii, 2), pairs(ii, 1), 16) = ifevl{ii}.kl.mean;
    
    emat(pairs(ii, 1), pairs(ii, 2), 16) = efevl{ii}.kl.mean;
    emat(pairs(ii, 2), pairs(ii, 1), 16) = efevl{ii}.kl.mean;
    
    % create matrix output for endpoint to non-endpoint connectivity
    nmat(pairs(ii, 1), pairs(ii, 2), 1) = length(pconn{ii}.nid1.fibers);
    nmat(pairs(ii, 2), pairs(ii, 1), 1) = length(pconn{ii}.nid2.fibers);
    
    % clear temporary objects
    clear ifc ifd ifl ifdl wifc wifd wifl wifdl wisum wiwgh widen wwifl
    clear efc efd efl efdl wefc wefd wefl wefdl wesum wewgh weden wwefl

end

% fix zeros / nan
for ii = 1:size(emat, 3)
    ilen = imat(:,:,ii);
    ilen(isnan(ilen)) = 0;
    imat(:,:,ii) = ilen;
    
    elen = emat(:,:,ii);
    elen(isnan(elen)) = 0;
    emat(:,:,ii) = elen;
    
    clear ilen elen
end

clear ilen elen

isoe = imat(:,:,13);
isoe(isoe < 0) = 0;
imat(:,:,13) = isoe;

esoe = emat(:,:,13);
esoe(esoe < 0) = 0;
emat(:,:,13) = esoe;

ijd = imat(:,:,15);
ijd(ijd < 0) = 0;
imat(:,:,15) = ijd;

ejd = emat(:,:,15);
ejd(ejd < 0) = 0;
emat(:,:,15) = ejd;

ikl = imat(:,:,16);
ikl(ikl < 0) = 0;
imat(:,:,16) = ikl;

ekl = emat(:,:,16);
ekl(ekl < 0) = 0;
emat(:,:,16) = ekl;

clear isoe esoe ijd ejd ikl ekl

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

%dataGroup = 'hcp';
% switch dataGroup
%     case 'stn'
%         saveNames = char(strcat('7T_', file.stn.name.dir(file.name.indx), '_', dataModel, '_', file.name.lmax, '_rep', file.name.iter{rep}, '.mat'));
%     case 'hcp'
%         saveNames = char(strcat('7T_', file.hc7.name.dir(file.name.indx), '_', dataModel, '_', file.name.lmax, '_rep', file.name.iter{rep}, '.mat'));
%     otherwise
%         error('No valid data group selected');
% end

saveNames = char(strcat('7T_', file.hc7.name.dir(file.name.indx), '_', dataModel, '_', file.name.lmax, '_rep', file.name.iter{rep}, '.mat'));

time = toc;

save(fullfile(file.out.dir, saveNames), 'roiNames', 'out', 'pconn', 'emat', 'imat', 'nmat', 'fns', 'pval', 'fdat', 'time');

end

