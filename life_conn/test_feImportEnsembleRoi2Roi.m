%% 
% Brent McPherson
% 20160101
% script to run / test the ensemble import
%

ens_tck = '/N/dc2/projects/lifebid/HCP/Brent/vss-2016/mrtrix/ensemble_tracks2';
fout = '/N/dc2/projects/lifebid/HCP/Brent/vss-2016/mrtrix/feImportWB2.mat';

%tic;
%[ fg, dict ] = feImportEnsembleRoi2Roi(ens_tck, fout);
%toc;

% create dict of fiber counts
%dict = makeDict(ens_tck);

[fg, dict] = feImportEnsembleDict(ens_tck, fout);

%% import and combine whole brain fibers

wbens = dir('/N/dc2/projects/lifebid/HCP/Brent/vss-2016/mrtrix/tck_SD_*.tck');
wbout = fullfile('/N/dc2/projects/lifebid/HCP/Brent/vss-2016/mrtrix', 'whole_brain_200K_ensemble_fibers');

% Strip out the file name.
[ ~, wb ] = fileparts(wbout);
wbOut = dtiNewFiberGroup(wb);

for iwbfg = 1:length(wbens)
    
    display([ 'Importing .tck file: ', wbens(iwbfg).name ]);
    
    % import whole brain fiber group .tck file to matlab
    tfg = dtiImportFibersMrtrix(fullfile('/N/dc2/projects/lifebid/HCP/Brent/vss-2016/mrtrix', wbens(iwbfg).name));
    
    % merge each fiber group into final one
    wbOut = dtiMergeFiberGroups(wbOut, tfg, wbout);
    
end

% save whole brain fiber group
fgWrite(wbOut, wbout, 'mat');

%% Build LiFE Model

% build paths
file.dir.base    = '/N/dc2/projects/lifebid/HCP/Brent/vss-2016/mrtrix';
file.dir.anatomy = fullfile(file.dir.base, 'anatomy');
file.dir.fe      = file.dir.base;

file.name.fe = 'fe_test_new_life_connection_matrix';
file.name.dwi = 'dwi_data_b2000_aligned_trilin_dwi.nii.gz';
file.name.anatomy = 'T1w_acpc_dc_restore_1p25.nii.gz';
file.name.subj = '105115';
file.name.fg = 'feImportWB2.mat';

file.path.dwi     = fullfile(file.dir.base, file.name.dwi);
file.path.anatomy = fullfile(file.dir.anatomy, file.name.anatomy);
file.path.fg      = fullfile(file.dir.base, file.name.fg);

%% initialize and fit LiFE model

% create merge whole brain fiber file
wbfile = fullfile(file.dir.base, 'whole_brain_200k_ensemble_with_ROI2ROI_ensemble_fibers.mat');

% merge whole brain and ROI to ROI fibers
wbfg = dtiMergeFiberGroups(fg, wbOut, wbfile);

% pull random sample of 500000 - breaks dictionary
fg_full = fg;
fg.fibers = randsample(fg.fibers, 500000);

% run just the ROI to ROI generated fibers
fe = feConnectomeInit(file.path.dwi, fg, file.name.fg, file.dir.fe, file.path.dwi, file.path.anatomy);
fe = feSet(fe,'fit',feFitModel(feGet(fe,'mfiber'),feGet(fe,'dsigdemeaned'),'bbnnls'));

% run the whole brain fiber group merged with the ROI to ROI fibers
wbfe = feConnectomeInit(file.path.dwi, wbfg, file.name.fg, file.dir.fe, file.path.dwi, file.path.anatomy);
wbfe = feSet(wbfe,'fit',feFitModel(feGet(wbfe,'mfiber'),feGet(wbfe,'dsigdemeaned'),'bbnnls'));

%% Perform on Virtual Lesion per ROI-pair

parpool;

% clear out params?
fe.fg.params = {};
wbfe.fg.params = {};

% for each ROI pair
for iRoi = 1:length(dict)
    
    % empty fiber index for each ROI
    edge_indices = double.empty();
    
    % for each tracking type
    for iFile = 1:10
        
        % use the size of ifib to determine tracking properties
        tmp = dict(iRoi).ifib{iFile};
        
        % if there's only 1 number, go to next iteration
        if size(tmp, 2) == 1
            
            continue
        
        else
            
            % combine indices for ROI files with fibers
            edge_indices = [ edge_indices tmp(1):tmp(2) ];
            
        end
        
    end
    
    %display(['Perform Virtual Lesion on: ', num2str(length(edge_indices)), ' fibers']);
    %display(['Edge Indices: ', num2str(edge_indices)]);
    
    if length(edge_indices) > 1
        
        try
            SE{iRoi} = feVirtualLesion(fe, edge_indices');
        catch
            warning('Virtual Lesion may have removed all the fibers in the path neightborhood');
            SE{iRoi}.em.mean = 0;
            SE{iRoi}.s.mean = 0;
            SE{iRoi}.brokeVL = 1;
        end
        
    else
        
        SE{iRoi}.em.mean = 0;
        SE{iRoi}.s.mean = 0;
        SE{iRoi}.brokeVL = 0;
    
    end

end

%% index dictionary to get inter- / intra-hemispheric fiber counts 

% internal counters to keep lengths of outputs right
interc = 1;
intrac = 1;
intraL = 1;
intraR = 1;

% split into inter / intra connections
for kk = 1:length(dict)
    
    % pull the temp rois
    roi1 = dict(kk).roi1;
    roi2 = dict(kk).roi2;
    
    % pull the first index
    roi1 = roi1(1);
    roi2 = roi2(1);
        
    % make the indices of within / between
    if strcmp(roi1, roi2)
        intra(intrac) = kk;
        intrac = intrac + 1;
    else
        inter(interc) = kk;
        interc = interc + 1;
    end
    
end
    
% split intra connections into left / right
for kk = 1:length(intra)
    
    % pull the temp rois
    roi1 = dict(intra(kk)).roi1;
    roi2 = dict(intra(kk)).roi2;
    
    % pull the first index
    roi1 = roi1(1);
    roi2 = roi2(1);
    
    % if both rois are left in within group, add to left
    if (strcmp(roi1, 'l') && strcmp(roi2, 'l'))
        intLft(intraL) = intra(kk);
        intraL = intraL + 1;
    end
    
    % if both rois are right in within group, add to right
    if (strcmp(roi1, 'r') && strcmp(roi2, 'r'))
        intRgt(intraR) = intra(kk);
        intraR = intraR + 1;
    end
    
end
    

% sum counts of inter hemispheric connections
ninter = 0;
for kk = 1:length(inter)
    
    tmp = sum([dict(inter(kk)).nfib{:}]);
    ninter = tmp + ninter;
    
end
    
% sum counts of intra hemispheric connections
nintra = 0;
for kk = 1:length(intra)
    
    tmp = sum([dict(intra(kk)).nfib{:}]);
    nintra = tmp + nintra;
    
end

% sum counts of left hemispheric connections
nintLt = 0;
for kk = 1:length(intLft)
    
    tmp = sum([dict(intLft(kk)).nfib{:}]);
    nintLt = tmp + nintLt;
    
end

% sum counts of right hemispheric connections
nintRt = 0;
for kk = 1:length(intRgt)
    
    tmp = sum([dict(intRgt(kk)).nfib{:}]);
    nintRt = tmp + nintRt;
    
end




