%% edit of life_scripts function stub to make CC mask
% This script makes the white-matter mask used to track the connectomes in
% Pestilli et al., LIFE paper.
%
% Copyright Brent McPherson & Franco Pestilli (c) Stanford University, 2016


%% just use fsl...

% first make a mask of the ventricle b/c inflation takes too much ventricle IMO
% - make left/right separately, merge filling holes 
% - invert so eveything but ventricles is left out
% then mask CC, inflate 7 mask out ventricles, binarize and mask

% ventricle mask
% fslmaths aseg.nii.gz -thr 4 -uthr 5 -bin lh.wbMask
% fslmaths aseg.nii.gz -thr 43 -uthr 44 -bin rh.wbMask
% fslmaths lh.wbMask -add rh.wbMask -fillh -binv wbMask

% smaller inflation excluding ventricles
% fslmaths aseg.nii.gz -thr 251 -uthr 255 -kernel sphere 7 -dilD -mas wbMask.nii.gz -bin cc_test

% then resample to subject space - generic file names, no paths
% flirt -in norm.nii.gz -ref dwi_data -omat fs2dwi -out fs2dwi
% flirt -in fs_cc.nii.gz -ref dwi_data -init fs2dwi -applyxfm -out dwi_cc -interp nearestneighbour

% convert to mif
% mrconvert dwi_cc.nii.gz dwi_data_b2000_aligned_trilin_cc.mif

% add w/ wm mask for inflated wm mask and convert to mif
% fslmaths dwi_cc.nii.gz -add wm_aseg.nii.gz -bin wm_cc_inf
% mrconvert wm_cc_inf.nii.gz dwi_data_b2000_aligned_trilin_wm_cc.mif


%%

outpath = '/N/dc2/projects/lifebid/HCP/Brent/vss-2016/mrtrix';
outwm = fullfile(outpath, 'wmMask_aseg.nii.gz');
outcc = fullfile(outpath, 'wmMask_cc.nii.gz');
outcw = fullfile(outpath, 'wmMask_cw.nii.gz');

% Get the base directory for the data
anatomypath = '/N/dc2/projects/lifebid/HCP/Brent/anatomy';
subject = '105115';

% load the labeled aseg file
wmAsegFile = fullfile(anatomypath, subject, 'mri', 'aseg.nii.gz');
wm = niftiRead(wmAsegFile);

% copy the wm mask for isolating the corpus callosum and change output names
wm.fname = outwm;

cc = wm;
cc.fname = outcc;

cw = wm;
cw.fname = outcw;

% index of wm values to keep
invals = [ 2 41 16 17 28 60 51 53 12 52 13 18 54 50 11 251 252 253 254 255 10 49 46 7 ];
ccindx = [ 251 252 253 254 255 ];

% pull the original values
origvals = unique(wm.data(:));

% initialize counters for variables
wmCounter = 0;
noWMCounter = 0;
wmCCcounter = 0;
nwmCCcounter = 0;

% for every wmMask value
for ii = 1:length(origvals);
    
    % if it's in invals, it's masked as 1
    if any(origvals(ii) == invals)
        wm.data( wm.data == origvals(ii) ) = 1;
        %wmCounter = wmCounter+1;
    
    % if it's not in invals, it's masked as 0
    else
        wm.data( wm.data == origvals(ii) ) = 0;
        %noWMCounter = noWMCounter + 1;
    end
    
    % create cc mask from labels
    if any(origvals(ii) == ccindx)
        cc.data(cc.data == origvals(ii)) = 1;
        %wmCCcounter = wmCCcounter + 1;
    else
        cc.data(cc.data == origvals(ii)) = 0;
        %nwmCCcounter = nwmCCounter + 1;
    end

end
   
% smooth to enlarge CC
smcc = smooth3(cc.data, 'gaussian', [ 3 3 3 ]);
cc.data = smcc;

% merge back into combined mask


% save the files
niftiWrite(wm);
niftiWrite(cc);
%niftiWrite(cw);

