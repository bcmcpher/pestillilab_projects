%% procedure to loop over and create individual ROIs for tracking
% will need to fix paths once I confirm it works

% subjects to loop
%subj = '102311';

subjs = {'102311' '109123' '158035' '352738' '770352' '108323' '131217' '200614' '573249' '910241'};

for subj = 1:length(subjs)
    
    % build paths for each subject
    fsdir = ['/N/dc2/projects/lifebid/HCP/Brent/7t_rerun/freesurfer/7t_' subjs{subj}];
    anatdir = ['/N/dc2/projects/lifebid/HCP7/' subjs{subj} '/anatomy'];
    t1ref = 'T1w_acpc_dc_restore_1.05.nii.gz';
    
    % build file output paths
    aparcPath = [fsdir '/mri/aparc+aseg.mgz'];
    anatPath = [anatdir '/' t1ref];
    
    % output is moved to different folder for sanity sake
    outDir = ['/N/dc2/projects/lifebid/HCP/Brent/7t_rerun/' subjs{subj} '/label'];
    
    % make nifti ROIs
    fs_roisFromAllLabels(aparcPath, outDir, 'nii', anatPath);
    
    % fix names
    cd(outDir);
    !for i in *.nii.gz; do mv $i ${i/*ctx-/}; done

end

