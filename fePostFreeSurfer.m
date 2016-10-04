function [ out ] = fePostFreeSurfer(subj, subjdir)
%fePostFreeSurfer - initial MATLAB run post-freesufer to create files for
% tracking / network construction
% 

%% build paths

% build paths for each subject
% if ~notDefined(subjdir)
%     subjdir = getenv('SUBJECTS_DIR');
% end

fsdir = fullfile(subjdir, subj);
%anatdir = ['/N/dc2/projects/lifebid/HCP7/' subjs{subj} '/anatomy'];
%t1ref = 'T1w_acpc_dc_restore_1.05.nii.gz'; % sshould this be T1.mgz?

% build file output paths
aparcPath = fullfile(fsdir, 'mri', 'aparc.mgz');
%anatPath = fullfile(fsdir, 'mri', 'T1.mgz');
%anatPath = '/N/dc2/projects/lifebid/glue/subjects/100132/T1_denoised_100132_acpc.nii.gz';

% output is moved to different folder for sanity sake
%outDir = fullfile(fsdir, 'netrois');

%% create network ROIs

% create labels from aparc image
fs_annotationToLabelFiles(subj, 'aparc', []);

% find all label files and create cell arrays with file names for all steps
labelFileNames = dir(fullfile(fsdir, 'label','*.label'));
labelRoisNames = cell(length(labelFileNames),1);
niiRoiFullPath = cell(length(labelFileNames),1);
matRoiFullPath = cell(length(labelFileNames),1);

for il = 1:length(labelFileNames)
    labelRoisNames{il} = labelFileNames(il).name;
    niftiRoiName = labelRoisNames{il};
    niftiRoiName(niftiRoiName=='.') = '_';
    niiRoiFullPath{il} = fullfile(fsdir, 'label', niftiRoiName);
    matRoiFullPath{il} = fullfile(fsdir, 'label', niftiRoiName);
end

% debugging the function that is supposed to work
fs_labelFileToNiftiRoi('100132', fullfile(fsdir, 'label', 'lh.bankssts.label'), 'lh.bankssts.label', ...
                       'lh', [], 0, '/N/dc2/projects/lifebid/glue/fstest')

% iterate over labels files, convert to .nii.gz and save, save to vistasoft
% roi .mat, save both files smoothed w/ 3mm^3 kernel.
for il = 1:length(labelFileNames)
    if ~(exist(niiRoiFullPath{il},'file')==2)
        fs_labelFileToNiftiRoi(subj,labelRoisNames{il}, niiRoiFullPath{il}, labelFileNames(il).name(1:2),[],[],fsdir);
    else
        fprintf('[%s] Found ROI, skipping: \n%s\n',mfilename,niiRoiFullPath{il})
    end
    if ~(exist(matRoiFullPath{il},'file')==2)
        dtiImportRoiFromNifti(niiRoiFullPath{il}, matRoiFullPath{il});
    else
        fprintf('[%s] Found ROI, skipping: \n%s\n',mfilename,niiRoiFullPath{il})
    end
end

% make nifti ROIs
%fs_roisFromAllLabels(aparcPath, outDir, 'nii', anatPath);

% fix names
%cd(outDir);
%!for i in outDir/*.nii.gz; do mv $i ${i/*ctx-/}; done

%% create wm tracking mask

% wmMaskFile = fullfile(fsdir, 'wm_mask.nii.gz');
% 
% fs_wm = aparcPath;
% eval(sprintf('!mri_convert  --out_orientation RAS %s %s', , wmMaskFile));
% wm = niftiRead(wmMaskFile);
% invals  = [2 41 16 17 28 60 51 53 12 52 13 18 54 50 11 251 252 253 254 255 10 49 46 7];
% origvals = unique(wm.data(:));
% fprintf('\n[%s] Converting voxels... ',mfilename);
% wmCounter=0;noWMCounter=0;
% for ii = 1:length(origvals);
%     if any(origvals(ii) == invals)
%         wm.data( wm.data == origvals(ii) ) = 1;
%         wmCounter=wmCounter+1;
%     else
%         wm.data( wm.data == origvals(ii) ) = 0;
%         noWMCounter = noWMCounter + 1;
%     end
% end
% fprintf('converted %i regions to White-matter (%i regions left outside of WM)\n\n',wmCounter,noWMCounter);
% niftiWrite(wm);

out = 'done';

end

