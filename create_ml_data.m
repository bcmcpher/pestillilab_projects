
% subject data structure to loop over
data_subj{1} = {'hcp', '105115'};
data_subj{2} = {'hcp', '110411'};
data_subj{3} = {'hcp', '111312'};
data_subj{4} = {'hcp', '113619'};
data_subj{5} = {'stn', 'FP'};
data_subj{6} = {'stn', 'HT'};
data_subj{7} = {'stn', 'KK'};
data_subj{8} = {'stn', 'MP'};
data_subj{9} = {'7T', '108323'};
data_subj{10} = {'7T', '109123'};
data_subj{11} = {'7T', '131217'};
data_subj{12} = {'7T', '910241'};

% write out every subjects data set
for ii = 1:length(data_subj)
    tmp1 = feMlExport(data_subj{ii}{1}, data_subj{ii}{2}, 'prob', '8', 1);
    tmp2 = feMlExport(data_subj{ii}{1}, data_subj{ii}{2}, 'prob', '8', 2);
    tmp3 = feMlExport(data_subj{ii}{1}, data_subj{ii}{2}, 'prob', '8', 14);
end

clear ii tmp1 tmp2 tmp3

%% load a pconn and create lobe nested region labels

% list of ROI names used 
roiNames = {'lh_superiorfrontal_label.nii.gz'; 
    'lh_rostralmiddlefrontal_label.nii.gz'; 
    'lh_caudalmiddlefrontal_label.nii.gz'; 
    'lh_parsopercularis_label.nii.gz'; 
    'lh_parsorbitalis_label.nii.gz'; 
    'lh_parstriangularis_label.nii.gz'; 
    'lh_lateralorbitofrontal_label.nii.gz'; 
    'lh_medialorbitofrontal_label.nii.gz'; 
    'lh_precentral_label.nii.gz'; 
    'lh_paracentral_label.nii.gz'; 
    'lh_frontalpole_label.nii.gz'; 
    'lh_rostralanteriorcingulate_label.nii.gz'; 
    'lh_caudalanteriorcingulate_label.nii.gz'; 
    'lh_insula_label.nii.gz'; 
    'lh_superiorparietal_label.nii.gz'; 
    'lh_inferiorparietal_label.nii.gz'; 
    'lh_supramarginal_label.nii.gz'; 
    'lh_postcentral_label.nii.gz'; 
    'lh_precuneus_label.nii.gz'; 
    'lh_posteriorcingulate_label.nii.gz'; 
    'lh_isthmuscingulate_label.nii.gz'; 
    'lh_inferiortemporal_label.nii.gz'; 
    'lh_middletemporal_label.nii.gz'; 
    'lh_superiortemporal_label.nii.gz'; 
    'lh_bankssts_label.nii.gz'; 
    'lh_fusiform_label.nii.gz'; 
    'lh_transversetemporal_label.nii.gz'; 
    'lh_entorhinal_label.nii.gz'; 
    'lh_temporalpole_label.nii.gz'; 
    'lh_parahippocampal_label.nii.gz'; 
    'lh_lateraloccipital_label.nii.gz'; 
    'lh_lingual_label.nii.gz'; 
    'lh_cuneus_label.nii.gz'; 
    'lh_pericalcarine_label.nii.gz'; 
    'rh_superiorfrontal_label.nii.gz'; 
    'rh_rostralmiddlefrontal_label.nii.gz'; 
    'rh_caudalmiddlefrontal_label.nii.gz'; 
    'rh_parsopercularis_label.nii.gz'; 
    'rh_parsorbitalis_label.nii.gz'; 
    'rh_parstriangularis_label.nii.gz'; 
    'rh_lateralorbitofrontal_label.nii.gz'; 
    'rh_medialorbitofrontal_label.nii.gz'; 
    'rh_precentral_label.nii.gz'; 
    'rh_paracentral_label.nii.gz'; 
    'rh_frontalpole_label.nii.gz'; 
    'rh_rostralanteriorcingulate_label.nii.gz'; 
    'rh_caudalanteriorcingulate_label.nii.gz'; 
    'rh_insula_label.nii.gz'; 
    'rh_superiorparietal_label.nii.gz'; 
    'rh_inferiorparietal_label.nii.gz'; 
    'rh_supramarginal_label.nii.gz'; 
    'rh_postcentral_label.nii.gz'; 
    'rh_precuneus_label.nii.gz'; 
    'rh_posteriorcingulate_label.nii.gz'; 
    'rh_isthmuscingulate_label.nii.gz'; 
    'rh_inferiortemporal_label.nii.gz'; 
    'rh_middletemporal_label.nii.gz'; 
    'rh_superiortemporal_label.nii.gz'; 
    'rh_bankssts_label.nii.gz'; 
    'rh_fusiform_label.nii.gz'; 
    'rh_transversetemporal_label.nii.gz'; 
    'rh_entorhinal_label.nii.gz'; 
    'rh_temporalpole_label.nii.gz'; 
    'rh_parahippocampal_label.nii.gz'; 
    'rh_lateraloccipital_label.nii.gz'; 
    'rh_lingual_label.nii.gz'; 
    'rh_cuneus_label.nii.gz'; 
    'rh_pericalcarine_label.nii.gz'};

% re-create paired order of ROIs
pairs = nchoosek(roiNames, 2);

% create label names
for ii = 1:length(pairs)
    label{ii} = [pairs{ii, 1} '_to_' pairs{ii, 2}];
end

% correct spelling for variables - use as nesting for repeats?
rois_label = strrep(label, '_label.nii.gz', '');

% correct labels to hemispheres
lobe_label = strrep(rois_label, 'superiorfrontal', 'frontal');
lobe_label = strrep(lobe_label, 'rostralmiddlefrontal', 'frontal');
lobe_label = strrep(lobe_label, 'caudalmiddlefrontal', 'frontal');
lobe_label = strrep(lobe_label, 'parsopercularis', 'frontal');
lobe_label = strrep(lobe_label, 'parsorbitalis', 'frontal');
lobe_label = strrep(lobe_label, 'parstriangularis', 'frontal');
lobe_label = strrep(lobe_label, 'lateralorbitofrontal', 'frontal');
lobe_label = strrep(lobe_label, 'medialorbitofrontal', 'frontal');
lobe_label = strrep(lobe_label, 'precentral', 'frontal');
lobe_label = strrep(lobe_label, 'paracentral', 'frontal');
lobe_label = strrep(lobe_label, 'frontalpole', 'frontal');
lobe_label = strrep(lobe_label, 'rostralanteriorcingulate', 'frontal');
lobe_label = strrep(lobe_label, 'caudalanteriorcingulate', 'frontal');
lobe_label = strrep(lobe_label, 'insula', 'frontal');
lobe_label = strrep(lobe_label, 'superiorparietal', 'parietal');
lobe_label = strrep(lobe_label, 'inferiorparietal', 'parietal');
lobe_label = strrep(lobe_label, 'supramarginal', 'parietal');
lobe_label = strrep(lobe_label, 'postcentral', 'parietal');
lobe_label = strrep(lobe_label, 'precuneus', 'parietal');
lobe_label = strrep(lobe_label, 'posteriorcingulate', 'parietal');
lobe_label = strrep(lobe_label, 'isthmuscingulate', 'parietal');
lobe_label = strrep(lobe_label, 'inferiortemporal', 'temporal');
lobe_label = strrep(lobe_label, 'middletemporal', 'temporal');
lobe_label = strrep(lobe_label, 'superiortemporal', 'temporal');
lobe_label = strrep(lobe_label, 'bankssts', 'temporal');
lobe_label = strrep(lobe_label, 'fusiform', 'temporal');
lobe_label = strrep(lobe_label, 'transversetemporal', 'temporal');
lobe_label = strrep(lobe_label, 'entorhinal', 'temporal');
lobe_label = strrep(lobe_label, 'temporalpole', 'temporal');
lobe_label = strrep(lobe_label, 'parahippocampal', 'temporal');
lobe_label = strrep(lobe_label, 'lateraloccipital', 'occipital');
lobe_label = strrep(lobe_label, 'lingual', 'occipital');
lobe_label = strrep(lobe_label, 'cuneus', 'occipital');
lobe_label = strrep(lobe_label, 'pericalcarine', 'occipital');

%unique(lobe_label)';

% pull left / right / between labels too?
hemi_label = regexprep(lobe_label, 'lh*_to_lh*', 'left');
hemi_label = strrep(hemi_label, 'rh*_to_rh*', 'right');
hemi_label = strrep(hemi_label, 'lh*_to_rh*', 'btwn');
hemi_label = strrep(hemi_label, 'rh*_to_lh*', 'btwn');


% save labels for clustering / groups
% because why would it be user friendly? value-added my ass...
%dlmwrite('ml_data/00_col_names.csv', rois_label, ',');
%dlmwrite('ml_data/00_col_lobes.csv', lobe_label, ',');

 fid = fopen('ml_data/00_col_names.csv','wt');
 if fid > 0
     for k = 1:size(rois_label, 1)
         fprintf(fid,'%s\n',rois_label{k,:});
     end
     fclose(fid);
 end

fid = fopen('ml_data/00_col_lobes.csv','wt');
 if fid > 0
     for k = 1:size(lobe_label, 1)
         fprintf(fid,'%s\n',lobe_label{k,:});
     end
     fclose(fid);
 end
