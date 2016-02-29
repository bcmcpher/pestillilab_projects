%% parse dictionary object to index output for Virtual Lesion values (em/st)
% 

% pull ROI labels
[ uLabel1, iLabel1, uIndex1 ] = unique({dict(:).roi1}, 'sorted');
[ uLabel2, iLabel2, uIndex2 ] = unique({dict(:).roi2}, 'sorted');

% merge labels into complete list
labs = unique([uLabel1 uLabel2], 'sorted');

% create empty output objects
em = zeros(length(labs), length(labs));
st = zeros(length(labs), length(labs));

% create list of uniqe pairs of indices and numeric indices
pairs = nchoosek(labs, 2);
pairn = nchoosek(1:length(labs), 2);

% for each pair of indices, select the ROI pair and add the value to the
% output matrix
for ii = 1:length(pairs)
    
    % index each pair of combinations and find overlapping entry
    apr1 = strmatch(pairs(ii, 1), {dict(:).roi1});
    apr2 = strmatch(pairs(ii, 2), {dict(:).roi2});
    aint = intersect(apr1, apr2);
    
    % if the ROI combination isn't found, flip the ROIs and search again
    if isempty(aint)
        
        % flip which is ROI1 / ROI2 because they are not necessarily run in that order
        bpr1 = strmatch(pairs(ii, 2), {dict(:).roi1});
        bpr2 = strmatch(pairs(ii, 1), {dict(:).roi2});
        bint = intersect(bpr1, bpr2);
        
        % create index in output object for ROI pair
        indx = bint;
    
    else
        
        % create index in output object for ROI pair
        indx = aint;
        
    end
    
    % fill in matrices
    em(pairn(ii, 1), pairn(ii, 2)) = SE{indx}.em.mean;
    em(pairn(ii, 2), pairn(ii, 1)) = SE{indx}.em.mean;
    
    st(pairn(ii, 1), pairn(ii, 2)) = SE{indx}.s.mean;
    st(pairn(ii, 2), pairn(ii, 1)) = SE{indx}.s.mean;

end

% save matrices
dlmwrite('../../mrtrix/em_20160118.csv', em);
dlmwrite('../../mrtrix/st_20160118.csv', st);

% write out necessary results
save ../../mrtrix/matrix_20160119.mat labs em st

%% create some more detailed matrices to make sure these are ok

% set proportional / absolute thresholding
em_out = threshold_absolute(em, 7);
em_out = threshold_proportional(em, 0.75);

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

% clustering coefficient
em_ccoef = clustering_coef_wu(em_tha);
em_clen = clustering_coef_wu(em_tha);

% rich club / transitivity coefficient
em_rc = rich_club_wu(em_tha);
em_tr = transitivity_wu(em_tha);

% global efficiency
em_glb = efficiency_wei(em_tha);
em_loc = efficiency_wei(em_tha, 1);



% plot figure
h = figure('name','ROI-to-ROI EMD Connectome Matrix', 'color','w');
colormap('hot');
hax = axes;
imagesc(em_out);
hold on;
title('Thresheld EMD Connectome Matrix');
set(gca, 'YTick', 1:1:length(labs), 'YTickLabel', labs,'tickdir','out');
xlabel('FreeSurfer Defined Anatomical Regions');
ylabel('FreeSurfer Defined Anatomical Regions');
axis('square'); axis('equal'); axis('tight');
y = colorbar;
ylabel(y, 'EMD')
SP = size(em_out, 1) / 2;
line([SP SP], get(hax, 'YLim'), 'Color', [0 0 1]);
line(get(hax, 'XLim'), [SP SP], 'Color', [0 0 1]);
line([SP*2 0], [SP*2 0], 'Color', [0 0 1]);

h = figure('name','ROI-to-ROI SOE Connectome Matrix', 'color','w');
colormap('hot');
hax = axes;
imagesc(st_out);
hold on;
title('Thresheld SOE Connectome Matrix');
set(gca, 'YTick', 1:1:length(labs), 'YTickLabel', labs,'tickdir','out');
xlabel('FreeSurfer Defined Anatomical Regions');
ylabel('FreeSurfer Defined Anatomical Regions');
axis('square'); axis('equal'); axis('tight');
y = colorbar;
ylabel(y, 'SOE')
SP = size(st_out, 1) / 2;
line([SP SP], get(hax, 'YLim'), 'Color', [0 0 1]);
line(get(hax, 'XLim'), [SP SP], 'Color', [0 0 1]);
line([SP*2 0], [SP*2 0], 'Color', [0 0 1]);
