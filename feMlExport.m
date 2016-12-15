function [ out ] = feMlExport(dgrp, subj, dmdl, lmax, indx)
%% load edge data into data frames for export into ML project

if indx == 1
    nam = 'count';
end

if indx == 2
    nam = 'density';
end

if indx == 14
    nam = 'emd';
end

outfile = ['ml_data/' dgrp '_' subj '_' dmdl '_' lmax '_' nam '.csv'];

% load subjects data
[ ~, emat ] = feMergeRepeats(dgrp, subj, dmdl, lmax);

% pull upper diagonal for all repeats 

for ii = 1:size(emat{indx}, 3)
    tmp = emat{indx}(:,:,ii);

    % linearly index in order of pconn the unique connections
    out(ii,:) = tmp(find(~triu(ones(size(tmp))))); 
    
end

clear ii tmp

dlmwrite(outfile, out, ',');

end