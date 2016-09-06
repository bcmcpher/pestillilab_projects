%% find what repeats weren't made
% Brent McPherson
% 20160829
%

% read in all files
files = dir('/N/dc2/projects/lifebid/HCP/Brent/cogs610/reps_data/*.mat');

% find the names
names = {files.name};

% build full list of file names
subj = {'hcp_105115', 'hcp_110411', 'hcp_111312', 'hcp_113619', 'stn_FP', 'stn_HT', 'stn_KK', 'stn_MP'};
fsbj = {'hcp_subj1', 'hcp_subj2', 'hcp_subj3', 'hcp_subj4', 'stn_subj1', 'stn_subj3', 'stn_subj4', 'stn_subj2'};
trck = {'detr', 'prob'};
lmax = {'2', '4', '6', '8', '10', '12'};
reps = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10'};

% run order
%file.stn.name.dir = {'FP', 'MP', 'HT', 'KK', 'KW'};
%file.hcp.name.dir = {'105115', '110411', '111312', '113619'};

iter = 1;
for ii = 1:8 % subjects / dataset
    for jj = 1:2 % detr / prob
        for kk = 1:6 % lmax
            for ll = 1:10 % reps
                all_names{iter} = [ subj{ii} '_' trck{jj} '_lmax' lmax{kk} '_rep' reps{ll} '.mat' ];
                moab_name{iter} = [ fsbj{ii} '_' trck{jj} '_lmax' lmax{kk} '_rep' reps{ll} '.moab' ];
                iter = iter + 1;
            end
        end
    end
end

clear ii jj kk ll iter

% jobs to run
[ nrun, indx ] = setdiff(all_names, names);

% replace names to find qsub scripts
nrun = strrep(nrun, 'hcp_105115', 'hcp_subj1');
nrun = strrep(nrun, 'hcp_110411', 'hcp_subj2');
nrun = strrep(nrun, 'hcp_111312', 'hcp_subj3');
nrun = strrep(nrun, 'hcp_113619', 'hcp_subj4');
nrun = strrep(nrun, 'stn_FP', 'stn_subj1');
nrun = strrep(nrun, 'stn_MP', 'stn_subj2');
nrun = strrep(nrun, 'stn_HT', 'stn_subj3');
nrun = strrep(nrun, 'stn_KK', 'stn_subj4');
nrun = strrep(nrun, 'lmax2', 'lmax02');
nrun = strrep(nrun, 'lmax4', 'lmax04');
nrun = strrep(nrun, 'lmax6', 'lmax06');
nrun = strrep(nrun, 'lmax8', 'lmax08');
nrun = strrep(nrun, '.mat', '.moab');

% nrun is now the .moab jobs to rerun on BRII
fid = fopen('job_rerun.txt', 'w') ;
for ii = 1:length(nrun);
    fprintf(fid, '%s\n', nrun{ii}) ;
end
fclose(fid);



dlmwrite('fix_jobs.txt', nrun);



