%% 
% Brent McPherson
% 20160101
% script to run / test the ensemble import
%

ens_tck = '/N/dc2/projects/lifebid/HCP/Brent/vss-2016/mrtrix/ensemble_tracks';
fout = '/N/dc2/projects/lifebid/HCP/Brent/vss-2016/mrtrix/feImportTest.mat';

tic;
[ fgOut, dict ] = feImportEnsembleRoi2Roi(ens_tck, fout);
toc;