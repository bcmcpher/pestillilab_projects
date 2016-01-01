%%
% Brent McPherson
% 20151231
% import series of .tck files and create a key of # of fibers for each file
%

%% read out necessary unique values

% read in all the files
files = dir('ensemble_tracks');

% parse info out of file name
tmp = files(15).name;

% find indices of underscores
und = strfind(tmp, '_');

% find recon algorithm
rc = tmp(und(2)+1:und(3)-1);

% find curvature
cr = tmp(und(3)+1:und(4)-1);

%% split roi names

% isolate roi descriptions and trim file exention
rois = tmp(und(4)+1:length(tmp)-4);

% find middle index and separate rois
splt = strfind(rois, '_to_');
roi1 = rois(1:splt-1);
roi2 = rois(splt+4:length(rois));

%% import fiber group

% import fiber group, save as .mat
fg = dtiImportFibersMrtrix(fullfile('ensemble_tracks',tmp));
fgWrite(fg, fg.name, 'mat');

%% create dictionary

dict.file = tmp;
dict.algo = rc;
dict.curv = cr; 
dict.roi1 = roi1;
dict.roi2 = roi2;

if isempty(fg.fibers{1})
    dict.nfib = 0;
else
    dict.nfib = length(fg.fibers);
end




