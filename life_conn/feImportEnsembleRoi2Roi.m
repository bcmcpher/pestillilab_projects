function [ fgOut, dict ] = feImportEnsembleRoi2Roi(ens_tck, fout)
%%
% Brent McPherson
% 20151231
% import series of .tck files and create a key of # of fibers for each file
% skips files that are empty
%
% Inputs:
%   ens_tck: string with path to ensemble track outputs
%
% Outputs:
%   fgOut: merged fiber group
%   dict: meta information for determining # of fibers during virtual
%   lesions
%

%% read out necessary unique values

display('indexing .tck files...');

% read in all the files
files = dir(ens_tck);

% remove dot folders
files = files(3:length(files));

% create fout if it doesn't exist
% copied from dtiImportFibersMrtrix.m
if ~exist(fout, 'file')
    
    % Strip out the file name.
    [~, f] = fileparts(fout);
    
    % Build an empty mrDiffusion fier group.
    fgOut = dtiNewFiberGroup(f);

end

display('Beginning Loop of .tck files...');

% loop over every file
for ii = 1:length(files)
    
    % if there are fibers in the file
    if feReadFiberCount(fullfile(ens_tck, files(ii).name)) > 0
        
        % parse info out of file name
        tmp = files(ii).name;
        
        % read in, add to dictionary
        fg = dtiImportFibersMrtrix(fullfile(ens_tck, tmp));
        
        % merge new fiber group w/ output fiber group
        fgOut = dtiMergeFiberGroups(fgOut, fg, fout);
        
        % write down fgOut
        fgWrite(fgOut, fout, 'mat');        
                              
        % find indices of underscores
        und = strfind(tmp, '_');
        
        % find recon algorithm
        algo = tmp(und(2)+1:und(3)-1);
        
        % find curvature
        curv = tmp(und(3)+1:und(4)-1);
        
        % split roi names
        % isolate roi descriptions and trim file exention
        rois = tmp(und(4)+1:length(tmp)-4);
        
        % find middle index and separate rois
        splt = strfind(rois, '_to_');
        roi1 = rois(1:splt-1);
        roi2 = rois(splt+4:length(rois));
        
        % create dictionary entry
        dict(ii).file = tmp;
        dict(ii).algo = algo;
        dict(ii).curv = curv;
        dict(ii).roi1 = roi1;
        dict(ii).roi2 = roi2;
        dict(ii).nfib = length(fg.fibers);
       
    else
        
        % skip to next iteration
        continue
        
    end
    
end

% loop index for fiber counts
findx = 1;

% create indice pairs for fibers for simplified subsets
for jj = 1:length(dict)
    
    dict(jj).ifib = [ findx, findx + dict(jj).nfib ];
    findx = findx + dict(jj).nfib;

end

end
% end of function

%% internal functions
function count = feReadFiberCount(filename)
% copied from dtiImportFibersMrtrix.m
% just enough to get fiber counts

% Strip out the file name.
[~,f] = fileparts(filename);

% Build an empty mrDiffusion fier group.
fg = dtiNewFiberGroup(f);

% Read a binary fiber tracking file (.tck) output from mrTrix. 
fid = fopen(filename ,'r','ieee-le'); % Note that we assume that the data 
                                      % always little-endian. 
if(fid==-1), error('Unable to access file %s\n', filename);end

% Read the .tck file just opened.
try
    % Read the text header, line-by-line, until the 'END' keyword. We'll
    % store all header fields in a cell array and then pull out the ones
    % that we need below.
    ln = fgetl(fid);
    ii = 1;
    while(~strcmp(ln,'END'))
        header{ii} = ln;
        ln = fgetl(fid);
        ii = ii+1;
    end
   
    % Get the datatype from the header cell array.
    dt = header{strmatch('datatype:',header)};
    if(isempty(findstr(dt,'Float32LE')))
        % *** FIXME: we should close the file and reopen in big-endian.
        error('Only Float32LE data supported!');
    end
    
    % Get the number of tracts from the header cell array. There seem to
    % be two possible keywords for this field.
    numIndx = strmatch('num_tracks:',header);
    if(isempty(numIndx))
        numIndx = strmatch('count:',header); 
        count = str2double(header{numIndx}(7:end));
    else
        count = str2double(header{numIndx}(12:end));
    end
    
    % fprintf(1,'Reading Fiber Data for %d fibers...\n', count);
    end
end