function [ dict ] = makeDict(ens_tck)
% makeDict Summary of this function goes here
%   Detailed explanation goes here

display('indexing .tck files...');

% read in all the files
files = dir(ens_tck);

% remove dot folders
files = files(3:length(files));

% testing short length
%files = files(1:300);

% more generalized dictionary structure
% dict(1).roi1 = 'ROI1';
% dict(1).roi2 = 'ROI2';
% dict(1).file = {'...tck', '...tck', ...};
% dict(1).algo = {'PROB', 'PROB', ...};
% dict(1).nfib = {[1 2], [3 7], ...};
% etc...

% create an empty structure object
%entry = struct('roi1', char, 'roi2', char, 'file', {}, 'algo', {}, 'curve', {}, 'nfib', {});
%dict = repmat(entry, 10);

% for all ROI pairs; include empty files
% define dictionary structure / size before entering loop
% don't write fiber group to disk until the end
% save -append when writing fiber

iter = 1;

display('Beginning Loop of .tck files...');

% for every ROI pair
for ii = 1:length(1:10:length(files))
    
    %display(['ROI pairs loop: ', num2str(ii)]);
    
    % for every combination of tracking parameters
    for jj = 1:10
        
        %display(['ROI file: ', num2str(jj)]);
         
        %display(['iter = ', num2str(iter)]);
        
        % parse info out of file name
        tmp = files(iter).name;
        
        % find indices of underscores
        und = strfind(tmp, '_');
        
        % split roi names
        rois = tmp(und(4)+1:length(tmp)-4);
        
        % find recon algorithm
        algo = tmp(und(2)+1:und(3)-1);
        
        % find curvature
        curv = str2double(tmp(und(3)+1:und(4)-1));
        
        % isolate roi descriptions and trim file exention
        % find middle index and separate rois
        splt = strfind(rois, '_to_');
        roi1 = rois(1:splt-1);
        roi2 = rois(splt+4:length(rois));
        
        % if there are fibers in the file
        count = feReadFiberCount(fullfile(ens_tck, tmp));

        % read in, add to dictionary
        % fg = dtiImportFibersMrtrix(fullfile(ens_tck, tmp));
        
        % create dictionary entry
        % assuming it ever gets intialized correctly...
        dict(ii).roi1 = roi1;
        dict(ii).roi2 = roi2;
        dict(ii).file{jj} = tmp;
        dict(ii).algo{jj} = algo;
        dict(ii).curv{jj} = curv;
        dict(ii).nfib{jj} = count;
        
        iter = iter + 1;
        
    end 
end

% loop index for fiber counts
%findx = 1;

% create indice pairs for fibers for simplified subsets

prev = 1;

% for every dictionary entry
for roi = 1:length(dict)
    
    % for every file in roi
    for file = 1:10
        
        % if there are no fibers in the file
        if dict(roi).nfib{file} == 0
            
            % mark as empty
            dict(roi).ifib{file} = -1;
        else
            
            % create the indices fibers for the file based on the indices
            % of the previous entry
            dict(roi).ifib{file} = [ prev (prev + dict(roi).nfib{file} - 1) ];
            
            % update previous indices
            prev = prev + dict(roi).nfib{file};
            
        end
        
    end
    
end

end

%% internal fxn
function count = feReadFiberCount(filename)
% copied from dtiImportFibersMrtrix.m
% just enough to get fiber counts

% Strip out the file name.
%[ ~, f ] = fileparts(filename);

% Build an empty mrDiffusion fier group.
%fg = dtiNewFiberGroup(f);

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
    
    fclose(fid);
    % fprintf(1,'Reading Fiber Data for %d fibers...\n', count);
    end
end

