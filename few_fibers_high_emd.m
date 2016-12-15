%% check how often the strongest connections are a really small number of fibers

% % index into pconn is always the same
% indx = nchoosek(1:68, 2);
% tmp = 1:length(indx);
% tmp = tmp';
% indx = [indx tmp];
% clear tmp

% the connection to compare across repeats
indx = 2199;

for ii = 1:10
    
    num = sprintf('%02d', ii);
    
    % load data
    load(['data/hcp_105115_prob_lmax10_rep' num '.mat']);
    
    % pull emd matrix
    remd(:,:,ii) = emat(:,:,14);
    %emdl = tril(emd(:,:,ii));
    
    % pull fiber indices
    fibs{ii} = pconn{indx}.end.fibers;
    nzfb{ii} = pconn{indx}.end.nzfibs;
    
    % sort and find x / y indices in matrix for strongest connections
    %[ vals, mInd ] = sort(emdl(:), 1, 'descend');
    %[x, y] = ind2sub(size(emdl), mInd);
    
    % pull the highest index
    %indx(indx(:,3) == vals(1), 1:2)
    
    % pull the same index
    %matcoord = indx(indx(:,3) == 2171, 1:2);

    clear num emat fdat fns imat nmat out pconn pval roiNames time
end

% collapse output
emd_val_fibs = [reshape(remd(66, 53, :), [10, 1]), cellfun('size', fibs, 1)', cellfun('size', nzfb, 1)']

% emd_val_fibs =
% 
%    28.6698    8.0000    1.0000
%    45.4964    4.0000         0
%    50.2971    7.0000         0
%    46.3656   10.0000         0
%    13.2140    7.0000         0
%    15.3618    2.0000         0
%    18.9040    3.0000         0
%    21.2896    5.0000         0
%    30.8231    6.0000         0
%    42.1248    6.0000         0
