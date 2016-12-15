%% group averaged network stats for SfN 2016 poster
%

% identify HCP3T
subjLoad{1} = {'hcp', '105115'};
subjLoad{2} = {'hcp', '110411'};
subjLoad{3} = {'hcp', '111312'};
subjLoad{4} = {'hcp', '113619'};

% identify Stanford
subjLoad{5} = {'stn', 'FP'}; 
subjLoad{6} = {'stn', 'HT'};
subjLoad{7} = {'stn', 'KK'};
subjLoad{8} = {'stn', 'MP'};

% identify 7T data
subjLoad{9}  = {'7T', '108323'};
subjLoad{10} = {'7T', '109123'};
subjLoad{11} = {'7T', '131217'};
subjLoad{12} = {'7T', '910241'};

%% loop over, load, and average data

% for every subject listed above
for ii = 1:length(subjLoad)
    
    % load data
    tmp = feMergeRepeats(subjLoad{ii}{1}, subjLoad{ii}{2}, 'detr', '2');

    % capture averages of interest
    densMat(:,:,ii) = tmp{2}.emat.mean;
    lemdMat(:,:,ii) = tmp{14}.emat.mean;
    
    % clean out tmp
    clear tmp
    
end

clear ii

% get group averaged mean matrices
densAvg = mean(densMat, 3, 'omitnan');
lemdAvg = mean(lemdMat, 3, 'omitnan');

% calculate strength / degree
densStr = strengths_und(densAvg);
lemdStr = strengths_und(lemdAvg);
densDeg = degrees_und(densAvg);
lemdDeg = degrees_und(lemdAvg);

% sort and index edges
[ ~, densStrIndx ] = sort(densStr);
[ ~, densDegIndx ] = sort(densDeg);
roiLabel(densStrIndx)
roiLabel(densDegIndx)

[ ~, lemdStrIndx ] = sort(lemdStr);
[ ~, lemdDegIndx ] = sort(lemdDeg);
roiLabel(lemdStrIndx)
roiLabel(lemdDegIndx)

% pull the top 10 from each 
% repeat for deterministic


%% Probabilistic Lmax 10

% Dens - Strength
%     'rh.entorhinal'
%     'rh.parahippocampal'
%     'rh.medialorbitofrontal'
%     'rh.temporalpole'
%     'lh.entorhinal'
%     'rh.isthmuscingulate'
%     'rh.transversetemporal'
%     'lh.parahippocampal'
%     'lh.lateralorbitofrontal'
%     'lh.parsorbitalis'

% Dens - Degree
%     'rh.entorhinal'
%     'rh.parahippocampal'
%     'lh.bankssts'
%     'rh.transversetemporal'
%     'lh.entorhinal'
%     'rh.rostralanteriorcingulate'
%     'lh.parahippocampal'
%     'rh.parsorbitalis'
%     'rh.paracentral'
%     'rh.bankssts'

% EMD - Strength
%     'rh.entorhinal'
%     'rh.parahippocampal'
%     'rh.temporalpole'
%     'rh.transversetemporal'
%     'rh.frontalpole'
%     'rh.rostralanteriorcingulate'
%     'lh.transversetemporal'
%     'rh.bankssts'
%     'lh.entorhinal'
%     'rh.medialorbitofrontal'

% EMD - Degree
%     'rh.entorhinal'
%     'rh.parahippocampal'
%     'lh.parahippocampal'
%     'lh.entorhinal'
%     'rh.rostralanteriorcingulate'
%     'rh.bankssts'
%     'rh.transversetemporal'
%     'rh.medialorbitofrontal'
%     'rh.temporalpole'
%     'lh.bankssts'   

%% Deterministic Lmax 2

% Dens - Strength
%     'rh.entorhinal'
%     'rh.parahippocampal'
%     'lh.entorhinal'
%     'lh.parahippocampal'
%     'rh.transversetemporal'
%     'rh.temporalpole'
%     'rh.medialorbitofrontal'
%     'rh.lateralorbitofrontal'
%     'lh.temporalpole'
%     'lh.lateralorbitofrontal'

% Dens - Degree
%     'rh.transversetemporal'
%     'lh.entorhinal'
%     'lh.parahippocampal'
%     'rh.paracentral'
%     'rh.posteriorcingulate'
%     'rh.bankssts'
%     'lh.supramarginal'
%     'rh.entorhinal'
%     'lh.caudalmiddlefrontal'
%     'rh.caudalmiddlefrontal'

% EMD - Strength
%     'rh.entorhinal'
%     'rh.temporalpole'
%     'lh.entorhinal'
%     'rh.transversetemporal'
%     'rh.rostralanteriorcingulate'
%     'rh.parahippocampal'
%     'lh.transversetemporal'
%     'lh.frontalpole'
%     'lh.parsorbitalis'
%     'lh.temporalpole'

% EMD - Degree
%     'rh.entorhinal'
%     'rh.transversetemporal'
%     'lh.parahippocampal'
%     'lh.bankssts'
%     'lh.caudalmiddlefrontal'
%     'rh.paracentral'
%     'rh.supramarginal'
%     'lh.paracentral'
%     'lh.entorhinal'
%     'rh.caudalanteriorcingulate'





    