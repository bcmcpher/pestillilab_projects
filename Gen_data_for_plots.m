function [ all_data ] = Gen_data_for_plots()

%% Edit this section to load different parameters for all subjects

% build the parameter set to load all subjects with

% essentially static - fewer lmax's would be fine, alg_set will always do
% tensor after the fact, so don't add - also makes removing one dumb
lparam_set = {'2', '4', '6', '8', '10', '12'}; % lmax
alg_set = {'detr', 'prob'};

% these can be changed to whatever as long as labels match
glbVars = {'smwrld', 'mcoef', 'trans', 'glEff', 'stat', 'mod'};
nodVars = {'str', 'deg', 'btw', 'eigv', 'pcoef', 'modz'};
connMod = {2, 6, 14};

% labels for output, condense to single object for shorter param list
glbLabel = {'Small Worldness', 'Mean Clustering Coefficient', 'Transitivity', ... 
           'Global Efficiency', 'Community Structure', 'Modularity'};
nodLabel = {'Strength', 'Degree', 'Betweenness', 'Eigenvector Centrality', ...
            'Partial Coefficiency', 'Modularity'};
modLabel = {'Fiber Density', 'NonZero Density', 'EMD'};

%% Read dataset HCP3T_90dir

sbj_set = {'105115', '110411', '111312', '113619'};

disp('Gen results HCP3T 90dir')
[data] = gen_results(sbj_set, glbVars, nodVars, lparam_set, alg_set, connMod);

all_data = data;

%% Read dataset STN_96dir

sbj_set = {'FP', 'HT', 'KK', 'MP'};

disp('Gen results STN 96dir')
[data] = gen_results(sbj_set, glbVars, nodVars, lparam_set, alg_set, connMod);

% all_data = [ all_data data ];
% 
% %% Read dataset HCP7T_60dir
% 
% sbj_set = {'108323','109123','131217','910241'};
% 
% disp('Gen results HCP7T 60dir')
% [data] = gen_resutls(sbj_set, glbVars, nodVars, lparam_set, alg_set, connMod);
% 
all_data = [ all_data data ];

end

function [data] = gen_results(sbj_set, glbVars, nodVars, lparam_set, alg_set, connMod)

% hard coded output labels
glbLabel = {'Small Worldness', 'Mean Clustering Coefficient', 'Transitivity', ... 
           'Global Efficiency', 'Community Structure', 'Modularity'};
nodLabel = {'Strength', 'Degree', 'Betweenness', 'Eigenvector Centrality', ...
            'Partial Coefficiency', 'Modularity'};
modLabel = {'Fiber Density', 'NonZero Density', 'EMD'};

% create empty data cell-array
data = {};

% do the thing with the stuff
mStem = 'avgSub{mInd}.fns.nrm.mean.';
sStem = 'avgSub{mInd}.fns.nrm.std.';
rStem = 'rawSub{mInd}.nrm.';

% for every subject
for isbj =1:length(sbj_set)
    subject = sbj_set{isbj};
    disp(['Retrieving Connectome Statistics for: ', subject, '...']);
    % for detr / prob
    for ialg = 1:length(alg_set)
        alg = alg_set{ialg};
        % for lmax
        for iparam = 1:length(lparam_set)
            lparam = lparam_set{iparam};
            % find and load subject repeats and create data observations
            %disp(['Retrieving Connectome Statistics for: ', subject, ' with ',alg,' tracking and Lmax ',lparam, '...']);
            switch subject
                % HCP 3T
                case {'105115','110411','111312','113619'} 
                    % load repeat data
                    [ avgSub, ~, rawSub ] = feMergeRepeats('hcp', subject, alg, lparam);
                    % for every connectome I want stats for
                    for imod = 1:length(connMod)
                        mInd = connMod{imod};
                        % for every stat I want
                        for iglb = 1:length(glbVars)
                            % pull stats by row, into a cell array of connectomes
                            dat{imod}(iglb,1) = eval([mStem glbVars{iglb}]);
                            % pull calculations of error
                            t1 = eval([sStem glbVars{iglb}]);
                            t2 = size(eval([rStem glbVars{iglb}]), 1);
                            t3 = t1 / t2;
                            % write lower / upper bounds of estimate error
                            % as 2nd / 3rd column of data
                            dat{imod}(iglb, 2) = dat{imod}(iglb, 1) - t3;
                            dat{imod}(iglb, 3) = dat{imod}(iglb, 1) + t3;
                            clear t1 t2 t3
                        end
                        % create a new temporary data set for nodes
                        for inod = 1:length(nodVars)
                            % pull stats by row, into a cell array of connectomes
                            ndat{inod, imod}(:, 1) = eval([mStem nodVars{inod}]);
                            % pull calculations of error
                            t1 = eval([sStem nodVars{inod}]);
                            t2 = size(eval([rStem nodVars{inod}]), 1);
                            t3 = t1 / t2;
                            % write lower / upper bounds of estimate error
                            % as 2nd / 3rd column of data
                            ndat{inod, imod}(:, 2) = ndat{inod, imod}(:, 1) - t3;
                            ndat{inod, imod}(:, 3) = ndat{inod, imod}(:, 1) + t3;
                            clear t1 t2 t3
                        end
                    end
                % STN                        
                case {'FP','HT','KK','MP'} 
                    % load repeat data
                    [ avgSub, ~, rawSub ] = feMergeRepeats('stn', subject, alg, lparam);
                    for imod = 1:length(connMod)
                        mInd = connMod{imod};
                        % for every stat I want
                        for iglb = 1:length(glbVars)
                            % pull stats by row, into a cell array of connectomes
                            dat{imod}(iglb,1) = eval([mStem glbVars{iglb}]);
                            % pull calculations of error
                            t1 = eval([sStem glbVars{iglb}]);
                            t2 = size(eval([rStem glbVars{iglb}]), 1);
                            t3 = t1 / t2;
                            % write lower / upper bounds of estimate error
                            % as 2nd / 3rd column of data
                            dat{imod}(iglb, 2) = dat{imod}(iglb, 1) - t3;
                            dat{imod}(iglb, 3) = dat{imod}(iglb, 1) + t3;
                            clear t1 t2 t3
                        end
                        % create a new temporary data set for nodes
                        for inod = 1:length(nodVars)
                            % pull stats by row, into a cell array of connectomes
                            ndat{inod, imod}(:, 1) = eval([mStem nodVars{inod}]);
                            % pull calculations of error
                            t1 = eval([sStem nodVars{inod}]);
                            t2 = size(eval([rStem nodVars{inod}]), 1);
                            t3 = t1 / t2;
                            % write lower / upper bounds of estimate error
                            % as 2nd / 3rd column of data
                            ndat{inod, imod}(:, 2) = ndat{inod, imod}(:, 1) - t3;
                            ndat{inod, imod}(:, 3) = ndat{inod, imod}(:, 1) + t3;
                            clear t1 t2 t3
                        end
                    end
                % HCP 7T    
                case {'108323','109123','131217','910241'} 
                    % load repeat data
                    [ avgSub, ~, rawSub ] = feMergeRepeats('hc7', subject, alg, lparam);
                    for imod = 1:length(connMod)
                        mInd = connMod{imod};
                        % for every stat I want
                        for iglb = 1:length(glbVars)
                            % pull stats by row, into a cell array of connectomes
                            dat{imod}(iglb,1) = eval([mStem glbVars{iglb}]);
                            % pull calculations of error
                            t1 = eval([sStem glbVars{iglb}]);
                            t2 = size(eval([rStem glbVars{iglb}]), 1);
                            t3 = t1 / t2;
                            % write lower / upper bounds of estimate error
                            % as 2nd / 3rd column of data
                            dat{imod}(iglb, 2) = dat{imod}(iglb, 1) - t3;
                            dat{imod}(iglb, 3) = dat{imod}(iglb, 1) + t3;
                            clear t1 t2 t3
                        end
                        % create a new temporary data set for nodes
                        for inod = 1:length(nodVars)
                            % pull stats by row, into a cell array of connectomes
                            ndat{inod, imod}(:, 1) = eval([mStem nodVars{inod}]);
                            % pull calculations of error
                            t1 = eval([sStem nodVars{inod}]);
                            t2 = size(eval([rStem nodVars{inod}]), 1);
                            t3 = t1 / t2;
                            % write lower / upper bounds of estimate error
                            % as 2nd / 3rd column of data
                            ndat{inod, imod}(:, 2) = ndat{inod, imod}(:, 1) - t3;
                            ndat{inod, imod}(:, 3) = ndat{inod, imod}(:, 1) + t3;
                            clear t1 t2 t3
                        end
                    end
            end
            % add specific data to loop
            dlab = [alg '_' lparam];
            data{isbj}.subj = subject;
            data{isbj}.mods_cols = modLabel;
            data{isbj}.glob_rows = glbLabel';
            data{isbj}.node_rows = nodLabel';
            eval(['data{isbj}.' dlab '.glob = dat;'])
            eval(['data{isbj}.' dlab '.node = ndat;'])
            %eval(['data{isbj}.' dlab ' = dat;'])
            %data{isbj}.ndat = ndat;
            clear dat ndat;
        end
    end
    
    % at the end of every subject,
    switch subject
        case {'105115','110411','111312','113619'}
            [ avgSub, ~, rawSub ] = feMergeRepeats('hcp', subject, 'tens', '2');
            for imod = 1:length(connMod)
                mInd = connMod{imod};
                % for every stat I want
                for iglb = 1:length(glbVars)
                    % pull stats by row, into a cell array of connectomes
                    dat{imod}(iglb,1) = eval([mStem glbVars{iglb}]);
                    % pull calculations of error
                    t1 = eval([sStem glbVars{iglb}]);
                    t2 = size(eval([rStem glbVars{iglb}]), 1);
                    t3 = t1 / t2;
                    % write lower / upper bounds of estimate error
                    % as 2nd / 3rd column of data
                    dat{imod}(iglb, 2) = dat{imod}(iglb, 1) - t3;
                    dat{imod}(iglb, 3) = dat{imod}(iglb, 1) + t3;
                    clear t1 t2 t3
                end
                % create a new temporary data set for nodes
                for inod = 1:length(nodVars)
                    % pull stats by row, into a cell array of connectomes
                    ndat{inod, imod}(:, 1) = eval([mStem nodVars{inod}]);
                    % pull calculations of error
                    t1 = eval([sStem nodVars{inod}]);
                    t2 = size(eval([rStem nodVars{inod}]), 1);
                    t3 = t1 / t2;
                    % write lower / upper bounds of estimate error
                    % as 2nd / 3rd column of data
                    ndat{inod, imod}(:, 2) = ndat{inod, imod}(:, 1) - t3;
                    ndat{inod, imod}(:, 3) = ndat{inod, imod}(:, 1) + t3;
                    clear t1 t2 t3
                end
            end
            data{isbj}.tens.glob = dat;
            data{isbj}.tens.node = ndat;
            clear dat ndat;
%         case {'FP','HT','KK','MP'}
%             [ avgSub, ~, rawSub ] = feMergeRepeats('stn', subject, 'tens', '2');
%             for imod = 1:length(connMod)
%                 mInd = connMod{imod};
%                 for iglb = 1:length(glbVars)
%                     dat{imod}(iglb,1) = eval([mStem glbVars{iglb}]);
%                     t1 = eval([sStem glbVars{iglb}]);
%                     t2 = size(eval([rStem glbVars{iglb}]), 1);
%                     t3 = t1 / t2;
%                     dat{imod}(iglb, 2) = dat{imod}(iglb, 1) - t3;
%                     dat{imod}(iglb, 3) = dat{imod}(iglb, 1) + t3;
%                     clear t1 t2 t3
%                 end
%             end
%                 % create a new temporary data set for nodes
%                 for inod = 1:length(nodVars)
%                     % pull stats by row, into a cell array of connectomes
%                     ndat{inod, imod}(:, 1) = eval([mStem nodVars{inod}]);
%                     % pull calculations of error
%                     t1 = eval([sStem nodVars{inod}]);
%                     t2 = size(eval([rStem nodVars{inod}]), 1);
%                     t3 = t1 / t2;
%                     % write lower / upper bounds of estimate error
%                     % as 2nd / 3rd column of data
%                     ndat{inod, imod}(:, 2) = ndat{inod, imod}(:, 1) - t3;
%                     ndat{inod, imod}(:, 3) = ndat{inod, imod}(:, 1) + t3;
%                     clear t1 t2 t3
%                 end
%             end
%             data{isbj}.tens.glob = dat;
%             data{isbj}.tens.node = ndat;
%             clear dat ndat;
%         case {'108323','109123','131217','910241'}
%             [ avgSub, ~, rawSub ] = feMergeRepeats('hc7', subject, 'tens', '2');
%             for imod = 1:length(connMod)
%                 mInd = connMod{imod};
%                 for iglb = 1:length(glbVars)
%                     dat{imod}(iglb,1) = eval([mStem glbVars{iglb}]);
%                     t1 = eval([sStem glbVars{iglb}]);
%                     t2 = size(eval([rStem glbVars{iglb}]), 1);
%                     t3 = t1 / t2;
%                     dat{imod}(iglb, 2) = dat{imod}(iglb, 1) - t3;
%                     dat{imod}(iglb, 3) = dat{imod}(iglb, 1) + t3;
%                     clear t1 t2 t3
%                 end
%             end
%                 % create a new temporary data set for nodes
%                 for inod = 1:length(nodVars)
%                     % pull stats by row, into a cell array of connectomes
%                     ndat{inod, imod}(:, 1) = eval([mStem nodVars{inod}]);
%                     % pull calculations of error
%                     t1 = eval([sStem nodVars{inod}]);
%                     t2 = size(eval([rStem nodVars{inod}]), 1);
%                     t3 = t1 / t2;
%                     % write lower / upper bounds of estimate error
%                     % as 2nd / 3rd column of data
%                     ndat{inod, imod}(:, 2) = ndat{inod, imod}(:, 1) - t3;
%                     ndat{inod, imod}(:, 3) = ndat{inod, imod}(:, 1) + t3;
%                     clear t1 t2 t3
%                 end
%             end
%             data{isbj}.tens.glob = dat;
%             data{isbj}.tens.node = ndat;
%             clear dat ndat;
    end
    
end

end
