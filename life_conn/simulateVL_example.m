
clear all
load wrk_vl_test_20160202.mat

num_tot_fascicles = numel(feGet(fe,'fiber weights'));
num_tract_fas = [ 10, 25, 50, 100, 150 ];
nTests = 100;

% for itest  = 1:nTests
%     test(itest).fascicles       = sort(randsample(1:num_tot_fascicles,num_tract_fas));
%     [se, feLesion,  feNoLesion] = feVirtualLesion(fe,test(itest).fascicles);
%     test(itest).weightsmatch    = all(ismember(feLesion.life.fit.weights,feNoLesion.life.fit.weights));
%     test(itest).s = se.s.mean;
% end
% horzcat(test(:).s)
% printf('I found %i negative S out of %i tests.',sum(S < 0),nTests)

for itest = 1:nTests
    for jtest = 1:length(num_tract_fas)
        test(itest).fascicles{jtest} = sort(randsample(1:num_tot_fascicles,num_tract_fas(jtest)));
        [ se, feLesion, feNoLesion, pn ] = feVirtualLesion(fe, test(itest).fascicles{jtest});
        test(itest).wghtsmtch{jtest} = all(ismember(feLesion.life.fit.weights, feNoLesion.life.fit.weights));
        test(itest).pn{jtest} = pn;
        test(itest).s{jtest} = se.s.mean;
    end
end

% combine means into rows of sample by nFibers
S = vertcat(test(:).s);

%sprintf('I found %i negative S out of %i tests.',sum(S < 0),nTests)

%% Test large / small weigted fibers
% ia = highest to lowest fiber weights
% ib = index of sorted fibers in tensor
[ia, ib] = sort(fe.life.fit.weights, 'descend');

% identify highest weights
hghWght = ib(1:150); % highest weighted 15 fibers

% pull a bunch of 0 weights
lowWght = ib(40000:40150); % all 0 weighted fibers

% compare high to low
for ii = 1:500
    htmp = feVirtualLesion(fe, hghWght);
    seHgh(ii) = htmp.s.mean;
    ltmp = feVirtualLesion(fe, lowWght);
    seLow(ii) = ltmp.s.mean;
end

% quick descriptives
minmax(seHgh)
minmax(seLow)

% plots
hist(seHgh, 40);
hist(seLow, 40);

% not consistent between re-runs of the same fibers...

% Same parameters, different answers
% seHgh.s.mean
% ans = 0.0080
% 
% seLow.s.mean
% ans = 0.0064
%    
% seHgh.s.mean
% ans = 0.0076
% 
% seLow.s.mean
% ans = -0.0069
   



%% tweak histogram iterations in feComputeEvidence

% later




%% Simulate the weight subtract and virtual lesion

numFas    = 100;
numTractFas = 10;

weights   = rand(numFas,1);
fas_preds = randn(numFas,1);
PredNoLes = fas_preds.*weights;
MesSig    = PredNoLes+ (randn(numFas,1).*0.01);
orig_rmse = mean(sqrt( (PredNoLes - MesSig).^2 ));

numOfMontacarloSimulations = 200;
for imc = 1:numOfMontacarloSimulations
    w_i = randsample(1:numFas,numTractFas);
    weights(w_i) = 0;
    PredLes = fas_preds.*weights;
    les_rmse(imc) = mean(sqrt( (PredLes - MesSig).^2 ));
end

% is any instance negative?
any( (les_rmse - orig_rmse) < 0)
