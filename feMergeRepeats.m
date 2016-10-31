function [ avgOut, ematOut, fnsOut, pvalOut, fdatOut ] = feMergeRepeats(dgrp, subj, dmdl, lmax)
% merge repeated data by subject into a single parsable object
%   
% dgrp = 'stn'; subj = 'FP'; dmdl = 'prob'; lmax = '10';
%

%keyboard; 

%% indentify files based on subject

% hard-coded reference to data directory
topdir = '/N/dc2/projects/lifebid/HCP/Brent/cogs610/matlab/data';

% format numeric lmax interface to 
%lmax = sprintf('%02d', lmax);

% create stem for dir to read in list of .mat objects
fnam = [dgrp, '_', subj, '_', dmdl, '_lmax', lmax, '_rep*'];

% read in specified file names
files = dir(fullfile(topdir, fnam));
files = {files.name};

nreps = size(files, 2);

% load data
for ii = 1:nreps
    emat{ii} = load(fullfile(topdir, files{ii}), 'emat');
    fns{ii} = load(fullfile(topdir, files{ii}), 'fns');
    pval{ii} = load(fullfile(topdir, files{ii}), 'pval');
    fdat{ii} = load(fullfile(topdir, files{ii}), 'fdat');
end

%% make subject data object

% create empty output object
rmat = zeros(68, 68, nreps);
for ii = 1:16
    ematOut{ii} = rmat; 
end

% fill output object;
% for every matrix type
for ii = 1:16
    % for every repeat
    for jj = 1:nreps
        % fill in values
        ematOut{ii}(:,:,jj) = emat{jj}.emat(:,:,ii);
    end
end

% fh = figure;
% for kk = 1:9
%     subplot(3, 3, kk);
%     colormap('hot');
%     imagesc(out{16}(:,:,kk));
%     xlabel('FS DK Regions');
%     ylabel('FS DK Regions');
%     axis('square'); axis('equal'); axis('tight');
%     y = colorbar;
%     ylabel(y, 'Strength of Connection');
% end

%% make subject statistics object

% make empty template for nrm
nrm.cnsns = [];
nrm.deg = [];
nrm.str = [];
nrm.btw = [];
nrm.lcEff = [];
nrm.ccoef = [];
nrm.bcv = [];
nrm.eigv = [];
nrm.asgn = [];
nrm.stat = [];
nrm.struc = [];
nrm.mod = [];
nrm.pcoef = [];
nrm.modz = [];
nrm.mcoef = [];
nrm.dens = [];
nrm.dvrt = [];
nrm.dedg = [];
nrm.trans = [];
nrm.glEff = [];
nrm.chpl = [];
nrm.rcc = {};
nrm.smwrld = [];

% make empty template for bin
bin.kn = [];
bin.plord = {};
bin.pllvl = {};
bin.core = [];
bin.ckn = [];

% create empty output
for ii = 1:16
    fnsOut{ii}.nrm = nrm;
    fnsOut{ii}.bin = bin;
end

% merge repeats by matrix
for ii = 1:16
    for jj = 1:nreps
        
        % merge weighted stats
        fnsOut{ii}.nrm.cnsns = cat(2, fnsOut{ii}.nrm.cnsns, fns{jj}.fns{ii}.nrm.cnsns);
        fnsOut{ii}.nrm.deg = cat(2, fnsOut{ii}.nrm.deg, fns{jj}.fns{ii}.nrm.deg);
        fnsOut{ii}.nrm.str = cat(2, fnsOut{ii}.nrm.str, fns{jj}.fns{ii}.nrm.str);
        fnsOut{ii}.nrm.btw = cat(2, fnsOut{ii}.nrm.btw, fns{jj}.fns{ii}.nrm.btw);
        fnsOut{ii}.nrm.lcEff = cat(2, fnsOut{ii}.nrm.lcEff, fns{jj}.fns{ii}.nrm.lcEff);
        fnsOut{ii}.nrm.ccoef = cat(2, fnsOut{ii}.nrm.ccoef, fns{jj}.fns{ii}.nrm.ccoef);
        fnsOut{ii}.nrm.bcv = cat(2, fnsOut{ii}.nrm.bcv, fns{jj}.fns{ii}.nrm.bcv);
        fnsOut{ii}.nrm.eigv = cat(2, fnsOut{ii}.nrm.eigv, fns{jj}.fns{ii}.nrm.eigv);
        fnsOut{ii}.nrm.asgn = cat(2, fnsOut{ii}.nrm.asgn, fns{jj}.fns{ii}.nrm.asgn);
        fnsOut{ii}.nrm.stat = cat(1, fnsOut{ii}.nrm.stat, fns{jj}.fns{ii}.nrm.stat);
        fnsOut{ii}.nrm.struc = cat(2, fnsOut{ii}.nrm.struc, fns{jj}.fns{ii}.nrm.struc);
        fnsOut{ii}.nrm.mod = cat(1, fnsOut{ii}.nrm.mod, fns{jj}.fns{ii}.nrm.mod);
        fnsOut{ii}.nrm.pcoef = cat(2, fnsOut{ii}.nrm.pcoef, fns{jj}.fns{ii}.nrm.pcoef);
        fnsOut{ii}.nrm.modz = cat(2, fnsOut{ii}.nrm.modz, fns{jj}.fns{ii}.nrm.modz);
        fnsOut{ii}.nrm.mcoef = cat(1, fnsOut{ii}.nrm.mcoef, fns{jj}.fns{ii}.nrm.mcoef);
        fnsOut{ii}.nrm.dens = cat(1, fnsOut{ii}.nrm.dens, fns{jj}.fns{ii}.nrm.dens);
        fnsOut{ii}.nrm.dvrt = cat(1, fnsOut{ii}.nrm.dvrt, fns{jj}.fns{ii}.nrm.dvrt);
        fnsOut{ii}.nrm.dedg = cat(1, fnsOut{ii}.nrm.dedg, fns{jj}.fns{ii}.nrm.dedg);
        fnsOut{ii}.nrm.trans = cat(1, fnsOut{ii}.nrm.trans, fns{jj}.fns{ii}.nrm.trans);
        fnsOut{ii}.nrm.glEff = cat(1, fnsOut{ii}.nrm.glEff, fns{jj}.fns{ii}.nrm.glEff);
        fnsOut{ii}.nrm.chpl = cat(1, fnsOut{ii}.nrm.chpl, fns{jj}.fns{ii}.nrm.chpl);
        fnsOut{ii}.nrm.rcc(jj,:) = {fns{jj}.fns{ii}.nrm.rcc};
        fnsOut{ii}.nrm.smwrld = cat(1, fnsOut{ii}.nrm.smwrld, fns{jj}.fns{ii}.nrm.smwrld);
        
        % merge binarized stats
        fnsOut{ii}.bin.kn = cat(1, fnsOut{ii}.bin.kn, fns{jj}.fns{ii}.bin.kn);
        fnsOut{ii}.bin.plord(jj, :) = {fns{jj}.fns{ii}.bin.plord};
        fnsOut{ii}.bin.pllvl(jj, :) = {fns{jj}.fns{ii}.bin.pllvl};
        fnsOut{ii}.bin.core = cat(2, fnsOut{ii}.bin.core, fns{jj}.fns{ii}.bin.core);
        fnsOut{ii}.bin.ckn = cat(2, fnsOut{ii}.bin.ckn, fns{jj}.fns{ii}.bin.ckn);
        
    end
end

%% make subject p-value object

clear nrm

nrm.deg = [];
nrm.str = [];
nrm.btw = [];
nrm.lcEff = [];
nrm.ccoef = [];
nrm.bcv = [];
nrm.eigv = [];
nrm.asgn = [];
nrm.stat = [];
nrm.struc = [];
nrm.mod = [];
nrm.pcoef = [];
nrm.modz = [];
nrm.mcoef = [];
nrm.dens = [];
nrm.dvrt = [];
nrm.dedg = [];
nrm.trans = [];
nrm.glEff = [];
nrm.chpl = [];
nrm.rcc = {};

% create empty p-value object
for ii = 1:16
    pvalOut{ii} = nrm;
end

% merge repeats by matrix
for ii = 1:16
    for jj = 1:nreps
        
        % merge weighted stats
        pvalOut{ii}.deg = cat(2, pvalOut{ii}.deg, pval{jj}.pval{ii}.deg);
        pvalOut{ii}.str = cat(2, pvalOut{ii}.str, pval{jj}.pval{ii}.str);
        pvalOut{ii}.btw = cat(2, pvalOut{ii}.btw, pval{jj}.pval{ii}.btw);
        pvalOut{ii}.lcEff = cat(2, pvalOut{ii}.lcEff, pval{jj}.pval{ii}.lcEff);
        pvalOut{ii}.ccoef = cat(2, pvalOut{ii}.ccoef, pval{jj}.pval{ii}.ccoef);
        pvalOut{ii}.bcv = cat(2, pvalOut{ii}.bcv, pval{jj}.pval{ii}.bcv);
        pvalOut{ii}.eigv = cat(2, pvalOut{ii}.eigv, pval{jj}.pval{ii}.eigv);
        pvalOut{ii}.asgn = cat(2, pvalOut{ii}.asgn, pval{jj}.pval{ii}.asgn);
        pvalOut{ii}.stat = cat(1, pvalOut{ii}.stat, pval{jj}.pval{ii}.stat);
        pvalOut{ii}.struc = cat(2, pvalOut{ii}.struc, pval{jj}.pval{ii}.struc);
        pvalOut{ii}.mod = cat(1, pvalOut{ii}.mod, pval{jj}.pval{ii}.mod);
        pvalOut{ii}.pcoef = cat(2, pvalOut{ii}.pcoef, pval{jj}.pval{ii}.pcoef);
        pvalOut{ii}.modz = cat(2, pvalOut{ii}.modz, pval{jj}.pval{ii}.modz);
        pvalOut{ii}.mcoef = cat(1, pvalOut{ii}.mcoef, pval{jj}.pval{ii}.mcoef);
        pvalOut{ii}.dens = cat(1, pvalOut{ii}.dens, pval{jj}.pval{ii}.dens);
        pvalOut{ii}.dvrt = cat(1, pvalOut{ii}.dvrt, pval{jj}.pval{ii}.dvrt);
        pvalOut{ii}.dedg = cat(1, pvalOut{ii}.dedg, pval{jj}.pval{ii}.dedg);
        pvalOut{ii}.trans = cat(1, pvalOut{ii}.trans, pval{jj}.pval{ii}.trans);
        pvalOut{ii}.glEff = cat(1, pvalOut{ii}.glEff, pval{jj}.pval{ii}.glEff);
        pvalOut{ii}.chpl = cat(1, pvalOut{ii}.chpl, pval{jj}.pval{ii}.chpl);
        pvalOut{ii}.rcc(jj,:) = {pval{jj}.pval{ii}.rcc};
        
    end
end

%% make subject fdat object 

clear nrm

nrm.prp = [];
nrm.kcore = [];
nrm.raw = [];
nrm.thr = [];
nrm.nrm = [];
nrm.len = [];
nrm.bin = [];
nrm.dist = [];
nrm.edge = [];
nrm.agree = [];
nrm.nsrta = [];
nrm.btcm = [];

for ii = 1:16
    fdatOut{ii} = nrm;
end

% merge fdat objects
for ii = 1:16
    for jj = 1:nreps
        
        fdatOut{ii}.prp = cat(1, fdatOut{ii}.prp, fdat{jj}.fdat{ii}.prp);
        fdatOut{ii}.kcore = cat(3, fdatOut{ii}.kcore, fdat{jj}.fdat{ii}.kcore);
        fdatOut{ii}.raw = cat(3, fdatOut{ii}.raw, fdat{jj}.fdat{ii}.raw);
        fdatOut{ii}.thr = cat(3, fdatOut{ii}.thr, fdat{jj}.fdat{ii}.thr);
        fdatOut{ii}.nrm = cat(3, fdatOut{ii}.nrm, fdat{jj}.fdat{ii}.nrm);
        fdatOut{ii}.len = cat(3, fdatOut{ii}.len, fdat{jj}.fdat{ii}.len);
        fdatOut{ii}.bin = cat(3, fdatOut{ii}.bin, fdat{jj}.fdat{ii}.bin);
        fdatOut{ii}.dist = cat(3, fdatOut{ii}.dist, fdat{jj}.fdat{ii}.dist);
        fdatOut{ii}.edge = cat(3, fdatOut{ii}.edge, fdat{jj}.fdat{ii}.edge);
        fdatOut{ii}.agree = cat(3, fdatOut{ii}.agree, fdat{jj}.fdat{ii}.agree);
        fdatOut{ii}.nsrta = cat(3, fdatOut{ii}.nsrta, fdat{jj}.fdat{ii}.nsrta);
        fdatOut{ii}.btcm = cat(3, fdatOut{ii}.btcm, fdat{jj}.fdat{ii}.btcm);
        
    end
end

%% average objects

% mean / std of all the things
for ii = 1:16
    
    % data matrix
    avgOut{ii}.emat.mean = mean(ematOut{ii}, 3);
    
    % average weighted network statistics output 
    avgOut{ii}.fns.nrm.mean.cnsns = mean(fnsOut{ii}.nrm.cnsns, 2);
    avgOut{ii}.fns.nrm.mean.deg = mean(fnsOut{ii}.nrm.deg, 2);
    avgOut{ii}.fns.nrm.mean.str = mean(fnsOut{ii}.nrm.str, 2);
    avgOut{ii}.fns.nrm.mean.btw = mean(fnsOut{ii}.nrm.btw, 2);
    avgOut{ii}.fns.nrm.mean.lcEff = mean(fnsOut{ii}.nrm.lcEff, 2);
    avgOut{ii}.fns.nrm.mean.ccoef = mean(fnsOut{ii}.nrm.ccoef, 2);
    avgOut{ii}.fns.nrm.mean.bcv = mean(fnsOut{ii}.nrm.bcv, 2);
    avgOut{ii}.fns.nrm.mean.eigv = mean(fnsOut{ii}.nrm.eigv, 2);
    avgOut{ii}.fns.nrm.mean.asgn = mean(fnsOut{ii}.nrm.asgn, 2);
    avgOut{ii}.fns.nrm.mean.stat = mean(fnsOut{ii}.nrm.stat, 1);
    avgOut{ii}.fns.nrm.mean.struc = mean(fnsOut{ii}.nrm.struc, 2);
    avgOut{ii}.fns.nrm.mean.mod = mean(fnsOut{ii}.nrm.mod, 1);
    avgOut{ii}.fns.nrm.mean.pcoef = mean(fnsOut{ii}.nrm.pcoef, 2);
    avgOut{ii}.fns.nrm.mean.modz = mean(fnsOut{ii}.nrm.modz, 2);
    avgOut{ii}.fns.nrm.mean.mcoef = mean(fnsOut{ii}.nrm.mcoef, 1);
    avgOut{ii}.fns.nrm.mean.dens = mean(fnsOut{ii}.nrm.dens, 1);
    avgOut{ii}.fns.nrm.mean.dvrt = mean(fnsOut{ii}.nrm.dvrt, 1);
    avgOut{ii}.fns.nrm.mean.dedg = mean(fnsOut{ii}.nrm.dedg, 1);
    avgOut{ii}.fns.nrm.mean.trans = mean(fnsOut{ii}.nrm.trans, 1);
    avgOut{ii}.fns.nrm.mean.glEff = mean(fnsOut{ii}.nrm.glEff, 1);
    avgOut{ii}.fns.nrm.mean.chpl = mean(fnsOut{ii}.nrm.chpl, 1);
    avgOut{ii}.fns.nrm.mean.smwrld = mean(fnsOut{ii}.nrm.smwrld, 1);

    % average binary statistics output
    avgOut{ii}.fns.bin.mean.kn = mean(fnsOut{ii}.bin.kn, 1);
    avgOut{ii}.fns.bin.mean.core = mean(fnsOut{ii}.bin.core, 2);
    avgOut{ii}.fns.bin.mean.ckn = mean(fnsOut{ii}.bin.ckn, 2);
    
    % average p-values
    avgOut{ii}.pval.mean.deg = mean(pvalOut{ii}.deg, 2);
    avgOut{ii}.pval.mean.str = mean(pvalOut{ii}.str, 2);
    avgOut{ii}.pval.mean.btw = mean(pvalOut{ii}.btw, 2);
    avgOut{ii}.pval.mean.lcEff = mean(pvalOut{ii}.lcEff, 2);
    avgOut{ii}.pval.mean.ccoef = mean(pvalOut{ii}.ccoef, 2);
    avgOut{ii}.pval.mean.bcv = mean(pvalOut{ii}.bcv, 2);
    avgOut{ii}.pval.mean.eigv = mean(pvalOut{ii}.eigv, 2);
    avgOut{ii}.pval.mean.asgn = mean(pvalOut{ii}.asgn, 2);
    avgOut{ii}.pval.mean.stat = mean(pvalOut{ii}.stat, 1);
    avgOut{ii}.pval.mean.struc = mean(pvalOut{ii}.struc, 2);
    avgOut{ii}.pval.mean.mod = mean(pvalOut{ii}.mod, 1);
    avgOut{ii}.pval.mean.pcoef = mean(pvalOut{ii}.pcoef, 2);
    avgOut{ii}.pval.mean.modz = mean(pvalOut{ii}.modz, 2);
    avgOut{ii}.pval.mean.mcoef = mean(pvalOut{ii}.mcoef, 1);
    avgOut{ii}.pval.mean.dens = mean(pvalOut{ii}.dens, 1);
    avgOut{ii}.pval.mean.dvrt = mean(pvalOut{ii}.dvrt, 1);
    avgOut{ii}.pval.mean.dedg = mean(pvalOut{ii}.dedg, 1);
    avgOut{ii}.pval.mean.trans = mean(pvalOut{ii}.trans, 1);
    avgOut{ii}.pval.mean.glEff = mean(pvalOut{ii}.glEff, 1);
    avgOut{ii}.pval.mean.chpl = mean(pvalOut{ii}.chpl, 1);

    % average descriptive matrices
    avgOut{ii}.fdat.mean.prp = mean(fdatOut{ii}.prp, 1);
    avgOut{ii}.fdat.mean.kcore = mean(fdatOut{ii}.kcore, 3);
    avgOut{ii}.fdat.mean.raw = mean(fdatOut{ii}.raw, 3);
    avgOut{ii}.fdat.mean.thr = mean(fdatOut{ii}.thr, 3);
    avgOut{ii}.fdat.mean.nrm = mean(fdatOut{ii}.nrm, 3);
    avgOut{ii}.fdat.mean.len = mean(fdatOut{ii}.len, 3);
    avgOut{ii}.fdat.mean.bin = mean(fdatOut{ii}.bin, 3);
    avgOut{ii}.fdat.mean.dist = mean(fdatOut{ii}.dist, 3);
    avgOut{ii}.fdat.mean.edge = mean(fdatOut{ii}.edge, 3);
    avgOut{ii}.fdat.mean.agree = mean(fdatOut{ii}.agree, 3);
    avgOut{ii}.fdat.mean.nsrta = mean(fdatOut{ii}.nsrta, 3);
    avgOut{ii}.fdat.mean.btcm = mean(fdatOut{ii}.btcm, 3);
 
    % standard deviation of output
    
    % data matrix
    avgOut{ii}.emat.std = std(ematOut{ii}, 0, 3);
    
    % average weighted network statistics output 
    avgOut{ii}.fns.nrm.std.cnsns = std(fnsOut{ii}.nrm.cnsns, 0, 2);
    avgOut{ii}.fns.nrm.std.deg = std(fnsOut{ii}.nrm.deg, 0, 2);
    avgOut{ii}.fns.nrm.std.str = std(fnsOut{ii}.nrm.str, 0, 2);
    avgOut{ii}.fns.nrm.std.btw = std(fnsOut{ii}.nrm.btw, 0, 2);
    avgOut{ii}.fns.nrm.std.lcEff = std(fnsOut{ii}.nrm.lcEff, 0, 2);
    avgOut{ii}.fns.nrm.std.ccoef = std(fnsOut{ii}.nrm.ccoef, 0, 2);
    avgOut{ii}.fns.nrm.std.bcv = std(fnsOut{ii}.nrm.bcv, 0, 2);
    avgOut{ii}.fns.nrm.std.eigv = std(fnsOut{ii}.nrm.eigv, 0, 2);
    avgOut{ii}.fns.nrm.std.asgn = std(fnsOut{ii}.nrm.asgn, 0, 2);
    avgOut{ii}.fns.nrm.std.stat = std(fnsOut{ii}.nrm.stat, 0, 1);
    avgOut{ii}.fns.nrm.std.struc = std(fnsOut{ii}.nrm.struc, 0, 2);
    avgOut{ii}.fns.nrm.std.mod = std(fnsOut{ii}.nrm.mod, 0, 1);
    avgOut{ii}.fns.nrm.std.pcoef = std(fnsOut{ii}.nrm.pcoef, 0, 2);
    avgOut{ii}.fns.nrm.std.modz = std(fnsOut{ii}.nrm.modz, 0, 2);
    avgOut{ii}.fns.nrm.std.mcoef = std(fnsOut{ii}.nrm.mcoef, 0, 1);
    avgOut{ii}.fns.nrm.std.dens = std(fnsOut{ii}.nrm.dens, 0, 1);
    avgOut{ii}.fns.nrm.std.dvrt = std(fnsOut{ii}.nrm.dvrt, 0, 1);
    avgOut{ii}.fns.nrm.std.dedg = std(fnsOut{ii}.nrm.dedg, 0, 1);
    avgOut{ii}.fns.nrm.std.trans = std(fnsOut{ii}.nrm.trans, 0, 1);
    avgOut{ii}.fns.nrm.std.glEff = std(fnsOut{ii}.nrm.glEff, 0, 1);
    avgOut{ii}.fns.nrm.std.chpl = std(fnsOut{ii}.nrm.chpl, 0, 1);
    avgOut{ii}.fns.nrm.std.smwrld = std(fnsOut{ii}.nrm.smwrld, 0, 1);

    % average binary statistics output
    avgOut{ii}.fns.bin.std.kn = std(fnsOut{ii}.bin.kn, 0, 1);
    avgOut{ii}.fns.bin.std.core = std(fnsOut{ii}.bin.core, 0, 2);
    avgOut{ii}.fns.bin.std.ckn = std(fnsOut{ii}.bin.ckn, 0, 2);
    
    % average p-values
    avgOut{ii}.pval.std.deg = std(pvalOut{ii}.deg, 0, 2);
    avgOut{ii}.pval.std.str = std(pvalOut{ii}.str, 0, 2);
    avgOut{ii}.pval.std.btw = std(pvalOut{ii}.btw, 0, 2);
    avgOut{ii}.pval.std.lcEff = std(pvalOut{ii}.lcEff, 0, 2);
    avgOut{ii}.pval.std.ccoef = std(pvalOut{ii}.ccoef, 0, 2);
    avgOut{ii}.pval.std.bcv = std(pvalOut{ii}.bcv, 0, 2);
    avgOut{ii}.pval.std.eigv = std(pvalOut{ii}.eigv, 0, 2);
    avgOut{ii}.pval.std.asgn = std(pvalOut{ii}.asgn, 0, 2);
    avgOut{ii}.pval.std.stat = std(pvalOut{ii}.stat, 0, 1);
    avgOut{ii}.pval.std.struc = std(pvalOut{ii}.struc, 0, 2);
    avgOut{ii}.pval.std.mod = std(pvalOut{ii}.mod, 0, 1);
    avgOut{ii}.pval.std.pcoef = std(pvalOut{ii}.pcoef, 0, 2);
    avgOut{ii}.pval.std.modz = std(pvalOut{ii}.modz, 0, 2);
    avgOut{ii}.pval.std.mcoef = std(pvalOut{ii}.mcoef, 0, 1);
    avgOut{ii}.pval.std.dens = std(pvalOut{ii}.dens, 0, 1);
    avgOut{ii}.pval.std.dvrt = std(pvalOut{ii}.dvrt, 0, 1);
    avgOut{ii}.pval.std.dedg = std(pvalOut{ii}.dedg, 0, 1);
    avgOut{ii}.pval.std.trans = std(pvalOut{ii}.trans, 0, 1);
    avgOut{ii}.pval.std.glEff = std(pvalOut{ii}.glEff, 0, 1);
    avgOut{ii}.pval.std.chpl = std(pvalOut{ii}.chpl, 0, 1);

    % average descriptive matrices
    avgOut{ii}.fdat.std.prp = std(fdatOut{ii}.prp, 0, 1);
    avgOut{ii}.fdat.std.kcore = std(fdatOut{ii}.kcore, 0, 3);
    avgOut{ii}.fdat.std.raw = std(fdatOut{ii}.raw, 0, 3);
    avgOut{ii}.fdat.std.thr = std(fdatOut{ii}.thr, 0, 3);
    avgOut{ii}.fdat.std.nrm = std(fdatOut{ii}.nrm, 0, 3);
    avgOut{ii}.fdat.std.len = std(fdatOut{ii}.len, 0, 3);
    avgOut{ii}.fdat.std.bin = std(fdatOut{ii}.bin, 0, 3);
    avgOut{ii}.fdat.std.dist = std(fdatOut{ii}.dist, 0, 3);
    avgOut{ii}.fdat.std.edge = std(fdatOut{ii}.edge, 0, 3);
    avgOut{ii}.fdat.std.agree = std(fdatOut{ii}.agree, 0, 3);
    avgOut{ii}.fdat.std.nsrta = std(fdatOut{ii}.nsrta, 0, 3);
    avgOut{ii}.fdat.std.btcm = std(fdatOut{ii}.btcm, 0, 3);
    
end