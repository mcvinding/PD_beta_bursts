% PD-BB: Cluster-based permutation tests on PSD from running window and FFT
% annalysis
% set paths
clear all
close all

addpath /home/mikkel/PD_motor/global_scripts
[dirs, subs, ~] = PD_proj_setup('betaburst');

cd(dirs.megDir);

subjs = find_subs(dirs.megDir);     %Find subjects in folder
sub_ptns = intersect(subs.PD, subjs);
sub_ctrl = intersect(subs.ctrl, subjs);

%% Load data
disp('Loading...');
load('/home/mikkel/PD_motor/rest_ec/groupanalysis/PSD_data.mat')
disp('done');

%% Difference
ptns_diff = cell(length(sub_ptns),1);
ctrl_diff = cell(length(sub_ctrl),1);

cfg           = [];
cfg.parameter = 'powspctrm';
cfg.operation = '(x1-x2)';
for ii = 1:length(ptns_dat1)
    ptns_diff{ii} = ft_math(cfg, ptns_dat2{ii},ptns_dat1{ii});
end

for ii = 1:length(ctrl_dat1)
    ctrl_diff{ii} = ft_math(cfg,ctrl_dat2{ii},ctrl_dat1{ii});
end

%% Cluster stats
cfg = [];
cfg.method              = 'montecarlo';
cfg.correctm            = 'cluster';
cfg.clustertail         = 0;
cfg.clusteralpha        = 0.05;
cfg.clusterstatistic    = 'maxsum';
cfg.computeprob         = 'yes';
cfg.computecritval      = 'yes';
cfg.alpha               = .025;
cfg.numrandomization    = 1000;
cfg.tail                = 0;

% Between groups and interaction
cfg.statistic = 'ft_statfun_indepsamplesT'; 
cfg.design = [ones(length(ptns_dat1),1); ones(length(ctrl_dat1),1)*2];
cfg.ivar = 1;

result1 = ft_freqstatistics(cfg, ptns_dat1{:}, ctrl_dat1{:});
result2 = ft_freqstatistics(cfg, ptns_dat2{:}, ctrl_dat2{:});
resultX = ft_freqstatistics(cfg, ptns_diff{:}, ctrl_diff{:});

% Within groups
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.design = [ones(length(ptns_dat1),1); ones(length(ptns_dat2),1)*2];
cfg.design = [cfg.design, repmat(1:length(ptns_dat1),1,2)'];
cfg.ivar = 1;
cfg.uvar = 2;

result_ptns = ft_freqstatistics(cfg, ptns_dat1{:}, ptns_dat2{:});

cfg.design = [ones(length(ctrl_dat1),1); ones(length(ctrl_dat1),1)*2];
cfg.design = [cfg.design, repmat(1:length(ctrl_dat1),1,2)'];

result_ctrl = ft_freqstatistics(cfg, ctrl_dat1{:}, ctrl_dat2{:});

%% Save results
disp('Saving...');
save('/home/mikkel/PD_motor/rest_ec/groupanalysis/PSD_stats.mat', ...
    'result1','result2','result_ctrl','result_ptns','resultX')
disp('done');

%END