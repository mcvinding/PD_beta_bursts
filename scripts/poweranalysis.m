% PD-BB: Traditional spectral analysis.
% set paths
clear all
close all

addpath /home/mikkel/PD_motor/global_scripts
[dirs, subs, ~] = PD_proj_setup('betaburst');

cd(dirs.megDir);

subjs = find_subs(dirs.megDir);     %Find subjects in folder
sub_ptns = intersect(subs.PD, subjs);
sub_ctrl = intersect(subs.ctrl, subjs);

%% Find PSD of "raw" ROI time-series
ptns_dat1 = cell(length(sub_ptns),1);
ptns_dat2 = cell(length(sub_ptns),1);
ctrl_dat1 = cell(length(sub_ptns),1);
ctrl_dat2 = cell(length(sub_ptns),1);

ptns_relpow1 = zeros(length(sub_ptns),1);
ptns_relpow2 = zeros(length(sub_ptns),1);
ctrl_relpow1 = zeros(length(sub_ptns),1);
ctrl_relpow2 = zeros(length(sub_ptns),1);

ptns_bpow1 = zeros(length(sub_ptns),1);
ptns_bpow2 = zeros(length(sub_ptns),1);
ctrl_bpow1 = zeros(length(sub_ptns),1);
ctrl_bpow2 = zeros(length(sub_ptns),1);

% Run
for ss = 1:length(subjs)
    subID = subjs{ss};  
    sub_dir = fullfile(dirs.megDir,subID);
    files = dir(sub_dir);
    files = {files.name};    
    file_idx = find(~cellfun(@isempty,strfind(files,'-rawtc.mat'))); % Name of imported cropped file
    infiles = files(file_idx);
    infiles = sort(infiles);
    
    for f = 1:length(infiles)
        fname = infiles{f};
        load(fullfile(sub_dir,fname))
        
        data.trial = {[raw_lh]};
        data.time = {[1:length(raw_lh)]};
        data.label = {'lh_roi'};
        data.fsample = 1000;
        
        % Make pseudo-epochs
        cfg = [];
        cfg.length  = 3;
        cfg.overlap = .5;       % 0% overlap [?]
        epo = ft_redefinetrial(cfg, data);
        
        % Get PSD
        cfg = [];
        cfg.method  = 'mtmfft';
        cfg.output  = 'pow';
        cfg.taper   = 'hanning';
        cfg.foilim  = [1 48];
        cfg.pad     = 'nextpow2';
        
        pow = ft_freqanalysis(cfg, epo);
        
        % Relative beta power
        b_pow = bandpower(pow.powspctrm, pow.freq, [13 30], 'psd');
        all_pow = bandpower(pow.powspctrm, pow.freq, 'psd');
        relpow = b_pow/all_pow;
        
        % Arrange data [you should preallocate]
        if any(strcmp(subID,  sub_ptns))
            kk = find(~cellfun(@isempty,strfind(sub_ptns,subID))); % Name of imported cropped file
            if any(strfind(fname,'RsEc1'))
                ptns_dat1{kk} = pow;
                ptns_relpow1(kk) = relpow;
                ptns_bpow1(kk) = b_pow;
            elseif any(strfind(fname,'RsEc2'))
                ptns_dat2{kk} = pow;
                ptns_relpow2(kk) = relpow;  
                ptns_bpow2(kk) = b_pow;
            end
        elseif any(strcmp(subID,  sub_ctrl))      
            kk = find(~cellfun(@isempty,strfind(sub_ctrl,subID))); % Name of imported cropped file
            if any(strfind(fname,'RsEc1'))
                ctrl_dat1{kk} = pow;
                ctrl_relpow1(kk) = relpow;
                ctrl_bpow1(kk) = b_pow;
            elseif any(strfind(fname,'RsEc2'))
                ctrl_dat2{kk} = pow;
                ctrl_relpow2(kk) = relpow;  
                ctrl_bpow2(kk) = b_pow;
            end
        end
    end
    fprintf('Done with sub %s\n',subID)
end

%% Save
disp('Saving...');
save('/home/mikkel/PD_motor/rest_ec/groupanalysis/PSD_data.mat', 'ptns_dat1','ptns_dat2','ctrl_dat1','ctrl_dat2')
save('/home/mikkel/PD_motor/rest_ec/groupanalysis/B_relpow.mat', 'ptns_relpow1','ptns_relpow2','ctrl_relpow1','ctrl_relpow2')
save('/home/mikkel/PD_motor/rest_ec/groupanalysis/B_pow.mat', 'ptns_bpow1','ptns_bpow2','ctrl_bpow1','ctrl_bpow2')
disp('done');

%% Reload
disp('Loading...');
load('/home/mikkel/PD_motor/rest_ec/groupanalysis/PSD_data.mat')
disp('done');

% %% Intermediate t-tests
% 
% [h,p] = ttest2(log(ptns_bpow1),log(ctrl_bpow1))
% [h,p] = ttest2(log(ptns_bpow2),log(ctrl_bpow2))
% [h,p] = ttest(log(ptns_bpow1),log(ptns_bpow2))
% [h,p] = ttest(log(ctrl_bpow1),log(ctrl_bpow2))
% 
[h,p] = ttest2(ptns_relpow1,ctrl_relpow1)
[h,p] = ttest2(ptns_relpow2,ctrl_relpow2)
[h,p] = ttest((ptns_relpow1),(ptns_relpow2))
[h,p] = ttest(ctrl_relpow1,(ctrl_relpow2))
% 
% % Plot histograms
% figure; hold on; title('Band power')
% subplot(1,2,1); h1 = histogram(ptns_bpow1,20); hold on
% subplot(1,2,1); h2 = histogram(ctrl_bpow1,20);
% subplot(1,2,2); h3 = histogram(ptns_bpow2,20); hold on 
% subplot(1,2,2); h4 = histogram(ctrl_bpow2,20);
% 
figure; hold on; title('Relative power')
subplot(1,2,1); h1 = histogram(ptns_relpow1,20); hold on
subplot(1,2,1); h2 = histogram(ctrl_relpow1,20);
subplot(1,2,2); h3 = histogram(ptns_relpow2,20); hold on 
subplot(1,2,2); h4 = histogram(ctrl_relpow2,20);
% 
% Plot log histograms
figure; hold on; title('log-relative power')
subplot(1,2,1); h1 = histogram(log(ptns_relpow1),20); hold on
subplot(1,2,1); h2 = histogram(log(ctrl_relpow1),20);
subplot(1,2,2); h3 = histogram(log(ptns_relpow2),20); hold on 
subplot(1,2,2); h4 = histogram(log(ctrl_relpow2),20);

figure; hold on; title('log-band power')
subplot(1,2,1); h1 = histogram(log(ptns_bpow1),20); hold on
subplot(1,2,1); h2 = histogram(log(ctrl_bpow1),20);
subplot(1,2,2); h3 = histogram(log(ptns_bpow2),20); hold on 
subplot(1,2,2); h4 = histogram(log(ctrl_bpow2),20);

cfg = [];
cfg.parameter = 'powspctrm';
subplot(2,2,1); ft_singleplotER(cfg, ptns_dat1{:})
subplot(2,2,2); ft_singleplotER(cfg, ptns_dat2{:})
subplot(2,2,3); ft_singleplotER(cfg, ctrl_dat1{:})
subplot(2,2,4); ft_singleplotER(cfg, ctrl_dat2{:})

%% grand average [move to plotting scripts]
cfg = [];
cfg.keepindividual = 'yes';
GA.ptns1 = ft_freqgrandaverage(cfg, ptns_dat1{:});
GA.ptns2 = ft_freqgrandaverage(cfg, ptns_dat2{:});
GA.ctrl1 = ft_freqgrandaverage(cfg, ctrl_dat1{:});
GA.ctrl2 = ft_freqgrandaverage(cfg, ctrl_dat2{:});

% Plot
figure; hold on
plot(GA.ptns1.freq,median(squeeze(GA.ptns1.powspctrm)),'b-')
plot(GA.ptns2.freq,median(squeeze(GA.ptns2.powspctrm)),'b:')
plot(GA.ctrl1.freq,median(squeeze(GA.ctrl1.powspctrm)),'r-')
plot(GA.ctrl2.freq,median(squeeze(GA.ctrl2.powspctrm)),'r:'); hold off

figure; hold on
plot(GA.ptns1.freq,median(squeeze(GA.ptns1.powspctrm)),'b-')
plot(GA.ptns1.freq,mean(squeeze(GA.ptns1.powspctrm)),'b:')
plot(GA.ptns1.freq,median(squeeze(GA.ptns2.powspctrm)),'r-')
plot(GA.ptns1.freq,mean(squeeze(GA.ptns2.powspctrm)),'r:')

figure; hold on
plot(GA.ptns1.freq,median(squeeze(GA.ctrl1.powspctrm)),'b-')
plot(GA.ptns1.freq,mean(squeeze(GA.ctrl1.powspctrm)),'b:')
plot(GA.ptns1.freq,median(squeeze(GA.ctrl2.powspctrm)),'r-')
plot(GA.ptns1.freq,mean(squeeze(GA.ctrl2.powspctrm)),'r:')


figure; hold on
plot(GA.ptns1.freq,(GA.ptns1.powspctrm),'b-')
plot(GA.ptns2.freq,(GA.ptns2.powspctrm),'b--')
plot(GA.ctrl1.freq,(GA.ctrl1.powspctrm),'r-')
plot(GA.ctrl2.freq,(GA.ctrl2.powspctrm),'r--'); hold off

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
cfg.computeprob         = 'yes';
cfg.computecritval      = 'yes';
cfg.alpha               = .025;
cfg.numrandomization    = 1000;
cfg.tail                = 0;

% Between groups
cfg.statistic = 'ft_statfun_indepsamplesT'; 
cfg.design = [ones(length(ptns_dat1),1); ones(length(ctrl_dat1),1)*2];
cfg.ivar = 1;

result1 = ft_freqstatistics(cfg, ptns_dat1{:}, ctrl_dat1{:});
result2 = ft_freqstatistics(cfg, ptns_dat2{:}, ctrl_dat2{:});

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

% Interaction
cfg.statistic = 'ft_statfun_indepsamplesT'; 
cfg.design = [ones(length(ptns_dat1),1); ones(length(ctrl_dat1),1)*2];
cfg.ivar = 1;

resultX = ft_freqstatistics(cfg, ptns_diff{:}, ctrl_diff{:});

%% Save results
disp('Saving...');
save('/home/mikkel/PD_motor/rest_ec/groupanalysis/PSD_stats.mat', ...
    'result1','result2','result_ctrl','result_ptns','resultX')
save('/home/mikkel/PD_motor/rest_ec/groupanalysis/PSD_GA.mat', 'GA')
disp('done');

%END