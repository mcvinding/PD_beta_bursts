% PD-BB: Traditional spectral analysis.
% set paths
clear all
close all

addpath /home/mikkel/PD_motor/global_scripts
[dirs, subs, ~] = PD_proj_setup('betaburst');

cd(dirs.megDir);

subsj = dir(dirs.megDir);                                %Find subjects in folder
subsj = {subsj([subsj.isdir]).name};                       %Make list
subsj = subsj(~(strcmp('.',subs)|strcmp('..',subsj)));     %Remove dots

%% Settings
overwrite = 1;   % Overwirte old files 0=false or 1=true

%% Find PSD of "raw" ROI time-series

for ss = 1:length(subsj)
    subID = subsj{ss};
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
        cfg.overlap = 0;
        epo = ft_redefinetrial(cfg, data);
        
        % Get PSD
        cfg = [];
        cfg.method  = 'mtmfft';
        cfg.output  = 'pow';
        cfg.taper   = 'hanning';
        cfg.foilim  = [1 45];
        cfg.pad     = 'nextpow2';
        
        pow = ft_freqanalysis(cfg, epo);
        rpband = bandpower(pow.powspctrm, pow.freq, [13 30], 'psd');

        % Arrange data
        if any(strcmp(subID,  subs.PD))
            if any(strfind(fname,'RsEc1'))
                ptns_dat1{ss} = pow;
                ptns_relpow1(ss) = rpband;
            elseif any(strfind(fname,'RsEc2'))
                ptns_dat2{ss} = pow;
                ptns_relpow2(ss) = rpband;  
            end
        elseif any(strcmp(subID,  subs.ctrl))              
            if any(strfind(fname,'RsEc1'))
                ctrl_dat1{ss} = pow;
                ctrl_relpow1(ss) = rpband;
            elseif any(strfind(fname,'RsEc2'))
                ctrl_dat2{ss} = pow;
                ctrl_relpow2(ss) = rpband;       
            end
        end
    end
    fprintf('Done with sub %s\n',subID)
end

% Fix stupid array
ptns_relpow1 = ptns_relpow1(ptns_relpow1 ~= 0);
ptns_relpow2 = ptns_relpow2(ptns_relpow2 ~= 0);
ctrl_relpow1 = ctrl_relpow1(ctrl_relpow1 ~= 0);
ctrl_relpow2 = ctrl_relpow2(ctrl_relpow2 ~= 0);

ptns_dat1 = ptns_dat1(~cellfun(@isempty,ptns_dat1));
ptns_dat2 = ptns_dat2(~cellfun(@isempty,ptns_dat2));
ctrl_dat1 = ctrl_dat1(~cellfun(@isempty,ctrl_dat1));
ctrl_dat2 = ctrl_dat2(~cellfun(@isempty,ctrl_dat2));

%% Intermediate t-tests

[h,p] = ttest2(log(ptns_relpow1),log(ctrl_relpow1))
[h,p] = ttest2(log(ptns_relpow2),log(ctrl_relpow2))
[h,p] = ttest(log(ptns_relpow1),log(ptns_relpow2))
[h,p] = ttest(log(ctrl_relpow1),log(ctrl_relpow2))

% Plot histograms
figure; hold on
subplot(1,2,1); h1 = histogram(log(ptns_relpow1),20); hold on
subplot(1,2,1); h2 = histogram(log(ctrl_relpow1),20);
subplot(1,2,2); h3 = histogram(log(ptns_relpow2),20); hold on 
subplot(1,2,2); h4 = histogram(log(ctrl_relpow2),20);

%% grand average
cfg = [];
GA_ptns1 = ft_freqgrandaverage(cfg, ptns_dat1{:});
GA_ptns2 = ft_freqgrandaverage(cfg, ptns_dat2{:});
GA_ctrl1 = ft_freqgrandaverage(cfg, ctrl_dat1{:});
GA_ctrl2 = ft_freqgrandaverage(cfg, ctrl_dat2{:});

% Plot
figure; hold on
plot(GA_ptns1.freq,log(GA_ptns1.powspctrm),'b-')
plot(GA_ptns2.freq,log(GA_ptns2.powspctrm),'b--')
plot(GA_ctrl1.freq,log(GA_ctrl1.powspctrm),'r-')
plot(GA_ctrl2.freq,log(GA_ctrl2.powspctrm),'r--'); hold off

% 
%     plot(pow.freq,pow.powspctrm); hold on;
%     area(pow.freq(51:119),pow.powspctrm(51:119))

%% Cluster stats
cfg = [];
cfg.method = 'montecarlo';
cfg.design = [ones(length(ptns_dat1),1); ones(length(ctrl_dat1),1)*2];
cfg.statistic = 'ft_statfun_indepsamplesT'; 
cfg.numrandomization = 1000;
cfg.tail = 0;
cfg.ivar = 1;
cfg.correctm            = 'cluster';
cfg.computeprob         = 'yes';
cfg.computecritval      = 'yes';
cfg.alpha               = .025;

result1 = ft_freqstatistics(cfg, ptns_dat1{:}, ctrl_dat1{:});
result2 = ft_freqstatistics(cfg, ptns_dat2{:}, ctrl_dat2{:});

cfg.statistic = 'ft_statfun_depsamplesT';
cfg.design = [ones(length(ptns_dat1),1); ones(length(ptns_dat2),1)*2];
cfg.design = [cfg.design, repmat(1:length(ptns_dat1),1,2)'];
cfg.ivar = 1;
cfg.uvar = 2;

result_ptns = ft_freqstatistics(cfg, ptns_dat1{:}, ptns_dat2{:});

cfg.design = [ones(length(ctrl_dat1),1); ones(length(ctrl_dat1),1)*2];
cfg.design = [cfg.design, repmat(1:length(ctrl_dat1),1,2)'];
result_ctrl = ft_freqstatistics(cfg, ctrl_dat1{:}, ctrl_dat2{:});




