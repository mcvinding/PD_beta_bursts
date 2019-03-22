% Cross-correlation between lh-rh hemispheres
addpath /home/mikkel/PD_motor/global_scripts
[dirs, ~, ~] = PD_proj_setup('betaburst');

subs = find_subs(dirs.megDir);
cd(dirs.megDir);

[PD_subs, PDidx] = intersect(subs,subjs.PD);
[ctrl_subs, ctrlidx] = intersect(subs, subjs.ctrl);

%% Load data
ptns_epo1 = cell(size(PD_subs));
ctrl_epo1 = cell(size(ctrl_subs));
ptns_epo2 = cell(size(PD_subs));
ctrl_epo2 = cell(size(ctrl_subs));

% Ptns
for ii = 1:length(PD_subs)
    % Load and prepare data
    load(fullfile(dirs.megDir,PD_subs{ii},'subvals.mat'))
        
    % Session #1
    load(fullfile(dirs.megDir,PD_subs{ii},[PD_subs{ii},'_RsEc1-rawft.mat']))
    [pkval_lh, pkidx_lh] = findpeaks(rwft_lh);
    [pkval_rh, pkidx_rh] = findpeaks(rwft_rh);
    
    plot(rwft_lh); hold on
    plot(pkidx_lh,'rx')
    
    pkdat.label     = {'lh_roi','rh_roi'};
    pkdat.timestamp = {[pkidx_lh],[pkidx_rh]};
    pkdat.dimord    = '{chan}_lead_time_spike';
    
    events = subvals{1}.bdat.event;
    relpk = -(subvals{1}.bdat.maxidx-subvals{1}.bdat.event(:,1));
    cfg = [];
    cfg.trl = [events,  relpk];
    cfg.timestampspersecond = 1000;
    pkdat_trl = ft_spike_maketrials(cfg,pkdat)

    
    cfg       = [];
    cfg.bins  = [0:0.0005:0.1]; % use bins of 0.5 milliseconds
    cfg.param = 'coeffvar'; % compute the coefficient of variation (sd/mn of isis)
    isih = ft_spike_isi(cfg,pkdat_trl);
    
    cfg             = [];
    cfg.binsize     =  0.02; % if cfgPsth.binsize = 'scott' or 'sqrt', we estimate the optimal bin size from the data itself
    cfg.outputunit  = 'rate'; % give as an output the firing rate
    cfg.latency     = [-1 3]; % between -1 and 3 sec.
    cfg.vartriallen = 'yes'; % variable trial lengths are accepted
    cfg.keeptrials  = 'yes'; % keep the psth per trial in the output
    psth = ft_spike_psth(cfg,pkdat_trl);

    cfg              = [];
    cfg.topplotfunc  = 'line'; % plot as a line
    cfg.spikechannel = pkdat_trl.label([1 2]);
    cfg.latency      = [-1 1];
    cfg.errorbars    = 'std'; % plot with the standard deviation
    cfg.interactive  = 'no'; % toggle off interactive mode
    figure, ft_spike_plot_raster(cfg,pkdat_trl, psth)

    
    cfg            = [];
    cfg.latency    = [-1 max(pkdat_trl.trialtime(:))]; % sustained response period
    cfg.keeptrials = 'yes';
    rate = ft_spike_rate(cfg,pkdat_trl);

    cfg             = [];
cfg.maxlag      = 0.2; % maximum 200 ms
cfg.binsize     = 0.01; % bins of 1 ms
cfg.outputunit  = 'proportion'; % make unit area
cfg.latency     = [-1 1];
cfg.vartriallen = 'no'; % do not allow variable trial lengths
cfg.method      = 'xcorr'; % compute the normal cross-correlogram
    Xc = ft_spike_xcorr(cfg,pkdat_trl);

cfg.method      = 'shiftpredictor'; % compute the shift predictor
Xshuff = ft_spike_xcorr(cfg,pkdat_trl);
    
iCmb = 1;
jCmb = 2;
figure
xcSmoothed = conv(squeeze(Xc.xcorr(iCmb,jCmb,:)),ones(1,5)./5,'same'); % do some smoothing
hd = plot(Xc.time(3:end-2),xcSmoothed(3:end-2),'k'); % leave out borders (because of smoothing)
hold on
% xcSmoothed = conv(squeeze(Xshuff.shiftpredictor(iCmb,jCmb,:)),ones(1,5)./5,'same');
% plot(Xc.time(3:end-2),xcSmoothed(3:end-2),'r')
% hold on
xlabel('delay')
ylabel('proportion of coincidences')
title([Xc.label{iCmb} Xc.label{jCmb}])
axis tight






    cfg             = [];
    cfg.maxlag      = 0.02; % maximum 200 ms
    cfg.binsize     = 0.001; % bins of 1 ms
    cfg.outputunit  = 'proportion'; % make unit area
    cfg.latency     = [-2.5 0];
    cfg.vartriallen = 'no'; % do not allow variable trial lengths
    cfg.method      = 'xcorr'; % compute the normal cross-correlogram
    Xc = ft_spike_xcorr(cfg,pkdat_trl);
end
    
    data.trial = {[rwft_lh;rwft_rh]};
    data.time = {[1:length(rwft_lh)]};
    data.label = {'lh_roi','rh_roi'};
    data.fsample = 1000;
    data.dimord = 'chan_time'
    
    
    events = subvals{1}.bdat.event;
    relpk = -(subvals{1}.bdat.maxidx-subvals{1}.bdat.event(:,1));
    cfgt = [];
    cfg.trl = [events,  relpk];
    
    dat_epo = ft_redefinetrial(cfg,data)
        dat_epo.dimord = 'chan_time'

    cfg = [];
    cfg.method = 'mtmfft';
    cfg.output = 'powandcsd';
    cfg.pad = 'nextpow2';
    cfg.taper = 'hanning';
    cfg.foilim = [13 30]
    cfg.tapsmofrq = 1;
    powcsd = ft_freqanalysis(cfg, dat_epo);
    
    cfg = [];
    cfg.method = 'plv';
    coh = ft_connectivityanalysis(cfg,powcsd)
    
    
    mxidx = subvals{1}.bdat.maxidx; 
    trl = [mxidx-150, mxidx+150,repmat(-150,size(mxidx))];
    cfg = [];
    cfg.trl = trl(2:end-2,:); % Skip first/last tirals or NaN
	mx_epo = ft_redefinetrial(cfg,data);
    mx_evo = ft_timelockanalysis([],mx_epo);
    
    ptns_epo1{ii} = mx_evo;
    clear rwft* lhdata* mx_*

    % Session #2
    load(fullfile(dirs.megDir,PD_subs{ii},[PD_subs{ii},'_RsEc2-rawft.mat']))
    data.trial = {[rwft_lh]};
    data.time = {[1:length(rwft_lh)]};
    data.label = {'lh_roi'};
    data.fsample = 1000;
    
    mxidx = subvals{2}.bdat.maxidx; 
    trl = [mxidx-150, mxidx+150,repmat(-150,size(mxidx))];
    cfg = [];
    cfg.trl = trl(2:end-2,:);
	mx_epo = ft_redefinetrial(cfg,data);
    mx_evo = ft_timelockanalysis([],mx_epo);
    ptns_epo2{ii} = mx_evo; 
end

% Ctrl
for ii = 1:length(ctrl_subs)
    % Load and prepare data
    load(fullfile(dirs.megDir,ctrl_subs{ii},'subvals.mat'))
        
    % Session #1
    load(fullfile(dirs.megDir,ctrl_subs{ii},[ctrl_subs{ii},'_RsEc1-rawft.mat']))
    data.trial = {[rwft_lh]};
    data.time = {[1:length(rwft_lh)]};
    data.label = {'lh_roi'};
    data.fsample = 1000;
    
    mxidx = subvals{1}.bdat.maxidx; 
    trl = [mxidx-150, mxidx+150,repmat(-150,size(mxidx))];
    cfg = [];
    cfg.trl = trl(2:end-2,:); % Skip first/last tirals or NaN
	mx_epo = ft_redefinetrial(cfg,data);
    mx_evo = ft_timelockanalysis([],mx_epo);
    
    ctrl_epo1{ii} = mx_evo;
    clear rwft* lhdata* mx_*

    % Session #2
    load(fullfile(dirs.megDir,ctrl_subs{ii},[ctrl_subs{ii},'_RsEc2-rawft.mat']))
    data.trial = {[rwft_lh]};
    data.time = {[1:length(rwft_lh)]};
    data.label = {'lh_roi'};
    data.fsample = 1000;
    
    mxidx = subvals{2}.bdat.maxidx; 
    trl = [mxidx-150, mxidx+150,repmat(-150,size(mxidx))];
    cfg = [];
    cfg.trl = trl(2:end-2,:);
	mx_epo = ft_redefinetrial(cfg,data);
    mx_evo = ft_timelockanalysis([],mx_epo);
    
    ctrl_epo2{ii} = mx_evo; 
end

% Grand average
cfg = [];
cfg.method = 'within';
ptns_ga1 = ft_timelockgrandaverage(cfg,ptns_epo1{:});
ptns_ga2 = ft_timelockgrandaverage(cfg,ptns_epo2{:});
ctrl_ga1 = ft_timelockgrandaverage(cfg,ctrl_epo1{:});
ctrl_ga2 = ft_timelockgrandaverage(cfg,ctrl_epo2{:});

% Plot stuff...
figure; 
subplot(1,2,1); hold on
for i = 1:length(ptns_epo1)
    plot(ptns_epo1{i}.time,ptns_epo1{i}.avg,'-','LineWidth',0.2)
    plot(ptns_epo2{i}.time,ptns_epo2{i}.avg,'--','LineWidth',0.2)

end
plot(ptns_ga1.time,ptns_ga1.avg, 'k-','LineWidth',3);
plot(ptns_ga2.time,ptns_ga2.avg, 'k--','LineWidth',3);
ylim([-0.7 1.1])
subplot(1,2,2); hold on
for i = 1:length(ctrl_epo1)
    plot(ctrl_epo1{i}.time,ctrl_epo1{i}.avg,'-','LineWidth',0.2)
    plot(ctrl_epo2{i}.time,ctrl_epo2{i}.avg,'--','LineWidth',0.2)

end
plot(ptns_ga1.time,ctrl_ga1.avg, 'k-','LineWidth',3);
plot(ptns_ga2.time,ctrl_ga2.avg, 'k--','LineWidth',3);
ylim([-0.7 1.1])

%END