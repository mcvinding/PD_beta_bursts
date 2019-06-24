% Cross-correlation between lh-rh hemispheres
addpath /home/mikkel/PD_motor/global_scripts
[dirs, subjs, ~] = PD_proj_setup('betaburst');

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
% for ii = 1:length(PD_subs)
ii = 1;

    % Load and prepare data
    load(fullfile(dirs.megDir,PD_subs{ii},'subvals.mat'))
        
    % Session #1
    load(fullfile(dirs.megDir,PD_subs{ii},[PD_subs{ii},'_RsEc1-rawft.mat']))
    [pkval_lh, pkidx_lh] = findpeaks(rwft_lh);
    [pkval_rh, pkidx_rh] = findpeaks(rwft_rh);
    
%     plot(rwft_lh); hold on
    
    x = rand(1,1000)*1000*180;
    
    pkdat.label     = {'lh_roi','rh_roi'};
    pkdat.timestamp = {[pkidx_lh],[pkidx_rh]};
    pkdat.dimord    = '{chan}_lead_time_spike';
    
    events = subvals{1}.bdat.event;
    nonevents = [events(1:end-1,2), events(2:end,1)];
    relpk = -(subvals{1}.bdat.maxidx-subvals{1}.bdat.event(:,1));
    cfg = [];
    cfg.trl = [events, zeros(length(events),1)];
    cfg.timestampspersecond = 1000;
    pkdat_trl = ft_spike_maketrials(cfg,pkdat)
    cfg.trl = [nonevents, zeros(length(nonevents),1)];
    pkdat_ntrl = ft_spike_maketrials(cfg,pkdat)
    cfg.trl = [0, 180000, 0];
    pkdat_all = ft_spike_maketrials(cfg,pkdat)


    cfg             = [];
    
%%    
cfg.maxlag      = 0.1;  % maximum 200 ms
cfg.binsize     = 0.02; % bins of 10 ms
cfg.outputunit  = 'proportion'; % make unit area
cfg.latency     = [-inf inf];
cfg.vartriallen = 'yes'; % do not allow variable trial lengths
cfg.method      = 'xcorr'; % compute the normal cross-correlogram
Xc = ft_spike_xcorr(cfg,pkdat_trl);
Xcn = ft_spike_xcorr(cfg,pkdat_ntrl);
Xca = ft_spike_xcorr(cfg,pkdat_all);

figure; hold on; 
plot(Xc.time,squeeze(Xc.xcorr(1,2,:)))
plot(Xcn.time,squeeze(Xcn.xcorr(1,2,:)))
plot(Xca.time,squeeze(Xca.xcorr(1,2,:)))

cfg.method      = 'shiftpredictor'; % compute the shift predictor
Xshuff = ft_spike_xcorr(cfg,pkdat_trl);
    
iCmb = 1;
jCmb = 2;
figure
xcSmoothed = conv(squeeze(Xc.xcorr(iCmb,jCmb,:)),ones(1,5)./5,'same'); % do some smoothing
hd = plot(Xc.time(3:end-2),xcSmoothed(3:end-2),'k'); % leave out borders (because of smoothing)
hold on
xcSmoothed = conv(squeeze(Xshuff.shiftpredictor(iCmb,jCmb,:)),ones(1,5)./5,'same');
plot(Xc.time(3:end-2),xcSmoothed(3:end-2),'r')
hold on
xlabel('delay')
ylabel('proportion of coincidences')
title([Xc.label{iCmb} Xc.label{jCmb}])
axis tight

%% Native xcorr


data.trial = {[rwft_lh;rwft_rh]};
data.time = {[1:length(rwft_lh)]};
data.label = {'lh_roi','rh_roi'};
data.fsample = 1000;
data.dimord = 'chan_time';

events = subvals{1}.bdat.event;
relpk = -(subvals{1}.bdat.maxidx-subvals{1}.bdat.event(:,1));
cfg = [];
cfg.trl = [events,  relpk];
data_epo = ft_redefinetrial(cfg,data)

tt = data_epo.trial{1};
x = tt(1,:);
y = tt(2,:);

[Xx, lag] = xcorr(x,y,'biased');



dt = linspace(min(data_epo.time{1}),max(data_epo.time{1}),2*length(data_epo.time{1})-1)

figure; plot(data_epo.time{1},tt)
figure; plot(lag/1000,Xx)

xtrl = cell(length(data_epo.trial),1);
xtim = cell(length(data_epo.trial),1);
figure; hold on
for i = 1:length(data_epo.trial)
    tt = data_epo.trial{i};
    x = tt(1,:);
    y = tt(2,:);
    [Xx, lag] = xcorr(x,y,'biased');
    xtrl{i} = Xx;
    xtim{i} = lag/1000;
    
%     dt = linspace(min(data_epo.time{1}),max(data_epo.time{1}),2*length(data_epo.time{1})-1)
    plot(lag/1000,Xx)
end

xdata.trial = xtrl;
xdata.time  = xtim;
xdata.label = {'xcor'};
xdata.fsample = 1000;
xdata.dimord = 'chan_time';

cfg = [];
cfg.vartrllength = 2;
xavg = ft_timelockanalysis(cfg,xdata)

figure; hold on
plot(xavg.time,xavg.avg)
plot(xavg.time,xavg.avg+xavg.var,'b--')
plot(xavg.time,xavg.avg-xavg.var,'b--')


    cfg             = [];
    cfg.maxlag      = 0.02; % maximum 200 ms
    cfg.binsize     = 0.001; % bins of 1 ms
    cfg.outputunit  = 'proportion'; % make unit area
    cfg.latency     = [-2.5 0];
    cfg.vartriallen = 'no'; % do not allow variable trial lengths
    cfg.method      = 'xcorr'; % compute the normal cross-correlogram
    Xc = ft_spike_xcorr(cfg,pkdat_trl);
end
    

    
    
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