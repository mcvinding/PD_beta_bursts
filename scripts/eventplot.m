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
    lhdata.trial = {[rwft_lh]};
    lhdata.time = {[1:length(rwft_lh)]};
    lhdata.label = {'lh_roi'};
    lhdata.fsample = 1000;
    
    mxidx = subvals{1}.bdat.maxidx; 
    trl = [mxidx-150, mxidx+150,repmat(-150,size(mxidx))];
    cfg = [];
    cfg.trl = trl(2:end-2,:); % Skip first/last tirals or NaN
	mx_epo = ft_redefinetrial(cfg,lhdata);
    mx_evo = ft_timelockanalysis([],mx_epo);
    
    ptns_epo1{ii} = mx_evo;
    clear rwft* lhdata* mx_*

    % Session #2
    load(fullfile(dirs.megDir,PD_subs{ii},[PD_subs{ii},'_RsEc2-rawft.mat']))
    lhdata.trial = {[rwft_lh]};
    lhdata.time = {[1:length(rwft_lh)]};
    lhdata.label = {'lh_roi'};
    lhdata.fsample = 1000;
    
    mxidx = subvals{2}.bdat.maxidx; 
    trl = [mxidx-150, mxidx+150,repmat(-150,size(mxidx))];
    cfg = [];
    cfg.trl = trl(2:end-2,:);
	mx_epo = ft_redefinetrial(cfg,lhdata);
    mx_evo = ft_timelockanalysis([],mx_epo);
    ptns_epo2{ii} = mx_evo; 
end

% Ctrl
for ii = 1:length(ctrl_subs)
    % Load and prepare data
    load(fullfile(dirs.megDir,ctrl_subs{ii},'subvals.mat'))
        
    % Session #1
    load(fullfile(dirs.megDir,ctrl_subs{ii},[ctrl_subs{ii},'_RsEc1-rawft.mat']))
    lhdata.trial = {[rwft_lh]};
    lhdata.time = {[1:length(rwft_lh)]};
    lhdata.label = {'lh_roi'};
    lhdata.fsample = 1000;
    
    mxidx = subvals{1}.bdat.maxidx; 
    trl = [mxidx-150, mxidx+150,repmat(-150,size(mxidx))];
    cfg = [];
    cfg.trl = trl(2:end-2,:); % Skip first/last tirals or NaN
	mx_epo = ft_redefinetrial(cfg,lhdata);
    mx_evo = ft_timelockanalysis([],mx_epo);
    
    ctrl_epo1{ii} = mx_evo;
    clear rwft* lhdata* mx_*

    % Session #2
    load(fullfile(dirs.megDir,ctrl_subs{ii},[ctrl_subs{ii},'_RsEc2-rawft.mat']))
    lhdata.trial = {[rwft_lh]};
    lhdata.time = {[1:length(rwft_lh)]};
    lhdata.label = {'lh_roi'};
    lhdata.fsample = 1000;
    
    mxidx = subvals{2}.bdat.maxidx; 
    trl = [mxidx-150, mxidx+150,repmat(-150,size(mxidx))];
    cfg = [];
    cfg.trl = trl(2:end-2,:);
	mx_epo = ft_redefinetrial(cfg,lhdata);
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