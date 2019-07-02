% Plot events timelocked to event max
addpath /home/mikkel/PD_motor/global_scripts
[dirs, subjs, ~] = PD_proj_setup('betaburst');
addpath /home/mikkel/matlab/export_fig/

subs = find_subs(dirs.megDir);
cd(dirs.megDir);

[PD_subs, PDidx] = intersect(subs,subjs.PD);
[ctrl_subs, ctrlidx] = intersect(subs, subjs.ctrl);

%% Load and arrange data
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

%% Plot settings
sub_lnwdt       = 0.5;
avg_lnwdt       = 2.5;
axis_lnwdt      = 2;
def_fontsize    = 9;
label_fontsize  = 10; 
ttl_fontsize    = 11;

yrange          = [-0.7 1.1];
xrange          = [-150,150];
xtck            = [-100:50:100];
ytck            = [-0.5:0.5:1];

%% Make plot
close all
fig = figure; hold on
set(fig,'Position', [0 0 800 600], 'color','w')

% Ptns 1
subplot(2,2,1); hold on
for i = 1:length(ptns_epo1)
    plot(ptns_epo1{i}.time*1000,ptns_epo1{i}.avg,'-','LineWidth',sub_lnwdt)
end
p1 = plot(ptns_ga1.time*1000,ptns_ga1.avg, 'k-','LineWidth',avg_lnwdt);
ylim(yrange); xlim(xrange);
set(gca, 'LineWidth', 1,'fontweight','bold','fontsize',def_fontsize, ...
        'XTick', xtck, 'xticklabel',{xtck},...
        'YTick', ytck, 'yticklabel',{ytck});
% xlabel('Time (ms)','fontsize',label_fontsize);
ylabel('Amp. (F-score)','fontsize',label_fontsize)
title('PD 1/OFF', 'fontsize', ttl_fontsize);

% Ptns 2
subplot(2,2,3); hold on
for i = 1:length(ptns_epo1)
    plot(ptns_epo2{i}.time*1000,ptns_epo2{i}.avg,'-','LineWidth',sub_lnwdt)
end
p2 = plot(ptns_ga2.time*1000,ptns_ga2.avg, 'k-','LineWidth',avg_lnwdt);
ylim(yrange); xlim(xrange);
set(gca, 'LineWidth', 1,'fontweight','bold','fontsize',def_fontsize, ...
        'XTick', xtck, 'xticklabel',{xtck},...
        'YTick', ytck, 'yticklabel',{ytck});
xlabel('Time (ms)','fontsize',label_fontsize);
ylabel('Amp. (F-score)','fontsize',label_fontsize)
title('PD 2/ON', 'fontsize', ttl_fontsize)

subplot(2,2,2); hold on
for i = 1:length(ctrl_epo1)
    plot(ctrl_epo1{i}.time*1000,ctrl_epo1{i}.avg,'-','LineWidth',sub_lnwdt)
end
p3 = plot(ctrl_ga1.time*1000,ctrl_ga1.avg, 'k-','LineWidth',avg_lnwdt);
ylim([-0.7 1.1]); xlim([-150,150]);
ylim(yrange); xlim(xrange);
set(gca, 'LineWidth', 1,'fontweight','bold','fontsize',def_fontsize, ...
        'XTick', xtck, 'xticklabel',{xtck},...
        'YTick', ytck, 'yticklabel',{ytck});
% xlabel('Time (ms)','fontsize',label_fontsize);
% ylabel('Amplitude (F-score)','fontsize',label_fontsize)
title('Controls 1', 'fontsize', ttl_fontsize);

subplot(2,2,4); hold on
for i = 1:length(ctrl_epo1)
    plot(ctrl_epo2{i}.time*1000,ctrl_epo2{i}.avg,'-','LineWidth',sub_lnwdt)
end
p4 = plot(ctrl_ga2.time*1000,ctrl_ga2.avg, 'k-','LineWidth',avg_lnwdt);
ylim(yrange); xlim(xrange);
set(gca, 'LineWidth', 1,'fontweight','bold','fontsize',def_fontsize, ...
        'XTick', xtck, 'xticklabel',{xtck},...
        'YTick', ytck, 'yticklabel',{ytck});
xlabel('Time (ms)','fontsize',label_fontsize);
% ylabel('Amplitude (F-score)','fontsize',label_fontsize)
title('Controls 2', 'fontsize', ttl_fontsize);

%% Export
export_fig(fullfile(dirs.figures,'timelockedEvent.png'), '-r500', '-p0.05', '-CMYK')

%END