% Plot events timelocked to event max
addpath /home/mikkel/PD_motor/global_scripts
[dirs, subjs, ~] = PD_proj_setup('betaburst');
addpath /home/mikkel/matlab/export_fig/
addpath /home/mikkel/matlab/subtightplot

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
    load(fullfile(dirs.megDir,PD_subs{ii},[PD_subs{ii},'_RsEc1-rawtc.mat']))
    lhdata.trial = {[raw_lh]};
    lhdata.time = {[1:length(raw_lh)]};
    lhdata.label = {'lh_roi'};
    lhdata.fsample = 1000;
    
    mxidx = subvals{1}.bdat.maxidx; 
    trl = [mxidx-150, mxidx+150,repmat(-150,size(mxidx))];
    cfg = [];
    cfg.trl = trl(2:end-2,:); % Skip first/last tirals or NaN
	mx_epo = ft_redefinetrial(cfg,lhdata);
    cfg = [];
    cfg.demean          = 'yes';
    cfg.baselinewindow  = [-0.150 -0.100];
    mx_epo_bs = ft_preprocessing(cfg, mx_epo);
    mx_evo = ft_timelockanalysis([],mx_epo_bs);
    
    ptns_epo1{ii} = mx_evo;
    clear rwft* lhdata* mx_*

    % Session #2
    load(fullfile(dirs.megDir,PD_subs{ii},[PD_subs{ii},'_RsEc2-rawtc.mat']))
    lhdata.trial = {[raw_lh]};
    lhdata.time = {[1:length(raw_lh)]};
    lhdata.label = {'lh_roi'};
    lhdata.fsample = 1000;
    
    mxidx = subvals{2}.bdat.maxidx; 
    trl = [mxidx-150, mxidx+150,repmat(-150,size(mxidx))];
    cfg = [];
    cfg.trl = trl(2:end-2,:);
	mx_epo = ft_redefinetrial(cfg,lhdata);
    cfg = [];
    cfg.demean          = 'yes';
    cfg.baselinewindow  = [-0.150 -0.100];
    mx_epo_bs = ft_preprocessing(cfg, mx_epo);
    mx_evo = ft_timelockanalysis([],mx_epo_bs);
    
    ptns_epo2{ii} = mx_evo; 
    clear rwft* lhdata* mx_*
end

% Ctrl
for ii = 1:length(ctrl_subs)
    % Load and prepare data
    load(fullfile(dirs.megDir,ctrl_subs{ii},'subvals.mat'))
        
    % Session #1
    load(fullfile(dirs.megDir,ctrl_subs{ii},[ctrl_subs{ii},'_RsEc1-rawtc.mat']))
    lhdata.trial = {[raw_lh]};
    lhdata.time = {[1:length(raw_lh)]};
    lhdata.label = {'lh_roi'};
    lhdata.fsample = 1000;
        
    mxidx = subvals{1}.bdat.maxidx; 
    trl = [mxidx-150, mxidx+150,repmat(-150,size(mxidx))];
    cfg = [];
    cfg.trl = trl(2:end-2,:); % Skip first/last tirals or NaN
	mx_epo = ft_redefinetrial(cfg,lhdata);
    cfg = [];
    cfg.demean          = 'yes';
    cfg.baselinewindow  = [-0.150 -0.100];
    mx_epo_bs = ft_preprocessing(cfg, mx_epo);
    mx_evo = ft_timelockanalysis([],mx_epo_bs);
    
    ctrl_epo1{ii} = mx_evo;
    clear rwft* lhdata* mx_*

    % Session #2
    load(fullfile(dirs.megDir,ctrl_subs{ii},[ctrl_subs{ii},'_RsEc2-rawtc.mat']))
    lhdata.trial = {[raw_lh]};
    lhdata.time = {[1:length(raw_lh)]};
    lhdata.label = {'lh_roi'};
    lhdata.fsample = 1000;
    
    mxidx = subvals{2}.bdat.maxidx; 
    trl = [mxidx-150, mxidx+150,repmat(-150,size(mxidx))];
    cfg = [];
    cfg.trl = trl(2:end-2,:);
	mx_epo = ft_redefinetrial(cfg,lhdata);
    cfg = [];
    cfg.demean          = 'yes';
    cfg.baselinewindow  = [-0.150 -0.100];
    mx_epo_bs = ft_preprocessing(cfg, mx_epo);
    mx_evo = ft_timelockanalysis([],mx_epo_bs);
    
    ctrl_epo2{ii} = mx_evo; 
    clear rwft* lhdata* mx_*

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
def_fontsize    = 12;
label_fontsize  = 13; 
ttl_fontsize    = 16;

yrange          = [-1.0 2.0];
xrange          = [-150,150];
xtck            = [-100:50:100];
ytck            = [-0.5:0.5:2.0];
submargin       = [0.1 0.05];
marg_h          = [.1 .05];
marg_w          = [.1 .05];

%% Make plot
close all
fig = figure; hold on
set(fig,'Position', [500 0 800 600], 'color','w')

% Ptns 1
subtightplot(2,2,1,submargin,marg_h,marg_w); hold on
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
subtightplot(2,2,3,submargin,marg_h,marg_w); hold on
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

subtightplot(2,2,2,submargin,marg_h,marg_w); hold on
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

subtightplot(2,2,4,submargin,marg_h,marg_w); hold on
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
export_fig(fullfile(dirs.figures,'timelockedEvent.png'), '-r600', '-p0.05', '-CMYK')

%END
