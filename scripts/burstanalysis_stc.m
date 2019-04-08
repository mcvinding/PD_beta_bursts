% Get beta bursts: 1) find threshold 2) use threshold to determine beta
% events.
% set paths
clear all
close all

addpath /home/mikkel/PD_motor/global_scripts
[dirs, ~, ~] = PD_proj_setup('betaburst');

subs = find_subs(dirs.megDir);                                %Find subjects in folder
cd(dirs.megDir);

%% Settings
overwrite = 0;   % Overwirte old files 0=false or 1=true

steps = 0:0.1:5;

%% Find peaks: testing multiple methods based on the litterature
rho_sdamp = nan(length(subs),length(steps),2);
rho_sdpow = nan(length(subs),length(steps),2);
rho_mdamp = nan(length(subs),length(steps),2);
rho_mdpow = nan(length(subs),length(steps),2);

for ss = 1:length(subs)
    subID = subs{ss};
    sub_dir = fullfile(dirs.megDir,subID);
    files = dir(sub_dir);
    files = {files.name};    
    file_idx = find(~cellfun(@isempty,strfind(files,'-hilbt.mat'))); % Name of imported cropped file
    infiles = files(file_idx);
    infiles = sort(infiles);
    
    for f = 1:length(infiles)
        fname = infiles{f};
        load(fullfile(sub_dir,fname))
            
        % Make pseudo data
        lhdata.trial = {[hilb_lh]};
        lhdata.time = {[1:length(hilb_lh)]};
        lhdata.label = {'lh_roi'};
        lhdata.fsample = 1000;
            
        % Find correlations across thresholds
        cfg = [];
        cfg.length      = 3;
        cfg.overlap     = 0;
        cfg.steps       = steps;
        
        cfg.corrtype    = 'amp';
        cfg.cutofftype  = 'sd';
        [~, rho_sdamp(ss,:,f)] = find_betaevents(cfg,lhdata);
        cfg.corrtype    = 'pow';
        [~, rho_sdpow(ss,:,f)] = find_betaevents(cfg,lhdata);
        cfg.cutofftype  = 'med';
        [~, rho_mdpow(ss,:,f)] = find_betaevents(cfg,lhdata);
        cfg.corrtype    = 'amp';
        [~, rho_mdamp(ss,:,f)] = find_betaevents(cfg,lhdata);

    end
end
save('/home/mikkel/PD_motor/rest_ec/groupanalysis/rhomats.mat','rho_sdamp','rho_sdpow','rho_mdamp','rho_mdpow')
disp('done')

%% Compare methods
% Load data
load('/home/mikkel/PD_motor/rest_ec/groupanalysis/rhomats.mat')
disp('done')

% Rho1: sd based threshold ~ mean amplitude of epoch
% Rho2: median based threshold ~ mean amplitude of epoch
% Rho3: sd based threshold ~ mean power (amp^2) of epoch
% Rho4: median based threshold ~ mean power (amp^2) of epoch

figure; plot(steps,nanmean(rho_sdamp,3));
figure; plot(steps,nanmean(rho_sdpow,3)); 
figure; plot(steps,nanmean(rho_mdamp,3)); 
figure; plot(steps,nanmean(rho_mdpow,3)); 

% cutoff_sdamp = find_threshold(rho_sdamp, steps, 1); title('sd amp')
% cutoff_sdpow = find_threshold(rho_sdpow, steps, 1); title('sd pow')
cutoff_mdamp = find_threshold(rho_mdamp, steps, 1); title('med amp')
% cutoff_mdpow = find_threshold(rho_mdpow, steps, 1); title('med pow')

%% Get beta summary
cfg = [];
cfg.length      = 3;
cfg.overlap     = 0;
cfg.steps       = cutoff_mdamp;
cfg.corrtype    = 'amp';
cfg.cutofftype  = 'med'; 
cfg.halfmax     = 'mixed';
cfg.makeplot    = 'yes';

for ss = 1:length(subs)
    subvals  = [];
    subID    = subs{ss};
    sub_dir  = fullfile(dirs.megDir,subID);
    files    = dir(sub_dir);
    files    = {files.name};
    file_idx = find(~cellfun(@isempty,strfind(files,'-hilbt.mat'))); % Name of imported cropped file
    infiles  = files(file_idx);
    infiles  = sort(infiles);
    outfname = fullfile(sub_dir,'subvals.mat');
    if exist(outfname,'file') && ~overwrite
        continue
    end
    
    for f = 1:length(infiles)
        fname = infiles{f};
        load(fullfile(sub_dir,fname))
        
        % Make ft style data
        lhdata.trial   = {hilb_lh};
        lhdata.time    = {1:length(hilb_lh)};
        lhdata.label   = {'lh_roi'};
        lhdata.fsample = 1000;

        [subvals{f}] = find_betaevents(cfg,lhdata);
    end
    save(outfname,'subvals')
end