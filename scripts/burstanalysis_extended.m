% Get beta bursts summary across all thresholds
% set paths
clear all
close all

addpath /home/mikkel/PD_motor/global_scripts
[dirs, subjs, ~] = PD_proj_setup('betaburst');

subs = find_subs(dirs.megDir);                                %Find subjects in folder
cd(dirs.megDir);

[PD_subs] = intersect(subs,subjs.PD);
[ctrl_subs] = intersect(subs, subjs.ctrl);

%% Settings
overwrite = 0;   % Overwirte old files 0=false or 1=true

steps = 0.1:0.1:5;

%% Get beta summary across thresholds
% Function specs
cfg = [];
cfg.length      = 3;
cfg.overlap     = 0;
cfg.steps       = steps;
cfg.corrtype    = 'amp';
cfg.cutofftype  = 'med'; 
cfg.halfmax     = 'mixed';
cfg.makeplot    = 'no';

% Init
neve_long = [];
step_long = [];
subs_long = [];
sesi_long = [];
grup_long = [];

% Run
for ss = 1:length(subs)
    subvals  = [];
    subID    = subs{ss};
    sub_dir  = fullfile(dirs.megDir,subID);
    files    = dir(sub_dir);
    files    = {files.name};
    file_idx = find(~cellfun(@isempty,strfind(files,'-hilbt.mat'))); % Name of imported cropped file
    infiles  = files(file_idx);
    infiles  = sort(infiles);
    outfname = fullfile(sub_dir,'subvals_ext.mat');
    if exist(outfname,'file') && ~overwrite
        fprintf('File %s already exists.\Pass\n',outfname)
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
        
        % Arrange for export
        neve_long = [neve_long, subvals{f}.n_events];
        step_long = [step_long, steps];
        subs_long = [subs_long, repmat({subID},1,length(steps))];
        sesi_long = [sesi_long, repmat(f,1,length(steps))];
        if any(strcmp(ctrl_subs,subID))
            grup_long = [grup_long, repmat({'ctrl'},1,length(steps))];
        else
            grup_long = [grup_long, repmat({'ptns'},1,length(steps))];
        end
    end
    
    save(outfname,'subvals')
end

save('/home/mikkel/PD_motor/rest_ec/groupanalysis/nevent_extended.mat', ...
    'neve_long', 'step_long', 'subs_long', 'sesi_long', 'grup_long');


% End