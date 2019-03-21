% Cross-correlation between lh-rh hemispheres
addpath /home/mikkel/PD_motor/global_scripts
[dirs, ~, ~] = PD_proj_setup('betaburst');

subs = find_subs(dirs.megDir);
cd(dirs.megDir);

[PD_subs, PDidx] = intersect(subs,subjs.PD);
[ctrl_subs, ctrlidx] = intersect(subs, subjs.ctrl);

%% Load data
ii=1; % Loop here

load(fullfile(dirs.megDir,subs{ii},[subs{ii},'_RsEc1-rawtc.mat']))  % raw_lh, raw_rh



