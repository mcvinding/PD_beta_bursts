% Extract data timelocked to event maximum and plot
clear all
close all

addpath /home/mikkel/PD_motor/global_scripts
[dirs, ~, ~] = PD_proj_setup('betaburst');

subs = find_subs(dirs.megDir);
cd(dirs.megDir);

[PD_subs, PDidx] = intersect(subs,subjs.PD);
[ctrl_subs, ctrlidx] = intersect(subs, subjs.ctrl);

%% Load data

for ii = 1:length(subs)
    % Session 1
    load(fullfile(dirs.megDir,subs{ii},'pkidx1.mat'))
    load(fullfile(dirs.megDir,subs{ii},[subs{ii},'_RsEc1-rawtc.mat']))  % raw_lh, raw_rh
    load(fullfile(dirs.megDir,subs{ii},[subs{ii},'_RsEc1-hilbt.mat']))  % raw_lh, raw_rh


    lh_filt = ft_preproc_bandpassfilter(raw_lh, 1000, [13 30]);
    
    plot(raw_lh)
    
    
    % session 2
    load(fullfile(dirs.megDir,subs{ii},'pkidx2.mat'))
    
      
    % Save pkidx for reading and plotting in MNE-Py
    maxidx1 = subvals{1}.bdat.maxidx;
    maxidx2 = subvals{2}.bdat.maxidx;

end




