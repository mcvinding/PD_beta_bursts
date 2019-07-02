% PD-BB: Traditional spectral analysis using running window and FFT.
% set paths
clear all
close all

addpath /home/mikkel/PD_motor/global_scripts
[dirs, subs, ~] = PD_proj_setup('betaburst');

cd(dirs.megDir);

subjs = find_subs(dirs.megDir);     %Find subjects in folder
sub_ptns = intersect(subs.PD, subjs);
sub_ctrl = intersect(subs.ctrl, subjs);

%% Find PSD of "raw" ROI time-series
ptns_dat1 = cell(length(sub_ptns),1);
ptns_dat2 = cell(length(sub_ptns),1);
ctrl_dat1 = cell(length(sub_ptns),1);
ctrl_dat2 = cell(length(sub_ptns),1);

ptns_relpow1 = zeros(length(sub_ptns),1);
ptns_relpow2 = zeros(length(sub_ptns),1);
ctrl_relpow1 = zeros(length(sub_ptns),1);
ctrl_relpow2 = zeros(length(sub_ptns),1);

ptns_bpow1 = zeros(length(sub_ptns),1);
ptns_bpow2 = zeros(length(sub_ptns),1);
ctrl_bpow1 = zeros(length(sub_ptns),1);
ctrl_bpow2 = zeros(length(sub_ptns),1);

% Run
for ss = 1:length(subjs)
    subID = subjs{ss};  
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
        cfg.overlap = .5;       % 0% overlap [?]
        epo = ft_redefinetrial(cfg, data);
        
        % Get PSD
        cfg = [];
        cfg.method  = 'mtmfft';
        cfg.output  = 'pow';
        cfg.taper   = 'hanning';
        cfg.foilim  = [1 48];
        cfg.pad     = 'nextpow2';
        
        pow = ft_freqanalysis(cfg, epo);
        
        % Relative beta power
        b_pow = bandpower(pow.powspctrm, pow.freq, [13 30], 'psd');
        all_pow = bandpower(pow.powspctrm, pow.freq, 'psd');
        relpow = b_pow/all_pow;
        
        % Arrange data [you should preallocate]
        if any(strcmp(subID,  sub_ptns))
            kk = find(~cellfun(@isempty,strfind(sub_ptns,subID))); % Name of imported cropped file
            if any(strfind(fname,'RsEc1'))
                ptns_dat1{kk} = pow;
                ptns_relpow1(kk) = relpow;
                ptns_bpow1(kk) = b_pow;
            elseif any(strfind(fname,'RsEc2'))
                ptns_dat2{kk} = pow;
                ptns_relpow2(kk) = relpow;  
                ptns_bpow2(kk) = b_pow;
            end
        elseif any(strcmp(subID,  sub_ctrl))      
            kk = find(~cellfun(@isempty,strfind(sub_ctrl,subID))); % Name of imported cropped file
            if any(strfind(fname,'RsEc1'))
                ctrl_dat1{kk} = pow;
                ctrl_relpow1(kk) = relpow;
                ctrl_bpow1(kk) = b_pow;
            elseif any(strfind(fname,'RsEc2'))
                ctrl_dat2{kk} = pow;
                ctrl_relpow2(kk) = relpow;  
                ctrl_bpow2(kk) = b_pow;
            end
        end
    end
    fprintf('Done with sub %s\n',subID)
end

%% grand average
cfg = [];
cfg.keepindividual = 'yes';
GA.ptns1 = ft_freqgrandaverage(cfg, ptns_dat1{:});
GA.ptns2 = ft_freqgrandaverage(cfg, ptns_dat2{:});
GA.ctrl1 = ft_freqgrandaverage(cfg, ctrl_dat1{:});
GA.ctrl2 = ft_freqgrandaverage(cfg, ctrl_dat2{:});

%% Save
disp('Saving...');
save('/home/mikkel/PD_motor/rest_ec/groupanalysis/PSD_data.mat', ...
     'ptns_dat1','ptns_dat2','ctrl_dat1','ctrl_dat2')
save('/home/mikkel/PD_motor/rest_ec/groupanalysis/B_relpow.mat', ...
     'ptns_relpow1','ptns_relpow2','ctrl_relpow1','ctrl_relpow2')
save('/home/mikkel/PD_motor/rest_ec/groupanalysis/B_pow.mat', ... 
     'ptns_bpow1','ptns_bpow2','ctrl_bpow1','ctrl_bpow2')
save('/home/mikkel/PD_motor/rest_ec/groupanalysis/PSD_GA.mat', 'GA')
disp('done');

