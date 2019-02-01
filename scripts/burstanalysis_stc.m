% Get beta bursts
% * Clean up: decide on cutoff algorithm.
% * save decision curves
% set paths
clear all
close all

addpath /home/mikkel/PD_motor/global_scripts
[dirs, ~, ~] = PD_proj_setup('betaburst');

cd(dirs.megDir);

subs = dir(dirs.megDir);                                %Find subjects in folder
subs = {subs([subs.isdir]).name};                       %Make list
subs = subs(~(strcmp('.',subs)|strcmp('..',subs)));     %Remove dots

%% Settings
overwrite = 1;   % Overwirte old files 0=false or 1=true

%% Find peaks: testing multiple methods based on the litterature
% * Correalte epoch pow with N-peaks in epoch base on whole data peak
%   threshold: both sd and median based cutoffs
% * Correalte epoch pow with N-peaks in epoch base on epoch based peak
%   threshold: both sd and median based cutoffs

steps = 0:0.1:5;
rho1mat = zeros(length(subs),length(steps),2);
rho2mat = zeros(length(subs),length(steps),2);
rho3mat = zeros(length(subs),length(steps),2);
rho4mat = zeros(length(subs),length(steps),2);

for ss = 1:length(subs)
    subID = subs{ss};
    sub_dir = fullfile(dirs.megDir,subID);
    files = dir(sub_dir);
    files = {files.name};
    
    outfname = fullfile(sub_dir,'subvals2.mat');
    if exist(outfname,'file') && ~overwrite
        continue
    end
   
    % Skip for now
    if strcmp('0327',subID)
        continue
    end
    % *Sensor level analysis*
%     file_idx = find(~cellfun(@isempty,strfind(files,'-hil.mat'))); % Name of imported cropped file
%     infiles = files(file_idx);
    
    % Load data... (LOOP HERE)
%     load(fullfile(sub_dir,infiles{1}));
    
    % Select data
%     cfg = [];
%     cfg.channel = {'MEG0232+0233', 'MEG0442+0443', 'MEG1622+1623', 'MEG1812+1813'};
%     cfg.avgoverchan = 'yes';
%     datachan = ft_selectdata(cfg, data_hilbt);

    % *Source level analysis*    
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
    
        % Make pseudo-tirals
        cfg = [];
        cfg.length  = 4;
        cfg.overlap = 0;
        dat = ft_redefinetrial(cfg,lhdata);

        trl = dat.sampleinfo;
        trl = trl-trl(1,1)+1; %Corrent for non-zero sample offset

        % Get PSD
        cfg = [];
        cfg.method      = 'mtmfft';
        cfg.taper       = 'hanning';
        cfg.foilim      = [10 30];
        cfg.tapsmofrq   = 2;
        cfg.keeptrials = 'yes';
        cfg.output = 'pow';
        pow = ft_freqanalysis(cfg,dat);
        pow.avgpowspctrm = squeeze(mean(pow.powspctrm,1));

        cfg = [];
        cfg.avgoverfreq = 'yes';
        powsum = ft_selectdata(cfg, pow);

        % Plot PSD
%         plot(pow.freq,pow.avgpowspctrm);

        % Cutoffs
        data = hilb_lh;
    %     tim = length(data)/1000;
        pkmat = zeros(length(trl),1);
        pkmat2 = zeros(length(trl),1);
        pkmat3 = zeros(length(trl),1);
        pkmat4 = zeros(length(trl),1);
        n_events = zeros(length(trl),1);
%     bdat = struct();

        med = median(data);
        sd = std(data);

        rho1 = zeros(length(steps),1);
        rho2 = zeros(length(steps),1);
        rho3 = zeros(length(steps),1);
        rho4 = zeros(length(steps),1);
    
        bdat = struct();
        for i = 1:length(steps)
            burst = data >= med+sd*steps(i);   
            burst2 = data >= med+med*steps(i);

            for k = 1:length(trl)
                tmp = burst(trl(k,1):trl(k,2));
                pkmat(k) = sum((diff(tmp)==1));

                % units of median (cf. Chin et al)
                tmp2 = burst2(trl(k,1):trl(k,2));
                pkmat3(k) = sum((diff(tmp2)==1)); 

                % Single tiral median
                tmpdat = data(trl(k,1):trl(k,2));
                tmpmed = median(tmpdat);
                tmpsd = std(tmpdat);
                tmpbst = tmpdat >= tmpmed+tmpsd*steps(i);
%                 tim = 1:length(tmpdat);

                pkmat2(k) = sum((diff(tmpbst)==1));

                tmpbst2 = tmpdat >= tmpmed+tmpmed*steps(i);
                pkmat4(k) = sum((diff(tmpbst2)==1));

    %             plot(tim,tmpdat); hold on
    %             plot(tim(tmpbst),tmpdat(tmpbst),'r');
    %             plot(tim(tmp),tmpdat(tmp),'g');
            end
            rho1(i) = corr(pkmat,powsum.powspctrm);
            rho2(i) = corr(pkmat2,powsum.powspctrm);
            rho3(i) = corr(pkmat3,powsum.powspctrm);
            rho4(i) = corr(pkmat4,powsum.powspctrm);
        
            % Get summaries
            dburst = diff([0 burst 0]);
            n_events(i) = sum(dburst==1);

            startb  = find(dburst==1);      % Start of burst
            endb    = find(dburst==-1);     % End of burst´
            blen    = endb-startb;          % Length of burst

            if any(blen<0)
                error('negative length')
            end

            evelen   = blen/dat.fsample;      % Length of burst is seconds  

            % Arrange data
            bdat(i).begsam = startb;
            bdat(i).endsam = endb;
            bdat(i).evelen = evelen;
        end
        
        rho1mat(ss,:,f) = rho1;
        rho2mat(ss,:,f) = rho2;
        rho3mat(ss,:,f) = rho3;
        rho4mat(ss,:,f) = rho4;
        
%         if any(strfind(fname,'RsEc1'))
%             ses = 'ses1'
%         elseif any(strfind(fname,'RsEc2'))
%             ses = 'ses2'
%         end
    
        % Combine summary values in struct as save
        subvals{f}.steps = steps;
        subvals{f}.bdat = bdat;
        subvals{f}.rho = rho1';
    end
%     save(fullfile(sub_dir,'subvals2.mat'),'subvals')
    clear bdat subvals
end
disp('done')

%% Compare methods

plot(steps,mean(rho1mat)); hold on
plot(steps,mean(rho2mat),'r'); 
plot(steps,mean(rho3mat),'k'); 
plot(steps,mean(rho4mat),'m'); hold off

plot(steps,rho1mat); 

rhoavg = mean(rho3mat)
rhosd = std(rho1mat)

plot(steps,rhoavg); hold on
plot(steps,rhoavg+rhosd, 'b--')
plot(steps,rhoavg-rhosd, 'b--')

[maxval,idx] = max(rhoavg);
cutoff = steps(idx);

%% Misc
figure; plot(steps,rho1,'-or'); hold on
plot(steps,rho2,'-xb')

[~,idx] = max(rho1);
burst = data > med+sd*steps(idx);   
burstdat = data(burst);
    
% Histograms and misc plots
hist(data,100)
tim = 1:length(data);

figure();
plot(tim,data); hold on
plot(tim(burst),burstdat, 'r');
mline = refline([0 med]);
mline.Color = 'r';
for i = 1:5
    sdline = refline([0 med+sd*i]);
    sdline.Color = 'g';
end

figure;
hist1 = histogram(data, 50, 'EdgeColor','b', 'FaceColor', 'b');
hold on;
hist2 = histogram(data, 50, 'EdgeColor','r', 'FaceColor', 'r', 'FaceAlpha',0.2);


figure;
plot(tim,data); hold on
plot(tim(burst),burstdat, 'r');
mline = refline([0 prctile(data, 97)]);

%% Plots


cfg = [];
cfg.lpfilter    = 'yes';
cfg.lpfreq      = 35;
cfg.hpfilter    = 'yes';
cfg.hpfreq      = 10;
data_bp = ft_preprocessing(cfg, raw_crop);

cfg = [];
cfg.method = 'svd';                         % NB! Must be SVD!
data_bpcmb = ft_combineplanar(cfg,data_bp);
  
cfg = [];
cfg.rectify = 'yes';
data_rect = ft_preprocessing(cfg, data_bpcmb);

cfg = [];
cfg.layout          = 'neuromag306cmb';
cfg.baseline        = 'no';
ft_multiplotER(cfg,data_bpcmb);


% Select data
cfg = [];
cfg.method = 'svd';                         % NB! Must be SVD!
raw_crop_cmb = ft_combineplanar(cfg,raw_crop);
cfg = [];
cfg.channel = {'MEG0232+0233', 'MEG0442+0443', 'MEG1622+1623', 'MEG1812+1813'};
cfg.avgoverchan = 'yes';
raw_datachan = ft_selectdata(cfg, raw_crop_cmb);




