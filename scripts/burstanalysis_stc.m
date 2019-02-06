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

steps = 0:0.1:5;

%% Find peaks: testing multiple methods based on the litterature
% * Correalte epoch pow with N-peaks in epoch base on whole data peak
%   threshold: both sd and median based cutoffs
% * Correalte epoch pow with N-peaks in epoch base on epoch based peak
%   threshold: both sd and median based cutoffs

rho1mat = nan(length(subs),length(steps),2);
rho2mat = nan(length(subs),length(steps),2);
rho3mat = nan(length(subs),length(steps),2);
rho4mat = nan(length(subs),length(steps),2);

for ss = 1:length(subs)
    subID = subs{ss};
    sub_dir = fullfile(dirs.megDir,subID);
    files = dir(sub_dir);
    files = {files.name};
    
    outfname = fullfile(sub_dir,'subvals2.mat');
    if exist(outfname,'file') && ~overwrite
        continue
    end

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
        cfg.length  = 3;
        cfg.overlap = 0;
        epo = ft_redefinetrial(cfg,lhdata);

        trl = epo.sampleinfo;
        trl = trl-trl(1,1)+1; %Corrent for non-zero sample offset

        % Get epoch power and amplitude
        epoamp = nan(length(epo.trial),1);
        epopow = nan(length(epo.trial),1);
        for k = 1:length(trl)
            epoamp(k) = mean(epo.trial{k});
            epopow(k) = mean(epo.trial{k}.^2);
        end

        % Cutoffs
        data = hilb_lh;
        pkmat1 = zeros(length(trl),1);
        pkmat2 = zeros(length(trl),1);

        n_events = zeros(length(trl),1);

        med = median(data);
        sd = std(data);

        rho1 = nan(length(steps),1);
        rho2 = nan(length(steps),1);
        rho3 = nan(length(steps),1);
        rho4 = nan(length(steps),1);
    
        bdat = struct();
        for ii = 1:length(steps)
            burst = data >= med+sd*steps(ii);   
            burst2 = data >= med+med*steps(ii);

            for k = 1:length(trl)
                % Little et al (2018)
                tmp = burst(trl(k,1):trl(k,2));
                pkmat1(k) = sum((diff(tmp)==1));

                % units of median (cf. Chin et al)
                tmp2 = burst2(trl(k,1):trl(k,2));
                pkmat2(k) = sum((diff(tmp2)==1));

            end
            rho1(ii) = corr(pkmat1, epoamp);
            rho2(ii) = corr(pkmat2, epoamp);
            rho3(ii) = corr(pkmat1, epopow);
            rho4(ii) = corr(pkmat2, epopow);
        
            % Get summaries
            dburst = diff([0 burst 0]);
            n_events(ii) = sum(dburst==1);

            startb  = find(dburst==1);      % Start of burst
            endb    = find(dburst==-1);     % End of burst´
            blen    = endb-startb;          % Length of burst

            if any(blen<0)
                error('negative length')
            end

            evelen   = blen/epo.fsample;      % Length of burst is seconds  

            % Arrange data
            bdat(ii).begsam = startb;
            bdat(ii).endsam = endb;
            bdat(ii).evelen = evelen;
        end
        
        rho1mat(ss,:,f) = rho1;
        rho2mat(ss,:,f) = rho2;
        rho3mat(ss,:,f) = rho3;
        rho4mat(ss,:,f) = rho4;
    
        % Combine summary values in struct as save
        subvals{f}.steps = steps;
        subvals{f}.bdat = bdat;
        subvals{f}.rho = rho1';
    end
%     save(fullfile(sub_dir,'subvals2.mat'),'subvals')
    clear bdat subvals
end
save('/home/mikkel/PD_motor/rest_ec/groupanalysis/rhomats.mat','rho1mat','rho2mat','rho3mat','rho4mat')
disp('done')

%% Compare methods
% Load data
load('/home/mikkel/PD_motor/rest_ec/groupanalysis/rhomats.mat')
disp('done')

% Rho1: sd based threshold ~ mean amplitude of epoch
% Rho2: median based threshold ~ mean amplitude of epoch
% Rho3: sd based threshold ~ mean power (amp^2) of epoch
% Rho4: median based threshold ~ mean power (amp^2) of epoch

figure; plot(steps,nanmean(rho1mat,3));
figure; plot(steps,nanmean(rho2mat,3)); 
figure; plot(steps,nanmean(rho3mat,3)); 
figure; plot(steps,nanmean(rho4mat,3)); 


rhoavg = nanmean(nanmean(rho1mat,3));
rhosd = nanstd(nanmean(rho1mat,3));

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
for ii = 1:5
    sdline = refline([0 med+sd*ii]);
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




