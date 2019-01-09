%%%%% Frequency analysis of beta-burst data.
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
overwrite_old_files = 1;   % Overwirte old files 0=false or 1=true

%% Loop through subjects
for ii = 1:length(subs)
    subID = subs{ii};
    sub_dir = [dirs.megDir,'/',subID];
    files = dir(sub_dir);
    files = {files.name};
    file_idx = find(~cellfun(@isempty,strfind(files,'-ft-raw.mat'))); % Name of imported cropped file
    infiles = files(file_idx);
        
    cd(sub_dir);
    
    fprintf('Now processing subject = %s.', subID);
    disp(['Found ',num2str(length(infiles)),' files for sub ', subID,''])


    for kk = 1:length(infiles)
        data_file = infiles{kk};
        outfname_wav = [data_file(1:10),'-wav.mat'];
        outfname_hil = [data_file(1:10),'-hil.mat'];
        outfname_bp = [data_file(1:10),'-bpcmb.mat'];

        if exist(fullfile(sub_dir,outfname_wav),'file') && exist(fullfile(sub_dir,outfname_hil),'file') && ~overwrite_old_files
            fprintf('Files %s and %s already exist. Will not overwrite!', outfname_wav,outfname_hil)
            continue
        end
        
        fprintf('Loading %s ...\n', data_file)
        load(data_file)
        
%         % Plot
%         cfg = [];
%         cfg.layout          = 'neuromag306all';
%         cfg.baseline        = 'no';
%         ft_multiplotER(cfg,raw_crop);        
        
%         if ~exist(fullfile(sub_dir,outfname_wav)) || overwrite_old_files
%             disp('Running wavelet...')
%             % WAVELET ANALYSIS
%             cfg = [];
%             cfg.method  = 'wavelet';
%             cfg.output  = 'pow';
%             cfg.foi     = 10:1:35;                              % frequency band (freqs in steps of one)
%             cfg.width   = 7;                                    % Wavelet cycles
%             cfg.toi     = min(raw_crop.time{:}):0.001:max(raw_crop.time{:});          % time bins. Steps of 100 ms
%             cfg.pad     ='nextpow2';
%             tic;
%             data_wavelet = ft_freqanalysis(cfg, raw_crop);
%             toc;                                                % ~43 s (dt=0.005, f [15,30]
%                                                                 % 123 s (dt=0.001, f [15,30]
%             cfg = [];
%             cfg.method = 'sum';
%             data_waveletcmb = ft_combineplanar([],data_wavelet);
%         
%             %Save wavelet
%             fprintf('Saving %s ...', outfname_wav)
%             save(fullfile(sub_dir,outfname_wav),'data_waveletcmb')
%             disp('done')
%             
%             clear data_wavelet data_waveletcmb
%         end
        
%         % Freq baseline
%         cfg = [];
%         cfg.baseline        = [-inf inf];
%         cfg.baselinetype    = 'normchange';
%         data_waveletbs = ft_freqbaseline(cfg,data_waveletcmb);
%         
%         % Plot
%         figure;
%         cfg = [];
%         cfg.layout          = 'neuromag306cmb';
%         cfg.baseline        = 'no';
%         ft_multiplotTFR(cfg,data_waveletbs);

        if ~exist(fullfile(sub_dir,outfname_hil)) || overwrite_old_files
            disp('Running Hilbert...')
            % HILBERT ANALYSIS
            cfg = [];
            cfg.lpfilter    = 'yes';
            cfg.lpfreq      = 30;
            cfg.hpfilter    = 'yes';
            cfg.hpfreq      = 10;
            data_bp = ft_preprocessing(cfg, raw_crop);

            cfg = [];
            cfg.method = 'svd';                         % NB! Must be SVD!
            data_bpcmb = ft_combineplanar(cfg,data_bp);

            cfg = [];
            cfg.hilbert = 'abs';
            data_hilbt = ft_preprocessing(cfg,data_bpcmb);
            
            % Save Hilbert and band-pass data
            fprintf('Saving %s ...', outfname_bp)
            save(fullfile(sub_dir,outfname_bp),'data_bpcmb')
            fprintf('Saving %s ...', outfname_hil)
            save(fullfile(sub_dir,outfname_hil),'data_hilbt')
            disp('done')
            
            clear data_bp data_bpcmb data_hilbt
        end
        
%         % Plot
%         cfg = [];
%         cfg.layout          = 'neuromag306cmb';
%         cfg.baseline        = 'no';
%         ft_multiplotER(cfg,data_hilbt);

    end
end