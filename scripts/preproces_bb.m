%%%%% Import data for beta-burst part.
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

%% Import
% Loop through subjects
for ii = 1:length(subs)
    subID = subs{ii};
    sub_dir = [dirs.megDir,'/',subID];
    files = dir(sub_dir);
    files = {files.name};
    fif_idx = find(~cellfun(@isempty,strfind(files,'ica_raw.fif')));
    infiles = files(fif_idx);
        
    cd(sub_dir);
    
    fprintf('Now processing subject = %s.', subID);
    disp(['Found ',num2str(length(infiles)),' FIF files for sub ', subID,''])


    for kk = 1:length(infiles)
        data_file = infiles{kk};
        
        outfname = fullfile(sub_dir,[data_file(1:10),'-ft-raw.mat']);

        if ~exist(outfname, 'file') || overwrite_old_files
            
            % Find start and stop
            eve = ft_read_event(data_file);
            eve = eve(strcmpi('STI101',{eve.type}));
            
            if strcmp(data_file(1:10),'0327_RsEc1') % Subject with missing triggers.
                begsam = 30000;
                endsam = 180102+begsam;
            elseif strcmp(data_file(1:10),'0333_RsEc2') % Subject with missing triggers.
                begsam = 20000;
                endsam = 180102+begsam;
            else
                begsam = eve([eve.value]==1).sample;
                endsam = eve([eve.value]==64).sample;
            end  
            
            % Add buffer
            begsam = begsam-1000; % Add one second before start trigger
            endsam = endsam+1000; % Add one secind after end trigger
            
            % Find grad channels
            hdr = ft_read_header(data_file);
%             chidx = find(strcmp('megplanar', hdr.chantype));
            begtim = begsam/hdr.Fs;
            endtim = endsam/hdr.Fs;
                       
%             raw = ft_read_data(data_file,'begsample',begsam,'endsample',endsam,'chanindx',chidx)
            
            % Read raw data
            cfg = [];
            cfg.dataset     = data_file;
            cfg.channel     = 'MEGGRAD';
            cfg.demean      = 'no';
            cfg.continuous  = 'yes';
   
            raw = ft_preprocessing(cfg);
                                    
            % Inspect            
%             ft_databrowser([],raw)
            
            % Crop data based on start/stop triggers
            cfg = [];
            cfg.latency = [begtim endtim];
            raw_crop = ft_selectdata(cfg, raw);
            
            % Inspect            
%             ft_databrowser([],raw_crop)
            
            % Save
            save(outfname,'raw_crop','-v7.3')

        end
    end
end
