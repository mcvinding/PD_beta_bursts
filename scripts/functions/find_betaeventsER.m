function [output, rhomat] = find_betaeventsER(cfg, data)
% Output is a data structure of N length, containing A) number of beta
% events, B) start sample of beta events in data, C) end sample of beta
% event in data, C) length of events. rhomat is a vector of length N
% contining correaltion between number of events and amplitude/power in
% epochs. N is the number of steps used to determine threshold. When
% thereshold is determined then a single scalar is used. [NB revwrite!]
%
% USE: [output, rhomat] = find_betaevents(cfg, data)
% INPUT:
% cfg.steps         = [Nx1] Steps to find correlation and threshold (required)
% cfg.cutofftype    = ['sd'/'med'] Should cutoff be determined relative to
%                     x*standard deviations away from median (Little et al, 2018)
%                     or x*medians away from median (Shin et al. 2016).
%                     (default='med').
% cfg.corrtype      = ['amp'/'pow'] What to correalte: amplitude of epochs
%                     vs. number of events or power (amp^2) of epoch vs.
%                     number of beta events (Little et al, 2018).
%                     (default='amp').
% cfg.halfmax       = ['yes'/'no'/'mixed'] redefine timing of events as half-
%                     maximum of the peak value. Overlapping peaks are counted 
%                     as one event. Otherwise, event length is defined by
%                     threshold crossing. 'Mixed' replace half-maximum
%                     values above threshold with the value of the
%                     threshold crossing (default='no').
% cfg.length        = [num] length of epoch window in seconds. Passed to
%                     FT_REDEFINETRIAL (default=3).
% cfg.overlap       = [num] overlap between epochs assed to FT_REDEFINETRIAL
%                     (default=0, i.e. no overlap)
% cfg.makeplot      = ['yes'/'no'] plot the time series with cutoff, peaks,
%                     and events marked (default='no').
% OUTPUT:
% ...

% ft_checkdata.
dat = ft_checkdata(data, 'datatype','timelock');

% opts
cfg = ft_checkconfig(cfg, 'required', 'steps');
steps = cfg.steps;

% cfg.getcutoff   = ft_getopt(cfg, 'getcutoff', 'no');
cfg.cutofftype  = ft_getopt(cfg, 'cutofftype', 'med');
cfg.corrtype    = ft_getopt(cfg, 'corrtype', 'amp');
cfg.makeplot    = ft_getopt(cfg, 'makeplot', 'no');
cfg.halfmax     = ft_getopt(cfg, 'halfmax', 'no');

% Check variables
if ~any(strcmp({'sd','med'},cfg.cutofftype))
    error('Unknown cutoff method \"%s\".',cfg.cutofftype)
end

if ~any(strcmp({'amp','pow'}, cfg.corrtype))
    error('Unknown cutoff method \"%s\".',cfg.cutofftype)
end

trl = dat.sampleinfo;
% trl = trl-trl(1,1)+1; %Corrent for non-zero sample offset

% Get mean epoch power and amplitude (per epoch)
epoamp = nan(length(trl),1);
epopow = nan(length(trl),1);
for n = 1:length(trl)
    epoamp(n) = mean(data.trial{n});
    epopow(n) = mean(data.trial{n}.^2);
end

% Cutoffs
% Initiate values
pkmat    = zeros(length(trl),1);
n_events = zeros(1,length(steps));
rhomat   = nan(length(steps),1);
cutoff   = zeros(1,length(steps));
bdat     = struct();

% Find values
tempdat = squeeze(dat.trial);
med = median(tempdat(:));         % Median all across trials
sd = std(tempdat(:));             % SD acriss all trials
fprintf('Median of all trials: %.3f. sd: %.3f.\n', med, sd)

for ii = 1:length(steps)
    if strcmp(cfg.cutofftype, 'sd')
        cutoff(ii) = med+sd*steps(ii);
        burst = tempdat >= cutoff(ii);   
    elseif strcmp(cfg.cutofftype, 'med')
        cutoff(ii) = med+med*steps(ii);
        burst = tempdat >= cutoff(ii);
    end
    burst = squeeze(burst);        % Remove empty dimension
    
    % Get correlation [MOVE TO END]
    if strcmp('amp', cfg.corrtype)
        rhomat(ii) = corr(n_events, epoamp);
    elseif strcmp('pow', cfg.corrtype)
        rhomat(ii) = corr(n_events, epopow);
    end
% ------------------------------------------------------------------- %
    % Get summaries
    burst(:,end) = 0;
    dburst = diff([zeros(size(burst,1),1) burst], [], 2);
    n_events = sum(dburst==1,2);    % Events per trial
    cburst = burst;
    
    for kk = 1:size(tempdata,1)
        maxarray = zeros(n_events(kk),1);
        maxidx   = zeros(n_events(kk),1);
        trldat = tempdat(kk,:);
        trldbs = dburst(kk,:);
        startb  = find(trldbs==1);      % Start of burst
        endb    = find(trldbs==-1);     % End of burst´
        
        for n = 1:n_events(kk)
            [maxarray(n), maxidx(n)] = max(trldat(startb(n):endb(n)));
            maxidx(n) = maxidx(n)+startb(n)-1;
        end
    
%     tmp = diff(burst,[],2);
%     pkmat = sum(tmp==1,2);

        % start-stop based on half-max width
        if strcmp(cfg.halfmax, 'yes') || strcmp(cfg.halfmax, 'mixed')
            begsam = zeros(1,n_events(kk));
            endsam = zeros(1,n_events(kk));
            hlfmx = maxarray/2;

            if  strcmp(cfg.halfmax, 'mixed')
                hlfmx(hlfmx>cutoff(ii)) = cutoff(ii);
            end

            for n = 1:n_events(kk)
                %start
                idx = maxidx(n);
                xval = trldat(idx);
                while xval > hlfmx(n)
                    if idx == 1
                        break
                    end
                    idx = idx-1;
                    xval = trldat(idx);
                end
                begsam(n) = idx;
                % end
                idx = maxidx(n);
                xval = trldat(idx);
                while xval > hlfmx(n)
                    if idx == length(dat)
                        break
                    end
                    idx = idx+1;
                    xval = trldat(idx);
                end
                endsam(n) = idx;
            end

            trlb = zeros(1,length(trldat));
            for n = 1:n_events(kk)
                trlb(begsam(n):endsam(n)) = 1;
            end
            
            cburst(kk,:) = trlb;   % Corrected matrix
        end
        
        % Get summaries (again)
        if trlb(end) == 1; trlb(end) = 0; end
        dcburst = diff([0 trlb']);
        neve = sum(dcburst==1);
        startb  = find(dcburst==1);      % Start of burst
        endb    = find(dcburst==-1);     % End of burst´

        maxarray = zeros(neve,1);
        maxidx   = zeros(neve,1);
        for n = 1:neve
            [maxarray(n), maxidx(n)] = max(dat(startb(n):endb(n)));
        end
        maxidx = maxidx+startb'-1;
        
        n_events(ii) = neve;
    end
    
    evemark = nan(length(dat.trial),1);
    for n = 1:n_events(ii)
        evemark(startb(n):endb(n)) = dat(startb(n):endb(n));
    end

    
    if strcmp(cfg.makeplot,'yes')
%         xidx = 1:length(dat)/data.fsample;
        if length(steps) >= 8 
            dim = ceil(length(steps)/8);
            subplot(dim,8,ii); hold on
        else
            figure; hold on
        end
        plot(dat);
        plot(repmat(cutoff(ii),length(dat),1),'r--');
        plot(maxidx,maxarray, 'ko');
        plot(1:length(dat),evemark, 'linewidth',2);
        xlim([0 length(dat)]);
        title(ii);
    end
    
    % Length of burst
    blen    = endb-startb;          % Length of burst
    evelen   = blen/data.fsample;      % Length of burst is seconds  
    if any(blen<0)
        error('negative length of event')
    end

    % Arrange data
    bdat(ii).event  = [startb; endb]';
    bdat(ii).evelen = evelen';
    bdat(ii).maxpk  = maxarray;
    bdat(ii).maxidx = maxidx;
end

% Make output
output.steps    = steps;
output.n_events = n_events;
output.cutoff   = cutoff;
output.bdat     = bdat;

% End