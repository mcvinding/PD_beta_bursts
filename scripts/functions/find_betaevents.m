function [output, rhomat] = find_betaevents(cfg, data)
% USE: [output, rhomat] = find_betaevents(cfg, data)
% Output is a data structure of N length, containing A) number of beta
% events, B) start sample of beta events in data, C) end sample of beta
% event in data, C) length of events. rhomat is a vector of length N
% contining correaltion between number of events and amplitude/power in
% epochs. N is the number of steps used to determine threshold. When
% thereshold is determined then a single scalar is used.
% INPUT:
% cfg.steps         = [Nx1] Steps to find correlation and threshold (required)
% cfg.cutofftype    = ['sd'/'med'] Should cutoff be determined relative to
%                     x*standard deviations away from median (Little et al, 2018)
%                     or x*medians away from median (Shin et al. 2016).
%                     Default = 'med'.
% cfg.corrtype      = ['amp'/'pow'] What to correalte: amplitude of epochs
%                     vs. number of events or power (amp^2) of epoch vs.
%                     number of beta events (Little et al, 2018).
%                     Default='amp'.
% cfg.length        = [num] length of epoch window in seconds (default=3)
% cfg.overlap       = [num] overlap between epochs (default=0, i.e. no overlap)
%
% OUTPUT:
% ...

% TO DO:
% * Need a ft_checkdata section.

% opts
cfg = ft_checkconfig(cfg, 'required', 'steps');
steps = cfg.steps;

cfg.cutofftype  = ft_getopt(cfg, 'cutofftype', 'med');
cfg.corrtype    = ft_getopt(cfg, 'corrtype', 'amp');

% Check variables
if ~any(strcmp({'sd','med'},cfg.cutofftype))
    error('Unknown cutoff method \"%s\".',cfg.cutofftype)
end

if ~any(strcmp({'amp','pow'}, cfg.corrtype))
    error('Unknown cutoff method \"%s\".',cfg.cutofftype)
end

% Make pseudo-tirals.
cfg.length      = ft_getopt(cfg, 'length', 3);
cfg.overlap     = ft_getopt(cfg, 'overlap', 0);
epo = ft_redefinetrial(cfg, data);

trl = epo.sampleinfo;
% trl = trl-trl(1,1)+1; %Corrent for non-zero sample offset

% Get epoch power and amplitude
epoamp = nan(length(epo.trial),1);
epopow = nan(length(epo.trial),1);
for k = 1:length(trl)
    epoamp(k) = mean(epo.trial{k});
    epopow(k) = mean(epo.trial{k}.^2);
end

% Cutoffs
% Initiate values
pkmat       = zeros(length(trl),1);
n_events    = zeros(1,length(steps));
rhomat      = nan(length(steps),1);
bdat = struct();

% Find values
dat = data.trial{:};
med = median(dat);
sd = std(dat);
fprintf('Median of time-series: %.3f. sd: %.3f.\n', med, sd)

for ii = 1:length(steps)
    if strcmp(cfg.cutofftype, 'sd')
        cutoff = med+sd*steps(ii);
        burst = dat >= cutoff;   
    elseif strcmp(cfg.cutofftype, 'med')
        cutoff = med+med*steps(ii);
        burst = dat >= cutoff;
    end

    for k = 1:length(trl)
        tmp = burst(trl(k,1):trl(k,2));
        pkmat(k) = sum((diff(tmp)==1));

    end
    
    if strcmp('amp', cfg.corrtype)
        rhomat(ii) = corr(pkmat, epoamp);
    elseif strcmp('pow', cfg.corrtype)
        rhomat(ii) = corr(pkmat, epopow);
    end
    
    % Get summaries
    dburst = diff([0 burst 0]);
    n_events(ii) = sum(dburst==1);

    startb  = find(dburst==1);      % Start of burst
    endb    = find(dburst==-1);     % End of burst´
    blen    = endb-startb;          % Length of burst

    if any(blen<0)
        error('negative length of event')
    end

    evelen   = blen/data.fsample;      % Length of burst is seconds  

    % Arrange data
    bdat(ii).begsam = startb;
    bdat(ii).endsam = endb;
    bdat(ii).evelen = evelen;
    bdat(ii).cutoff = cutoff;
end

% Make output
output.steps    = steps;
output.n_events = n_events;
output.bdat     = bdat;

% End