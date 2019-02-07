function [cutoff] = find_threshold(rhomat, steps, plots)
% USE: [cutoff] = find_threshold(rhomat, steps, plots);
% Args:
% rhomat = N subjects x M steps (x J repetitions). Repetitions are averaged
%          per subject before averagining across all subjects.
% steps  = range of steps used for testing thresholds [Mx1].
% plots  = [logical] Plot curves steps x rho (default=0)

if nargin < 3
    plots = 0;
end

dim = size(rhomat);
if length(dim) == 3
    fprintf('Assuming repeated measures across subjects:\n   Averaging %i repetitions per unit\n', dim(3))
    rhomat = nanmean(rhomat,3);
elseif length(dim) > 3
    error('Cannot handle input with %i dimensions.', length(dim))
end

if dim(2) ~= length(steps)
    error('steps and size of rhomat does not fit');
end

% Do stuff
rhoavg = nanmean(rhomat);
rhosd = nanstd(rhomat);

[maxval,idx] = max(rhoavg);
cutoff = steps(idx);
fprintf('Max correlation value is RHO=%.3f at cutoff=%.2f.\n', maxval, cutoff)

if plots==1
    figure; hold on
    plot(steps,rhoavg, 'k', 'linewidth', 2);
    plot(steps,rhoavg+rhosd, 'r--')
    plot(steps,rhoavg-rhosd, 'r--')
    line([cutoff,cutoff],[0 1])
    hold off
end

%end

