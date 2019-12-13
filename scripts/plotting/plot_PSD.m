% PD-BB: Plot PSD from running window and FFT.
% set paths
addpath /home/mikkel/PD_motor/global_scripts
[dirs, ~, ~] = PD_proj_setup('betaburst');
addpath /home/mikkel/matlab/export_fig/

%% Load data
disp('Loading...');
load('/home/mikkel/PD_motor/rest_ec/groupanalysis/PSD_GA.mat') % GA struct
disp('done');

%% plot settings
plot_lnwdt          = 1.5;
axis_lnwdt          = 2;
def_fontsize        = 12;
label_fontsize      = 14; 
xrange              = [0 45];
% yrange              = [0 0.01];       % native range
yrange              = [-8.5 -4.5];      % log range

%% Plot
close all
fig = figure; hold on
set(fig,'Position', [0 0 600 500], 'color','w');

% patch([12 30 30 12],[min(yrange) min(yrange) max(yrange) max(yrange)], ...
%     'b','FaceAlpha',0.075,'EdgeColor','none')
patch([0 12 12 0],[min(yrange) min(yrange) max(yrange) max(yrange)], ...
    'b','FaceAlpha',0.05,'EdgeColor','none')
patch([30 max(xrange) max(xrange) 30],[min(yrange) min(yrange) max(yrange) max(yrange)], ...
    'b','FaceAlpha',0.05,'EdgeColor','none')

p1 = plot(GA.ptns1.freq, mean(log(squeeze(GA.ptns1.powspctrm))),'b-', ... 
    'LineWidth', plot_lnwdt);
p2 = plot(GA.ptns2.freq, mean(log(squeeze(GA.ptns2.powspctrm))),'b:', ... 
    'LineWidth', plot_lnwdt*2);
p3 = plot(GA.ctrl1.freq, mean(log(squeeze(GA.ctrl1.powspctrm))),'r-', ... 
    'LineWidth', plot_lnwdt);
p4 = plot(GA.ctrl2.freq, mean(log(squeeze(GA.ctrl2.powspctrm))),'r:', ... 
    'LineWidth', plot_lnwdt*2); 

legend([p1,p2,p3,p4],{'PD 1/OFF','PD 2/ON','Ctrl 1','Ctrl 2'}, ...
    'Location','NorthEast'); legend BOXOFF

set(gca, 'LineWidth', axis_lnwdt,'fontweight','bold','fontsize',def_fontsize, ...
        'XTick', [0:5:45], 'xticklabel',{0:5:45});
%         'YTick', [-0.5:0.5:1], 'yticklabel',{-0.5:0.5:1});
xlabel('Frequency (Hz)','fontsize', label_fontsize);
ylabel('log(power)','fontsize', label_fontsize)
xlim(xrange); ylim(yrange);
hold off

%% Export
export_fig(fullfile(dirs.figures,'PSD.png'), ...
    '-r600', '-p0.05', '-CMYK')

%END