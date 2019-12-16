% Make nice figure of threhold cut-off (Figure 1B)
addpath /home/mikkel/PD_motor/global_scripts
[dirs, ~, ~] = PD_proj_setup('betaburst');
addpath /home/mikkel/matlab/export_fig/

% Load data
load('/home/mikkel/PD_motor/rest_ec/groupanalysis/rhomats.mat')
disp('done')

%% plot settings
axis_lnwdt          = 2;
def_fontsize        = 12;
label_fontsize      = 14; 

%% Plot 
close all
find_threshold(rho_mdamp, steps, 1);        % This function has a build in plot option
fig = gcf; hold on
set(fig,'Position', [0 0 500 400], 'color','w');

set(gca, 'LineWidth', axis_lnwdt,'fontweight','bold','fontsize',def_fontsize, ...
        'XTick', [0:1:5], 'xticklabel',{0:1:5}, ...
         'YTick', [0:0.2:1], 'yticklabel',{0:0.2:1});
xlabel('Threshold (median+median*x)','fontsize', label_fontsize);
ylabel('Perason''s \rho', 'fontsize', label_fontsize)
title('')

%% Export
export_fig(fullfile(dirs.figures, 'threshold_figure.png'), ...
    '-r500', '-p0.05', '-CMYK')

%END