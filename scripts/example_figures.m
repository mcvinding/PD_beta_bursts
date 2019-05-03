%% Example plots
addpath /home/mikkel/matlab/export_fig/
addpath /home/mikkel/PD_motor/global_scripts
[dirs, ~, ~] = PD_proj_setup('betaburst');

subs = find_subs(dirs.megDir);                                %Find subjects in folder
cd(dirs.megDir);

%% "raw" ROI beta band with envelope
subID = subs{1};
sub_dir = fullfile(dirs.megDir,subID);
files = dir(sub_dir);
files = {files.name};    
fileH_idx = find(~cellfun(@isempty,strfind(files,'-hilbt.mat'))); % Name of imported cropped file
infiles = sort(files(fileH_idx));
load(fullfile(sub_dir,infiles{1})) % hilb_xx
fileR_idx = find(~cellfun(@isempty,strfind(files,'-rawft.mat'))); % Name of imported cropped file
infiles = sort(files(fileR_idx));
load(fullfile(sub_dir,infiles{1})) % rwft_xx

tim = 0:0.001:5;

% Plot
fig = figure; hold on
set(fig,'Position', [0 0 800 300], 'color','w');
plot(tim, rwft_lh(9000:14000), 'b', 'lineWidth',1)
plot(tim, hilb_lh(9000:14000), 'r', 'lineWidth',1.5)
set(gca, 'LineWidth', 1,'fontweight','bold','fontsize',14, ...
    'YTick', [], 'yticklabel','', 'YColor','none', ...
    'XTick', [0,1,2,3,4,5], 'xticklabel',{0,1,2,3,4,5});
xlabel('Time (s)','fontsize',16);
ylabel('','fontsize',16)

export_fig(fullfile(dirs.figures,'hilbt.png'), '-r500', '-p0.05', '-CMYK')

%% "raw" ROI beta
snip = rwft_lh(134000:139000);
fig = figure; hold on
set(fig,'Position', [0 0 800 400], 'color','w');
plot(tim,snip, 'b', 'lineWidth',1)
xlim([min(tim),max(tim)])
ylim([min(snip)+0.2*min(snip), max(snip)+0.5*max(snip)]);
set(gca, 'LineWidth', 1,'fontweight','bold','fontsize',14, ...
    'YTick', [], 'yticklabel','', 'YColor','none', ...
    'XTick', [0,1,2,3,4,5], 'xticklabel',{0,1,2,3,4,5});
xlabel('Time (s)','fontsize',16);

export_fig(fullfile(dirs.figures,'beta5s.png'), '-r500', '-p0.05', '-CMYK')


plot(rwft_lh, 'b', 'lineWidth',1)
