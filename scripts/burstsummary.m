% Get various summaries of beta burst (explorative)
clear all
close all

addpath /home/mikkel/PD_motor/global_scripts
[dirs, subjs, ~] = PD_proj_setup('betaburst');

cd(dirs.megDir);
subs = find_subs(dirs.megDir);                                %Find subjects in folder

[PD_subs, PDidx] = intersect(subs,subjs.PD);
[ctrl_subs, ctrlidx] = intersect(subs, subjs.ctrl);

%% N events
nevent1 = nan(length(subs), 1);
nevent2 = nan(length(subs), 1);
b1_1 = nan(length(subs), 1);
b2_1 = nan(length(subs), 1);
b3_1 = nan(length(subs), 1);
b1_2 = nan(length(subs), 1);
b2_2 = nan(length(subs), 1);
b3_2 = nan(length(subs), 1);

for ii = 1:length(subs)
    load(fullfile(dirs.megDir,subs{ii},'subvals.mat'))
    
    nevent1(ii) = subvals{1}.n_events;
    nevent2(ii) = subvals{2}.n_events;
    
    b1_1(ii) = sum(subvals{1}.bdat.maxidx < 60000);
    b2_1(ii) = sum(subvals{1}.bdat.maxidx > 60000 & subvals{1}.bdat.maxidx < 120000);
    b3_1(ii) = sum(subvals{1}.bdat.maxidx > 120000);
    b1_2(ii) = sum(subvals{2}.bdat.maxidx < 60000);
    b2_2(ii) = sum(subvals{2}.bdat.maxidx > 60000 & subvals{2}.bdat.maxidx < 120000);
    b3_2(ii) = sum(subvals{2}.bdat.maxidx > 120000);
end

PDn1 = nevent1(PDidx);
ctrln1 = nevent1(ctrlidx);
PDn2 = nevent2(PDidx);
ctrln2 = nevent2(ctrlidx);

PDn1_1 = b1_1(PDidx);
PDn2_1 = b2_1(PDidx);
PDn3_1 = b3_1(PDidx);
PDn1_2 = b1_2(PDidx);
PDn2_2 = b2_2(PDidx);
PDn3_2 = b3_2(PDidx);

ctrln1_1 = b1_1(ctrlidx);
ctrln2_1 = b2_1(ctrlidx);
ctrln3_1 = b3_1(ctrlidx);
ctrln1_2 = b1_2(ctrlidx);
ctrln2_2 = b2_2(ctrlidx);
ctrln3_2 = b3_2(ctrlidx);

% PDnavg1 = nanmean(PDn1);
% ctrlnavg1 = nanmean(ctrln1);
% PDnavg2 = nanmean(PDn2);
% ctrlnavg2 = nanmean(ctrln2);
% 
% PDnsd1 = nanstd(PDn1);
% ctrlnsd1 = nanstd(ctrln1);
% PDnsd2 = nanstd(PDn2);
% ctrlnsd2 = nanstd(ctrln2);

% [~, pt1,~,t1] = ttest2(PDn1,ctrln1);
% [~, pt2,~,t2] = ttest2(PDn2,ctrln2);
% [~, ptPt,~,tPt] = ttest(PDn1,PDn2);
% [~, ptCt,~,tCt] = ttest(ctrln1,ctrln2);

figure;
subplot(1,2,1); histogram(PDn1,10); hold on
subplot(1,2,1); histogram(ctrln1,10); hold off
subplot(1,2,2); histogram(PDn2,10); hold on
subplot(1,2,2); histogram(ctrln2,10); hold off

% mns = [PDnavg1, ctrlnavg1; PDnavg2, ctrlnavg2];
% sds = [PDnsd1, ctrlnsd1; PDnsd2, ctrlnsd2];

save('/home/mikkel/PD_motor/rest_ec/groupanalysis/nevent.mat', ...
    'PDn1', 'ctrln1', 'PDn2', 'ctrln2','PD_subs','ctrl_subs');

save('/home/mikkel/PD_motor/rest_ec/groupanalysis/nevent_min.mat', ...
    'PDn1_1', 'PDn2_1', 'PDn3_1', 'PDn1_2', 'PDn2_2', 'PDn3_2', ...
    'ctrln1_1', 'ctrln2_1', 'ctrln3_1', 'ctrln1_2', 'ctrln2_2', 'ctrln3_2', ...
    'PD_subs','ctrl_subs');

%% Event duration
lenmean1 = zeros(length(subs), 1);   % subjects x number of steps
lenmedn1 = zeros(length(subs), 1);   % subjects x number of steps
lenmean2 = zeros(length(subs), 1);   % subjects x number of steps
lenmedn2 = zeros(length(subs), 1);   % subjects x number of steps

len1 = [];
len2 = [];
sub1 = [];
sub2 = [];

for ii = 1:length(subs)
    load(fullfile(dirs.megDir,subs{ii},'subvals.mat'))
    
    len1 = [len1; subvals{1}.bdat.evelen];
    len2 = [len2; subvals{2}.bdat.evelen];
    
    sub1 = [sub1; repmat(subs{ii},length(subvals{1}.bdat.evelen),1)];
    sub2 = [sub2; repmat(subs{ii},length(subvals{2}.bdat.evelen),1)];

    lenmedn1(ii) = median(subvals{1}.bdat.evelen);
    lenmean1(ii) = mean(subvals{1}.bdat.evelen);
    lenmedn2(ii) = median(subvals{2}.bdat.evelen);
    lenmean2(ii) = mean(subvals{2}.bdat.evelen);
end

PDlenmn1 = lenmedn1(PDidx);
ctrllenmn1 = lenmedn1(ctrlidx);
PDlenmn2 = lenmedn2(PDidx);
ctrllenmn2 = lenmedn2(ctrlidx);

PDlenavg1 = mean(PDlenmn1);
ctrllenavg1 = mean(ctrllenmn1);
PDlenavg2 = mean(PDlenmn2);
ctrllenavg2 = mean(ctrllenmn2);

PDlensd1 = std(PDlenmn1);
ctrllensd1 = std(ctrllenmn1);
PDlensd2 = std(PDlenmn2);
ctrllensd2 = std(ctrllenmn2);

figure; 
subplot(1,2,1); histogram(PDlenmn1,20); hold on
subplot(1,2,1); histogram(ctrllenmn1,20); hold off
subplot(1,2,2); histogram(PDlenmn2,20); hold on
subplot(1,2,2); histogram(ctrllenmn2,20); hold off

% Pooled errors t-test
[~, pt1,~,t1] = ttest2(PDlenmn1,ctrllenmn1);
[~, pt2,~,t2] = ttest2(PDlenmn2,ctrllenmn2);
[~, ptPt,~,tPt] = ttest(PDlenmn1,PDlenmn2);
[~, ptCt,~,tCt] = ttest(ctrllenmn1,ctrllenmn2);

mns = [PDlenavg1, ctrllenavg1; PDlenavg2, ctrllenavg2];
sds = [PDlensd1, ctrllensd1; PDlensd2, ctrllensd2];

figure; hold on
bar(1:2,mns,0.6,'b','grouped')
set(gca,'xtick',[])

% Save for export
save('/home/mikkel/PD_motor/rest_ec/groupanalysis/lenevent.mat', ...
    'len1', 'len2', 'sub1', 'sub2');

%% Time to next event
toemedn1 = zeros(length(subs), 1);   % subjects x number of steps
toemedn2 = zeros(length(subs), 1);   % subjects x number of steps

alltoe1 = [];
alltoe2 = [];
sub1 = [];
sub2 = [];

for ii = 1:length(subs)
    load(fullfile(dirs.megDir,subs{ii},'subvals.mat'))
    
    toe1 = zeros(length(subvals{1}.bdat.event-1)-1,1);
    toe2 = zeros(length(subvals{2}.bdat.event-1)-1,1);

    for k = 1:length(subvals{1}.bdat.event)-1
        toe1(k) = (subvals{1}.bdat.event(k+1,1)-subvals{1}.bdat.event(k,2))/1000;
    end
    
    for k = 1:length(subvals{2}.bdat.event)-1
        toe2(k) = (subvals{2}.bdat.event(k+1,1)-subvals{2}.bdat.event(k,2))/1000;
    end
    
    toemedn1(ii) = median(toe1);
    toemedn2(ii) = median(toe2);
    
    alltoe1 = [alltoe1; toe1];
    alltoe2 = [alltoe2; toe2];
    
    sub1 = [sub1; repmat(subs{ii},length(toe1),1)];
    sub2 = [sub2; repmat(subs{ii},length(toe2),1)];
end

PDtoemn1 = toemedn1(PDidx);
ctrltoemn1 = toemedn1(ctrlidx);
PDtoemn2 = toemedn2(PDidx);
ctrltoemn2 = toemedn2(ctrlidx);

figure
subplot(1,2,1); histogram(PDtoemn1,10); hold on
subplot(1,2,1); histogram(ctrltoemn1,10); hold off
subplot(1,2,2); histogram(PDtoemn2,10); hold on
subplot(1,2,2); histogram(ctrltoemn2,10); hold off

PDtoeavg1 = mean(PDtoemn1);
ctrltoeavg1 = mean(ctrltoemn1);
PDtoesd1 = std(PDtoemn1);
ctrltoesd1 = std(ctrltoemn1);

[~, pt1,~,t1] = ttest2(PDtoemn1,ctrltoemn1);
[~, pt2,~,t2] = ttest2(PDtoemn2,ctrltoemn2);
[~, ptPt,~,tPt] = ttest(PDtoemn1,PDtoemn2);
[~, ptCt,~,tCt] = ttest(ctrltoemn1,ctrltoemn2);

% Save for export
save('/home/mikkel/PD_motor/rest_ec/groupanalysis/toevent.mat', ...
    'alltoe1', 'alltoe2', 'sub1', 'sub2');

%% Max peak in events
maxmedn1 = zeros(length(subs), 1);   % subjects x number of steps
maxmedn2 = zeros(length(subs), 1);   % subjects x number of steps

max1 = [];
max2 = [];
sub1 = [];
sub2 = [];

for ii = 1:length(subs)
    load(fullfile(dirs.megDir,subs{ii},'subvals.mat'))
    
    max1 = [max1; subvals{1}.bdat.maxpk];
    max2 = [max2; subvals{2}.bdat.maxpk];
    
    sub1 = [sub1; repmat(subs{ii},length(subvals{1}.bdat.maxpk),1)];
    sub2 = [sub2; repmat(subs{ii},length(subvals{2}.bdat.maxpk),1)];

    maxmedn1(ii) = median(subvals{1}.bdat.maxpk);
    maxmedn2(ii) = median(subvals{2}.bdat.maxpk);
    
    % Save pkidx for reading and plotting in MNE-Py
    maxidx1 = subvals{1}.bdat.maxidx;
    maxidx2 = subvals{2}.bdat.maxidx;
    save(fullfile(dirs.megDir,subs{ii},'pkidx1.mat'),'maxidx1')
    save(fullfile(dirs.megDir,subs{ii},'pkidx2.mat'),'maxidx2')
end

PDmaxmn1 = maxmedn1(PDidx);
ctrlmaxmn1 = maxmedn1(ctrlidx);
PDmaxmn2 = maxmedn2(PDidx);
ctrlmaxmn2 = maxmedn2(ctrlidx);

PDlenavg1 = mean(PDmaxmn1);
ctrllenavg1 = mean(ctrlmaxmn1);
PDlenavg2 = mean(PDmaxmn2);
ctrllenavg2 = mean(ctrlmaxmn2);

PDlensd1 = std(PDmaxmn1);
ctrllensd1 = std(ctrlmaxmn1);
PDlensd2 = std(PDmaxmn2);
ctrllensd2 = std(ctrlmaxmn2);

figure; 
subplot(1,2,1); histogram(PDmaxmn1,10); hold on
subplot(1,2,1); histogram(ctrlmaxmn1,10); hold off
subplot(1,2,2); histogram(PDmaxmn2,10); hold on
subplot(1,2,2); histogram(ctrlmaxmn2,10); hold off

[~, pt1,~,t1] = ttest2(PDmaxmn1,ctrlmaxmn1);
[~, pt2,~,t2] = ttest2(PDmaxmn2,ctrlmaxmn2);
[~, ptPt,~,tPt] = ttest(PDmaxmn1,PDmaxmn2);
[~, ptCt,~,tCt] = ttest(ctrlmaxmn1,ctrlmaxmn2);

save('/home/mikkel/PD_motor/rest_ec/groupanalysis/pkmaxevent.mat', ...
    'max1', 'max2', 'sub1', 'sub2');


% END