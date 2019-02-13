% Get various summaries of beta burst (explorative)
clear all
close all

addpath /home/mikkel/PD_motor/global_scripts
[dirs, subjs, ~] = PD_proj_setup('betaburst');

cd(dirs.megDir);

subs = dir(dirs.megDir);                                %Find subjects in folder
subs = {subs([subs.isdir]).name};                       %Make list
subs = subs(~(strcmp('.',subs)|strcmp('..',subs)));     %Remove dots

[PD_subs, PDidx] = intersect(subs,subjs.PD);
[ctrl_subs, ctrlidx] = intersect(subs, subjs.ctrl);

%% Load threshold 

%% N events
nevent1 = nan(length(subs), 1);
nevent2 = nan(length(subs), 1);

for ii = 1:length(subs)
    load(fullfile(dirs.megDir,subs{ii},'subvals.mat'))
    
    nevent1(ii) = subvals{1}.n_events;
    nevent2(ii) = subvals{2}.n_events;
end

PDn1 = nevent1(PDidx);
ctrln1 = nevent1(ctrlidx);
PDn2 = nevent2(PDidx);
ctrln2 = nevent2(ctrlidx);

PDnavg1 = nanmean(PDn1);
ctrlnavg1 = nanmean(ctrln1);
PDnavg2 = nanmean(PDn2);
ctrlnavg2 = nanmean(ctrln2);

PDnsd1 = nanstd(PDn1);
ctrlnsd1 = nanstd(ctrln1);
PDnsd2 = nanstd(PDn2);
ctrlnsd2 = nanstd(ctrln2);

[~, pt1,~,t1] = ttest2(PDn1,ctrln1);
[~, pt2,~,t2] = ttest2(PDn2,ctrln2);
[~, ptPt,~,tPt] = ttest(PDn1,PDn2);
[~, ptCt,~,tCt] = ttest(ctrln1,ctrln2);


figure
h1 = histogram(PDn1,10); hold on
h2 = histogram(ctrln1,10); hold off

mns = [PDnavg1, ctrlnavg1; PDnavg2, ctrlnavg2];
sds = [PDnsd1, ctrlnsd1; PDnsd2, ctrlnsd2];

%% Event duration
lenmean1 = zeros(length(subs), 1);   % subjects x number of steps
lenmedn1 = zeros(length(subs), 1);   % subjects x number of steps
lenmean2 = zeros(length(subs), 1);   % subjects x number of steps
lenmedn2 = zeros(length(subs), 1);   % subjects x number of steps

for ii = 1:length(subs)
    load(fullfile(dirs.megDir,subs{ii},'subvals.mat'))
    
    len1 = subvals{1}.bdat.evelen;
    len2 = subvals{2}.bdat.evelen;

    lenmedn1(ii) = median(subvals{1}.bdat.evelen);
    lenmean1(ii) = mean(subvals{1}.bdat.evelen);
    lenmedn2(ii) = median(subvals{2}.bdat.evelen);
    lenmean2(ii) = mean(subvals{2}.bdat.evelen);
end


mean(lenmean1)
mean(lenmedn1)

PDlenmn1 = lenmedn1(PDidx);
ctrllenmn1 = lenmedn1(ctrlidx);
PDlenmn2 = lenmedn2(PDidx);
ctrllenmn2 = lenmedn2(ctrlidx);

PDlenavg1 = mean(PDlenmn1);
ctrllenavg1 = mean(ctrllenmn1);
PDlenavg2 = mean(PDlenmn2);
ctrllenavg2 = mean(ctrllenmn2);

figure; 
subplot(1,2,1); histogram(PDlenmn1,20); hold on
subplot(1,2,1); histogram(ctrllenmn1,20); hold off
subplot(1,2,2); histogram(PDlenmn2,20); hold on
subplot(1,2,2); histogram(ctrllenmn2,20); hold off


[~, pt1,~,t1] = ttest2(PDlenmn1,ctrllenmn1);
[~, pt2,~,t2] = ttest2(PDlenmn2,ctrllenmn2);
[~, ptPt,~,tPt] = ttest2(PDlenmn1,PDlenmn2);
[~, ptCt,~,tCt] = ttest2(ctrllenmn1,ctrllenmn2);



mns = [PDlenavg1, ctrllenavg1; ];
sds = [PDlensd(17), ctrllensd(17)];
sem = [PDlensd(17)/sqrt(length(PD_subs)), ctrllensd(17)/sqrt(length(ctrl_subs))];

figure; hold on
bar(1:2,mns,0.6,'b','grouped')
errorbar(1:2,mns,sem*2,'r.')
set(gca,'xtick',[])

%% Time to next event
toemedn = zeros(length(subs), 51);   % subjects x number of steps
toemean = zeros(length(subs), 51);   % subjects x number of steps

for ii = 1:length(subs)
    load(fullfile(dirs.megDir,subs{ii},'subvals.mat'))
    
    for j = 1:length(subvals.steps)
        toe = zeros(length(subvals.bdat{j}.endsam)-1,1);
        for k = 1:length(subvals.bdat{j}.endsam)-1
            toe(k) = (subvals.bdat{j}.begsam(k+1)-subvals.bdat{j}.endsam(k))/1000;
        end
        toemedn(ii,j) = median(toe);
        toemean(ii,j) = mean(toe);
    end
%     clear subvals
end

mean(toemedn)
mean(toemean)

figure; hold on
plot(steps,toemean, 'bo-');
plot(steps,toemedn, 'rx-'); hold off

PDlenmn1 = toemedn(PDidx,:);
ctrllenmn1 = toemedn(ctrlidx,:);

PDlenavg = median(PDlenmn1);
ctrllenavg = median(ctrllenmn1);
PDlensd = std(PDlenmn1);
ctrllensd = std(ctrllenmn1);

figure; hold on
plot(subvals.steps,PDlenavg,'or-');
plot(subvals.steps,ctrllenavg,'ob-'); hold off

for ii = 1:length(subvals.steps)
    [~, pt(ii)] = ttest2(PDlenmn1(:,ii),ctrllenmn1(:,ii));
end
figure; plot(subvals.steps,pt)

mns = [PDlenavg(17), ctrllenavg(17)];
sds = [PDlensd(17), ctrllensd(17)];
sem = [PDlensd(17)/sqrt(length(PD_subs)), ctrllensd(17)/sqrt(length(ctrl_subs))];

figure; hold on
bar(1:2,mns,0.6,'b','grouped')
errorbar(1:2,mns,sem*2,'r.')
set(gca,'xtick',[])




hist()

%% Max peak in events




%% 