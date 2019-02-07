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

%% N tirals
steps = subvals{1}.steps;
nevent1 = nan(length(subs), 1);
nevent2 = nan(length(subs), 1);

for ii = 1:length(subs)
    load(fullfile(dirs.megDir,subs{ii},'subvals2.mat'))
    
    nevent1(ii) = subvals{1}.n_events;
    nevent2(ii) = subvals{1}.n_events;
end

lala = nanmean(nevent1)
PDn = nevent1(PDidx,:);
ctrln = nevent1(ctrlidx,:);

PDnavg = nanmean(PDn)
ctrlnavg = nanmean(ctrln)
PDnsd = nanstd(PDn)
ctrlnsd = nanstd(ctrln)

plot(steps,lala,'o-')

plot(steps,PDnavg,'or-'); hold on
plot(steps,ctrlnavg,'ob-'); hold off

% for i = 1:length(steps)
%     [~, pt(i),~,t] = ttest2(PDn(:,i),ctrln(:,i));
% end
% 
% plot(steps,pt)

h1 = histogram(PDn,20); hold on
h2 = histogram(ctrln,20); hold off

mns = [PDnavg, ctrlnavg];
sds = [PDnsd, ctrlnsd];
sem = [PDnsd/sqrt(length(PD_subs)), ctrlnsd/sqrt(length(ctrl_subs))];

figure; hold on
bar(1:2,mns,0.6,'b','grouped')
errorbar(1:2,mns,sem*2,'r.')
set(gca,'xtick',[])

%% Duration
lenmean = zeros(length(subs), 1);   % subjects x number of steps
lenmedn = zeros(length(subs), 1);   % subjects x number of steps

for ii = 1:length(subs)
    load(fullfile(dirs.megDir,subs{ii},'subvals.mat'))
    
    for j = 1:length(subvals.steps)
        lenmedn(ii,j) = median(subvals.bdat{j}.evelen);
        lenmean(ii,j) = mean(subvals.bdat{j}.evelen);
    end
    
%     clear subvals
end


mean(lenmean)
mean(lenmedn)

plot(steps,lenmean, 'bo-'); hold on
plot(steps,lenmedn, 'rx-'); hold off


PDlenmn = lenmean(PDidx,:);
ctrllenmn = lenmean(ctrlidx,:);

PDlenavg = mean(PDlenmn);
ctrllenavg = mean(ctrllenmn);
PDlensd = std(PDlenmn);
ctrllensd = std(ctrllenmn);

figure;
plot(subvals.steps,PDlenavg,'or-'); hold on
plot(subvals.steps,ctrllenavg,'ob-'); hold off

for ii = 1:length(subvals.steps)
    [~, pt(ii)] = ttest2(PDlenmn(:,ii),ctrllenmn(:,ii));
end
plot(subvals.steps,pt)

% Median
PDlenmd = lenmedn(PDidx,:);
ctrllenmd = lenmedn(ctrlidx,:);

PDlenavg = mean(PDlenmd);
ctrllenavg = mean(ctrllenmd);
PDlensd = std(PDlenmd);
ctrllensd = std(ctrllenmd);

figure; hold on
plot(subvals.steps,PDlenavg,'or-');
plot(subvals.steps,ctrllenavg,'ob-'); hold off

for ii = 1:length(subvals.steps)
    [~, pt(ii)] = ttest2(PDlenmd(:,ii),ctrllenmd(:,ii));
end

figure; plot(subvals.steps,pt)



mns = [PDlenavg(17), ctrllenavg(17)];
sds = [PDlensd(17), ctrllensd(17)];
sem = [PDlensd(17)/sqrt(length(PD_subs)), ctrllensd(17)/sqrt(length(ctrl_subs))];

figure; hold on
bar(1:2,mns,0.6,'b','grouped')
errorbar(1:2,mns,sem*2,'r.')
set(gca,'xtick',[])

%% Time to next event
%% Duration
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

PDlenmn = toemedn(PDidx,:);
ctrllenmn = toemedn(ctrlidx,:);

PDlenavg = median(PDlenmn);
ctrllenavg = median(ctrllenmn);
PDlensd = std(PDlenmn);
ctrllensd = std(ctrllenmn);

figure; hold on
plot(subvals.steps,PDlenavg,'or-');
plot(subvals.steps,ctrllenavg,'ob-'); hold off

for ii = 1:length(subvals.steps)
    [~, pt(ii)] = ttest2(PDlenmn(:,ii),ctrllenmn(:,ii));
end
figure; plot(subvals.steps,pt)

mns = [PDlenavg(17), ctrllenavg(17)];
sds = [PDlensd(17), ctrllensd(17)];
sem = [PDlensd(17)/sqrt(length(PD_subs)), ctrllensd(17)/sqrt(length(ctrl_subs))];

figure; hold on
bar(1:2,mns,0.6,'b','grouped')
errorbar(1:2,mns,sem*2,'r.')
set(gca,'xtick',[])




hist(