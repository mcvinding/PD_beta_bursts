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

%% N tirals
steps = subvals{1}.steps;
nevent = nan(length(subs), 51);   % subjects x number of steps

for i = 1:length(subs)
    % Skip for now
    if strcmp('0327',subs{i})
        continue
    end
    load(fullfile(dirs.megDir,subs{i},'subvals2.mat'))
    
    for j = 1:length(subvals{1}.steps)
        nevent(i,j) = length(subvals{1}.bdat(j).evelen);
    end
%     clear subvals
end

lala = nanmean(nevent)
PDn = nevent(PDidx,:);
ctrln = nevent(ctrlidx,:);

PDnavg = nanmean(PDn)
ctrlnavg = nanmean(ctrln)
PDnsd = nanstd(PDn)
ctrlnsd = nanstd(ctrln)

plot(steps,lala,'o-')

plot(steps,PDnavg,'or-'); hold on
plot(steps,ctrlnavg,'ob-'); hold off

for i = 1:length(steps)
    [~, pt(i),~,t] = ttest2(PDn(:,i),ctrln(:,i));
end

plot(steps,pt)

h1 = histogram(PDn(:,17),20); hold on
h2 = histogram(ctrln(:,17),20);


mns = [PDnavg(17), ctrlnavg(17)];
sds = [PDnsd(17), ctrlnsd(17)];
sem = [PDnsd(17)/sqrt(length(PD_subs)), ctrlnsd(17)/sqrt(length(ctrl_subs))];

figure
hold on
bar(1:2,mns,0.6,'b','grouped')
errorbar(1:2,mns,sem*2,'r.')
set(gca,'xtick',[])

%% Duration
lenmean = zeros(length(subs), 51);   % subjects x number of steps
lenmedn = zeros(length(subs), 51);   % subjects x number of steps

for i = 1:length(subs)
    load(fullfile(dirs.megDir,subs{i},'subvals.mat'))
    
    for j = 1:length(subvals.steps)
        lenmedn(i,j) = median(subvals.bdat{j}.evelen);
        lenmean(i,j) = mean(subvals.bdat{j}.evelen);
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

plot(subvals.steps,PDlenavg,'or-'); hold on
plot(subvals.steps,ctrllenavg,'ob-'); hold off

for i = 1:length(subvals.steps)
    [~, pt(i)] = ttest2(PDlenmn(:,i),ctrllenmn(:,i));
end
plot(subvals.steps,pt)

% Median
PDlenmd = lenmedn(PDidx,:);
ctrllenmd = lenmedn(ctrlidx,:);

PDlenavg = mean(PDlenmd);
ctrllenavg = mean(ctrllenmd);
PDlensd = std(PDlenmd);
ctrllensd = std(ctrllenmd);

plot(subvals.steps,PDlenavg,'or-'); hold on
plot(subvals.steps,ctrllenavg,'ob-'); hold off

for i = 1:length(subvals.steps)
    [~, pt(i)] = ttest2(PDlenmd(:,i),ctrllenmd(:,i));
end

plot(subvals.steps,pt)



mns = [PDlenavg(12), ctrllenavg(12)];
sds = [PDlensd(12), ctrllensd(12)];
sem = [PDlensd(12)/sqrt(length(PD_subs)), ctrllensd(12)/sqrt(length(ctrl_subs))];

figure
hold on
bar(1:2,mns,0.6,'b','grouped')
errorbar(1:2,mns,sem,'r.')
set(gca,'xtick',[])

%% Time to next event
%% Duration
toemedn = zeros(length(subs), 51);   % subjects x number of steps
toemean = zeros(length(subs), 51);   % subjects x number of steps

for i = 1:length(subs)
    load(fullfile(dirs.megDir,subs{i},'subvals.mat'))
    
    for j = 1:length(subvals.steps)
        toe = zeros(length(subvals.bdat{j}.endsam)-1,1);
        for k = 1:length(subvals.bdat{j}.endsam)-1
            toe(k) = (subvals.bdat{j}.begsam(k+1)-subvals.bdat{j}.endsam(k))/1000;
        end
        toemedn(i,j) = median(toe);
        toemean(i,j) = mean(toe);
    end
%     clear subvals
end

mean(toemedn)
mean(toemean)

plot(steps,toemean, 'bo-'); hold on
plot(steps,toemedn, 'rx-'); hold off

PDlenmn = toemedn(PDidx,:);
ctrllenmn = toemedn(ctrlidx,:);

PDlenavg = mean(PDlenmn);
ctrllenavg = mean(ctrllenmn);
PDlensd = std(PDlenmn);
ctrllensd = std(ctrllenmn);

plot(subvals.steps,PDlenavg,'or-'); hold on
plot(subvals.steps,ctrllenavg,'ob-'); hold off

for i = 1:length(subvals.steps)
    [~, pt(i)] = ttest2(PDlenmn(:,i),ctrllenmn(:,i));
end
plot(subvals.steps,pt)

mns = [PDlenavg(12), ctrllenavg(12)];
sds = [PDlensd(12), ctrllensd(12)];
sem = [PDlensd(12)/sqrt(length(PD_subs)), ctrllensd(12)/sqrt(length(ctrl_subs))];

figure
hold on
bar(1:2,mns,0.6,'b','grouped')
errorbar(1:2,mns,sem,'r.')
set(gca,'xtick',[])