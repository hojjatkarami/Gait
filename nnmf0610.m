clear
close all
clc

nnmf_init

% all changes of this script will be applied to Synergy field of 'P' Database
%% Settings
load('P.mat')


%% NNMF stage 1
tic
subject =[5];
% test = {[4],[4],[4],[4],[4]};
% trial = {{100 100 100},{100 100 100},{100 100 100},{100 100 100},{100 100 100}};
SYN = 6;
option.type = 'EMG';    % EMG or kin2 or EmgKin
% state='Pre';
option.side = 'Right';  % Right: right side, Left: left side; 
option.rep = 1;

for i_sub = subject
    
%     for i_test = test{i_sub} 
if strcmp(option.side, 'Right')
        TrialNum =length(P(i_sub).Trajectory.Right);
else
        TrialNum =length(P(i_sub).Trajectory.Left);
end
%         if TrialNum==0
%             continue;
%         end
% TrialNum=1;
        for i_trial = 1:TrialNum
            for i_syn = SYN
                disp(['sub:',num2str(i_sub),', trial:',num2str(i_trial),', syn:',num2str(i_syn),', ' ])
%                 M = Pf(i_sub).Test(i_test).Synergy(i_trial).(option.type).mat;
                if strcmp(option.type, 'EmgKin')
                   M=P(i_sub);
                else
                   M=P(i_sub).(option.type);
                end
                switch option.type
                    case 'EMG'                        
                        [W_best, S_best,W_avg,S_all] = nnmfEMG0524(M,i_syn,option,i_trial);
                    case 'kin2'                        
                        [W_best, S_best,W_avg,S_all] = nnmfKin0524(M,i_syn,option,i_trial);
                    case 'EmgKin'
                        [W_best, S_best,W_avg,S_all] = nnmfEMGKin0524(M,i_syn,option,i_trial);
                        
                end
                P(i_sub).Synergy(i_trial).(option.type).(option.side).syn(i_syn).W_avg = W_avg;
                P(i_sub).Synergy(i_trial).(option.type).(option.side).syn(i_syn).W_best = W_best;
                P(i_sub).Synergy(i_trial).(option.type).(option.side).syn(i_syn).S_best = S_best;
            end
        end    
end

toc
save('P.mat','P')
disp('nnmf done!')





%% calculate goodness of fits

subject = [1:10];
test=1;
trial=1;
% test = {[4],[4],[4],[4],[4]};
% trial = {{100 100 100},{100 100 100},{100 100 100},{100 100 100},{100 100 100}};
SYN = 4:6;
option.type = 'EMG';    % EMG or kin2 or EmgKin
% state='Pre';
option.side = 'Right';  % Right: right side, Left: left side; 
W={};S={};M={};rec={}; name={};
k=1;
for i_sub = subject
    W={};S={};M={};rec={}; name={};
    k=1;
    for i_test = test 
if strcmp(option.side, 'Right')
        TrialNum =length(P(i_sub).Trajectory.Right);
else
        TrialNum =length(P(i_sub).Trajectory.Left);
end
        if TrialNum==0
            continue;
        end
%         TrialNum=1;
        for i_trial = 1:TrialNum
            for i_syn = SYN
                disp(['sub:',num2str(i_sub),', trial:',num2str(i_trial),', syn:',num2str(i_syn),', ' ])
%                 M = Pf(i_sub).Test(i_test).Synergy(i_trial).(option.type).mat;
                if strcmp(option.type, 'EmgKin')
                   M{k}=P(i_sub);
                else
                   M{k}=P(i_sub).(option.type).(option.side)(i_trial).M_R;
                end
                if strcmp(option.type,'kin2')
                   M{k}=kinProcess(M{k}); 
                end
                name{k} = [num2str(i_sub),num2str(i_test),num2str(i_trial),num2str(i_syn)];
                W{k} = P(i_sub).Synergy(i_trial).(option.type).(option.side).syn(i_syn).W_best;
                S{k} = P(i_sub).Synergy(i_trial).(option.type).(option.side).syn(i_syn).S_best;
                rec{k} = W{k}*S{k};
                
                
                %
                    [r2_bootstat_lb(k), r2_bootstat_ub(k)] = myBootStrap(100,95,'@rsq1',M{k},rec{k});
                    [vaf_bootstat_lb(k), vaf_bootstat_ub(k)] = myBootStrap(100,95,'@vaf1',M{k},rec{k});
                    VAF(k) = vaf1(rec{k},M{k},0);
                    RSQ(k) = rsq1(rec{k},M{k},0);
                    vafLocal = vaf1(rec{k},M{k},1);
                    [minVal, id] = min(vafLocal);
                    rsqLocal = rsq1(rec{k},M{k},1);
                    [minVal, id] = min(rsqLocal);
                    P(i_sub).Synergy(i_trial).(option.type).Right.syn(i_syn).r2_bootstat_lb = r2_bootstat_lb(k);
                    P(i_sub).Synergy(i_trial).(option.type).Right.syn(i_syn).VAF = VAF(k);
                    P(i_sub).Synergy(i_trial).(option.type).Right.syn(i_syn).RSQ = RSQ(k);
                    P(i_sub).Synergy(i_trial).(option.type).Right.syn(i_syn).vafLocal = vafLocal;
                    P(i_sub).Synergy(i_trial).(option.type).Right.syn(i_syn).minVafLocal = min(vafLocal);

                    P(i_sub).Synergy(i_trial).(option.type).Right.syn(i_syn).rsqLocal = rsqLocal;
                    P(i_sub).Synergy(i_trial).(option.type).Right.syn(i_syn).minRsqLocal = min(rsqLocal);

                    k=k+1;
            end
            
        end    
    end
end
%% find number of synergies
subject = [1 2 3 5 6 8];
option.type = 'EMG';    % EMG or kin2 or EmgKin

figure; hold on;    grid on
xlabel('Number of Subject')
ylabel('Number of Synergy')
title('number of synergy for each trial of subjects')
n={};
meanOfSubject=[];
stdOfSubject=[];
for i_sub=subject
    a=[];
   for i_trial=1: length(P(i_sub).Trajectory.Right)
       
      n{i_sub}(i_trial)=3+min(find(extractfield(P(i_sub).Synergy(i_trial).(option.type).Right.syn(:),'r2_bootstat_lb')>.8));

       plot(i_sub+i_trial*.02,n{i_sub}(i_trial),'k*')
       
       
    end  
%     plot(i_sub,mean(a),'b*')
meanOfSubject(i_sub) = mean(n{i_sub});
stdOfSubject(i_sub) = std(n{i_sub});

end
% figure; hold on;    grid on
s = bar(meanOfSubject);
alpha(s,.5)
errorbar(meanOfSubject,stdOfSubject,'LineStyle','none')
for i=1:length(n)
end

%% plot goodness of fits
% close all
figure;
hold on;
label = {'RTA','RPL','RSOL','RGC',...
        'RRF','RMH','RVL',...
        'RIP','RGMAX','RGMED','RAD','RTFL',...
        'RIC','RLG','RRA','REO'};
vaf_th=0.75;
rsq_th=0.6;
lineStyle={'--','-'};
VAF=[];
RSQ=[];


for i=1:length(M)
    
    [r2_bootstat_lb(i), r2_bootstat_ub(i)] = myBootStrap(100,95,'@rsq1',M{i},rec{i});
    [vaf_bootstat_lb(i), vaf_bootstat_ub(i)] = myBootStrap(100,95,'@vaf1',M{i},rec{i});
    VAF(i) = vaf1(rec{i},M{i},0);
    RSQ(i) = rsq1(rec{i},M{i},0);
    
    subplot(2,2,3); hold on;    title('Local VAF'); grid on
    set(gca,'ButtonDownFcn',{@plot_GoF_Click,M,rec,label}) 

    vafLocal = vaf1(rec{i},M{i},1);
    [minVal, id] = min(vafLocal);
    
    plot(vafLocal,'linestyle',lineStyle{(minVal>vaf_th)+1},'DisplayName',name{i})
    plot(id,minVal,'r*','handlevisibility','off')
    legend
    xticks(1:length(label))
    xtickangle(45)
    xticklabels(label)

    subplot(2,2,4); hold on;    title('Local R-sqaure'); grid on
    set(gca,'ButtonDownFcn',{@plot_GoF_Click,M,rec,label}) 

    vafLocal = rsq1(rec{i},M{i},1);
    [minVal, id] = min(vafLocal);
    plot(vafLocal,'linestyle',lineStyle{(minVal>rsq_th)+1},'DisplayName',name{i})
    plot(id,minVal,'r*','handlevisibility','off')
    legend
    xticks(1:length(label))
    xtickangle(45)
    xticklabels(label)

end
subplot(2,2,1); hold on;    title('Total VAF');grid on
meanVal = (vaf_bootstat_lb+vaf_bootstat_ub)/2;
err = vaf_bootstat_ub-meanVal;
errorbar(meanVal,err)
plot(VAF)
ylim([0.6 1])
% xlim([0 8])
subplot(2,2,2); hold on;    title('Total R-sqaure'); grid on
meanVal = (r2_bootstat_lb+r2_bootstat_ub)/2;
err = r2_bootstat_ub-meanVal;
errorbar(meanVal,err)
plot(RSQ)
ylim([0 1])
% xlim([0 8])

%% each muscle activity
close all
subject=[1:10];
musNo=14;
option.type = 'EMG';    % EMG or kin2 or EmgKin
option.side = 'Right';  % Right: right side, Left: left side; 
for i_sub = subject
    
if strcmp(option.side, 'Right')
        TrialNum =length(P(i_sub).Trajectory.Right);
else
        TrialNum =length(P(i_sub).Trajectory.Left);
end
        if TrialNum==0
            continue;
        end
        figure; hold on
        for i_mus=1:musNo
            subplot(4,4,i_mus); hold on;
            title(P(i_sub).(option.type).(option.side)(1).groupName{i_mus});
            sig=[];
            for i_trial = 1:TrialNum
                
                temp = P(i_sub).(option.type).(option.side)(i_trial).M_R(i_mus,:);
                temp2 = P(i_sub).Events.(option.side)(i_trial);
                [t,sig(i_trial,:)] = NormS(temp,double(temp2.PurtFrame));
%                 plot(t,sig(i_trial,:))

            end
            meanSig = mean(sig);
            stdSig = std(sig);
               plot([0:1:99],meanSig,'Color','k')
    x=[[0:1:99],fliplr([0:1:99])];
    y=[meanSig-stdSig,...
    fliplr(meanSig+stdSig)];
    s=fill(x,y,'k','EdgeColor','none');
    alpha(s,.1)
        end
    
end
