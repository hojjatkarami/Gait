clear
close all
clc

nnmf_init

% all changes of this script will be applied to Synergy field of 'P' Database
%% Settings
load('P.mat')
subject = [11];

%% EMG>syn : extract n1:n2 number of synergy
SYN = 1:8;
option.side = 'Right';  % Right: right side, Left: left side; 
option.NormalizeToUnitVariance = 1;     % 0:do not divide, 1:divide and then multiply, 3:divide and do not multiply
option.rep = 10;
%% NNMF
for i_sub = subject
    TrialNum =length( P(i_sub).EMG.(option.side));
    M=[];
    for i_trial = 1:TrialNum
        M = [M, P(i_sub).EMG.(option.side)(i_trial).Mn];
    end
        for i_syn = SYN
            disp(['sub:',num2str(i_sub), 'syn:',num2str(i_syn),', ' ])
            [W_best, S_best, M_rec] = nnmfEMG0720(M,i_syn,option);
%             P(i_sub).Synergy.EMG.(option.side)(i_trial).syn2(i_syn).W_best = W_best;
%             S_best = P(i_sub).Synergy.EMG.(option.side)(i_trial).syn2(i_syn).S_best;
%             M_rec = P(i_sub).Synergy.EMG.(option.side)(i_trial).syn2(i_syn).M_rec;
            for i_trial = 1:TrialNum
                P(i_sub).Synergy.EMG_con.(option.side)(i_trial).syn(i_syn).W_best = W_best;
                P(i_sub).Synergy.EMG_con.(option.side)(i_trial).syn(i_syn).S_best = S_best(:,((i_trial-1)*100+1):(i_trial*100));
                P(i_sub).Synergy.EMG_con.(option.side)(i_trial).syn(i_syn).M_rec = M_rec(:,((i_trial-1)*100+1):(i_trial*100));
            end
        end

end
save('P.mat','P');
%% calculate gof
for i_sub = subject
    TrialNum =length( P(i_sub).EMG.(option.side));
    M=[];
    M_rec=[];
    for i_trial = 1:TrialNum
        M = [M, P(i_sub).EMG.(option.side)(i_trial).Mn];
    end
        for i_syn = SYN
            disp(['sub:',num2str(i_sub), 'syn:',num2str(i_syn),', ' ])
            M_rec=[];

            for i_trial = 1:TrialNum
                M_rec = [M_rec, P(i_sub).Synergy.EMG_con.(option.side)(i_trial).syn(i_syn).M_rec];
            end
            P(i_sub).Synergy.EMG_con.(option.side)(1).gof.VAF(i_syn) = vaf1(M_rec, M, 0);
            P(i_sub).Synergy.EMG_con.(option.side)(1).gof.RSQ(i_syn) = rsq1(M_rec, M, 0);
            P(i_sub).Synergy.EMG_con.(option.side)(1).gof.vaf(i_syn,:) = vaf1(M_rec, M, 1);
            P(i_sub).Synergy.EMG_con.(option.side)(1).gof.rsq(i_syn,:) = rsq1(M_rec, M, 1);

            [r2_bootstat_lb, r2_bootstat_ub] = myBootStrap(100,95,'@rsq1',M,M_rec);
            [vaf_bootstat_lb, vaf_bootstat_ub] = myBootStrap(100,95,'@vaf1',M,M_rec);
            P(i_sub).Synergy.EMG_con.(option.side)(1).gof.BootStrap.r2_bootstat_lb(i_syn) = r2_bootstat_lb;
            P(i_sub).Synergy.EMG_con.(option.side)(1).gof.BootStrap.r2_bootstat_ub(i_syn) = r2_bootstat_ub;
            P(i_sub).Synergy.EMG_con.(option.side)(1).gof.BootStrap.vaf_bootstat_lb(i_syn) = vaf_bootstat_lb;
            P(i_sub).Synergy.EMG_con.(option.side)(1).gof.BootStrap.vaf_bootstat_ub(i_syn) = vaf_bootstat_ub;
            
        end

end
save('P.mat','P');
%% EMG>number
SYN = 1:8;
option.side = 'Right';  % Right: right side, Left: left side;
option.type = 'EMG_con';
vaf_th=0.85;    VAF_th = 0.95;
rsq_th=0.5;     RSQ_th = 0.8;
rsq_boot_th = 0.5;     vaf_boot_th=0.95;

for i_sub = subject
    i_sub
    TrialNum =length(P(i_sub).EMG.Right);
    i_trial = 1;
        gof = P(i_sub).Synergy.(option.type).Right(i_trial).gof  ;
        

        P(i_sub).Synergy.(option.type).(option.side)(i_trial).number.VAF = min(find(gof.VAF>=VAF_th));
        P(i_sub).Synergy.(option.type).(option.side)(i_trial).number.RSQ = min(find(gof.RSQ>=RSQ_th));
        P(i_sub).Synergy.(option.type).(option.side)(i_trial).number.vaf = min(find(min(gof.vaf')>=vaf_th));
        P(i_sub).Synergy.(option.type).(option.side)(i_trial).number.rsq = min(find(min(gof.rsq')>=rsq_th));
        if isempty(P(i_sub).Synergy.(option.type).(option.side)(i_trial).number.rsq)
            P(i_sub).Synergy.(option.type).(option.side)(i_trial).number.rsq=8;
        end

        P(i_sub).Synergy.(option.type).(option.side)(i_trial).number.BootStrap_r2 = min(find(gof.BootStrap.r2_bootstat_lb>=rsq_boot_th));
        P(i_sub).Synergy.(option.type).(option.side)(i_trial).number.BootStrap_vaf = min(find(gof.BootStrap.vaf_bootstat_lb>=vaf_boot_th));

        % Factor Analysis
        
    
end
save('P.mat','P');

%% plot recos of each muscle
SYN=6;
option.side = 'Right';  % Right: right side, Left: left side; 

for i_sub = [11]
for i_syn=SYN
    figure;
    
    M=[];
    M_rec=[];
    TrialNum =length( P(i_sub).EMG.(option.side));
    for i_trial = 1:TrialNum
        M = [M, P(i_sub).EMG.(option.side)(i_trial).Mn];
        M_rec = [M_rec, P(i_sub).Synergy.EMG_con.(option.side)(i_trial).syn(i_syn).M_rec];
    end
    for i=1:16
        subplot(16,1,i);  hold on
        set(gca,'YTick',[],'XTick',[])
        ylabel(muscleName(P(i_sub).EMG.Right(1).muscleOrder(i)))
        plot(M(i,:))
        plot(M_rec(i,:))
   end
end
end
%% plot GoF - Local rsq
subject = [1:11];
th=0.5;
SYN=3:8;
option.side = 'Right';  % Right: right side, Left: left side; 
figure;
% suptitle('local vaf')
for i_sub = subject
     M=[];
    M_rec=[];
    TrialNum =length( P(i_sub).EMG.(option.side));
    label = muscleName(P(i_sub).EMG.Right(1).muscleOrder);
        subplot(2,6,i_sub);  hold on
        title(['Subject ',num2str(i_sub),', N=',...
            num2str(P(i_sub).Synergy.EMG_con.(option.side)(1).number.rsq)])
        set(gca,'FontSize',6)
        grid on
        xlim([0 17])
%         ylim([0.65 1])
        xticks(1:length(label))
        xtickangle(45)
        xticklabels(label)
        line([0 17],[th th],'Color','red','DisplayName',[num2str(th*100) ,'% threshold'])
for i_syn=SYN
   
%         ylabel(muscleName()
        plot(P(i_sub).Synergy.EMG_con.Right(1).gof.rsq(i_syn,:),'DisplayName',[num2str(i_syn),' synergies'])
end
end
lgd = legend('Location','southoutside')
lgd.NumColumns = 1;
%% plot GoF - Global rsq
subject = [1:11];
th=0.8;
SYN=1:8;
option.side = 'Right';  % Right: right side, Left: left side; 
figure;
% suptitle('local vaf')
for i_sub = subject
    N = P(i_sub).Synergy.EMG_con.(option.side)(1).number.RSQ;
     M=[];
    M_rec=[];
    TrialNum =length( P(i_sub).EMG.(option.side));
    label = muscleName(P(i_sub).EMG.Right(1).muscleOrder);
        subplot(2,6,i_sub);  hold on
        title(['Subject ',num2str(i_sub),', N=',...
            num2str(N)])
%         set(gca,'FontSize',6)
        grid on
%         xlim([0 17])
        ylim([0 1])
        xticks(SYN(1):SYN(end))
%         xtickangle(45)
%         xticklabels(label)
        line([0 SYN(end)],[th th],'Color','red','DisplayName',[num2str(th*100) ,'% threshold'])
   
%         ylabel(muscleName()
        plot(SYN,P(i_sub).Synergy.EMG_con.Right(1).gof.RSQ(SYN),'*-','DisplayName','Global R2')
%         x = [N N];
%         y = [P(i_sub).Synergy.EMG_con.Right(1).gof.RSQ(N)+.05 P(i_sub).Synergy.EMG_con.Right(1).gof.RSQ(N)];
%         annotation(gca,'textarrow',x/max(x),y/max(y),'String',['N=',num2str(N)])
        plot(N,P(i_sub).Synergy.EMG_con.Right(1).gof.RSQ(N),'o','MarkerSize',6,'handlevisibility','off')
%         y = (P(i_sub).Synergy.EMG_con.Right(1).gof.BootStrap.r2_bootstat_ub + P(i_sub).Synergy.EMG_con.Right(1).gof.BootStrap.r2_bootstat_lb)/2;
%         y_pos = P(i_sub).Synergy.EMG_con.Right(1).gof.BootStrap.r2_bootstat_ub   - y;
%         y_neg = y - P(i_sub).Synergy.EMG_con.Right(1).gof.BootStrap.r2_bootstat_ub;
%         errorbar(SYN,y(SYN),-y_neg(SYN),y_pos(SYN))

end
legend
% lgd = legend('Location','eastoutside')
%% plot S intra subject

subject = [3:11];
option.side = 'Right';
% SYN = 1:6;
N=6;
load('intra')

for i_sub = subject
    if isempty(P(i_sub).Synergy.EMG_con)
        continue;
    end
    W={}; S={}; name={}; k=1;    
    TrialNum =length(P(i_sub).EMG.Right);
    for i_trial = 1:TrialNum
        W{k} = P(i_sub).Synergy.EMG_con...
                        .(option.side)(i_trial).syn(N).W_best;
        S{k} = P(i_sub).Synergy.EMG_con...
                        .(option.side)(i_trial).syn(N).S_best';
        R2=1;
        name{k} = {['s:',num2str(i_sub),'_te:',...
                            '_tr:',num2str(i_trial),'_syn:',...
                            num2str(N)];...
                            ['R2:',num2str(R2,'%.2f')]};
        k=k+1;
    end   
    
    [id,Stype,S] = kmeanW(S, N, 0.6);
    for i=1:length(S)
        S{i} = transpose(S{i});
    end
    for i=1:length(Stype)
        Stype{i} = transpose(Stype{i});
    end
   
    for j=1:N
        temp=[];
%         for i=1:size(Stype{j},1)
        for i=1:length(S)
%             temp = [temp; Stype{j}(i,:)];
            temp = [temp; S{i}(j,:)];

        end
        intra(i_sub).emg_con.syn(N).meanOfTypeS(j,:)=mean(temp);
        intra(i_sub).emg_con.syn(N).stdOfTypeS(j,:)=std(temp);

    end
    intra(i_sub).emg_con.syn(N).meanOfTypeW = W{1};
%     inter(i_sub).emg_con.syn(N).meanOfTypeS = cell2mat(meanOfTypeS');
                % plot W
%                     h1=figure(500+i_sub);
%                     supTitle = 'Intra Subject clustering';
%                     label = muscleName(P(i_sub).EMG.Right(1).muscleOrder);
%                     plotWtype(h1,W,Wtype,id,name,label,supTitle)
figure
for i=1:N
    subplot(N,4,(i-1)*4+1); hold on;
        bar(intra(i_sub).emg_con.syn(N).meanOfTypeW(:,i))
    subplot(N,4,(i-1)*4+2); hold on;
        bar(intra(i_sub).emg.syn(N).meanOfTypeW(:,i))
        errorbar(intra(i_sub).emg.syn(N).meanOfTypeW(:,i),2*intra(i_sub).emg.syn(N).stdOfTypeW(:,i),'LineStyle','None')
    subplot(N,4,(i-1)*4+3); hold on;
        plot([0:1:99],intra(i_sub).emg_con.syn(N).meanOfTypeS(i,:))
        x=[[0:1:99],fliplr([0:1:99])];
        y=[intra(i_sub).emg_con.syn(N).meanOfTypeS(i,:)-intra(i_sub).emg_con.syn(N).stdOfTypeS(i,:),...
        fliplr(intra(i_sub).emg_con.syn(N).meanOfTypeS(i,:)+intra(i_sub).emg_con.syn(N).stdOfTypeS(i,:))];
        s=fill(x,y,'k','EdgeColor','none');
        alpha(s,.1)
        ylim([0 6])

    subplot(N,4,(i-1)*4+4); hold on;
%         plot(intra(i_sub).emg.syn(N).meanOfTypeS(i,:))
        plot([0:1:99],intra(i_sub).emg.syn(N).meanOfTypeS(i,:))
        x=[[0:1:99],fliplr([0:1:99])];
        y=[intra(i_sub).emg.syn(N).meanOfTypeS(i,:)-intra(i_sub).emg.syn(N).stdOfTypeS(i,:),...
        fliplr(intra(i_sub).emg.syn(N).meanOfTypeS(i,:)+intra(i_sub).emg.syn(N).stdOfTypeS(i,:))];
        s=fill(x,y,'k','EdgeColor','none');
        alpha(s,.1)
        ylim([0 6])

end
                % plot S
%                     h2=figure(5500+i_sub);
%                     plotStype0720(h2,S,Stype,id,name,label,supTitle)
end
save('intra.mat','intra')


%% 900 - EMG SYNERGY concat inter
allN=6;
N=6;
W={}; S={}; name={}; k=1;
load('intra')
%find all muscles of subjects
allMuscles=[];
sameMuscles = [];
for i_sub = subject
    allMuscles = [allMuscles,intra(i_sub).emg.muscleOrder];
end
allMuscles = sort(allMuscles);
% allMuscles(find(diff(allMuscles)==0))=[];   % remove repeating elements
for i=1:18
   if length(find(allMuscles==i) )==length(subject)
       sameMuscles = [sameMuscles i];
   end
end
for i_N = allN

for i_sub = subject
    if isempty(intra(i_sub).emg)
        continue;
    end
%     for i_test = test{i_sub}
%         for i=1:length(sameMuscles)
%             x = find(intra(i_sub).emg.muscleOrder==sameMuscles(i));
%             
%                 W{k}(i,:) = intra(i_sub).emg.syn(N).meanOfTypeW(x,:);
%             
%         end
        W{k} = intra(i_sub).emg_con.syn(i_N).meanOfTypeW;
        
        S{k} = intra(i_sub).emg_con.syn(i_N).meanOfTypeS;
        
        name{k} = ['s:',num2str(i_sub),'N:',num2str(i_N)];
        k=k+1;
        
%     end
end

end
[id,Wtype,~] = kmeanW(W,N,0.6);
n=length(id);

Stype = cell(1,N);
for i=1:N    
    for j=1:n
        x = find(id{j}==i);
        if isempty(x) || id{j}(x)==0
            continue;
        end
         Stype{i}=[Stype{i} ; S{j}(x,:)];
    end
end
    meanOfTypeW = cell(1,N);
    stdOfTypeW = cell(1,N);
    meanOfTypeS = cell(1,N);
    stdOfTypeS = cell(1,N);
    for i=1:N
        meanOfTypeW{i} =  mean(Wtype{i},2);
        stdOfTypeW{i} = std(Wtype{i},0,2);
        meanOfTypeS{i} =  mean(Stype{i},1);
        stdOfTypeS{i} = std(Stype{i},0,1);
    end
    

        % plot W
            h1=figure();
            supTitle = 'Inter Subject clustering';
            label = muscleName(sameMuscles);
            plotWtype(h1,W,Wtype,id,name,label,supTitle)
            
        
        % plot S
            h2=figure();
            plotStype0720(h2,S,Stype,id,name,label,supTitle)
        
        
badW={}; badS={}; badN=2;
for i=1:n
    badW{i} = W{i}(:,find(id{i}==0));
    badS{i} = S{i}(find(id{i}==0),:);
end
[badId,badWtype,~] = kmeanW(badW,badN,0.1);
badStype = cell(1,badN);
for i=1:badN    
    for j=1:n
        x = find(badId{j}==i);
        if isempty(x) || id{j}(x)==0
            continue;
        end
         badStype{i}=[badStype{i} ; badS{j}(x,:)];
    end
end
    meanOfTypeW = cell(1,badN);
    stdOfTypeW = cell(1,badN);
    meanOfTypeS = cell(1,badN);
    stdOfTypeS = cell(1,badN);
    for i=1:badN
        meanOfTypeW{i} =  mean(badWtype{i},2);
        stdOfTypeW{i} = std(badWtype{i},0,2);
        meanOfTypeS{i} =  mean(badStype{i},1);
        stdOfTypeS{i} = std(badStype{i},0,1);
    end
    
% plot bad W
    h1=figure();
    supTitle = 'Inter Subject clustering';
    label = muscleName(sameMuscles);
    plotWtype(h1,badW,badWtype,badId,name,label,supTitle)
% plot bad S
    h2=figure();
    plotStype0720(h2,badS,badStype,badId,name,label,supTitle)
        
