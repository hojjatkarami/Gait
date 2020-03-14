clc
clear
close all
nnmf_init
%% Settings
load('P.mat')
subject = [3:11];
SYN = 1:8;
option.side = 'Right';  % Right: right side, Left: left side;
% vaf_th=0.9;    VAF_th = 0.95;
% rsq_th=0.7;     RSQ_th = 0.9;
% rsq_boot_th = 0.8;     vaf_boot_th=0.90;
% lineStyle={'--','-'};
% VAF=[];
% RSQ=[];
% 
% label2 = {'Ank 1','Ank 2','Knee','Hip 1','Hip 2','Hip 3','Pelv 1','Pelv2'};
%% for trials of a subject
subject=5;
i_syn = 6;
for i_sub = subject
    if isempty(P(i_sub).Synergy.EMG)
        continue;
    end
    TrialNum =length(P(i_sub).EMG.Right);
    for i_trial = 1:1
        W_EMG = P(i_sub).Synergy.EMG.Right(i_trial).syn(i_syn).W_best;
        S_EMG = P(i_sub).Synergy.EMG.Right(i_trial).syn(i_syn).S_best;

        W_kin = P(i_sub).Synergy.kin.Right(i_trial).syn(i_syn).W_best;
        S_kin = P(i_sub).Synergy.kin.Right(i_trial).syn(i_syn).S_best;
        
        W_EMGkin = P(i_sub).Synergy.EMGkin.Right(i_trial).syn(i_syn).W_best;
        S_EMGkin = P(i_sub).Synergy.EMGkin.Right(i_trial).syn(i_syn).S_best;
        
        lbl_EMG = muscleName(P(i_sub).EMG.Right(1).muscleOrder);
        lbl_kin = angleName(1:18);
        plot_EMGkin(W_EMG, S_EMG, W_kin, S_kin,W_EMGkin,S_EMGkin,lbl_EMG, lbl_kin)
    
    end
    
    
end

%% Intra-subject
subject=11;
i_syn = 6;
load('intra')

for i_sub = subject
    if isempty(P(i_sub).Synergy.EMG)
        continue;
    end
    
        W_EMG = intra(i_sub).emg.syn(i_syn).meanOfTypeW;
        S_EMG = intra(i_sub).emg.syn(i_syn).meanOfTypeS;
        W_kin = intra(i_sub).kin.syn(i_syn).meanOfTypeW;
        S_kin = intra(i_sub).kin.syn(i_syn).meanOfTypeS;
        
        W_EMG_std = intra(i_sub).emg.syn(i_syn).stdOfTypeW;
        S_EMG_std = intra(i_sub).emg.syn(i_syn).stdOfTypeS;
        W_kin_std = intra(i_sub).kin.syn(i_syn).stdOfTypeW;
        S_kin_std = intra(i_sub).kin.syn(i_syn).stdOfTypeS;
        
        W_EMGkin = intra(i_sub).emgkin.syn(i_syn).meanOfTypeW;
        S_EMGkin = intra(i_sub).emgkin.syn(i_syn).meanOfTypeS;
        
        lbl_EMG = muscleName(P(i_sub).EMG.Right(1).muscleOrder);
        lbl_kin = angleName(1:18);
        plot_EMGkin(W_EMG, S_EMG, W_kin, S_kin,W_EMG_std, S_EMG_std, W_kin_std, S_kin_std,...
                                        W_EMGkin,S_EMGkin,lbl_EMG, lbl_kin)
    
    
    
    
end




%% Inter-subject

i_syn = 6;
load('inter')

    
        W_EMG = inter.emg_con.syn(i_syn).meanOfTypeW;
        S_EMG = inter.emg_con.syn(i_syn).meanOfTypeS;
        W_kin = inter.kin_con.syn(i_syn).meanOfTypeW;
        S_kin = inter.kin_con.syn(i_syn).meanOfTypeS;
        
        W_EMG_std = inter.emg_con.syn(i_syn).stdOfTypeW;
        S_EMG_std = inter.emg_con.syn(i_syn).stdOfTypeS;
        W_kin_std = inter.kin_con.syn(i_syn).stdOfTypeW;
        S_kin_std = inter.kin_con.syn(i_syn).stdOfTypeS;
        
        W_EMGkin = inter.emgkin.syn(i_syn).meanOfTypeW;
        S_EMGkin = inter.emgkin.syn(i_syn).meanOfTypeS;
        
%         W_EMG = max(inter.emg.syn(i_syn).meanOfTypeW-2*inter.emg.syn(i_syn).stdOfTypeW,0);
%         S_EMG = inter.emg.syn(i_syn).meanOfTypeS;
% 
%         W_kin = max(inter.kin.syn(i_syn).meanOfTypeW-2*inter.kin.syn(i_syn).stdOfTypeW,0);
%         S_kin = inter.kin.syn(i_syn).meanOfTypeS;
%         
%         W_EMGkin = max(inter.emgkin.syn(i_syn).meanOfTypeW-2*inter.emgkin.syn(i_syn).stdOfTypeW,0);
%         S_EMGkin = inter.emgkin.syn(i_syn).meanOfTypeS;
        
        lbl_EMG = muscleName(P(5).EMG.Right(1).muscleOrder);
        lbl_kin = angleName(1:18);
       plot_EMGkin(W_EMG, S_EMG, W_kin, S_kin,W_EMG_std, S_EMG_std, W_kin_std, S_kin_std,...
                                        W_EMGkin,S_EMGkin,lbl_EMG, lbl_kin)
    
    




%% functions

function plot_EMGkin(W_EMG, S_EMG, W_kin, S_kin,W_EMG_std, S_EMG_std, W_kin_std, S_kin_std,...
                                        W_EMGkin,S_EMGkin,lbl_EMG, lbl_kin)
    [S_EMG, order] = mySort(S_EMG);
    W_EMG = W_EMG(:,order); W_EMG_std = W_EMG_std(:,order); S_EMG_std = S_EMG_std(order,:);
    [S_kin, order] = mySort(S_kin);
    W_kin = W_kin(:,order); W_kin_std = W_kin_std(:,order); S_kin_std = S_kin_std(order,:);
    [S_EMGkin, order] = mySort(S_EMGkin);
    W_EMGkin = W_EMGkin(:,order);
    
    N=size(W_EMG,2);
    h1=figure;
    set(h1,'units','normalized','outerposition',[0 0 1 1]);% suptitle(supTitle)

    for i=1:N
        subplot(N,3,(i-1)*3+1);  hold on;   set(gca,'YTick',[],'XTick',[],'ButtonDownFcn',{@Click_showXtick,lbl_EMG});
%             h=bar([W_EMG(:,i) W_EMGkin(1:16,i)]);
             h=bar([W_EMG(:,i)]);
             errorbar(W_EMG(:,i),W_EMG_std(:,i),'linestyle','none','color','black')
            set(h(1), 'FaceColor','red');
%             set(h(2), 'FaceColor','green');
                xticks(1:length(lbl_EMG))
                xtickangle(45)
                xticklabels(lbl_EMG)  
                ylabel(['synergy ',num2str(i)])
        subplot(N,3,(i-1)*3+2);  hold on;   set(gca,'YTick',[],'XTick',[],'ButtonDownFcn',{@Click_showXtick,lbl_kin});
%             h=bar([W_kin(:,i) W_EMGkin(17:end,i)]);
            h=bar([W_kin(:,i) ]);
         errorbar(W_kin(:,i),W_kin_std(:,i),'linestyle','none','color','black')

            set(h(1), 'FaceColor','blue');
%             set(h(2), 'FaceColor','green');
                xticks(1:length(lbl_kin))
                xtickangle(45)
                xticklabels(lbl_kin)
        subplot(N,3,(i-1)*3+3);  hold on;   set(gca,'YTick',[],'XTick',[]);
            plot(S_EMG(i,:),'color','red');
                x=[[0:1:99],fliplr([0:1:99])];
                y=[S_EMG(i,:)-S_EMG_std(i,:),...
                fliplr(S_EMG(i,:)+S_EMG_std(i,:))];
                s=fill(x,y,'k','EdgeColor','none');
                alpha(s,.1)
            plot(S_kin(i,:),'color','blue');
                x=[[0:1:99],fliplr([0:1:99])];
                y=[S_kin(i,:)-S_kin_std(i,:),...
                fliplr(S_kin(i,:)+S_kin_std(i,:))];
                s=fill(x,y,'k','EdgeColor','none');
                alpha(s,.1)
            
%             plot(S_EMGkin(i,:),'color','green');
            

    end
%     subplot(N,3,(i-1)*3+1);  hold on;    
%         xticks(1:length(lbl_EMG))
%         xtickangle(45)
%         xticklabels(lbl_EMG)    
%     subplot(N,3,(i-1)*3+2);  hold on;
%         xticks(1:length(lbl_kin))
%         xtickangle(45)
%         xticklabels(lbl_kin)
%     subplot(N,3,(i-1)*3+3);  hold on;
%         xticks(1:length(lbl_EMG))
%         xtickangle(45)
%         xticklabels(lbl_EMG)
           

end

function [S2, order] = mySort(S1)
    [m,id] = max(S1');
    [id_sorted, order] = sort(id);
    S2 = S1(order,:);
    
%     figure;
%     for i=1:size(S2,1)
%         subplot(6,1,i)
%         plot(S2(i,:))
%     end

end



