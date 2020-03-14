clc
clear
close all
nnmf_init
load('P.mat')

%% Settings
subject = [3:11];
SYN = 1:8;
option.side = 'Right';  % Right: right side, Left: left side;
vaf_th=0.9;    VAF_th = 0.95;
rsq_th=0.75;     RSQ_th = 0.8;
rsq_boot_th = 0.8;     vaf_boot_th=0.90;
lineStyle={'--','-'};
VAF=[];
RSQ=[];

label2 = {'Ank 1','Ank 2','Knee','Hip 1','Hip 2','Hip 3','Pelv 1','Pelv2'};
%% 100 - Goodness of fit for each criteria vs number of synergy

for i_sub = subject
    figure(100+i_sub)
    labelEMG = muscleName(P(i_sub).EMG.Right(1).muscleOrder);
    labelkin =angleName([19:26]);
    TrialNum =length(P(i_sub).EMG.Right);
    for i_trial = 1:TrialNum
        gof = P(i_sub).Synergy.EMG.(option.side)(i_trial).gof;
        
        subplot(2,2,1); hold on;    title('Global VAF');    grid on
%             set(gca,'ButtonDownFcn',{@plot_GoF_Click,M,rec,label});
            meanVal = (gof.BootStrap.vaf_bootstat_lb + gof.BootStrap.vaf_bootstat_ub)/2;
            err = gof.BootStrap.vaf_bootstat_ub-meanVal;
            errorbar(meanVal,err,'linestyle','none')
            plot(gof.VAF)
            line([SYN(1)-.5,SYN(end)+.5],[1 1]*VAF_th,'color','red')
            xlim([[SYN(1)-.5,SYN(end)+.5]]);     ylim([0.6 1]);
            xticks(SYN);
            xlabel('Number of Synergy');
            
        subplot(2,2,2); hold on;    title('Global R-square');    grid on
%             set(gca,'ButtonDownFcn',{@plot_GoF_Click,M,rec,label});
            meanVal = (gof.BootStrap.r2_bootstat_lb + gof.BootStrap.r2_bootstat_ub)/2;
            err = gof.BootStrap.r2_bootstat_ub-meanVal;
            errorbar(meanVal,err)
            plot(gof.RSQ)
            line([SYN(1)-.5,SYN(end)+.5],[1 1]*RSQ_th,'color','red')
            xlim([[SYN(1)-.5,SYN(end)+.5]]);     ylim([0.6 1]);
            xticks(SYN); 
       
            
            
    end 
    subplot(2,2,3); hold on;    title('Local VAF');    grid on
        set(gca,'ButtonDownFcn',{@plot_open});        temp=[];
        for i_syn = SYN
            for i_trial=1:TrialNum
                gof = P(i_sub).Synergy.EMG.(option.side)(i_trial).gof;
                temp(i_trial,:) = gof.vaf(i_syn,:);
            end
            errorbar(mean(temp),std(temp))
            [minVal, id] = min(gof.vaf(i_syn,:));    
            plot(gof.vaf(i_syn,:),'linestyle','-')
        end
        
        line([.5 16.5],[1 1]*vaf_th,'color','red')
        xticks(1:length(labelEMG))
        xtickangle(45)
        xticklabels(labelEMG)
    subplot(2,2,4); hold on;    title('Local R-square');    grid on
        set(gca,'ButtonDownFcn',{@plot_open});
        temp=[];
        for i_syn = SYN
            for i_trial=1:TrialNum
                gof = P(i_sub).Synergy.EMG.(option.side)(i_trial).gof;
                temp(i_trial,:) = gof.rsq(i_syn,:);
            end
            errorbar(mean(temp),std(temp),'DisplayName',['syn',num2str(i_syn)])
            [minVal, id] = min(gof.rsq(i_syn,:));    

        end
        
        line([.5 16.5],[1 1]*rsq_th,'color','red','DisplayName','Thereshold')
        lgd = legend('location','southoutside');
        lgd.NumColumns=3;
        xticks(1:length(labelEMG))
        xtickangle(45)
        xticklabels(labelEMG)
end
%% 200 - number of synergy based on each criterion vs subject
figure()
option.type = 'EMG';
for i_sub = subject
    TrialNum =length(P(i_sub).EMG.Right);
    for i_trial = 1:TrialNum
        number = P(i_sub).Synergy.(option.type).(option.side)(i_trial).number;
        
        subplot(2,3,1); hold on;    title(['Global VAF/',num2str(VAF_th)]);    grid on
            xlim([0 11]);   ylim([0 9]);    xticks([1:10]); yticks([1:8]);
            plot(i_sub+(rand-.5)/2,number.VAF,'*')
        subplot(2,3,2); hold on;    title(['Global R2/',num2str(RSQ_th)]);    grid on
            xlim([0 11]);   ylim([0 9]);    xticks([1:10]); yticks([1:8]);
            plot(i_sub+(rand-.5)/2,number.RSQ,'*')
        subplot(2,3,3); hold on;    title(['Local VAF/',num2str(vaf_th)]);    grid on
            xlim([0 11]);   ylim([0 10]);    xticks([1:10]); yticks([1:8]);
            plot(i_sub+(rand-.5)/2,number.vaf,'*')
        subplot(2,3,4); hold on;    title(['Local R2/',num2str(rsq_th)]);    grid on
            xlim([0 11]);   ylim([0 9]);    xticks([1:10]); yticks([1:8]);
            plot(i_sub+(rand-.5)/2,number.rsq,'*')
        subplot(2,3,5); hold on;    title(['BootStrap VAF/',num2str(vaf_boot_th)]);    grid on
            xlim([0 11]);   ylim([0 9]);    xticks([1:10]); yticks([1:8]);
            plot(i_sub+(rand-.5)/2,number.BootStrap_vaf,'*')
        subplot(2,3,6); hold on;    title(['BootStrap R2/',num2str(rsq_boot_th)]);    grid on
            xlim([0 11]);   ylim([0 9]);    xticks([1:10]); yticks([1:8]);
            plot(i_sub+(rand-.5)/2,number.BootStrap_r2,'*')
    end

    
end
%% local R2 and global VAF for all subjects
option.type = 'EMG';
temp1=zeros(10,10);
temp2 = zeros(10,10);
for i_sub = subject
    TrialNum =length(P(i_sub).EMG.Right);
    for i_trial = 1:TrialNum
        number = P(i_sub).Synergy.(option.type).(option.side)(i_trial).number;
        temp1(i_sub,i_trial) = number.rsq;
        temp2(i_sub,i_trial) = number.VAF;
    end
    
end
temp1 = sort(temp1,2);
temp2 = sort(temp2,2);
for i=2:10
   row1 = temp1(i,:);
   x1 = find(row1==0);
   if isempty(x1)
       x1=0;
   end
   row1 = row1(x1(end)+1:end);
   row1_m(i) = mean(row1,2);
   row1_std(i) = std(row1,0,2);
   
   row2 = temp2(i,:);
   x2 = find(row2==0);
   if isempty(x2)
       x2=0;
   end
   row2 = row2(x2(end)+1:end);
   row2_m(i) = mean(row2,2);
   row2_std(i) = std(row2,0,2);
end
figure(); 
subplot(1,2,1); hold on; 
    title(['Local criterion: r2 over ',num2str(rsq_th)]);
    grid on
    xlim([0 11]);   ylim([0 8]);    xticks([1:9]); yticks([1:7]);
    xlabel('Subject')
    ylabel('Number of synergy')
    bar(temp1(2:end,:),'b')
    errorbar(row1_m(2:end),row1_std(2:end),'r*','linestyle','none')
subplot(1,2,2); hold on; 
    title(['Global criterion: VAF over ',num2str(VAF_th)]);
    grid on
    xlim([0 11]);   ylim([0 8]);    xticks([1:9]); yticks([1:7]);
    xlabel('Subject')
    ylabel('Number of synergy')
    bar(temp2(2:end,:),'g')
    errorbar(row2_m(2:end),row2_std(2:end),'r*','linestyle','none')
    

%% Reconstruction
option.type = 'EMG';

subject = 9;
SYN=6;
for i_sub = subject
    TrialNum =length(P(i_sub).EMG.Right);
    for i_syn = SYN
        
        for i_trial = 4
            figure;
        for i_mus = 1:size(P(i_sub).(option.type).Right(i_trial).Mn,1)
            M = P(i_sub).(option.type).Right(i_trial).Mn;
            h = subplot(4,4,i_mus); hold on;
            h.Tag = num2str(i_mus);
                
                plot(M(i_mus,:))
%                 for i_syn=SYN
                   plot(P(i_sub).Synergy.(option.type).Right(i_trial).syn(i_syn).M_rec(i_mus,:)) 
%                 end
%                 gof = P(i_sub).Synergy.(option.type).Right(i_trial).gof.vaf(i_syn,i_mus);
%                 title([cell2mat(muscleName(i_mus)),' ',num2str(gof)]);
                set(gca,'ButtonDownFcn',{@Click_recon,P,i_sub,i_trial,i_mus,SYN})
        end
    
        end    
    end
end


%% Dyamic Plots