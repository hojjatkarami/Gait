function [T, W, name, Wtype, id] = inter_period(T, type, subjects, period, trials, SYN, N, toPlot)


for i_sub = subjects
    if trials == 0
        trials = 1:length(T(i_sub).Synergy);
    end
    for i_trial = trials
        W={}; S={}; name={}; k=1;  Wstd={}; 
        for i_period = 1:length(period)
            W{k} = T(i_sub).Synergy(i_trial).(period{i_period}).Wmean;
            Wstd{k} = T(i_sub).Synergy(i_trial).(period{i_period}).Wstd;
            name{k} = ['s', num2str(i_sub)...
                '/t', num2str(i_trial)...
                '/', period{i_period}...
                '/r'];
            k=k+1;
        end
        
        [id,Wtype,Wmean, ~] = kmeanW1122(W,N,0.6);
        T(i_sub).Synergy(i_trial).inter_period.Wtype = Wtype;
        T(i_sub).Synergy(i_trial).inter_period.Wmean = Wmean;
        T(i_sub).Synergy(i_trial).inter_period.id = id;
        
        disp('done')
        if toPlot
            % plot W
            h1=figure();
            supTitle = 'Inter Subject clustering';
            label = muscleName(1:16);
            plotWtype1122(h1,W,Wstd,Wtype,id,name,label,supTitle)
        end
        
        
    end
end
disp('done')
%     if isempty(T(i_sub).Synergy.EMG_con)
%         continue;
%     end
%     W={}; S={}; name={}; k=1;
%     TrialNum =length(P(i_sub).EMG.Right);
%     for i_trial = 1:TrialNum
%         W{k} = P(i_sub).Synergy.EMG_con...
%             .(option.side)(i_trial).syn(N).W_best;
%         S{k} = P(i_sub).Synergy.EMG_con...
%             .(option.side)(i_trial).syn(N).S_best';
%         R2=1;
%         name{k} = {['s:',num2str(i_sub),'_te:',...
%             '_tr:',num2str(i_trial),'_syn:',...
%             num2str(N)];...
%             ['R2:',num2str(R2,'%.2f')]};
%         k=k+1;
%     end
%
%     [id,Stype,S] = kmeanW(S, N, 0.6);
%     for i=1:length(S)
%         S{i} = transpose(S{i});
%     end
%     for i=1:length(Stype)
%         Stype{i} = transpose(Stype{i});
%     end
%
%     for j=1:N
%         temp=[];
%         %         for i=1:size(Stype{j},1)
%         for i=1:length(S)
%             %             temp = [temp; Stype{j}(i,:)];
%             temp = [temp; S{i}(j,:)];
%
%         end
%         intra(i_sub).emg_con.syn(N).meanOfTypeS(j,:)=mean(temp);
%         intra(i_sub).emg_con.syn(N).stdOfTypeS(j,:)=std(temp);
%
%     end
%     intra(i_sub).emg_con.syn(N).meanOfTypeW = W{1};
%     %     inter(i_sub).emg_con.syn(N).meanOfTypeS = cell2mat(meanOfTypeS');
%     % plot W
%     %                     h1=figure(500+i_sub);
%     %                     supTitle = 'Intra Subject clustering';
%     %                     label = muscleName(P(i_sub).EMG.Right(1).muscleOrder);
%     %                     plotWtype(h1,W,Wtype,id,name,label,supTitle)
%     figure
%     for i=1:N
%         subplot(N,4,(i-1)*4+1); hold on;
%         bar(intra(i_sub).emg_con.syn(N).meanOfTypeW(:,i))
%         subplot(N,4,(i-1)*4+2); hold on;
%         bar(intra(i_sub).emg.syn(N).meanOfTypeW(:,i))
%         errorbar(intra(i_sub).emg.syn(N).meanOfTypeW(:,i),2*intra(i_sub).emg.syn(N).stdOfTypeW(:,i),'LineStyle','None')
%         subplot(N,4,(i-1)*4+3); hold on;
%         plot([0:1:99],intra(i_sub).emg_con.syn(N).meanOfTypeS(i,:))
%         x=[[0:1:99],fliplr([0:1:99])];
%         y=[intra(i_sub).emg_con.syn(N).meanOfTypeS(i,:)-intra(i_sub).emg_con.syn(N).stdOfTypeS(i,:),...
%             fliplr(intra(i_sub).emg_con.syn(N).meanOfTypeS(i,:)+intra(i_sub).emg_con.syn(N).stdOfTypeS(i,:))];
%         s=fill(x,y,'k','EdgeColor','none');
%         alpha(s,.1)
%         ylim([0 6])
%
%         subplot(N,4,(i-1)*4+4); hold on;
%         %         plot(intra(i_sub).emg.syn(N).meanOfTypeS(i,:))
%         plot([0:1:99],intra(i_sub).emg.syn(N).meanOfTypeS(i,:))
%         x=[[0:1:99],fliplr([0:1:99])];
%         y=[intra(i_sub).emg.syn(N).meanOfTypeS(i,:)-intra(i_sub).emg.syn(N).stdOfTypeS(i,:),...
%             fliplr(intra(i_sub).emg.syn(N).meanOfTypeS(i,:)+intra(i_sub).emg.syn(N).stdOfTypeS(i,:))];
%         s=fill(x,y,'k','EdgeColor','none');
%         alpha(s,.1)
%         ylim([0 6])
%
%     end
%     % plot S
%     %                     h2=figure(5500+i_sub);
%     %                     plotStype0720(h2,S,Stype,id,name,label,supTitle)
%
%
%
%
%
