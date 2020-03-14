
% subject = [5];
label2 = angleName(19:27);
toSave = 0;
N=6;
option.side = 'Right';  % Right: right side, Left: left side;
%% 300 - EMG Channels
load('intra')


for i_sub = subject
    if isempty(T(i_sub).EMG)
        continue;
    end
    fig = figure(300+i_sub);
    TrialNum =length(T(i_sub).EMG.(option.side));
    
    musOrd = T(i_sub).EMG.(option.side)(1).muscleOrder;
    musNum = length(musOrd);
    musLbl = muscleName(musOrd);
    intra(i_sub).emg.muscleOrder = musOrd;
    for i_mus=1:musNum
        subplot(8,2,i_mus); hold on;    ylabel(musLbl{i_mus})
        set(gca,'ButtonDownFcn',{@plot_open});

        temp = [];
        temp_purtNorm = [];
        for i_trial = 1:TrialNum
            temp(i_trial,:) = T(i_sub).EMG.(option.side)(i_trial).Mn(i_mus,:);
            temp_purtNorm(i_trial,:) = T(i_sub).Events.(option.side)(i_trial).PurtFrameNorm; 
%             tNorm=T(i_sub).Events.Right(i_trial).tNorm;
%             plot(tNorm,temp(i_trial,:))
        end   
            intra(i_sub).emg.avg(i_mus,:) = mean(temp);
            intra(i_sub).emg.std(i_mus,:) = std(temp);
            plot(intra(i_sub).emg.avg(i_mus,:),'linewidth',1.5,'linestyle','-.')
            m = mean(temp_purtNorm); myLine('v',m(2:end-1),[0 1])
            xticks([])
            
            plot(temp')
            x=[[1:100],fliplr([1:1:100])];
            y=[intra(i_sub).emg.avg(i_mus,:)-intra(i_sub).emg.std(i_mus,:),...
            fliplr(intra(i_sub).emg.avg(i_mus,:)+intra(i_sub).emg.std(i_mus,:))];
            s=fill(x,y,'k','EdgeColor','none');
            alpha(s,.1)
    end
    xticks(m)
    xtickangle(45)
    xticklabels({'RFC1','LTO','LFC','RTO','RFC2'})
    if toSave
        savefig(fig,['output/intra-EMG',num2str(i_sub)]);
        saveas(fig,['output/intra-EMG',num2str(i_sub)],'jpg');
        close all
    end
end

save('intra.mat','intra')


%% 400 - Kinematic Angles
load('intra')
for i_sub = subject
    if isempty(T(i_sub).KIN)
        continue;
    end
    fig = figure(400+i_sub);
    TrialNum =length(T(i_sub).KIN.(option.side));
    for i_kin=1:9
        subplot(3,3,i_kin); hold on;    ylabel(label2{i_kin})
        set(gca,'ButtonDownFcn',{@plot_open});

        temp = [];
        temp_purtNorm = [];
        for i_trial = 1:TrialNum
            temp(i_trial,:) = T(i_sub).KIN.(option.side)(i_trial).Mn(i_kin,:);
            temp_purtNorm(i_trial,:) = T(i_sub).Events.(option.side)(i_trial).PurtFrameNorm; 
%             tNorm=T(i_sub).Events.Right(i_trial).tNorm;
%             plot(tNorm,temp(i_trial,:))
        end   
            intra(i_sub).kin.avg(i_kin,:) = mean(temp);
            intra(i_sub).kin.std(i_kin,:) = std(temp);
            plot(intra(i_sub).kin.avg(i_kin,:),'linewidth',1.5,'linestyle','-.')
            m = mean(temp_purtNorm); %myLine('v',m(2:end-1),[0 1])
            xticks([])
            
            plot(temp')
            x=[[1:100],fliplr([1:1:100])];
            y=[intra(i_sub).kin.avg(i_kin,:)-intra(i_sub).kin.std(i_kin,:),...
            fliplr(intra(i_sub).kin.avg(i_kin,:)+intra(i_sub).kin.std(i_kin,:))];
            s=fill(x,y,'k','EdgeColor','none');
            alpha(s,.1)
    end
    xticks(m)
    xtickangle(45)
    xticklabels({'RFC1','LTO','LFC','RTO','RFC2'})
    if toSave
        savefig(fig,['output/intra-kin',num2str(i_sub)]);
        saveas(fig,['output/intra-kin',num2str(i_sub)],'jpg');
        close all
    end
end
save('intra.mat','intra')

%% 500 - EMG SYNERGY

load('intra')
for i_sub = subject
    if isempty(T(i_sub).Synergy.EMG)
        continue;
    end
    W={}; S={}; name={}; k=1;    
    TrialNum =length(T(i_sub).EMG.Right);
    for i_trial = 1:TrialNum
        W{k} = T(i_sub).Synergy.EMG...
                        .(option.side)(i_trial).syn(N).W_best;
        S{k} = T(i_sub).Synergy.EMG...
                        .(option.side)(i_trial).syn(N).S_best;
        R2=1;
        name{k} = {['s',num2str(i_sub),...
                            '/t',num2str(i_trial),'/syn',...
                            num2str(N)];...
                            ['R2:',num2str(R2,'%.2f')]};
        k=k+1;
    end   
    
    [id,Wtype,W] = kmeanW(W, N,0.6);
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
    
    intra(i_sub).emg.syn(N).meanOfTypeW = cell2mat(meanOfTypeW);
    intra(i_sub).emg.syn(N).meanOfTypeS = cell2mat(meanOfTypeS');
    intra(i_sub).emg.syn(N).stdOfTypeW = cell2mat(stdOfTypeW);
    intra(i_sub).emg.syn(N).stdOfTypeS = cell2mat(stdOfTypeS');

        
                % plot W
%                     h1=figure();
%                     supTitle = 'Intra Subject clustering';
%                     label = muscleName(T(i_sub).EMG.Right(1).muscleOrder);
%                     plotWtype(h1,W,Wtype,id,name,label,supTitle)
                % plot S
%                     h2=figure();
%                     plotStype0720(h2,S,Stype,id,name,label,supTitle)
%                 if toSave
%                     savefig(h1,['output/intra-EMG-syn-W',num2str(i_sub)]);
%                     saveas(h1,['output/intra-EMG-syn-W',num2str(i_sub)],'jpg');
%                     savefig(h2,['output/intra-EMG-syn-S',num2str(i_sub)]);
%                     saveas(h2,['output/intra-EMG-syn-S',num2str(i_sub)],'jpg');
%                     close all
%                 end
end
save('intra.mat','intra')

        
%% 600 - Kinematic SYNERGY
load('intra')
for i_sub = subject
    if isempty(T(i_sub).Synergy.kin)
        continue;
    end
    W={}; S={}; Name={}; k=1;    
    TrialNum =length(T(i_sub).EMG.Right);
    for i_trial = 1:TrialNum
        W{k} = T(i_sub).Synergy.kin...
                        .(option.side)(i_trial).syn(N).W_best;
        S{k} = T(i_sub).Synergy.kin...
                        .(option.side)(i_trial).syn(N).S_best;
        R2=1;
        name{k} = {['s:',num2str(i_sub),'_te:',...
                            '_tr:',num2str(i_trial),'_syn:',...
                            num2str(N)];...
                            ['R2:',num2str(R2,'%.2f')]};
        k=k+1;
    end   
    
    [id,Wtype,W] = kmeanW(W, N,0.6);
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
    
    intra(i_sub).kin.syn(N).meanOfTypeW = cell2mat(meanOfTypeW);
    intra(i_sub).kin.syn(N).meanOfTypeS = cell2mat(meanOfTypeS');
    intra(i_sub).kin.syn(N).stdOfTypeW = cell2mat(stdOfTypeW);
    intra(i_sub).kin.syn(N).stdOfTypeS = cell2mat(stdOfTypeS');
    
                % plot W
                    h1=figure(600+i_sub);
                    supTitle = 'Intra Subject clustering';
                    label = angleName(1:16);
                    plotWtype(h1,W,Wtype,id,name,label,supTitle)
                % plot S
                    h2=figure(650+i_sub);
                    plotStype0720(h2,S,Stype,id,name,label,supTitle)
%                 if toSave
%                     savefig(h1,['output/intra-kin-syn-W',num2str(i_sub)]);
%                     saveas(h1,['output/intra-kin-syn-W',num2str(i_sub)],'jpg');
%                     savefig(h2,['output/intra-kin-syn-S',num2str(i_sub)]);
%                     saveas(h2,['output/intra-kin-syn-S',num2str(i_sub)],'jpg');
%                     close all
%                 end
end
save('intra.mat','intra')
%% EMGkin

load('intra')
for i_sub = subject
    if isempty(T(i_sub).Synergy.EMG)
        continue;
    end
    W={}; S={}; name={}; k=1;    
    TrialNum =length(T(i_sub).EMG.Right);
    for i_trial = 1:TrialNum
        W{k} = T(i_sub).Synergy.EMGkin...
                        .(option.side)(i_trial).syn(N).W_best;
        S{k} = T(i_sub).Synergy.EMGkin...
                        .(option.side)(i_trial).syn(N).S_best;
        R2=1;
        name{k} = {['s',num2str(i_sub),...
                            '/t',num2str(i_trial),'/syn',...
                            num2str(N)];...
                            ['R2:',num2str(R2,'%.2f')]};
        k=k+1;
    end   
    
    [id,Wtype,W] = kmeanW(W, N,0.6);
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
    
    intra(i_sub).emgkin.syn(N).meanOfTypeW = cell2mat(meanOfTypeW);
    intra(i_sub).emgkin.syn(N).meanOfTypeS = cell2mat(meanOfTypeS');
    intra(i_sub).emgkin.syn(N).stdOfTypeW = cell2mat(stdOfTypeW);
    intra(i_sub).emgkin.syn(N).stdOfTypeS = cell2mat(stdOfTypeS');
                % plot W
%                     h1=figure();
%                     supTitle = 'Intra Subject clustering';
%                     label = muscleName(T(i_sub).EMG.Right(1).muscleOrder);
%                     plotWtype(h1,W,Wtype,id,name,label,supTitle)
                % plot S
%                     h2=figure();
%                     plotStype0720(h2,S,Stype,id,name,label,supTitle)
%                 if toSave
%                     savefig(h1,['output/intra-EMG-syn-W',num2str(i_sub)]);
%                     saveas(h1,['output/intra-EMG-syn-W',num2str(i_sub)],'jpg');
%                     savefig(h2,['output/intra-EMG-syn-S',num2str(i_sub)]);
%                     saveas(h2,['output/intra-EMG-syn-S',num2str(i_sub)],'jpg');
%                     close all
%                 end
end
save('intra.mat','intra')

%% Classify EMG and Kinematic
load('intra')
for i_sub = subject
    % sort EMG
    [m,i] = max(intra(i_sub).emg.syn(6).meanOfTypeS');    
    [i_sorted, order] = sort(i);
    intra(i_sub).emg.syn(6).meanOfTypeS= intra(i_sub).emg.syn(6).meanOfTypeS(order,:);
    intra(i_sub).emg.syn(6).meanOfTypeW= intra(i_sub).emg.syn(6).meanOfTypeW(:,order);

    % sort kin
    [m,i] = max(intra(i_sub).kin.syn(6).meanOfTypeS');    
    [i_sorted, order] = sort(i);
    intra(i_sub).kin.syn(6).meanOfTypeS= intra(i_sub).kin.syn(6).meanOfTypeS(order,:);
    intra(i_sub).kin.syn(6).meanOfTypeW= intra(i_sub).kin.syn(6).meanOfTypeW(:,order);

    figure(i_sub)
    for i_syn=1:6
        
        subplot(6,3,3*i_syn-2); hold on;
            bar(intra(i_sub).emg.syn(6).meanOfTypeW(:,i_syn));
        subplot(6,3,3*i_syn-1); hold on;
            bar(intra(i_sub).kin.syn(6).meanOfTypeW(:,i_syn));
        subplot(6,3,3*i_syn); hold on;
            plot(intra(i_sub).emg.syn(6).meanOfTypeS(i_syn,:));
            plot(intra(i_sub).kin.syn(6).meanOfTypeS(i_syn,:));
    end
end
save('intra.mat','intra')
%% 300 - clustering all synergies(1:N) for one subject
% close all
n_common = 6;
N=8;
i_trial = 1;
subject = 5;
option.type = 'kin';
load('intra')
for i_sub = subject
    if isempty(T(i_sub).Synergy.(option.type))
        continue;
    end
    TrialNum =length(T(i_sub).EMG.Right);

    for i_trial=1
    W={}; S={}; name={}; k=1;    
    for i_syn = 1:8
        W{k} = T(i_sub).Synergy.(option.type)...
                        .(option.side)(i_trial).syn(i_syn).W_best;
        S{k} = T(i_sub).Synergy.(option.type)...
                        .(option.side)(i_trial).syn(i_syn).S_best;
        R2=1;
        name{k} = {['s:',num2str(i_sub),'_te:',...
                            '_tr:',num2str(i_trial),'_syn:',...
                            num2str(i_syn)];...
                            ['R2:',num2str(R2,'%.2f')]};
        k=k+1;
    end   
    
    [id,Wtype,W] = kmeanW(W, N);
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

    meanOfTypeW = {1,N};
    stdOfTypeW = {1,N};
    meanOfTypeS = {1,N};
    stdOfTypeS = {1,N};    
    for i=1:N
        meanOfTypeW{i} =  mean(Wtype{i},2);
        stdOfTypeW{i} = std(Wtype{i},0,2);
        meanOfTypeS{i} =  mean(Stype{i},1);
        stdOfTypeS{i} = std(Stype{i},0,1);
    end
    
    intra(i_sub).emg.syn(N).meanOfTypeW = cell2mat(meanOfTypeW);
    intra(i_sub).emg.syn(N).meanOfTypeS = cell2mat(meanOfTypeS');
                % plot W
%                     h1=figure();
                    supTitle = 'Intra Subject clustering';
                    supTitle = num2str(T(i_sub).Synergy.(option.type).Right(i_trial).number.rsq);
                    if strcmp(option.type,'EMG')
                        label = muscleName(T(i_sub).EMG.Right(1).muscleOrder);
                    else
                        label = angleName(1:16);
                    end
%                     plotWtype(h1,W,Wtype,id,name,label,supTitle)
                % plot S
                    h2=figure();
                    plotStype0720(h2,S,Stype,id,name,label,supTitle)
                if toSave
                    savefig(h1,['output/intra-EMG-syn-W',num2str(i_sub)]);
                    saveas(h1,['output/intra-EMG-syn-W',num2str(i_sub)],'jpg');
                    savefig(h2,['output/intra-EMG-syn-S',num2str(i_sub)]);
                    saveas(h2,['output/intra-EMG-syn-S',num2str(i_sub)],'jpg');
                    close all
                end
    end
end
save('intra.mat','intra')
