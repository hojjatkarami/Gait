subject = [2 3 4 5 6 7 8 9 10];

%% 700 - EMG Channels

load('intra')
fig = figure(700);
k=1;
all = intra(1).emg.muscleOrder;
for i_mus=1:16
    subplot(4,4,i_mus); hold on;    ylabel(muscleName(all(i_mus)))
    set(gca,'ButtonDownFcn',{@plot_open});

    temp = [];
    dispName = {};
    for i_sub = subject
        if isempty(intra(i_sub).emg)
            continue;
        end
%         x = find(intra(i_sub).emg.muscleOrder==i_mus);
x = i_mus;
        if isempty(x)
            continue;
        end
        temp(i_sub,:) = intra(i_sub).emg.avg(x,:);
        dispName{i_sub} = num2str(i_sub);
        plot(temp(i_sub,:),'DisplayName',dispName{i_sub})
    end   
    if isempty(temp)
        continue;
    end
        nonZeroRows = any(temp,2);
        temp = temp(nonZeroRows,:);
%         intra.emg.avg = mean(temp);
%         intra.emg.std = std(temp);
        plot(mean(temp),'linewidth',1.5,'linestyle','-.')
%             m = mean(temp_purtNorm); myLine('v',m(2:end-1),[0 1])
        xticks([])
        
%         plot(temp','DisplayName',dispName)
        x=[[1:100],fliplr([1:1:100])];
        y=[mean(temp)-std(temp),...
        fliplr( mean(temp)+std(temp))];
        s=fill(x,y,'k','EdgeColor','none');
        alpha(s,.1)
end
savefig(fig,'output/intra-EMG');
saveas(fig,'output/intra-EMG','svg');
save('intra.mat','intra')
%% 800 - Kinematic Angles
load('intra')
fig=figure(701);
for i_kin=1:9
    subplot(3,3,i_kin); hold on;    ylabel(label2{i_kin})
    set(gca,'ButtonDownFcn',{@plot_open});

    temp = [];

    for i_sub = subject
        if isempty(intra(i_sub).kin)
            continue;
        end
        temp(i_sub,:) = intra(i_sub).kin.avg(i_kin,:);
        dispName{i_sub} = num2str(i_sub);
        plot(temp(i_sub,:),'DisplayName',dispName{i_sub})
    end   
        nonZeroRows = any(temp,2);
        temp = temp(nonZeroRows,:);
%         intra.kin.avg = mean(temp);
%         intra.kin.std = std(temp);
        plot(mean(temp),'linewidth',1.5,'linestyle','-.')
%             m = mean(temp_purtNorm); myLine('v',m(2:end-1),[0 1])
        xticks([])
%         plot((temp'))
        x=[[1:100],fliplr([1:1:100])];
        y=[mean(temp)-std(temp),...
        fliplr(mean(temp)+std(temp))];
        s=fill(x,y,'k','EdgeColor','none');
        alpha(s,.1)
end
savefig(fig,'output/intra-kin');
saveas(fig,'output/intra-kin','svg')
save('intra.mat','intra')

%% 900 - EMG SYNERGY 
N=6;
W={}; S={}; name={}; k=1;
load('intra')
load('inter')
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

for i_sub = subject
    if isempty(intra(i_sub).emg)
        continue;
    end
%     for i_test = test{i_sub}
        for i=1:length(sameMuscles)
            x = find(intra(i_sub).emg.muscleOrder==sameMuscles(i));
            
                W{k}(i,:) = intra(i_sub).emg.syn(N).meanOfTypeW(x,:);
            
        end
%         W{k} = intra(i_sub).emg.syn(N).meanOfTypeW;
        
        S{k} = intra(i_sub).emg.syn(N).meanOfTypeS;
        
        name{k} = [['s:',num2str(i_sub)],'_te:',...
                        '_syn:',num2str(N)];
        k=k+1;
        
%     end
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
        
save('inter.mat','inter')

%% 9000 - Kinematic SYNERGY 
N=6;
W={}; S={}; name={}; k=1;
load('intra')
load('inter')

for i_sub = subject
    if isempty(intra(i_sub).kin)
        continue;
    end

    
%     for i_test = test{i_sub}
%         for i=1:length(allMuscles)
%             x = find(intra(i_sub).emg.muscleOrder==allMuscles(i));
%             if isempty(x)
%                 W{k}(i,:)=0;
%             else
%                 W{k}(i,:) = intra(i_sub).emg.syn(N).meanOfTypeW(x,:);
%             end
%         end
        W{k} = intra(i_sub).kin.syn(N).meanOfTypeW;
        
        S{k} = intra(i_sub).kin.syn(N).meanOfTypeS;
        
        name{k} = [['s:',num2str(i_sub)],'_te:',...
                        '_syn:',num2str(N)];
        k=k+1;
        
%     end
end
[id,Wtype,W] = kmeanW(W,N);
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
    
    inter.kin.syn(N).meanOfTypeW = cell2mat(meanOfTypeW);
    inter.kin.syn(N).meanOfTypeS = cell2mat(meanOfTypeS');
    inter.kin.syn(N).stdOfTypeW = cell2mat(stdOfTypeW);
    inter.kin.syn(N).stdOfTypeS = cell2mat(stdOfTypeS');
        % plot W
        h1=figure();
        supTitle = 'Inter Subject clustering';
        label = angleName(1:16);
        plotWtype(h1,W,Wtype,id,name,label,supTitle)
        % plot S
        h2=figure();
        plotStype0720(h2,S,Stype,id,name,label,supTitle)
        
%         savefig(h1,'output/intra-kin-syn-W');
%         saveas(h1,'output/intra-kin-syn-W','svg');
%         savefig(h2,'output/intra-kin-syn-S');
%         saveas(h2,'output/intra-kin-syn-S','svg');
save('inter.mat','inter')
%% EMGkin SYNERGY
N=6;
W={}; S={}; name={}; k=1;
load('intra')
load('inter')

for i_sub = subject
    if isempty(intra(i_sub).kin)
        continue;
    end

    
%     for i_test = test{i_sub}
%         for i=1:length(allMuscles)
%             x = find(intra(i_sub).emg.muscleOrder==allMuscles(i));
%             if isempty(x)
%                 W{k}(i,:)=0;
%             else
%                 W{k}(i,:) = intra(i_sub).emg.syn(N).meanOfTypeW(x,:);
%             end
%         end
        W{k} = intra(i_sub).emgkin.syn(N).meanOfTypeW;
        
        S{k} = intra(i_sub).emgkin.syn(N).meanOfTypeS;
        
        name{k} = [['s:',num2str(i_sub)],'_te:',...
                        '_syn:',num2str(N)];
        k=k+1;
        
%     end
end
[id,Wtype,W] = kmeanW(W,N);
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
    
    inter.emgkin.syn(N).meanOfTypeW = cell2mat(meanOfTypeW);
    inter.emgkin.syn(N).meanOfTypeS = cell2mat(meanOfTypeS');
    inter.emgkin.syn(N).stdOfTypeW = cell2mat(stdOfTypeW);
    inter.emgkin.syn(N).stdOfTypeS = cell2mat(stdOfTypeS');
        % plot W
        h1=figure();
        supTitle = 'Inter Subject clustering';
        label = angleName(1:16);
        plotWtype(h1,W,Wtype,id,name,label,supTitle)
        % plot S
        h2=figure();
        plotStype0720(h2,S,Stype,id,name,label,supTitle)
        
%         savefig(h1,'output/intra-kin-syn-W');
%         saveas(h1,'output/intra-kin-syn-W','svg');
%         savefig(h2,'output/intra-kin-syn-S');
%         saveas(h2,'output/intra-kin-syn-S','svg');

save('inter.mat','inter')

%% Classify EMG and Kinematic
load('intra')
    figure(1000)
    
    for i_syn=1:6
        tempW_EMG=[];   tempS_EMG=[];
        tempW_kin=[];   tempS_kin=[];
        for i_sub = subject
            if isempty(intra(i_sub).kin)
                continue;
            end
            tempW_EMG = [tempW_EMG, intra(i_sub).emg.syn(6).meanOfTypeW(:,i_syn)]; 
            tempS_EMG = [tempS_EMG; intra(i_sub).emg.syn(6).meanOfTypeS(i_syn,:)]; 
            tempW_kin = [tempW_kin, intra(i_sub).kin.syn(6).meanOfTypeW(:,i_syn)]; 
            tempS_kin = [tempS_kin; intra(i_sub).kin.syn(6).meanOfTypeS(i_syn,:)]; 
        end
        subplot(6,3,3*i_syn-2); hold on;
            bar(tempW_EMG)
%             bar(mean(tempW_EMG'));
%             errorbar( mean(tempW_EMG'),  2*std(tempW_EMG'),'LineStyle','none','Color','black')

        subplot(6,3,3*i_syn-1); hold on;
            bar(mean(tempW_kin'));
        subplot(6,3,3*i_syn); hold on;
%             plot(mean(tempS_EMG));
%             plot(mean(tempS_kin));
            plot((tempS_EMG'));

        
    end
        
