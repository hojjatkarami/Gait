

%% intra subject clustering
% subject = {1,,3,4,5};
nnmf_init;
currentFigure_S=1006;
subject = 5;
SYN = 6;
option.type = 'EMG'; %EMG, kin2, EmgKin
option.side = 'Right';  % Right: right side, Left: left side;
%%
 i_sub = subject  ;
    
        if strcmp(option.side, 'Right')
           TrialNum =length(P(i_sub).Trajectory.Right);
        else
           TrialNum =length(P(i_sub).Trajectory.Left);
        end

        for i_syn=SYN  
            W={};        S={}; tempS={}; t={};       Event={};        Norm={};        name = {};        k=1;
            if isempty(P(i_sub).Synergy(1).(option.type).(option.side).syn(i_syn).W_best)
                continue;                    
            end
            for i_trial = 1:TrialNum
                
                    
        %             label = Pf(i_sub).Test(i_test).(option.type)(1).sigName;
                W{k} = P(i_sub).Synergy(i_trial).(option.type)...
                        .(option.side).syn(i_syn).W_best;
                S{k} = P(i_sub).Synergy(i_trial).(option.type)...
                        .(option.side).syn(i_syn).S_best;
                 if strcmp(option.side, 'Right') 
                    f1 = P(i_sub).Events.Right(i_trial).PurtFrame;
                 else
                    f1 = P(i_sub).Events.Left(i_trial).PurtFrame;
                 end

                Event=f1;
                Event=double(Event);
                [t{k} S{k}] = NormS(S{k},Event);
               
%                 R2 = Pf(i_sub).Test(i_test).Synergy(i_trial).(option.type).(option.side).syn(i_syn).R2;
                    R2=1;
                name{k} = {['s:',num2str(subject),'_te:',...
                            '_tr:',num2str(i_trial),'_syn:',...
                            num2str(i_syn)];...
                            ['R2:',num2str(R2,'%.2f')]};
                k=k+1;
                
             
            end
                N = i_syn;
                [id,Wtype,W] = kmeanW(W, N);
                n=length(id);
                
                Stype = cell(1,N);
                for i=1:N    
                    for j=1:n
                        x = find(id{j}==i);
                        if isempty(x) || id{j}(x)==0
                            continue;
                        end

%                         dim=length(S{j}(x,:));
%                         [t2, y2] = rescale111(1:dim,S{j}(x,:),Event{j},Norm{j},1000);
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
                %     bad = find((meanOfType{i}-2*stdOfType{i})<0);
                %     meanOfType{i}(bad) = 0;
                end
                inter.(['s',num2str(subject)]).(option.type).syn(i_syn).meanOfTypeW = cell2mat(meanOfTypeW);
                inter.(['s',num2str(subject)]).(option.type).syn(i_syn).meanOfTypeS = cell2mat(meanOfTypeS');
                inter.(['s',num2str(subject)]).(option.type).syn(i_syn).t = t{1};

           end
%     end
    
  save('inter.mat','inter')
disp('Done')
 %% plot W
h1=figure(currentFigure_S-1);
supTitle = 'Intra Subject clustering';
label = P(1).(option.type).Right(1).groupName;
plotWtype(h1,W,Wtype,id,name,label,supTitle)

%% plot S
h2=figure(currentFigure_S);
plotStype(h2,t,S,Stype,id,name,label,supTitle)
%% Kinematics  and EMG matrix selection

groupNameKM = {'RTA','RPL','RSOL','RGC',...
                    'RRF','RMH','RVL',...
                    'RIP','RGMAX','RGMED','RAD','RTFL',...
                    'RIC','RLG','RRA','REO'...
                    'Ankle Dorsi Flex','Ankle Plant Flex',...
                   'Ankle Inv','Ankle Evr','Knee Flex','Knee Ext',...
                '   Hip Flex','Hip Ext','Hip Add','Hip Abd','Hip Int Rot','Hip Ext Rot'};
groupPartitionLineKM = [4 7 12 16 20 22 28];
P(1).EmgKin.Right(1).groupName=groupNameKM;

