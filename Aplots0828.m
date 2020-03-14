

%% select and compare from subjects

IDtoSave = 1;
showFigure = 0;

file = 'T1';
subject = [3];
% trial = 1:3;
period = {'period2'};
SYN=5;
type = 'W';
minCorr = 0.6;
N=SYN;

dataW={}; dataS={}; name={}; k=0;
eval(['T = ',file,';']);
for i_sub = subject
    trialNum = length(T(i_sub).ForcePlate);
    for i_trial = 1:trialNum
        for i_period = 1:length(period)
            if length(fieldnames(T(i_sub).Synergy(i_trial).(period{i_period})))==2
                continue;
            end
            for i_interval = 1:length(T(i_sub).Synergy(i_trial).(period{i_period}).interval)
                k=k+1;
                dataW{k} = T(i_sub).Synergy(i_trial).(period{i_period}).interval(i_interval).syn(SYN).W_best;
                dataS{k} = T(i_sub).Synergy(i_trial).(period{i_period}).interval(i_interval).syn(SYN).S_best;
                name{k} = 'r';
            end
        end
    end
end

[id,Wtype,dataW,Stype,dataS] = kmeanWS0828(dataW,dataS,type, N,minCorr);
n=length(id);



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

% saving
load('selectedID.mat')
selectedID(IDtoSave).Wtype = Wtype;
selectedID(IDtoSave).Stype = Stype;
selectedID(IDtoSave).meanOfTypeW =  cell2mat(meanOfTypeW);
selectedID(IDtoSave).meanOfTypeS =  cell2mat(meanOfTypeS');
selectedID(IDtoSave).stdOfTypeW =  cell2mat(stdOfTypeW);
selectedID(IDtoSave).stdOfTypeS =  cell2mat(stdOfTypeS');
save('selectedID.mat','selectedID')


if showFigure
    % plot W
    h1=figure();
    supTitle = 'Intra Subject clustering';
    label = muscleName(T(i_sub).EMG.Right(1).muscleOrder);
    plotmeanWtype(h1,dataW,Wtype,id,name,label,supTitle)
    % plot S
    %                     h2=figure();
    %                     plotStype0720(h2,dataS,Stype,id,name,label,supTitle)
    
end
disp('Done')
%% select and compare from id
ID = [14 19];
IDtoSave = 21;
showFigure = 1;

dataW={}; dataS={}; name={}; k=0;
load('selectedID.mat')

for i_ID = ID
    
    k=k+1;
    dataW{k} = selectedID(i_ID).meanOfTypeW;
    dataS{k} = selectedID(i_ID).meanOfTypeS;
    name{k} = 'r';
    
end

[id,Wtype,dataW,Stype,dataS] = kmeanWS0828(dataW,dataS,type, N,minCorr);
n=length(id);



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

% saving
load('selectedID.mat')
selectedID(IDtoSave).Wtype = Wtype;
selectedID(IDtoSave).Stype = Stype;
selectedID(IDtoSave).meanOfTypeW =  cell2mat(meanOfTypeW);
selectedID(IDtoSave).meanOfTypeS =  cell2mat(meanOfTypeS');
selectedID(IDtoSave).stdOfTypeW =  cell2mat(stdOfTypeW);
selectedID(IDtoSave).stdOfTypeS =  cell2mat(stdOfTypeS');
save('selectedID.mat','selectedID')


if showFigure
    % plot W
    h1=figure();
    supTitle = 'Intra Subject clustering';
    label = muscleName(T(i_sub).EMG.Right(1).muscleOrder);
    plotWtype(h1,dataW,Wtype,id,name,label,supTitle)
    % plot S
    %                     h2=figure();
    %                     plotStype0720(h2,dataS,Stype,id,name,label,supTitle)
    
end
disp('Done')


