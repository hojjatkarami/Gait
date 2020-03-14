
% subject = {1,2,3,4,5};
% test = {[4 5 6],[4 5 6],[4 5 6],[4 5 6],[4 5 6]};
% SYN = 2:6;
% subject = {1,2};
% test = {[4],[4 5 6],[4 5 6],[4 5 6],[4 5 6]};
% SYN = 6;
% option.type = 'EMG';
% state='pre';
% option.side = 'Right';  % r: right side, l: left side; b: both side
compareFile = 'D:\University\MSc\Thesis project\Gait Data\Afternnmf\Compare\inter.mat';
load(compareFile);
option.type='kin2';
currentFigure_S=1006;
W={};        S={};        Event={};        Norm={};        name = {};        k=1;
%%
subject=length(inter);
SYN=6;
for i_sub = [2:8]
    
%     for i_test = test{i_sub}
        W{k} = (inter.(['s',num2str(i_sub)]).(option.type).syn(SYN).meanOfTypeW);
        S{k} = (inter.(['s',num2str(i_sub)]).(option.type).syn(SYN).meanOfTypeS);
        t{k} = (inter.(['s',num2str(i_sub)]).(option.type).syn(SYN).t);
        name{k} = [['s:',num2str(i_sub)],'_te:',...
                        '_syn:',num2str(SYN)];
        k=k+1;
        
%     end
end
[id,Wtype,W] = kmeanW(W, SYN);
n=length(id);

Stype = cell(1,SYN);
for i=1:SYN    
    for j=1:n
        x = find(id{j}==i);
        if isempty(x) || id{j}(x)==0
            continue;
        end

%         dim=length(S{j}(x,:));
%         [t2, y2] = rescale111(1:dim,S{j}(x,:),Event{j},Norm{j},1000);
         Stype{i}=[Stype{i} ; S{j}(x,:)];
    end
end
                

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
                '   Hip Flex','Hip Ext','Hip Add','Hip Abd','Hip Int Rot','Hip Ext Rot','Pelv Upwrd Obliq', 'Pelv Dw Obliq', 'Pelv Int Rot','Pelv Ext Rot'};
groupPartitionLineKM = [4 7 12 16 20 22 28];
P(1).EmgKin.Right(1).groupName=groupNameKM;

    