%%
clear
close all
clc
set(0, 'defaulttextinterpreter', 'latex');
set(0, 'defaultlegendinterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex')


%% Settings
path='D:\University\MSc\Thesis project\Gait Data\ProcessdataMat\';

sub=5;
% trial=1:5;

syn=5;

matName = 'P';
load([path,matName]);
% eval(['P',matName]);

option.devideByStdFlag = 1;
option.VAF_th = 0.9;  % total VAF thereshold
option.vaf_th = 0.85;  % VAF thereshold for each muscle
option.corr_th = 0.8;
option.rsq_th = 0.6;
option.RSQ_th = 0.9;
option.rsqCI_th = 0.99;
option.vafCI_th = 0.95;
option.n_bootstrap = 1000;
option.condition = 'none';  % VAF, vaf, corr, rsq, vafCI, rsqCI, none
option.type = 'emg';    % emg, kin or emgkin
option.side = 'r';  % r: right side, l: left side; b: both side
option.rep1 = 1;
option.rep2=1;


%% NNMF stage 1
tic
clc
% sub=5;
if option.side=='r'
    trial=1:length(P(sub).Trajectory.Right);
else
    trial=1:length(P(sub).Trajectory.Left);
end
for i_trial=trial
        h_trial =P(sub);
        [M, std_val] = preNNMF(h_trial,i_trial,option);
        for i_syn = syn
            disp([num2str(sub),', ',num2str(i_trial),', ',num2str(i_syn),', ' ])
            [W_best, S_best] = MyNNMF(M, std_val, i_syn, option);
            if option.side == 'r'
               P(sub).EMG.Right(i_trial).Syn(i_syn).W_best = W_best;
               P(sub).EMG.Right(i_trial).Syn(i_syn).S_best = S_best;
            else
               P(sub).EMG.Left(i_trial).Syn(i_syn).W_best = W_best;
               P(sub).EMG.Left(i_trial).Syn(i_syn).S_best = S_best; 
            end
        end
end
toc
% save all data to mat file
% eval([matName,'=data;']);
% save([path,matName],matName);
disp('nnmf done!')

%% structure clustering and selection for EMG
% close all
path='D:\University\MSc\Thesis project\Gait Data\Compair';
matName = 'P';
% sub=5;
% trial=1:5;
syn=5;

 name=[matName,'_Sub',num2str(sub),'_syn',num2str(syn)];

% load([path,matName]);
% eval(['data=',matName]);
% trial = 1:size(data.test(test).trial,2)


cellData1 = {};
cellData2 = {};
%% Right Side
clc
% sub=6;
    for i_trial = 1:length(P(sub).Trajectory.Right)
        for i_syn = syn
            purtFrame = P(sub).Events.Right(i_trial).PurtFrame;
            interval = length(purtFrame)-1;
            x1 = purtFrame(1):purtFrame(end);
            x1=double(x1);
            Y=[];
            X=[];
            for i=1:i_syn
                y1 = P(sub).EMG.Right(i_trial).Syn(i_syn).S_best(i,:);
                x1=double(x1);
                y1=double(y1);
                [x2, y2] = rescale11(x1,y1,purtFrame,[0:100/interval:100],1000);
                Y=[Y;y2];
                X=[X;x2];
            end
            cellData1{i_syn,i_trial} = P(sub).EMG.Right(i_trial).Syn(i_syn).W_best;        
            cellData2{i_syn,i_trial}.Y = Y;
            cellData2{i_syn,i_trial}.X = X;
        end
    end
 
    Rnew1 = cellData1(~cellfun(@isempty, cellData1));  % remove empty cells    
    Rnew1 = reshape(Rnew1,1,[]);  %convert to 1-D cell
    Rnew2 = cellData2(~cellfun(@isempty, cellData2));  % remove empty cells    
    Rnew2 = reshape(Rnew2,1,[]);  %convert to 1-D cell
    label =P(sub).EMG.Right(i_trial).groupName;
    partitionLine = P(sub).EMG.Right(i_trial).groupPartitionLine;
%     group = plotSyn(Rnew1,Rnew2,label,partitionLine);
        N=max(syn);
        group = kmean1(N,Rnew1,Rnew2,label,partitionLine,name);
    
    eval([name,'=group']);
     save(['D:\University\MSc\Thesis project\Gait Data\Compair\',name],name)
    
 %% Left Side   
 clc
sub=5;
    for i_trial = 1:length(P(sub).Trajectory.Left)
        for i_syn = syn
            purtFrame = P(sub).Events.Left(i_trial).PurtFrame;
            interval = length(purtFrame)-1;
            x1 = purtFrame(1):purtFrame(end);
            x1=double(x1);
            Y=[];
            X=[];
            for i=1:i_syn
                y1 = P(sub).EMG.Left(i_trial).Syn(i_syn).S_best(i,:);
                x1=double(x1);
                y1=double(y1);
                [x2, y2] = rescale11(x1,y1,purtFrame,[0:100/interval:100],1000);
                Y=[Y;y2];
                X=[X;x2];
            end
            cellData1{i_syn,i_trial} = P(sub).EMG.Left(i_trial).Syn(i_syn).W_best;        
            cellData2{i_syn,i_trial}.Y = Y;
            cellData2{i_syn,i_trial}.X = X;
        end
    end
 
    Rnew1 = cellData1(~cellfun(@isempty, cellData1));  % remove empty cells    
    Rnew1 = reshape(Rnew1,1,[]);  %convert to 1-D cell
    Rnew2 = cellData2(~cellfun(@isempty, cellData2));  % remove empty cells    
    Rnew2 = reshape(Rnew2,1,[]);  %convert to 1-D cell
    label =P(sub).EMG.Left(i_trial).groupName;
    partitionLine = P(sub).EMG.Left(i_trial).groupPartitionLine;
%     group = plotSyn(Rnew1,Rnew2,label,partitionLine);
        N=max(syn);
        group = kmean1(N,Rnew1,Rnew2,label,partitionLine,name);
    
%     eval([name,'=group']);
%     save(['Data\compare\',name],name)
 
%% compare tool
clc
subject = {'P_Sub5_syn5','P_Sub6_syn5','P_Sub7_syn5'};
n = length(subject);
% rgb = maxdistcolor(2*n,@srgb_to_Lab);
cellData = {};
rep=[];
for i=1:n
    matName = subject{i};
    load(['D:\University\MSc\Thesis project\Gait Data\Compair\',matName]);
    eval(['dataCmp=',matName]);
    
    cellData{i}(:,1) = mean(dataCmp.Wtype{1},2);
    rep(i,1) = size(dataCmp.Wtype{1},2);
    for j=2:dataCmp.WtypeNo
        cellData{i} = [cellData{i}, mean(dataCmp.Wtype{j},2)];
        rep(i,j) = size(dataCmp.Wtype{j},2);
        
    end


%     for j=1:dataCmp.typeNo
%         subplot(dataCmp.typeNo,n,(j-1)*n+i); hold on;
%         rep = size(dataCmp.type{j},2);
%         title(['rep:',num2str(rep)]);
%         bar(mean(dataCmp.type{j},2),'FaceColor',rgb(i,:))
%         errorbar(mean(dataCmp.type{j},2), std(dataCmp.type{j},0,2),'LineStyle','none','Color',rgb(n+i,:))
%         set(gca,'XTick',[],'YTick',[],'ButtonDownFcn',@createnew_fig)
%     end
end
    label = P(sub).EMG.Left(i_trial).groupName;
    partitionLine = P(sub).EMG.Left(i_trial).groupPartitionLine;
    cmpSyn(cellData,rep,label,partitionLine);

%% plot W and S
test=5;
trial=2;
syn=4;
figure;
for i_test = test    
    i_test
    cellData = {};
    for i_trial = trial
        for i_syn=1:syn
            % plot W
            subplot(syn,2,(i_syn-1)*2+1); hold on;
            bar(data.test(i_test).trial(i_trial).emg.syn(syn).W_best(:,i_syn))
            label = data.test(i_test).trial(i_trial).emg.groupName;
            partitionLine = data.test(i_test).trial(i_trial).emg.groupPartitionLine;
            for ii = 1:length(partitionLine)
                plot(partitionLine(ii)*[1 1]+0.5 , [0 1],'color','black') 
            end
            xticks(1:16)
            xtickangle(45)
            xticklabels(label)
            % plot S
            subplot(syn,2,(i_syn-1)*2+2); hold on;
%             set(gca,'XTick',[],'YTick',[],'ButtonDownFcn',@createnew_fig)        
            purtFrame = floor(data.test(i_test).trial(i_trial).purtFrame/2)+1;
            interval = length(purtFrame)-1;
            x1 = purtFrame(1):purtFrame(end);
            y1 = data.test(i_test).trial(i_trial).emg.syn(syn).S_best(i_syn,x1);
            [x2, y2] = rescale11(x1,y1,purtFrame,[0:100/interval:100],1000);
            plot(x2,y2)
            l=[0:100/interval:100];
            for i=1:length(purtFrame)
                line([1 1]*l(i),[0 1],'color','black')
                
            end
        end
    end
end
%% plot goodness of fit of trials
test=3;
trial=1:5;
syn=2:6;
for i_test = test    
    i_test
    cellData = {};
    for i_trial = trial
        for i_syn=syn
            ref = data.test(i_test).trial(i_trial).emg.M_R;
            rec = data.test(i_test).trial(i_trial).emg.syn(i_syn).W_best * data.test(i_test).trial(i_trial).emg.syn(i_syn).S_best;
            cellData{1,i_trial,i_syn} = ref;
            cellData{2,i_trial,i_syn} = rec; 
        end
    end
    plotGoF(cellData,syn)
end


%% plot Reconstructed vs Orginal signal
i_test=6;
i_trial=1;
syn=2:5;
figure('units','normalized','outerposition',[0 0 1 1]);
for j=1:16
    subplot(4,4,j); hold on;
    set(gca,'XTick',[],'YTick',[],'ButtonDownFcn',@createnew_fig)
    title(data.test(i_test).trial(i_trial).emg.chName{j})
    ref = data.test(i_test).trial(i_trial).emg.M_R(j,:);
    plot(ref,'DisplayName','Original')
    for i_syn=syn
        rec = data.test(i_test).trial(i_trial).emg.syn(i_syn).W_best * data.test(i_test).trial(i_trial).emg.syn(i_syn).S_best;
        plot(rec(j,:),'DisplayName',['SynNo. ',num2str(i_syn)])   
    end
    
end
legend

%% plot kinematic angles
test=[8];
eval(horzcat('data=',filename{i_test},';'));
dz = data.kinematics(:,3)-data.kinematics(:,9);
dxy = sqrt((data.kinematics(:,1)-data.kinematics(:,7)).^2 + (data.kinematics(:,2)-data.kinematics(:,8)).^2) ;
angFront = atand(-dz ./ dxy);
angFront = angFront(data.purt_frame(1):end);

dz = data.kinematics(:,6)-data.kinematics(:,12);
dxy = sqrt((data.kinematics(:,4)-data.kinematics(:,10)).^2 + (data.kinematics(:,5)-data.kinematics(:,11)).^2) ;
angRight = atand(dz ./ dxy);
angRight = angRight(data.purt_frame(1):end);

for i_syn=1:syn_no
   subplot(syn_no,2,2*i_syn); hold on;
   plot(mean1((angFront ./ max(angFront))',2,1), 'DisplayName', 'Front angle');

    plot(mean1((angRight ./ max(angRight))',2,1), 'DisplayName', 'Right angle');
    ylim([-1 1])
%     legend
end





%% 
