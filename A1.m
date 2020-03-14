cd('C:\DATA\COURSES\Msc Thesis\scripts\Balance0822')
clc
clear
close all
nnmf_init


%% Main Header
subject = [3 4 9 11];
file='T1';
folderToSave = 'raw Data5';
mkdir(folderToSave);
description = [];

%% checking dimensions
load(['raw Data\raw\' file '.mat'])
eval(['T = ',file,';']);

clc

for i_sub = subject
    trialNum = length(T(i_sub).ForcePlate);
    for i_trial = 1:trialNum
        com = size(T(i_sub).ForcePlate(i_trial).COM,1);
        fp = size(T(i_sub).ForcePlate(i_trial).Force,1)/10;
        cop = size(T(i_sub).ForcePlate(i_trial).COP,1)/10;
        % trajectory right
        allMarkerNames = fieldnames(T(i_sub).Trajectory.Right(i_trial));
        kinR = inf;
        for i_marker = 1:length(allMarkerNames)
            dim = size(T(i_sub).Trajectory.Right(i_trial).(allMarkerNames{i_marker}),1);
            if dim<kinR
                kinR = dim;
            end
        end
        % trajectory left
        allMarkerNames = fieldnames(T(i_sub).Trajectory.Left(i_trial));
        kinL = inf;
        for i_marker = 1:length(allMarkerNames)
            dim = size(T(i_sub).Trajectory.Left(i_trial).(allMarkerNames{i_marker}),1);
            if dim<kinL
                kinL = dim;
            end
        end
        % EMG right
        allmusNames = fieldnames(T(i_sub).EMG.Right(i_trial));
        emgR = inf;
        for i_mus = 1:length(allmusNames)
            dim = size(T(i_sub).EMG.Right(i_trial).(allmusNames{i_mus}),1);
            if dim<emgR
                emgR = dim/10;
            end
        end
        % EMG left
        allmusNames = fieldnames(T(i_sub).EMG.Left(i_trial));
        emgL = inf;
        for i_mus = 1:length(allmusNames)
            dim = size(T(i_sub).EMG.Left(i_trial).(allmusNames{i_mus}),1);
            if dim<emgL
                emgL = dim/10;
            end
        end
        
        % choosing min dimension
        if isequal(com,fp,cop,kinR,kinL,emgR,emgL)==0
            fprintf('%s sub:%d trial:%d com:%d fp:%d cop:%d kin:%d emg:%d\n',file,i_sub,i_trial,com,fp,cop,kinR,emgR)
            %             fprintf('COM:%d FP:%d COM:%d ',)
        end
        dimMin = floor(min([com , fp,cop,kinR,kinL,emgR,emgL]));
        
        % bad end of trial(from COP and COM)
        if (i_sub==10) && (i_trial ==1)
            dimMin = 2770;
        end
        
        
        % cutting data...
        
        T(i_sub).ForcePlate(i_trial).COM = T(i_sub).ForcePlate(i_trial).COM(1:dimMin,:);
        T(i_sub).ForcePlate(i_trial).Force = mean1(T(i_sub).ForcePlate(i_trial).Force(1:10*dimMin,:),10,2);
        T(i_sub).ForcePlate(i_trial).COP = mean1(T(i_sub).ForcePlate(i_trial).COP(1:10*dimMin,:),10,2);
        
        % trajectory right
        allMarkerNames = fieldnames(T(i_sub).Trajectory.Right(i_trial));
        for i_marker = 1:length(allMarkerNames)
            T(i_sub).Trajectory.Right(i_trial).(allMarkerNames{i_marker}) = T(i_sub).Trajectory.Right(i_trial).(allMarkerNames{i_marker})(1:dimMin,:);
        end
        % trajectory left
        allMarkerNames = fieldnames(T(i_sub).Trajectory.Left(i_trial));
        for i_marker = 1:length(allMarkerNames)
            T(i_sub).Trajectory.Left(i_trial).(allMarkerNames{i_marker}) = T(i_sub).Trajectory.Left(i_trial).(allMarkerNames{i_marker})(1:dimMin,:);
        end
        % EMG right
        allmusNames = fieldnames(T(i_sub).EMG.Right(i_trial));
        for i_mus = 1:length(allmusNames)
            T(i_sub).EMG.Right(i_trial).(allmusNames{i_mus}) = T(i_sub).EMG.Right(i_trial).(allmusNames{i_mus})(1:10*dimMin);
        end
        % EMG left
        allmusNames = fieldnames(T(i_sub).EMG.Left(i_trial));
        for i_mus = 1:length(allmusNames)
            T(i_sub).EMG.Left(i_trial).(allmusNames{i_mus}) = T(i_sub).EMG.Left(i_trial).(allmusNames{i_mus})(1:10*dimMin);
        end
        
    end
    
end
eval([file,' = T;']);
save([folderToSave '\' file '.mat'],file)
disp('checking dimensions done.')
%% ISB
% just for right side
load([folderToSave '\' file '.mat'])
eval(['T = ',file,';']);

T = ISB(T, subject);

eval([file,' = T;']);
save([folderToSave '\' file '.mat'],file)
%% KIN2mat
% just for right side

load([folderToSave '\' file '.mat'])
eval(['T = ',file,';']);

angles = {[1 2],[1],[1 2 3],[1 2 3]};   % {ankle, knee, hip, pelvis}
T = KIN2mat(T, subject);

eval([file,' = T;']);
save([folderToSave '\' file '.mat'],file)

%% EMG2mat -> filtering EMG
hp1 = 35;
lp1 = 0;
lp2=20;
timeBin = 10;

load([folderToSave '\' file '.mat'])
eval(['T = ',file,';']);

T = EMG2mat(T, subject, hp1, lp1,lp2, timeBin);

eval([file,' = T;']);
save([folderToSave '\' file '.mat'],file)
%% Tilt plate angle and removing bias of COM->COM2 and COP->COP2
load([folderToSave '\' file '.mat'])
eval(['T = ',file,';']);

for i_sub = subject
    fprintf('subject:%d \n',i_sub)
    trialNum = length(T(i_sub).ForcePlate);
    for i_trial = 1:trialNum
        %         T(i_sub).Para(i_trial).COPx = mean1((T(i_sub).ForcePlate(i_trial).COP(:,1) - mean(T(i_sub).ForcePlate(i_trial).COP(:,1)))',10,1);
        %         T(i_sub).Para(i_trial).COPy = mean1((T(i_sub).ForcePlate(i_trial).COP(:,2) - mean(T(i_sub).ForcePlate(i_trial).COP(:,2)))',10,1);
        
        tilt1 = T(i_sub).Trajectory.Right(i_trial).TILT1;
        tilt2 = T(i_sub).Trajectory.Right(i_trial).TILT2;
        tilt3 = T(i_sub).Trajectory.Right(i_trial).TILT3;
        tilt4 = T(i_sub).Trajectory.Right(i_trial).TILT4;
        
        [tiltAP, tiltML] = calTiltAngle(tilt1, tilt2, tilt3, tilt4);
        T(i_sub).KIN.tilt(i_trial).AP = tiltAP;
        T(i_sub).KIN.tilt(i_trial).ML = tiltML;
        
        % removing COM bias by midpoint of tilt plate
        midPoint = mean((tilt1+tilt2+tilt3+tilt4)/4);
        midPoint = midPoint(1:3);
        COM = T(i_sub).ForcePlate(i_trial).COM;
        
        T(i_sub).ForcePlate(i_trial).COM2 = COM - midPoint;
        
        % removing COM bias by substracting the mean
        T(i_sub).ForcePlate(i_trial).COM2 = T(i_sub).ForcePlate(i_trial).COM - mean(T(i_sub).ForcePlate(i_trial).COM);

        % removing COP bias by substracting the mean
        T(i_sub).ForcePlate(i_trial).COP2 = T(i_sub).ForcePlate(i_trial).COP - mean(T(i_sub).ForcePlate(i_trial).COP);

    end
    
    
end
eval([file,' = T;']);
save([folderToSave '\' file '.mat'],file)

disp('calTilt done')


%% Phasic Normalization
load(['raw Data\' file '.mat'])
eval(['T = ',file,';']);

T = PhasicNorm(T, subject);

eval([file,' = T;']);
save(['raw Data\' file '.mat'],file)
%%
plot0720_IntraSubject
plot0720_InterSubject
%% calculate Synergies
calSynergy
%% Visualizations
vis
%% calculate balance index
BI = {'COPx_RMS','COPy_RMS','COPv_RMS','std_vx','std_vy'};
file='T1';
eval(['T = ',file,';']);
for i_sub = subject
    
    trialNum = length(T(i_sub).ForcePlate);
    if trialNum==0
        continue;
    end
    for i_trial = 1:trialNum
        
        T(i_sub).BI(i_trial) = calIndex(T(i_sub).ForcePlate(i_trial).COP, T(i_sub).ForcePlate(i_trial).COM);
        
    end
    
    
end
eval([file,' = T;']);
disp('done calculating balance index')

%% balance index improvement
BI = fieldnames(T1(3).BI);
k_sub=0;
figure
for i_BI = 1:length(BI)
    subplot(3,3,i_BI);    hold on;
    xticks(subject)
    before=cell(length(subject),10);
    
    after=cell(length(subject),10);
    for i_sub = subject
        
        
        
        k_sub=k_sub+1;
        
        title(BI{i_BI})
        trialNum = length(T1(i_sub).ForcePlate);
        for i_trial = 1:trialNum
            before{i_sub}(i_trial) = T1(i_sub).BI(i_trial).(BI{i_BI});
        end
        
        trialNum = length(T2(i_sub).ForcePlate);
        
        for i_trial = 1:trialNum
            after{i_sub}(i_trial) = T2(i_sub).BI(i_trial).(BI{i_BI});
        end
        meanBefore(i_sub,:) = mean(before{i_sub});
        meanAfter(i_sub,:) = mean(after{i_sub});
        stdBefore(i_sub,:) = std(before{i_sub});
        stdAfter(i_sub,:) = std(after{i_sub});
        
    end
    
    barwitherr([stdBefore stdAfter],[meanBefore meanAfter])
    %     bar([meanBefore meanAfter])
    %         plot(ones(1,length(before)),before,'k.')
    %         plot(2*ones(1,length(after)),after,'k.')
    
    %         errorbar([1 2],[mean(before), mean(after)],[std(before) std(after)],...
    %             'linestyle','None')
    
end







%% function
function BI = calIndex(COP, COM)
COPx = COP(:,1) - mean(COP(:,1));
COPy = COP(:,2) - mean(COP(:,2));
COMx = COM(:,1) - mean(COM(:,1));
COMy = COM(:,2) - mean(COM(:,2));
COMz = COM(:,3) - mean(COM(:,3));

COPr = sqrt(COPx .^2 + COPy .^2);
r90CI = sort(COPr);
BI.r90CI = r90CI(floor(.9 * length(COPr)));
BI.COPx_std = std(COPx);
BI.COPy_std = std(COPy);
BI.COMx_std = std(COMx);
BI.COMy_std = std(COMy);
BI.COMz_std = std(COMz);


BI.COPrv_std = std(diff(COPr)*1200);


end
