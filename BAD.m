%%
load('T1T2backup\T1.mat')
load('T1T2backup\T2.mat')

%% EMG raw visualization
close all
N=5;   % window size
subjects = [4];
file='T2';
side='Right';
eval(['T = ',file,';']);
% notice to right or left side
for i_sub = subjects
    for i_trial=1:length(T(i_sub).EMG.(side))
%         figure
        idEMG=[];
        idKIN=[];
        for i=1:16
            sig = T(i_sub).EMG.(side)(i_trial).all(:,i);
%             subplot(4,4,i)
%             hold on
%             plot(sig)
            
            %             [sigF ~] = calFeature(sig,N,{'RMS'});
            
            badFrames = abs(sig)<8*std(sig);
            idEMG = [idEMG; find(badFrames==0)];
            
            newSig = sig(badFrames);
            plot(newSig)
        end
        
        idKIN = unique(floor((idEMG-1)/10));
        
        fprintf('\n %d',length(idEMG)/length(sig)*100);
        
        
        % now removing bad parts
        % emg right
        allFields = fieldnames(T(i_sub).EMG.Right(1));
        for i_mus = 1:length(allFields)
            if muscleNo(allFields{i_mus})~=0
                T(i_sub).EMG.Right(i_trial).(allFields{i_mus})(idEMG) = [];
            end
        end
        % emg left
        allFields = fieldnames(T(i_sub).EMG.Left(1));
        for i_mus = 1:length(allFields)
            if muscleNo(allFields{i_mus})~=0
                T(i_sub).EMG.Left(i_trial).(allFields{i_mus})(idEMG) = [];
            end
        end
        % cop moment force
        T(i_sub).ForcePlate(i_trial).COP(idEMG,:) = [];
        T(i_sub).ForcePlate(i_trial).Moment(idEMG,:) = [];
        T(i_sub).ForcePlate(i_trial).Force(idEMG,:) = [];
      
        % kin right
        allFields = fieldnames(T(i_sub).Trajectory.Right(1));
        for i_mus = 1:length(allFields)
            T(i_sub).Trajectory.Right(i_trial).(allFields{i_mus})(idKIN,:) = [];
        end
        
        % kin left
        allFields = fieldnames(T(i_sub).Trajectory.Left(1));
        for i_mus = 1:length(allFields)
            T(i_sub).Trajectory.Left(i_trial).(allFields{i_mus})(idKIN,:) = [];
        end
        % com
        T(i_sub).ForcePlate(i_trial).COM(idKIN,:) = [];
                       
        % end of removing

    end
end
%% test which feature is the best
N=5;   % window size
fea = {'MAV','RMS','WL','ZC','LD','MAX','VAR' };
fea = {'RMS'};
data = T1(4).EMG.Right(2).RPL;
[dataF ~] = calFeature(data,N,fea);

figure; hold on
subplot(2,2,4)
plot(data)
title('RAW')
r=[min(data) max(data)];
ylim(r)
for i=1:length(fea)
    subplot(2,2,i); hold on
    plot(data)
    newSig = data(abs(dataF(:,i))<mean(dataF(:,i))+5*std(dataF(:,i)));
    % newSig = data(abs(data)<4*std(data));
    
    plot(newSig)
    %     ylim(r)
    title([fea{i} ' ' num2str(max(newSig)-min(newSig))])
    plot(dataF(:,i)-4)
end
% plot(dataF(abs(dataF)<3*std(dataF))+2)
%%
load('T1T2backup\T1.mat')
load('T1T2backup\T2.mat')
%% T1 bad
% subject 3 all(11) left side is very bad

bad(4).trial(3).val = [15741 15810];    %right
bad(4).trial(2).val = [7811 8080; 22141 22630; 29801 30300; 45511,45610];   %right
bad(4).trial(1).val = [20481 20660];    %right

%% T2 bad
bad(4).trial(4).val = [32721 33480];    %left


%% T1 or T2?
T=T1;
%% remove
% just run once

for i_sub=4:11
    if isempty(T(i_sub).EMG)
        continue
    end
    trialNum = length(T(i_sub).EMG.Right);
    
    for i_trial=1:trialNum
        for i=1:size(bad(i_sub).trial(i_trial).val,1)
            range1 = [bad(i_sub).trial(i_trial).val(i,1):bad(i_sub).trial(i_trial).val(i,2)];
            
            
            % emg right
            allFields = fieldnames(T(i_sub).EMG.Right(1));
            for i_mus = 1:length(allFields)
                if muscleNo(allFields{i_mus})~=0
                    T(i_sub).EMG.Right(i_trial).(allFields{i_mus})(range1) = [];
                end
            end
            % emg left
            allFields = fieldnames(T(i_sub).EMG.Left(1));
            for i_mus = 1:length(allFields)
                if muscleNo(allFields{i_mus})~=0
                    T(i_sub).EMG.Left(i_trial).(allFields{i_mus})(range1) = [];
                end
            end
            % cop moment force
            T(i_sub).ForcePlate(i_trial).COP(range1,:) = [];
            T(i_sub).ForcePlate(i_trial).Moment(range1,:) = [];
            T(i_sub).ForcePlate(i_trial).Force(range1,:) = [];
            
            
            
            range2 = [(bad(i_sub).trial(i_trial).val(i,1)-1)/10+1:bad(i_sub).trial(i_trial).val(i,2)/10];
            
            % kin right
            allFields = fieldnames(T(i_sub).Trajectory.Right(1));
            for i_mus = 1:length(allFields)
                T(i_sub).Trajectory.Right(i_trial).(allFields{i_mus})(range2,:) = [];
            end
            
            % kin left
            allFields = fieldnames(T(i_sub).Trajectory.Left(1));
            for i_mus = 1:length(allFields)
                T(i_sub).Trajectory.Left(i_trial).(allFields{i_mus})(range2,:) = [];
            end
            % com
            T(i_sub).ForcePlate(i_trial).COM(range2,:) = [];
            
        end
    end
    
end

%% save T1
T1=T;
save('T1.mat','T1')
%% save T2
T2=T;
save('T1.mat','T1')
