% header
eval(['T = ',file,';']);

LR = {'Right' 'Left'};

highPassFreq = 40;
lowPassFreq = 15;
timeBin = 10;


LR={'Left' 'Right'};
Co={'PelCo' 'ThiCo' 'TibCo' 'UNTibCo' 'FootCo'};
xyz={'x' 'y' 'z'};

%% calculate tilt plate angle

for i_sub = subject
    
    trialNum = length(T(i_sub).ForcePlate);
    for i_trial = 1:trialNum
        T(i_sub).Para(i_trial).COPx = mean1((T(i_sub).ForcePlate(i_trial).COP(:,1) - mean(T(i_sub).ForcePlate(i_trial).COP(:,1)))',10,1);
        T(i_sub).Para(i_trial).COPy = mean1((T(i_sub).ForcePlate(i_trial).COP(:,2) - mean(T(i_sub).ForcePlate(i_trial).COP(:,2)))',10,1);
        
        tilt1 = T(i_sub).Trajectory.Right(i_trial).TILT1;
        tilt2 = T(i_sub).Trajectory.Right(i_trial).TILT2;
        tilt3 = T(i_sub).Trajectory.Right(i_trial).TILT3;
        tilt4 = T(i_sub).Trajectory.Right(i_trial).TILT4;
        
        [tiltAP, tiltML] = calTiltAngle(tilt1, tilt2, tilt3, tilt4);
        T(i_sub).Para(i_trial).tiltAP = tiltAP;
        T(i_sub).Para(i_trial).tiltML = tiltML;
    end
    
    
end
eval([file,' = T;']);
%%

%% concatenating EMG data into all matrix
for i_sub = subject
    
    trialNum = length(T(i_sub).ForcePlate);
    for i_trial = 1:trialNum
        for side=[1 2]
            T(i_sub).EMG.(LR{side})(i_trial).all = [];
            muscleOrder = [];
            T(i_sub).EMG.(LR{side})(i_trial).muscleOrder = [];
            allFields = fieldnames(T(i_sub).EMG.(LR{side})(1));
            for i_mus = 1:length(allFields)
                if muscleNo(allFields{i_mus})~=0
                    T(i_sub).EMG.(LR{side})(i_trial).all=[T(i_sub).EMG.(LR{side})(i_trial).all, ...
                        T(i_sub).EMG.(LR{side})(i_trial).(allFields{i_mus})];
                    muscleOrder = [muscleOrder muscleNo(allFields{i_mus})];
                end
            end
            [B,I]=sort(muscleOrder);
            T(i_sub).EMG.(LR{side})(i_trial).muscleOrder = B;
            T(i_sub).EMG.(LR{side})(i_trial).all = T(i_sub).EMG.(LR{side})(i_trial).all(:,I);
            
            
        end
    end
    
    
    
end

disp('EMG concatenatoin done!')

%% EMG remove bad parts
N=5;   % window size
side='Right';
eval(['T = ',file,';']);
% notice to right or left side
for i_sub = subject
    for i_trial=1:length(T(i_sub).EMG.(side))
%         figure
        idEMG=[];
        idKIN=[];
        for i=1:16
            sigR = T(i_sub).EMG.Right(i_trial).all(:,i);
            sigL = T(i_sub).EMG.Left(i_trial).all(:,i);

%             subplot(4,4,i)
%             hold on
%             plot(sig)
            
            %             [sigF ~] = calFeature(sig,N,{'RMS'});
            
            badFramesR = abs(sigR)<8*std(sigR);
            badFramesL = abs(sigR)<8*std(sigL);
            
            idEMG = [idEMG; find(badFramesR==0)];
            idEMG = [idEMG; find(badFramesL==0)];
            
%             newSig = sig(badFrames);
%             plot(newSig)
        end
        
        idKIN = unique(floor((idEMG-1)/10));
        
        fprintf('\n %d',length(idEMG)/length(sigR)*100);
        
        
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
%%
for side=[1 2]
    for i_trial = 1:trialNum
        EMG_raw_Mat = T(i_sub).EMG.(LR{side})(i_trial).all;
        EMG_filt_Mat = filterEMGmat0828(EMG_raw_Mat, highPassFreq, lowPassFreq, timeBin);
        T(i_sub).EMG.(LR{side})(i_trial).M = EMG_filt_Mat';
    end
end
eval([file,' = T;']);

%% Kinematic angle extraction
for i_sub=subject
    disp(['subject:',num2str(i_sub)])
    for side=[1 2]
        
        Static.Trajectory=T(i_sub).Static.(LR{side});
        Static=Coordinate(Static);  % Segments Coordinate Defenition
        Rotation Matrix for Static Offset
        for i=1:length(Co)
            for k=1:length(xyz)
                if i==1
                    Static0.(Co{i}).(xyz{k})=mean(Static.(Co{i}).(xyz{k}),1);
                else
                    for h=1:length(LR)     %Left,Right
                        Static0.(Co{i}).(LR{h}).(xyz{k})=mean(Static.(Co{i}).(LR{h}).(xyz{k}),1);
                    end
                end
            end
        end
        Static2=EulerAngle(Static,Static0);     % Euler Angle Calculation
        T(i_sub).Static.Static0.(LR{side})=Static0;
        T(i_sub).Static.Static2.(LR{side})=Static2;
        
        for i_trial=1:length(T(i_sub).Trajectory.(LR{side}))
            kin.Trajectory = T(i_sub).Trajectory.(LR{side})(i_trial);
            kin=Coordinate(kin);  % Segments Coordinate Defenition
            kin2=EulerAngle(kin,Static0);   % Euler Angle Calculation
            T(i_sub).KIN.(LR{side})(i_trial).value = kin2;
            
            clear kin kin2
        end
        clear Static Static0 Static2
    end
    
    
    
end
%% Kinematic matrix selection
for i_sub=subject
    disp(['subject:',num2str(i_sub)])
    for side=[1 2]
        for i_trial=1:length(T(i_sub).Trajectory.(LR{side}))
            T(i_sub).KIN.(LR{side})(i_trial).M=[];
            T(i_sub).KIN.(LR{side})(i_trial).M(1,:)= T(i_sub).KIN.(LR{side})(i_trial).value.EulAngAnk.(LR{side})(:,1);
            T(i_sub).KIN.(LR{side})(i_trial).M(2,:)=-1.*T(i_sub).KIN.(LR{side})(i_trial).value.EulAngAnk.(LR{side})(:,2);
            T(i_sub).KIN.(LR{side})(i_trial).M(3,:)=T(i_sub).KIN.(LR{side})(i_trial).value.EulAngKne.(LR{side})(:,1);
            T(i_sub).KIN.(LR{side})(i_trial).M(4,:)=T(i_sub).KIN.(LR{side})(i_trial).value.EulAngHip.(LR{side})(:,1);
            T(i_sub).KIN.(LR{side})(i_trial).M(5,:)=T(i_sub).KIN.(LR{side})(i_trial).value.EulAngHip.(LR{side})(:,2);
            T(i_sub).KIN.(LR{side})(i_trial).M(6,:)=T(i_sub).KIN.(LR{side})(i_trial).value.EulAngHip.(LR{side})(:,3);
            T(i_sub).KIN.(LR{side})(i_trial).M(7,:)=T(i_sub).KIN.(LR{side})(i_trial).value.EulAngPel(:,1);
            T(i_sub).KIN.(LR{side})(i_trial).M(8,:)=T(i_sub).KIN.(LR{side})(i_trial).value.EulAngPel(:,2);
            T(i_sub).KIN.(LR{side})(i_trial).M(9,:)=T(i_sub).KIN.(LR{side})(i_trial).value.EulAngPel(:,3);
            groupName = {'Ank Dorsi Flex','Ank Plant Flex',...
                'Ank Inv','Ank Evr','Kne Flex','Kne Ext',...
                'Hip Flex','Hip Ext','Hip Add','Hip Abd','Hip Int Rot','Hip Ext Rot','Pelv Upwrd Obliq', 'Pelv Dw Obliq', 'Pelv Int Rot','Pelv Ext Rot'};
            groupPartitionLine = [4 6 12];  % 2 ankle, 1 knee, 3 hip, 4 pelvis
            T(i_sub).KIN.(LR{side})(i_trial).groupPartitionLine = groupPartitionLine;
            
            T(i_sub).KIN.(LR{side})(i_trial).angleOrder = 1:18;
            T(i_sub).KIN.(LR{side})(i_trial).groupName = groupName;
        end
    end
end
