clc
clear
close all
%% Settings
load('raw Data/T2.mat');  % read Database P(created by Vicon2matlab) and modify it

subject = [1:11];  % subjects id to be proccessed here
LR={'Left' 'Right'};
Co={'PelCo' 'ThiCo' 'TibCo' 'UNTibCo' 'FootCo'};
xyz={'x' 'y' 'z'};

%% Kinematic angle extraction
for i_sub=subject
    disp(['subject:',num2str(i_sub)])
    for side=[1 2]
        
        Static.Trajectory=P(i_sub).Static.(LR{side});
        Static=Coordinate(Static);  % Segments Coordinate Defenition
        % Rotation Matrix for Static Offset _______________________________________
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
        P(i_sub).Static.Static0.(LR{side})=Static0;
        P(i_sub).Static.Static2.(LR{side})=Static2;
        
        for i_trial=1:length(P(i_sub).Trajectory.(LR{side}))
            kin.Trajectory = P(i_sub).Trajectory.(LR{side})(i_trial);            
            kin=Coordinate(kin);  % Segments Coordinate Defenition         
            kin2=EulerAngle(kin,Static0);   % Euler Angle Calculation             
            P(i_sub).kin2.(LR{side})(i_trial).value = kin2;

        clear kin kin2
        end
        clear Static Static0 Static2
    end
end
%% Kinematic matrix selection
for i_sub=subject
    disp(['subject:',num2str(i_sub)])
    for side=[1 2]
        for i_trial=1:length(P(i_sub).Trajectory.(LR{side}))
            P(i_sub).kin2.(LR{side})(i_trial).M=[];
            P(i_sub).kin2.(LR{side})(i_trial).M(1,:)= P(i_sub).kin2.(LR{side})(i_trial).value.EulAngAnk.(LR{side})(:,1);
            P(i_sub).kin2.(LR{side})(i_trial).M(2,:)=-P(i_sub).kin2.(LR{side})(i_trial).value.EulAngAnk.(LR{side})(:,2);
            P(i_sub).kin2.(LR{side})(i_trial).M(3,:)=P(i_sub).kin2.(LR{side})(i_trial).value.EulAngKne.(LR{side})(:,1);
            P(i_sub).kin2.(LR{side})(i_trial).M(4,:)=P(i_sub).kin2.(LR{side})(i_trial).value.EulAngHip.(LR{side})(:,1);
            P(i_sub).kin2.(LR{side})(i_trial).M(5,:)=P(i_sub).kin2.(LR{side})(i_trial).value.EulAngHip.(LR{side})(:,2);
            P(i_sub).kin2.(LR{side})(i_trial).M(6,:)=-P(i_sub).kin2.(LR{side})(i_trial).value.EulAngHip.(LR{side})(:,3);
            P(i_sub).kin2.(LR{side})(i_trial).M(7,:)=P(i_sub).kin2.(LR{side})(i_trial).value.EulAngPel(:,1);
            
            if strcmp(LR{side},'Right')
                P(i_sub).kin2.(LR{side})(i_trial).M(8,:)=-P(i_sub).kin2.(LR{side})(i_trial).value.EulAngPel(:,2);
                P(i_sub).kin2.(LR{side})(i_trial).M(9,:)=P(i_sub).kin2.(LR{side})(i_trial).value.EulAngPel(:,3);
            else
                P(i_sub).kin2.(LR{side})(i_trial).M(8,:)=P(i_sub).kin2.(LR{side})(i_trial).value.EulAngPel(:,2);
                P(i_sub).kin2.(LR{side})(i_trial).M(9,:)=-P(i_sub).kin2.(LR{side})(i_trial).value.EulAngPel(:,3);
            
            end
            %%%%%% ipsi lateral pos first half
            %%%%%% contral lateral neg second half
            %%%%%% *-1 for right side
            

%             groupName = {'Ank Dorsi Flex','Ank Plant Flex',...
%                          'Ank Inv','Ank Evr','Kne Flex','Kne Ext',...
%                          'Hip Flex','Hip Ext','Hip Add','Hip Abd','Hip Int Rot','Hip Ext Rot','Pelv Upwrd Obliq', 'Pelv Dw Obliq', 'Pelv Int Rot','Pelv Ext Rot'};
%             groupPartitionLine = [4 6 12];  % 2 ankle, 1 knee, 3 hip, 4 pelvis
%             P(i_sub).kin2.(LR{side})(i_trial).groupPartitionLine = groupPartitionLine;

            P(i_sub).kin2.(LR{side})(i_trial).angleOrder = 1:18;
%             P(i_sub).kin2.(LR{side})(i_trial).groupName = groupName;
        end
    end
end
%% EMG filteration and matrix extraction
for i_sub=subject
    disp(['subject:',num2str(i_sub)])
% EMG filteration
    % add all field to EMG field(contining all muscles)-each muscle in each
    % column for both Side
    for side=[1 2]
        for i_trial=1:length(P(i_sub).Trajectory.(LR{side}))
            P(i_sub).EMG2.(LR{side})(i_trial).all = [];
            muscleOrder = [];
            P(i_sub).EMG2.(LR{side})(i_trial).muscleOrder = [];
            allFields = fieldnames(P(i_sub).EMG.(LR{side})(i_trial));
            for i_mus = 1:length(allFields)
               if muscleNo(allFields{i_mus})~=0
                   P(i_sub).EMG2.(LR{side})(i_trial).all=[P(i_sub).EMG2.(LR{side})(i_trial).all, ...
                                        P(i_sub).EMG.(LR{side})(i_trial).(allFields{i_mus})];
                    muscleOrder = [muscleOrder muscleNo(allFields{i_mus})];
               end            
            end
            [B,I]=sort(muscleOrder);
            P(i_sub).EMG2.(LR{side})(i_trial).muscleOrder = B;
            P(i_sub).EMG2.(LR{side})(i_trial).all = P(i_sub).EMG2.(LR{side})(i_trial).all(:,I);

        end
    end
    
    % filter for each side and each trial
    for side=[1 2] 
        for i_trial = 1:length(P(i_sub).Trajectory.(LR{side}))
            EMG_raw_Mat = P(i_sub).EMG2.(LR{side})(i_trial).all;
            EMG_filt_Mat = filterEMGmat0828(EMG_raw_Mat,30,10,10);
%             P(i_sub).EMG2.(LR{side})(i_trial).Filtered = EMG_filt_Mat;
            P(i_sub).EMG2.(LR{side})(i_trial).M = EMG_filt_Mat';            
        end
%         P(i_sub).EMG2.(LR{side}) = rmfield(P(i_sub).EMG2.(LR{side}),'all');
    end    
 
end
%% phasic time normalization
for i_sub=subject
    disp(['subject:',num2str(i_sub)])
%     %% Remove Static Offset and Calculate Kinematic angles
% 
    
    % Calculate angles of plate and find purt frames
    % add all and purtFrame field to Event
    for i_trial= 1:length(P(i_sub).Trajectory.Left)
        P(i_sub).Events.Left(i_trial).all=double([P(i_sub).Events.Left(i_trial).LFC1 P(i_sub).Events.Left(i_trial).RTO...
        P(i_sub).Events.Left(i_trial).RFC P(i_sub).Events.Left(i_trial).LTO P(i_sub).Events.Left(i_trial).LFC2]);
        P(i_sub).Events.Left(i_trial).PurtFrame = P(i_sub).Events.Left(i_trial).all - P(i_sub).Events.Left(i_trial).all(1)+1;
    end
    
    for i_trial= 1:length(P(i_sub).Trajectory.Right)
        P(i_sub).Events.Right(i_trial).all=double([P(i_sub).Events.Right(i_trial).RFC1 P(i_sub).Events.Right(i_trial).LTO...
        P(i_sub).Events.Right(i_trial).LFC P(i_sub).Events.Right(i_trial).RTO P(i_sub).Events.Right(i_trial).RFC2]);
        P(i_sub).Events.Right(i_trial).PurtFrame = P(i_sub).Events.Right(i_trial).all - P(i_sub).Events.Right(i_trial).all(1)+1;

    end

    % phasic time normalization for EMG and kinematics
    for side=[1 2]
        for i_trial= 1:length(P(i_sub).Trajectory.(LR{side}))

            PurtFrame = P(i_sub).Events.(LR{side})(i_trial).PurtFrame;
            PurtFrameNorm = round(PurtFrame/PurtFrame(end)*100);
            if PurtFrameNorm(1)==0
                PurtFrameNorm(1)=1;
            end
            P(i_sub).Events.(LR{side})(i_trial).PurtFrameNorm = PurtFrameNorm;

            Mat1 = P(i_sub).EMG2.(LR{side})(i_trial).M;        
            [Mat2,t] = rescale0720(Mat1, PurtFrame, PurtFrameNorm);
            P(i_sub).EMG2.(LR{side})(i_trial).Mn = Mat2;
            P(i_sub).Events.(LR{side})(i_trial).t = t;
            P(i_sub).Events.(LR{side})(i_trial).tNorm = (t/t(end)*100);
            Mat1 = P(i_sub).kin2.(LR{side})(i_trial).M;        
            [Mat2,t] = rescale0720(Mat1, PurtFrame, PurtFrameNorm);
            P(i_sub).kin2.(LR{side})(i_trial).Mn = Mat2;
        end
    end

end

%% Save mat file
if strcmp(questdlg('Save P file?','Save file','Yes','No','Yes'),'Yes')
    
    save('raw Data/P.mat','P');
    G = P;
    clear P
    load ('P.mat')
    G = rmfield(G,'ForcePlate');
    G = rmfield(G,'Kinetics');
    G = rmfield(G,'Kinematics');
    G = rmfield(G,'JC');
    G = rmfield(G,'Trajectory');    
    G = rmfield(G,'EMG');
    [G.EMG] = G.EMG2; G = orderfields(G,[1:5,7,6:6]); G = rmfield(G,'EMG2');
    [G.kin] = G.kin2; G = orderfields(G,[1:3,7,4:6]); G = rmfield(G,'kin2');
    for i_sub=[1:11]
        G(i_sub).Synergy = P(i_sub).Synergy
    end
    disp('Done saving!')
    P=G;
    save('P.mat','P');
end