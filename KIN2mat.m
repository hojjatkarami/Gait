function T = KIN2mat(T, subject)

LR={'Left' 'Right'};

% Kinematic matrix selection
for i_sub=subject
    disp(['subject:',num2str(i_sub)])
    for side=[2]
        for i_trial=1:length(T(i_sub).Trajectory.(LR{side}))
            T(i_sub).KIN.(LR{side})(i_trial).M=[];
            T(i_sub).KIN.(LR{side})(i_trial).M(1,:)= T(i_sub).KIN.(LR{side})(i_trial).value.EulAngAnk.(LR{side})(:,1);
            T(i_sub).KIN.(LR{side})(i_trial).M(2,:)=-T(i_sub).KIN.(LR{side})(i_trial).value.EulAngAnk.(LR{side})(:,2);
            T(i_sub).KIN.(LR{side})(i_trial).M(3,:)=T(i_sub).KIN.(LR{side})(i_trial).value.EulAngKne.(LR{side})(:,1);
            T(i_sub).KIN.(LR{side})(i_trial).M(4,:)=T(i_sub).KIN.(LR{side})(i_trial).value.EulAngHip.(LR{side})(:,1);
            T(i_sub).KIN.(LR{side})(i_trial).M(5,:)=T(i_sub).KIN.(LR{side})(i_trial).value.EulAngHip.(LR{side})(:,2);
            T(i_sub).KIN.(LR{side})(i_trial).M(6,:)=-T(i_sub).KIN.(LR{side})(i_trial).value.EulAngHip.(LR{side})(:,3);
            T(i_sub).KIN.(LR{side})(i_trial).M(7,:)=T(i_sub).KIN.(LR{side})(i_trial).value.EulAngPel(:,1);
            
            if strcmp(LR{side},'Right')
                T(i_sub).KIN.(LR{side})(i_trial).M(8,:)=-T(i_sub).KIN.(LR{side})(i_trial).value.EulAngPel(:,2);
                T(i_sub).KIN.(LR{side})(i_trial).M(9,:)=T(i_sub).KIN.(LR{side})(i_trial).value.EulAngPel(:,3);
            else
                T(i_sub).KIN.(LR{side})(i_trial).M(8,:)=T(i_sub).KIN.(LR{side})(i_trial).value.EulAngPel(:,2);
                T(i_sub).KIN.(LR{side})(i_trial).M(9,:)=-T(i_sub).KIN.(LR{side})(i_trial).value.EulAngPel(:,3);            
            end
            %%%%%% ipsi lateral pos first half
            %%%%%% contral lateral neg second half
            %%%%%% *-1 for right side
            

%             groupName = {'Ank Dorsi Flex','Ank Plant Flex',...
%                          'Ank Inv','Ank Evr','Kne Flex','Kne Ext',...
%                          'Hip Flex','Hip Ext','Hip Add','Hip Abd','Hip Int Rot','Hip Ext Rot','Pelv Upwrd Obliq', 'Pelv Dw Obliq', 'Pelv Int Rot','Pelv Ext Rot'};
%             groupPartitionLine = [4 6 12];  % 2 ankle, 1 knee, 3 hip, 4 pelvis
%             T(i_sub).KIN.(LR{side})(i_trial).groupPartitionLine = groupPartitionLine;

            T(i_sub).KIN.(LR{side})(i_trial).angleOrder = 1:18;
%             T(i_sub).KIN.(LR{side})(i_trial).groupName = groupName;
        end
    end
end
disp('KIN2mat done')