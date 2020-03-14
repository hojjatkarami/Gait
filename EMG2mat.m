function T = EMG2mat(T, subject, hp1, lp1,lp2, timeBin)

LR={'Left' 'Right'};
% EMG filteration and matrix extraction
for i_sub=subject
    disp(['subject:',num2str(i_sub)])
    % EMG filteration
    % add all field to EMG field(contining all muscles)-each muscle in each
    % column for both Side
    for side=[1 2]
        for i_trial=1:length(T(i_sub).EMG.(LR{side}))
            
            T(i_sub).EMG.(LR{side})(i_trial).all = [];
            muscleOrder = [];
            T(i_sub).EMG.(LR{side})(i_trial).muscleOrder = [];
            allFields = fieldnames(T(i_sub).EMG.(LR{side})(i_trial));
            for i_mus = 1:length(allFields)
                if muscleNo(allFields{i_mus})~=0
                    
                    % cheacking if length of muscle sig is not consistent
%                     if mod(length(T(i_sub).EMG.(LR{side})(i_trial).(allFields{i_mus})), 10)~=0
%                         warning('wrong muscle dimensions/ sub:%d trial:%d side:%s mus:%s',i_sub, i_trial, LR{side},allFields{i_mus})
%                         T(i_sub).EMG.(LR{side})(i_trial).(allFields{i_mus})(end+1) = T(i_sub).EMG.(LR{side})(i_trial).(allFields{i_mus})(end);
%                         
%                     end
                    
                    T(i_sub).EMG.(LR{side})(i_trial).all=[T(i_sub).EMG.(LR{side})(i_trial).all, ...
                        T(i_sub).EMG.(LR{side})(i_trial).(allFields{i_mus})];
                    muscleOrder = [muscleOrder muscleNo(allFields{i_mus})];
                end
            end
            [B,I]=sort(muscleOrder);
            T(i_sub).EMG.(LR{side})(i_trial).muscleOrder = B;
            T(i_sub).EMG.(LR{side})(i_trial).muscleName = muscleName(B);
            T(i_sub).EMG.(LR{side})(i_trial).all = T(i_sub).EMG.(LR{side})(i_trial).all(:,I);
            
        end
    end
    
    % filter for each side and each trial
    for side=[1 2]
        for i_trial = 1:length(T(i_sub).Trajectory.(LR{side}))
            EMG_raw_Mat = T(i_sub).EMG.(LR{side})(i_trial).all;
            EMG_filt_Mat = filterEMGmat1222(EMG_raw_Mat,hp1,lp1,lp2,timeBin);
            T(i_sub).EMG.(LR{side})(i_trial).M = EMG_filt_Mat';
        end
    end
    T(i_sub).EMG.filter.hp1 = hp1;
    T(i_sub).EMG.filter.lp1 = lp1;
    T(i_sub).EMG.filter.lp2 = lp2;
end

disp('EMG2mat done')
