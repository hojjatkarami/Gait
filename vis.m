
%% EMG filtered visualization
close all
subjects = [4];
file='T1';
eval(['T = ',file,';']);

for i_sub = subjects
    for i_trial=1:length(T(i_sub).EMG.Right)
        figure
        for i=1:16
            subplot(4,4,i)
            plot(T(i_sub).EMG.Right(i_trial).M(i,:))
        end
    end
end

%% KIN visualization
close all
subjects = [4];
file='T1';
eval(['T = ',file,';']);

for i_sub = subjects
    for i_trial=1:length(T(i_sub).KIN.Right)
        figure
        for i=1:9
            subplot(3,3,i)
            plot(T(i_sub).KIN.Right(i_trial).M(i,:))
        end
    end
end
%% COM visulaization
subjects = [3];
file='T1';
eval(['T = ',file,';']);
figure
for i_sub = subjects
    numTrials = length(T(i_sub).KIN.Right);
    for i_trial=1:numTrials
        
        for i=1:3
            sig = T(i_sub).ForcePlate(i_trial).COM(:,i);
            subplot(numTrials,4,4*(i_trial-1)+i)
            plot(sig)
            rmsVal = rms(sig-mean(sig));
            stdVal = std(sig);
            title(['rms: ',num2str(rmsVal),' /std: ',num2str(stdVal)])
        end
        subplot(numTrials,4,4*(i_trial-1)+4)
        plot(T(i_sub).ForcePlate(i_trial).COM(:,1), T(i_sub).ForcePlate(i_trial).COM(:,2))
        axis equal

    end
end

%% COP visulaization
subjects = [3];
file='T1';
eval(['T = ',file,';']);
figure
for i_sub = subjects
    numTrials = length(T(i_sub).KIN.Right);
    for i_trial=1:numTrials
        
        for i=1:2
            sig = T(i_sub).ForcePlate(i_trial).COP(:,i);
            subplot(numTrials,3,3*(i_trial-1)+i)
            plot(sig)
            rmsVal = rms(sig-mean(sig));
            stdVal = std(sig);
            title(['rms: ',num2str(rmsVal),' /std: ',num2str(stdVal)])
        end
        subplot(numTrials,3,3*(i_trial-1)+3)
        plot(T(i_sub).ForcePlate(i_trial).COP(:,1), T(i_sub).ForcePlate(i_trial).COP(:,2))
        axis equal

    end
end