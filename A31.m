
% Non-synergy comparison

subject = [3 4 9 10 11];
file='T1';
load(['raw Data\T1.mat'])
load(['raw Data\T2.mat'])

LR={'Left' 'Right'};
allMuscleNames = T1(3).EMG.Right(1).muscleName;

disp('T1 T2 successfully loaded.')
%% calculate balance Indices
for i_sub=subject
    disp(['subject:',num2str(i_sub)])
    
    % before
    trialNum = length(T1(i_sub).Trajectory.Right);
    for i_trial=1:trialNum
        T1 = calIndex(T1, i_sub, i_trial);
        
    end
    
    % after
    trialNum = length(T2(i_sub).Trajectory.Right);
    for i_trial=1:trialNum
        T2 = calIndex(T2, i_sub, i_trial);
        
    end
    
end
%% compare BI for entire trial
index = {'r90CI','COPx_std','COPy_std','COMx_std','COMy_std','COMz_std','COPrv_std','AP_std','ML_std'};
alpha = 0.05;

figure(1)
figure(2)
for i = 1:length(index)
    before = {};
    after = {};
    
    for i_sub = subject
        
        before = [before transpose( [T1(i_sub).BI(:).(index{i})] )];
        
        
        after = [after transpose( [T2(i_sub).BI(:).(index{i})] )];
        
        
        
    end
    figure(1);
    ax = subplot(3,3,i);
    allMean = [cellfun(@mean,before);...
        cellfun(@mean,after)]';
    allStd = [cellfun(@std,before);...
        cellfun(@std,after)]';
    
    plotMultipleBarsWithError(ax, allMean, allStd);
    
    title(index{i})
    xticks([1:length(subject)]);
    xticklabels(num2cell(subject));
    
    figure(2)
    [h p] = cellfun(@(x,y) ttest2(x,y,'Alpha',alpha),before,after);
    subplot(3,3,i); hold on;
    plot(find(h==1),p(find(h==1)),'b*','linestyle','none')
    plot(find(h==0),p(find(h==0)),'r*','linestyle','none')
    
    %     ylim([0 alpha])
    xlim([0 length(subject)+1])
    title(index{i})
    xticks([1:length(subject)]);
    xticklabels(num2cell(subject));
end
%% plot EMG variations
index = {'EMG_C_duty' 'EMG_C_area' 'EMG_C_rmse'};
alpha = 0.05;


for i = 1:length(index)
    before = {};
    after = {};
    
    for i_sub = subject
        
        before = [before transpose( [T1(i_sub).BI(:).(index{i})] )];
        
        
        after = [after transpose( [T2(i_sub).BI(:).(index{i})] )];
        
        
        
    end
    figure(10+i);
    allMean=cell(1,16);
    allStd = cell(1,16);
    for i_mus=1:16
        
        for j=1:length(before)
            allMean{i_mus} = [allMean{i_mus}; [mean(before{j}(:,i_mus)) mean(after{j}(:,i_mus))]];
            allStd{i_mus} = [allStd{i_mus}; [std(before{j}(:,i_mus)) std(after{j}(:,i_mus))]];
            
        end
        subplot(4,4,i_mus);
        
        plotMultipleBarsWithError(ax, allMean{i_mus}, allStd{i_mus});
        
        title(allMuscleNames{i_mus})
        xticks([1:length(subject)]);
        xticklabels(num2cell(subject));
        
    end
    suptitle(index{i})
    figure(20+i)
    
    for i_mus=1:16
        subplot(4,4,i_mus); hold on;
        for i_sub=1:length(before)
            [h p] = ttest2(before{i_sub}(:,i_mus),after{i_sub}(:,i_mus),'Alpha',alpha);
            if h==1
                plot(i_sub,p,'b*','linestyle','none')
                
                
            else
                plot(i_sub,p,'r*','linestyle','none')
                
            end
            
            
        end
        
        
        %     ylim([0 alpha])
        xlim([0 length(subject)+1])
        title(allMuscleNames{i_mus})
        xticks([1:length(subject)]);
        xticklabels(num2cell(subject));
        
    end
    suptitle(index{i})
    
end

%% function
function T = calIndex(T, i_sub, i_trial)
COP = T(i_sub).ForcePlate(i_trial).COP2;
COM = T(i_sub).ForcePlate(i_trial).COM2;
angAP = T(i_sub).KIN.tilt(i_trial).AP;
angML = T(i_sub).KIN.tilt(i_trial).ML;

COPx = COP(:,1) - mean(COP(:,1));
COPy = COP(:,2) - mean(COP(:,2));
COMx = COM(:,1) - mean(COM(:,1));
COMy = COM(:,2) - mean(COM(:,2));
COMz = COM(:,3) - mean(COM(:,3));

COPr = sqrt(COPx .^2 + COPy .^2);
r90CI = sort(COPr);
T(i_sub).BI(i_trial).r90CI = r90CI(floor(.9 * length(COPr)));
T(i_sub).BI(i_trial).COPx_std = std(COPx);
T(i_sub).BI(i_trial).COPy_std = std(COPy);
T(i_sub).BI(i_trial).COMx_std = std(COMx);
T(i_sub).BI(i_trial).COMy_std = std(COMy);
T(i_sub).BI(i_trial).COMz_std = std(COMz);
T(i_sub).BI(i_trial).AP_std = std(angAP);
T(i_sub).BI(i_trial).ML_std = std(angML);



T(i_sub).BI(i_trial).COPrv_std = std(diff(COPr)*1200);

M = T(i_sub).EMG.Right(i_trial).M';
T(i_sub).BI(i_trial).EMG_C_duty = transpose(sum(M>.15) / size(M,1));
T(i_sub).BI(i_trial).EMG_C_area = transpose(trapz( max(0, M-0.15) ) / size(M,1));
T(i_sub).BI(i_trial).EMG_C_rmse = transpose(rms(M));

end

