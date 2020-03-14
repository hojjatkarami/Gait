%% visula inspection of signals

subject = [3 4 9 10 11];
file='T2';
load([folderToSave '\' file '.mat'])
eval(['T = ',file,';']);
LR={'Left' 'Right'};
allMuscleNames = T(3).EMG.Right(1).muscleName;
%% plot KIN angles by ISB

for i_sub=subject
    disp(['subject:',num2str(i_sub)])
    figure
    
    for i_trial=1:length(T(i_sub).Trajectory.Right)
        
        
        MR = T(i_sub).KIN.Right(i_trial).M;
        ML = T(i_sub).KIN.Left(i_trial).M;
        
        n=size(MR,1);
        
        for i=1:n
            subplot(n,2,2*i-1);  hold on;
            plot(ML(i,:))
            subplot(n,2,2*i) ;   hold on;
            plot(MR(i,:))
            title(angleName(18+i))
        end
        %             suptitle(['sub:' num2str(i_sub) ' trial:' num2str(i_trial)])
    end
    suptitle(['sub:' num2str(i_sub) ' all trials'])
    
end

%% plot force plate and tilt plate angles

for i_sub=subject
    disp(['subject:',num2str(i_sub)])
    figure
    trialNum = length(T(i_sub).Trajectory.Right);
    for i_trial=1:length(T(i_sub).Trajectory.Right)
        
        
        AP = T(i_sub).KIN.tilt(i_trial).AP;
        ML = T(i_sub).KIN.tilt(i_trial).ML;
        COMx = T(i_sub).ForcePlate(i_trial).COM2(:,1);
        COMy = T(i_sub).ForcePlate(i_trial).COM2(:,2);
        COPx = T(i_sub).ForcePlate(i_trial).COP2(:,1);
        COPy = T(i_sub).ForcePlate(i_trial).COP2(:,2);
        subplot(trialNum,3,3*i_trial-2);  hold on;
        plot(AP,ML)
        axis equal
        ylim([-20 20])
        xlim([-20 20])
        title('AP/ML')
        subplot(trialNum,3,3*i_trial-1) ;   hold on;
        plot(COMx,COMy)
        axis equal
                ylim([-100 100])
                xlim([-100 100])
        title('COMx/COMy')
        subplot(trialNum,3,3*i_trial-0) ;   hold on;
        plot(COPx,COPy)
        axis equal
                ylim([-150 150])
                xlim([-150 150])
        title('COPx/COPy')
        %             suptitle(['sub:' num2str(i_sub) ' trial:' num2str(i_trial)])
    end
    suptitle(['sub:' num2str(i_sub) ' all trials'])
    
end
%% EMG visualization raw and filtered
LR={'Left' 'Right'};
allMuscleNames = T(3).EMG.Right(1).muscleName;
i_mus = 4;
side = 1;

k=1;
h1=figure;
set(h1,'units','normalized','outerposition',[0 0 1 1]);

for i_sub=subject
    disp(['subject:',num2str(i_sub)])
    
    trialNum = length(T(i_sub).Trajectory.Right);
    
    for i_trial=1:length(T(i_sub).Trajectory.Right)
        
        
        
        M_raw = mean1(T(i_sub).EMG.(LR{side})(i_trial).all(:,i_mus)',10,1);
        M_filt = T(i_sub).EMG.(LR{side})(i_trial).M(i_mus,:);
        %         M_filt = 0;
        subplot(5,5,k)
        hold on;
        plot(M_raw)
        plot(M_filt)
        title(['s:' num2str(i_sub) ' t:' num2str(i_trial)])
        
        
        k=k+1;
        %         suptitle('')
    end
    
end

suptitle([file ' - ' allMuscleNames{i_mus} ' - ' (LR{side})])
%% fft of all muscle for single trial
i_sub = 3;
i_trial = 1;
side=1;



figure
Fs = 1200;            % Sampling frequency                    
TT = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*TT;        % Time vector
S = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);
for i_mus = 1:16
S = T(i_sub).EMG.(LR{side})(i_trial).all(:,i_mus);
L = length(S);
X = S + 2*randn(size(t));

Y = fft(S);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
subplot(4,4,i_mus)
plot(f,P1) 
title (allMuscleNames{i_mus})
ylim([0 .01])
xlim([0 250])
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
end
suptitle('Single-Sided Amplitude Spectrum of S(t)')





%% plot EMG by single feature
fea = 'RMS';
window = 100;

k=1;
h1=figure;
set(h1,'units','normalized','outerposition',[0 0 1 1]);

for i_sub=subject
    disp(['subject:',num2str(i_sub)])
    
    trialNum = length(T(i_sub).Trajectory.Right);
    
    for i_trial=1:length(T(i_sub).Trajectory.Right)
        
        
        M_fea = calFeature(T(i_sub).EMG.(LR{side})(i_trial).all,window,fea);
        M_raw = mean1(T(i_sub).EMG.(LR{side})(i_trial).all(:,i_mus)',10,1);
        M_filt = T(i_sub).EMG.(LR{side})(i_trial).M(i_mus,:);
        %         M_filt = 0;
        subplot(5,5,k)
        hold on;
%         plot(M_raw)
        plot(M_fea(:,i_mus))
        title(['s:' num2str(i_sub) ' t:' num2str(i_trial)])
        
        
        k=k+1;
        %         suptitle('')
    end
    
end

suptitle([file ' - ' allMuscleNames{i_mus} ' - ' (LR{side})])

%% plot mean of feature for all muscles
k=1;
h1=figure;
set(h1,'units','normalized','outerposition',[0 0 1 1]);
hold on;
fea = 'RMS';
window = 100;

for i_sub=subject
    disp(['subject:',num2str(i_sub)])
    
    trialNum = length(T(i_sub).Trajectory.Right);
    
    for i_trial=1:length(T(i_sub).Trajectory.Right)
        
        
        M_fea = calFeature(T(i_sub).EMG.(LR{side})(i_trial).all,window,fea);
%         M_raw = mean1(T(i_sub).EMG.(LR{side})(i_trial).all(:,i_mus)',10,1);
%         M_filt = T(i_sub).EMG.(LR{side})(i_trial).M(i_mus,:);
        %         M_filt = 0;
%         subplot(5,5,k)
%         hold on;
%         plot(M_raw)
for i_mus =1:16
    str = [file ' s:' i_sub ' t:' i_trial ' m:' allMuscleNames{i_mus} ' ' LR{side}];
    plot(k,mean(M_fea(:,i_mus)),'k*','DisplayName',str)
%     errorbar(k,mean(M_fea(:,i_mus)),std(M_fea(:,i_mus)),'DisplayName',num2str(i_mus))
%         title(['s:' num2str(i_sub) ' t:' num2str(i_trial)])
end
        
        k=k+1;
        %         suptitle('')
    end
    
end

%% clustering features !!!
features = {'ZC' 'RMS' 'WL' 'SSC'};
n = length(features);

k=1;
mat = [];
TYPE=[];
type_sub = [];
type_mus = [];
for i_sub=subject
    disp(['subject:',num2str(i_sub)])
    
    trialNum = length(T(i_sub).Trajectory.Right);
    
    for i_trial=1:length(T(i_sub).Trajectory.Right)
        temp = [];
        temp1 = [];
        temp2 = [];
        for i_n = 1:n
             M_fea = calFeature(T(i_sub).EMG.(LR{side})(i_trial).all,window,features{i_n});
            temp = [temp; mean(M_fea); std(M_fea)];
        end
        TYPE(k).subject = repmat(i_sub,16,1);
        TYPE(k).id = repmat(1,16,1);
        TYPE(k).muscle = allMuscleNames;
        mat = [mat temp];
        k=k+1;
        %         suptitle('')
    end
    
end
%% tSNE plot
Y = tsne(mat','Algorithm','exact','Distance','seuclidean');
tSNEPlot(Y, TYPE, 'muscle')


%% min mean std of signals
for i_mus=i_mus
    figure
    k=1;
    for i_sub=subject
        disp(['subject:',num2str(i_sub)])
        
        trialNum = length(T(i_sub).Trajectory.Right);
        
        for i_trial=1:length(T(i_sub).Trajectory.Right)
            
            
            
            M_raw = T(i_sub).EMG.(LR{side})(i_trial).all(:,i_mus);
            M_filt = T(i_sub).EMG.(LR{side})(i_trial).M(i_mus,:);
            %         M_filt = 0;
            
            
            subplot(2,1,1); hold on;
            plot(k,min(M_filt),'r*')
            plot(k,mean(M_filt),'bs')
            errorbar(k,mean(M_filt),std(M_filt),'k')
            plot(k,min(M_filt),'ko')
            plot(k,max(M_filt),'ko')
            
            subplot(2,1,2); hold on;
            plot(k,mean(M_raw),'bs')
            errorbar(k,mean(M_raw),std(M_raw),'k')
            plot(k,min(M_raw),'ko')
            plot(k,max(M_raw),'ko')
            
            k=k+1;
            %         suptitle('')
        end
        
    end
    suptitle([file ' - ' allMuscleNames{i_mus} ' - ' (LR{side})])
end