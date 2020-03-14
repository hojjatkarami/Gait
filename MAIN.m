clc
clear
close all
nnmf_init

%% save


%%  header
subject = [3 4];
nnmfField = {'period2' 'period5' 'period10' 'period20'};
nnmfID = [1 2 3 4];
SYN = 6;
file='T2';
f = 120;
L =[1 2 5 10 20];


%% ISB
load(['raw Data\' file '.mat'])
eval(['T = ',file,';']);

T = ISB(T, subject);

eval([file,' = T;']);
save(['raw Data\' file '.mat'],file)

%% KIN2mat

load(['raw Data\' file '.mat'])
eval(['T = ',file,';']);

angles = {[1 2],[1],[1 2 3],[1 2 3]};   % {ankle, knee, hip, pelvis}
T = KIN2mat(T, subject);

eval([file,' = T;']);
save(['raw Data\' file '.mat'],file)

%% EMG2mat
highPassFreq = 40;
lowPassFreq = 15;
timeBin = 10;

load(['raw Data\' file '.mat'])
eval(['T = ',file,';']);

T = EMG2mat(T, subject, highPassFreq, lowPassFreq, timeBin);

eval([file,' = T;']);
save(['raw Data\' file '.mat'],file)

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

%% data visualization

data1={}; k=1;
data2={};
for i_sub = subject
    
    trialNum = length(T1(i_sub).ForcePlate);
    for i_trial = 1:trialNum
        data1{i_trial,1} = T1(i_sub).Para(i_trial).COPx;
        data1{i_trial,2} = T1(i_sub).Para(i_trial).COPy;
        data1{i_trial,3} = T1(i_sub).Para(i_trial).tiltAP;
        data1{i_trial,4} = T1(i_sub).Para(i_trial).tiltML;
        
        
    end
    %     k=k+1;
    %     for i_trial = 1:trialNum
    %         data1{k,i_trial} = T2(i_sub).ForcePlate(i_trial).COP;
    %
    %
    %     end
    %     k=k+1;
end
for i_trial=1:size(data1,1)
    
    %     for i_trial = 1:size(data1,2)
    figure(1)
    for i_para = 1:size(data1,2)
        if isempty(data1{i_trial,i_para})
            continue;
        end
        subplot(size(data1,1),size(data1,2),(i_trial-1)*size(data1,2)+i_para); hold on;
        plot(data1{i_trial,i_para})
        
        
    end
    figure(2)
    subplot(size(data1,1),1,i_trial); hold on;
    plot(data1{i_trial,1}, data1{i_trial,2})
    axis equal
    xlim([-100 100])
    ylim([-100 100])
end



%% plot W S
file='T2';
eval(['T = ',file,';']);
option.side = 'Right';  % Right: right side, Left: left side;
X=[];
TYPE={};
Z=[];
nnmfField = {'period2' 'period5' 'period10' 'period20'};
nnmfID = [2];
i_syn=4;

for i_sub=subject
    trialNum = length(T(i_sub).ForcePlate);
    for i_trial = 1:trialNum
        for i_period=nnmfID
            M_t = T(i_sub).EMG.(option.side)(i_trial).M;
            start = T(i_sub).Synergy(i_trial).(nnmfField{i_period}).start;
            finish = T(i_sub).Synergy(i_trial).(nnmfField{i_period}).finish;
            for i_interval = 1:length(start)
                W = T(i_sub).Synergy(i_trial).(nnmfField{i_period}).interval(i_interval).syn(i_syn).W_best;
                S = T(i_sub).Synergy(i_trial).(nnmfField{i_period}).interval(i_interval).syn(i_syn).S_best;
                X=[X W];
                TYPE = [TYPE;repmat({file},size(W,2),1)];
                Z=[Z sum(mean(S'))];
            end
            
        end
        
    end
    
end
%% t-SNE
h=figure;
set(h,'units','normalized','outerposition',[0 0 1 1]);

distance = {'euclidean'  'seuclidean'  'cityblock'  'chebychev'...
    'minkowski'  'mahalanobis'  'cosine'  'correlation'...
    'spearman'  'hamming'  'jaccard'};
distance = {'euclidean'   'cityblock' ...
    'minkowski' 'cosine'  'correlation'...
    'spearman'};
for i_dis = 1:length(distance)
    subplot(3,4,i_dis)
    Y = tsne(X','Algorithm','exact','Distance',distance{i_dis});
    gscatter(Y(:,1),Y(:,2),TYPE)
    title(distance{i_dis})
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
