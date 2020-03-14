




%% EMG
period = [1 3]
%% extract period trials

eval(['T = ',file,';']);
for i_sub = subject
    trialNum = length(T(i_sub).ForcePlate);
    for i_trial = 1:trialNum
        
        n = size(T(i_sub).EMG2.Right(i_trial).M, 2);
        temp = [1:2*f:n];
        T(i_sub).Synergy(i_trial).period2.start = temp(1:end-1);
        T(i_sub).Synergy(i_trial).period2.finish = temp(2:end)-1;
        
        temp = [1:5*f:n];
        T(i_sub).Synergy(i_trial).period5.start = temp(1:end-2);
        T(i_sub).Synergy(i_trial).period5.finish = temp(2:end)-1;
        
        temp = [1:10*f:n];
        T(i_sub).Synergy(i_trial).period10.start = temp(1:end-1);
        T(i_sub).Synergy(i_trial).period10.finish = temp(2:end)-1;
        
        temp = [1:20*f:n];
        if temp==1
            temp=[1,n];
        end
        T(i_sub).Synergy(i_trial).period20.start = temp(end-1);
        T(i_sub).Synergy(i_trial).period20.finish = temp(end)-1;
        
        
    end
end
eval([file,' = T;']);


%% NNMF EMG

eval(['T = ',file,';']);
option.side = 'Right';  % Right: right side, Left: left side;
option.NormalizeToUnitVariance = 1;     % 0:do not divide, 1:divide and then multiply, 3:divide and do not multiply
option.rep = 1;
option.mus = [];
SYN = 5;

for i_sub = subject
    trialNum = length(T(i_sub).ForcePlate);
    for i_trial = 1:trialNum
        for i_period=nnmfID
            
            
            M_t = T(i_sub).EMG2.(option.side)(i_trial).M;
            start = T(i_sub).Synergy(i_trial).(nnmfField{i_period}).start;
            finish = T(i_sub).Synergy(i_trial).(nnmfField{i_period}).finish;
            for i_interval = 1:length(start)
%                 M = M_t(:,start(i_interval):min(finish(i_interval),size(M_t,2)));
                M = M_t(:,start(i_interval):finish(i_interval));

                %                 if finish(i_interval)< size(M_t,2)
                %                     M = M_t(:,start(i_interval):finish(i_interval));
                %                 else
                %
                %                     continue;
                %                 end
                for i_syn = SYN
                    disp(['sub:',num2str(i_sub),', trial:',num2str(i_trial),', syn:',num2str(i_syn),', ' ])
                    [W_best, S_best, M_rec] = nnmfEMG1122(M,i_syn,option);
                    T(i_sub).Synergy(i_trial).(nnmfField{i_period}).interval(i_interval).syn(i_syn).W_best = W_best;
                    T(i_sub).Synergy(i_trial).(nnmfField{i_period}).interval(i_interval).syn(i_syn).S_best = S_best;
                    T(i_sub).Synergy(i_trial).(nnmfField{i_period}).interval(i_interval).syn(i_syn).M_rec = M_rec;
                end
            end
            
        end
        
    end
    
end
eval([file,' = T;']);
save(['raw Data\' file '.mat'],file)

%%
file='T2'
load(['raw Data\' file '.mat'])
eval(['T = ',file,';']);
%% intra_period
type = 'EMG2';
subjects = [3];
period = {'period2' 'period5' 'period10' 'period20'};
% period = {'period2' 'period5'}
trial=0;
toPlot=0;
SYN = 5;
N=5;
[T, W, name, Wtype, id] = intra_period(T, type, subjects, period, trial, SYN, N, toPlot);

%% inter_period
type = 'EMG2';
subjects = [3 4];
period = {'period2' 'period5'};
trial=0;
toPlot=1;
SYN = 5;
N=5;
[T, W, name, Wtype, id] = inter_period(T, type, subjects, period, trial, SYN, N, toPlot);



%% inter_trial
type = 'EMG2';
subjects = [3 4];
period = {'period2' 'period5' 'period10' 'period20'};
trial=0;
toPlot=1;
SYN = 5;
N=5;
[T, W, name, Wtype, id] = inter_trial(T, type, subjects, period, trial, SYN, N, toPlot);


%% t-SNE Comp
type = 'EMG2';
subjects = [3];
files = {'T1' 'T2'};
period = {'period2' 'period5' 'period10' 'period20'};
period = {'period2'};
trial=0;
toPlot=1;
SYN = 5;
N=5;
[Y, TYPE] = tSNEComp(files, type, subjects, period, trial, SYN, N, toPlot);

%% tSNE plot

labels = {'subject' 'period' 'id'};
for i=1:length(labels)
    tSNEPlot(Y, TYPE, labels{i})
end

%%
for i=1:5
    subplot(5,2,2*i-1)
    bar(T2(3).Synergy(1).period2.interval(1).syn(5).W_best(:,i))
    subplot(5,2,2*i)
    bar(T3(3).Synergy(1).period2.interval(1).syn(5).W_best(:,i))
end