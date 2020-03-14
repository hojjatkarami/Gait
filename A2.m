
subject = [3];
file='T1';
sf = 120;
folderToSave = 'raw Data5';

load([folderToSave '\' file '.mat'])
eval(['T = ',file,';']);





%% extract period trials [2 5 10 20]

eval(['T = ',file,';']);
for i_sub = subject
    trialNum = length(T(i_sub).ForcePlate);
    for i_trial = 1:trialNum
        
        n = size(T(i_sub).EMG.Right(i_trial).M, 2);
        temp = [1:2*sf:n];
        T(i_sub).Synergy(i_trial).period2.start = temp(1:end-1);
        T(i_sub).Synergy(i_trial).period2.finish = temp(2:end)-1;
        
        temp = [1:5*sf:n];
        T(i_sub).Synergy(i_trial).period5.start = temp(1:end-1);
        T(i_sub).Synergy(i_trial).period5.finish = temp(2:end)-1;
        
        temp = [1:10*sf:n];
        T(i_sub).Synergy(i_trial).period10.start = temp(1:end-1);
        T(i_sub).Synergy(i_trial).period10.finish = temp(2:end)-1;
        
        temp = [1:20*sf:n];
        if temp==1
            temp=[1,n];
        end
        T(i_sub).Synergy(i_trial).period20.start = temp(end-1);
        T(i_sub).Synergy(i_trial).period20.finish = temp(end)-1;
        
        
    end
end
eval([file,' = T;']);


%% NNMF EMG
nnmfField = {'period2' 'period5' 'period10' 'period20'};
nnmfID = [1 2 3 4];
eval(['T = ',file,';']);
option.side = 'Right';  % Right: right side, Left: left side;
option.NormalizeToUnitVariance = 1;     % 0:do not divide, 1:divide and then multiply, 3:divide and do not multiply
option.rep = 1;
option.mus = [];
SYN = [2 3 4 5];

for i_sub = subject
    trialNum = length(T(i_sub).ForcePlate);
    for i_trial = 1:trialNum
        for i_period=nnmfID
            
            
            M_t = T(i_sub).EMG.(option.side)(i_trial).M;
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
save([folderToSave '\' file '.mat'],file)


%% intra_period
type = 'EMG';
subjects = [3];
period = {'period2' 'period5' 'period10' 'period20'};
% period = {'period2' 'period5'}
trial=1;
toPlot=0;
SYN = 5;
N=5;
[T, W, name, Wtype, id] = intra_period(T, type, subjects, period, trial, SYN, N, toPlot);

%% inter_period
type = 'EMG';
subjects = [3 4];
period = {'period2' 'period5'};
trial=0;
toPlot=1;
SYN = 5;
N=5;
[T, W, name, Wtype, id] = inter_period(T, type, subjects, period, trial, SYN, N, toPlot);



%% inter_trial
type = 'EMG';
subjects = [3 4];
period = {'period2' 'period5' 'period10' 'period20'};
trial=0;
toPlot=1;
SYN = 5;
N=5;
[T, W, name, Wtype, id] = inter_trial(T, type, subjects, period, trial, SYN, N, toPlot);


