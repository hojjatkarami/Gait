





%% t-SNE Comp
type = 'EMG';
subjects = [3];
files = {'T1'};
period = {'period2' 'period5' 'period10' 'period20'};
% period = {'period2'};
trial=0;
toPlot=0;
SYN = 5;
N=5;
[Y, TYPE] = tSNEComp(files, type, subjects, period, trial, SYN, N, toPlot);

%% tSNE plot

labels = {'subject' 'period' 'id'};
labels = {'subject'};

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