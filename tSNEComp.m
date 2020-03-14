function [Y, TYPE] = tSNEComp(files, type, subjects, period, trials, SYN, N, toPlot)

W={}; S={}; name={}; k=1;  Wstd={}; TYPE={};
for i_file = 1:length(files)
    load(['raw Data\' files{i_file} '.mat'])
    eval(['T = ',files{i_file},';']);
    f_label = [files{i_file} '/'];
    
    for i_sub = subjects
        if trials == 0
            allTrials = 1:length(T(i_sub).Synergy);
        end
        
        for i_trial = allTrials
            
            for i_period = 1:length(period)
                
                for i_int = 1:length(T(i_sub).Synergy(i_trial).(period{i_period}).interval)
                    for i_syn = SYN
                        W{k} = T(i_sub).Synergy(i_trial).(period{i_period}).interval(i_int).syn(i_syn).W_best;
                        TYPE(k).subject = repmat({[f_label,num2str(i_sub)]},1,size(W{k},2));
                        TYPE(k).period = repmat({[f_label,period{i_period}]},1,size(W{k},2));
                        
                        %             Wstd{k} = T(i_sub).Synergy(i_trial).(period{i_period}).Wstd;
                        name{k} = ['s', num2str(i_sub)...
                            '/t', num2str(i_trial)...
                            '/', period{i_period}...
                            '/r'];
                        k=k+1;
                    end
                end
            end
        end
        
        %     t=length(T(i_sub).Synergy);
        %     T(i_sub).Synergy(t+1).(period{i_period}).Wtype = Wtype;
        %
        %     T(i_sub).Synergy(t+1).(period{i_period}).Wmean = Wmean;
        %     T(i_sub).Synergy(t+1).(period{i_period}).id = id;
        
        
    end    
end
[id,Wtype,Wmean, ~] = kmeanW1122(W,N,0.1);
for i=1:length(W)
    TYPE(i).id = id{i};
    
end

disp('done')

%%
distance = {'euclidean'  'seuclidean'  'cityblock'  'chebychev'...
    'minkowski'  'mahalanobis'  'cosine'  'correlation'...
    'spearman'  'hamming'  'jaccard'};
distance = {'euclidean'   'cityblock' ...
    'minkowski' 'cosine'  'correlation'...
    'spearman'};
distance = {'seuclidean'};
%% 13) plot t-SNE
allColor = 'mcrgby';
X=[];
% N=6;
rgb = maxdistcolor(N+3,@srgb_to_Lab);    % create maximum dixtinct color
plotColor = rgb(1:N,:);
% type by predifined clustering
for i=1:length(W)
    X= [X W{i}];
    
    %                 TYPE = [TYPE;repmat(1,size(W{i},2),1)];
    
end
t= sum(X,2);
t= find(t==0);
X(t,:)=[];
Y = tsne(X','Algorithm','exact','Distance',distance{end});
disp('tSNECmp done')
