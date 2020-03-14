function [W_best, S_best,W_avg,S_all] = nnmfEMGKin0720(M_EMG,M_kin,syn,option)

%process EMG
std_val1 =  std(M_EMG,0,2)+eps;
M_EMG = M_EMG ./ std_val1;


% process kinematic
M_kin = kinProcess(M_kin);
std_val2 =  std(M_kin,0,2)+eps;  
M_kin = M_kin ./ std_val2;

% put together
M = [M_EMG; M_kin];
std_val = [std_val1;std_val2];
%%
clear vaf rsq VAF RSQ r2_bootstat_lb r2_bootstat_ub vaf_bootstat_lb vaf_bootstat_ub
W_all = {};
S_all = {};   
name={};

u2=1;
while u2 <= option.rep

    opt = statset('MaxIter',1000,'Display','off');
    [W0,S0] = nnmf(M,syn,'replicates',10,'options',opt,'algorithm','mult');
    opt = statset('MaxIter',1000,'Display','off');
    [W,S] = nnmf(M,syn,'w0',W0,'h0',S0,'options',opt,'replicates',20,'algorithm','als');
    rec = (W*S);
    if min(any(W,1)) == 0
        continue;
    end           
    

% Saving the results ...



    c = max(W);
    W = W ./ c;
    S = S .* repmat(c',1,size(S,2));   

%     if u2==1
%         continue;        
%     end
       W_all{u2} = W;
    S_all{u2} = S;
    name{u2}=num2str(u2);
    u2=u2+1;
end


%% cluster W
[id,Wtype] = kmeanW(W_all, syn,0.6);

% calculate mean of each type
for i=1:syn
   meanOfType{i} = mean(Wtype{i},2) ;                
end
W_avg = cell2mat(meanOfType);
meanRepCorr=[];
% choose the best corr with W_avg
for i=1:option.rep
   
   if min(id{i})==0
       continue;
   end
   for k=1:syn
        out{i}(k) = corr(W_all{i}(:,k), meanOfType{ id{i}(k) }); 

   end
   meanRepCorr(i) = mean(out{i});
end

if isempty(meanRepCorr)
    W_best = W_all{1};
    S_best = S_all{1};
else    
    [~,bestRepID] = max(meanRepCorr);
    W_best = W_all{bestRepID};
    S_best = S_all{bestRepID};
end
    
    M_rec = (W_best .* std_val) * S_best;


% plots

% for i=1:syn
%     figure; 
%     subplot(2,1,1);hold on;
%         bar(Wtype{1,i})
%     subplot(2,1,2);hold on;
%     bar([meanOfType{i}, W_best(:,i)])
%     
% 
%     
% end


