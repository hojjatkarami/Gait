function [W_best, S_best,W_avg,S_all] = nnmfEMG0524(Mat,syn,option,i_trial)

switch option.side
    case 'Right'
            M1 = Mat.EMG.Right(i_trial).M_R;
            M2 = kinProcess(Mat.kin2.Right(i_trial).M_R);
    case 'Left'
         M1 = Mat.EMG.Left(i_trial).M_L;
         M2 = kinProcess(Mat.kin2.Left(i_trial).M_L);
end

std_val1 =  std(M1,0,2)+eps;
for i=1:size(M2,1)
    a=M2(i,:);
    b=a(a~=0);
    std_val2(i)=std(b)+eps;
    if isnan(std_val2(i))
        std_val2(i)=1;
    end
end
% std_val2 =  std(M2,0,2)+eps;
std_val = [std_val1; std_val2'];        
M = [M1; M2];
M = M ./ std_val;

%%
clear vaf rsq VAF RSQ r2_bootstat_lb r2_bootstat_ub vaf_bootstat_lb vaf_bootstat_ub
W_all = {};
S_all = {};   
        


for u2 = 1:option.rep

    opt = statset('MaxIter',1000,'Display','off');
    [W0,S0] = nnmf(M,syn,'replicates',10,'options',opt,'algorithm','mult');
    opt = statset('MaxIter',1000,'Display','off');
    [W,S] = nnmf(M,syn,'w0',W0,'h0',S0,'options',opt,'replicates',20,'algorithm','als');
    rec = (W*S);
          
    VAF(u2) =  vaf1(rec, M, 0);
    RSQ(u2) =  rsq1(rec, M, 0);             
    vaf(u2,:) = vaf1(rec, M, 1)';   % local vaf
    rsq(u2,:) = rsq1(rec, M, 1)';   % local rsq     
               
  W_best = W;
  S_best = S;
(u2);
% bootstrap

% rec = W_best * S_best;
% [r2_bootstat_lb, r2_bootstat_ub] = myBootStrap(option.n_bootstrap,'@rsq1',M,rec);
% [vaf_bootstat_lb, vaf_bootstat_ub] = myBootStrap(option.n_bootstrap,'@vaf1',M,rec);

% criteria = best_VAF(i_syn)>=option.VAF_th || min(best_vaf(i_syn,:)) >= option.vaf_th || min(best_rsq(i_syn,:))>=option.rsq_th ||vaf_bootstat_lb(end)> option.vafCI_th || r2_bootstat_lb(end)> option.rsqCI_th || option.condition == 'none';




% Saving the results ...


if strcmp(option.type, 'EmgKin')
%     std_val=1;
W_best = W_best .* std_val;
c1 = max(W_best(1:16,:));
c2 = max(W_best(17:32,:));
W_best(1:16,:)= W_best(1:16,:)./c1;
W_best(17:32,:)= W_best(17:32,:)./c2;

else
W_best = W_best .* std_val;
c = max(W_best);
W_best = W_best ./ c;
S_best = S_best .* repmat(c',1,size(S_best,2));

end


%         res.W_best = W_best;
%         res.S_best = S_best;
%         res.W_all = W_all;
%         res.S_all = S_all;
%         res.SynNum = syn; 
%         res.VAF = VAF_best;
%         res.RSQ = RSQ_best;
%         res.vaf = vaf_best;
%         res.rsq = rsq_best;
% 
%         res.r2_bootstat_ub = r2_bootstat_ub;
%         res.r2_bootstat_lb = r2_bootstat_lb;
%         res.vaf_bootstat_ub = vaf_bootstat_ub;
%         res.vaf_bootstat_lb = vaf_bootstat_lb; 
    W_all{u2} = W_best;
    S_all{u2} = S_best;
    
end


% cluster W
[id,Wtype] = kmeanW(W_all, syn);

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
   meanRepCorr(i) = mean(out{i}(k));
end
if isempty(meanRepCorr)
    W_best = W_all{1};
    S_best = S_all{1};
else
    meanRepCorr;
    [~,bestRepID] = max(meanRepCorr);
    W_best = W_all{bestRepID};
    S_best = S_all{bestRepID};
end












