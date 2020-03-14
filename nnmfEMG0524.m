function [W_best, S_best,W_avg,S_all] = nnmfEMG0524(Mat,syn,option,i_trial)

switch option.side
    case 'Right'
        M = Mat.Right(i_trial).M_R;                 
    case 'Left'
        M =Mat.Left(i_trial).M_L;
%     case 'Both'
%         M = [h_trial.emg.M_R ;...
%              h_trial.emg.M_L];
end
% M=[];
% for i=1:10
%     M=[M Mat.Right(i).M_R ];
% end
std_val =  std(M,0,2)+eps;  
M = M ./ std_val;

%%
clear vaf rsq VAF RSQ r2_bootstat_lb r2_bootstat_ub vaf_bootstat_lb vaf_bootstat_ub
W_all = {};
S_all = {};   
        name={};


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

% Saving the results ...


if strcmp(option.type, 'kinemg')
W_best = W_best .* std_val;
c1 = max(W_best(1:16,:));
c2 = max(W_best(17:28,:));
W_best(1:16,:)= W_best(1:16,:)./c1;
W_best(17:28,:)= W_best(17:28,:)./c2;

else
%     std_val=1;
W_best = W_best .* std_val;
c = max(W_best);
W_best = W_best ./ c;
S_best = S_best .* repmat(c',1,size(S_best,2));

end

    W_all{u2} = W_best;
    S_all{u2} = S_best;
    name{u2}=num2str(u2);
end


% cluster W
[id,Wtype] = kmeanW(W_all, syn);
% h1=figure;
% label = 'Ankle Dorsi Flex','Ankle Plant Flex',...
%                    'Ankle Inv','Ankle Evr','Knee Flex','Knee Ext',...
%                 '   Hip Flex','Hip Ext','Hip Add','Hip Abd','Hip Int Rot','Hip Ext Rot'};
% plotWtype(h1,W_all,Wtype,id,name,label,'nnmf')
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












