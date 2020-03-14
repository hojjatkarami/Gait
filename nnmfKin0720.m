function [W_best, S_best,M_rec] = nnmfKin0720(Mat,syn,option)
M1 = kinProcess(Mat);

% for i=1:size(M,1)
%     a=M(i,:);
%     b=a(a~=0);
%     std_val(i)=std(b)+eps;
%     if isnan(std_val(i))
%         std_val(i)=1;
%     end
% end

std_val =  std(M1,0,2)+eps;  
%  std_val =  std(b,0,2)+eps;  
M = M1 ./ std_val;
%% Normalize to Unit Variance
% if option.NormalizeToUnitVariance
%     std_val =  std(M,0,2)+eps;
% else
%     std_val = ones(size(M,1),1);
% end
% M = M ./ std_val;

%%
clear vaf rsq VAF RSQ r2_bootstat_lb r2_bootstat_ub vaf_bootstat_lb vaf_bootstat_ub
W_all = {};
S_all = {};   
        

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

    c = max(W)+eps;
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


% cluster W
[id,Wtype] = kmeanW(W_all, syn,0.6);
% h1=figure;
% label = {'Ank Dorsi Flex','Ank Plant Flex',...
%                      'Ank Inv','Ank Evr','Kne Flex','Kne Ext',...
%                      'Hip Flex','Hip Ext','Hip Add','Hip Abd','Hip Int Rot','Hip Ext Rot','Pelv Upwrd Obliq', 'Pelv Dw Obliq', 'Pelv Int Rot','Pelv Ext Rot'};
%                  
% plotWtype(h1,W_all,Wtype,id,name,label,'nnmf')

% calculate mean of each type
for i=1:syn
   meanOfType{i} = mean(Wtype{i},2) ;                
end
W_avg = cell2mat(meanOfType);
meanRepCorr=[];
% for i=1:option.rep
%  if isempty(W_all{i})
%        W_all{i}=W_all{i+1}
%  end 
% end
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

M_temp = (W_best .* std_val) * S_best;
% M_rec = M_rec .* std_val;
M_rec = InvkinProcess(M_temp);

% figure(1)
% clf
% subplot(2,2,1); hold on;
%     plot(Mat(1,:))
%     plot(M_rec(1,:))
% subplot(2,2,2); hold on;
%     plot(M1(1,:))
%     plot(M_temp(1,:))
% subplot(2,2,3); hold on;
%     plot(M1(2,:))
%     plot(M_temp(2,:))
% 


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







