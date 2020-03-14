function [W_best, S_best] = MyNNMF0510(M, std_val, syn, option)



clear vaf rsq VAF RSQ r2_bootstat_lb r2_bootstat_ub vaf_bootstat_lb vaf_bootstat_ub
   
        
        W_all=[];
        S_all=[];

%         for u1=1:option.rep1

            
            maxmin_rsq = 0;

%             for u2 = 1:option.rep2

                opt = statset('MaxIter',1000,'Display','off');
                [W0,S0] = nnmf(M,syn,'replicates',10,'options',opt,'algorithm','mult');
                opt = statset('MaxIter',1000,'Display','off');
                [W,S] = nnmf(M,syn,'w0',W0,'h0',S0,'options',opt,'replicates',20,'algorithm','als');
                rec = (W*S);
                W_all = [W_all, W];
                
%                 VAF(u2) =  vaf1(rec, M, 0);
%                 RSQ(u2) =  rsq1(rec, M, 0);             
%                 vaf(u2,:) = vaf1(rec, M, 1)';   % local vaf
%                 rsq(u2,:) = rsq1(rec, M, 1)';   % local rsq            

%                 if min(rsq(u2,:)) > maxmin_rsq || (u2==option.rep2)
%                       maxmin_rsq = min(rsq(u2,:));                  
                      W_best = W;
                      S_best = S;
%                       vaf_best = vaf(u2,:);
%                       rsq_best = rsq(u2,:);
%                       VAF_best = VAF(u2);
%                       RSQ_best = RSQ(u2);
            
%                 end
%             end

% bootstrap

% rec = W_best * S_best;
% [r2_bootstat_lb, r2_bootstat_ub] = myBootStrap(option.n_bootstrap,'@rsq1',M,rec);
% [vaf_bootstat_lb, vaf_bootstat_ub] = myBootStrap(option.n_bootstrap,'@vaf1',M,rec);

% criteria = best_VAF(i_syn)>=option.VAF_th || min(best_vaf(i_syn,:)) >= option.vaf_th || min(best_rsq(i_syn,:))>=option.rsq_th ||vaf_bootstat_lb(end)> option.vafCI_th || r2_bootstat_lb(end)> option.rsqCI_th || option.condition == 'none';

    


% Saving the results ...

        
        if strcmp(option.type, 'emgkin')
            W_best = W_best .* std_val;
            c1 = max(W_best(1:16,:));
            c2 = max(W_best(17:28,:));
            W_best(1:16,:)= W_best(1:16,:)./c1;
            W_best(17:28,:)= W_best(17:28,:)./c2;
            
        else
            W_best = W_best .* std_val;
            c = max(W_best);
            W_best = W_best ./ c;
            S_best = S_best .* repmat(c',1,size(S_best,2));

        end
        
        
        res.W_best = W_best;
        res.S_best = S_best;
        res.W_all = W_all;
        res.S_all = S_all;
        res.SynNum = syn; 
%         res.VAF = VAF_best;
%         res.RSQ = RSQ_best;
%         res.vaf = vaf_best;
%         res.rsq = rsq_best;
% 
%         res.r2_bootstat_ub = r2_bootstat_ub;
%         res.r2_bootstat_lb = r2_bootstat_lb;
%         res.vaf_bootstat_ub = vaf_bootstat_ub;
%         res.vaf_bootstat_lb = vaf_bootstat_lb; 





