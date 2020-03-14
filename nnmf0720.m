clear
close all
clc

nnmf_init

% all changes of this script will be applied to Synergy field of 'P' Database
%% Settings
load('P.mat')
subject = [11];
SYN = 1:8;
option.side = 'Right';  % Right: right side, Left: left side; 
option.NormalizeToUnitVariance = 1;     % 0:do not divide, 1:divide and then multiply, 3:divide and do not multiply
option.rep = 10;
%% EMG>syn : extract n1:n2 number of synergy


for i_sub = subject
    TrialNum =length( P(i_sub).EMG.(option.side) );
    for i_trial = 1:TrialNum
        M = P(i_sub).EMG.(option.side)(i_trial).Mn;   
        for i_syn = SYN
            disp(['sub:',num2str(i_sub),', trial:',num2str(i_trial),', syn:',num2str(i_syn),', ' ])
            [W_best, S_best, M_rec] = nnmfEMG0720(M,i_syn,option);
            P(i_sub).Synergy.EMG.(option.side)(i_trial).syn(i_syn).W_best = W_best;
            P(i_sub).Synergy.EMG.(option.side)(i_trial).syn(i_syn).S_best = S_best;
            P(i_sub).Synergy.EMG.(option.side)(i_trial).syn(i_syn).M_rec = M_rec;
        end

    end
    
end
save('P.mat','P');

%% kin>syn

for i_sub = subject
    TrialNum =length( P(i_sub).EMG.(option.side) );
    for i_trial = 1:TrialNum
        M = P(i_sub).kin.(option.side)(i_trial).Mn; 

        for i_syn = SYN
            disp(['sub:',num2str(i_sub),', trial:',num2str(i_trial),', syn:',num2str(i_syn),', ' ])
            [W_best, S_best, M_rec] = nnmfKin0720(M,i_syn,option);
            P(i_sub).Synergy.kin.(option.side)(i_trial).syn(i_syn).W_best = W_best;
            P(i_sub).Synergy.kin.(option.side)(i_trial).syn(i_syn).S_best = S_best;
            P(i_sub).Synergy.kin.(option.side)(i_trial).syn(i_syn).M_rec = M_rec;
        end

    end
end
save('P.mat','P');
%% EMGkin>syn

for i_sub = subject
    TrialNum =length( P(i_sub).EMG.(option.side));
    for i_trial = 1:TrialNum

        M_EMG = P(i_sub).EMG.(option.side)(i_trial).Mn; 
        M_kin = P(i_sub).kin.(option.side)(i_trial).Mn; 

        for i_syn = SYN
            
            disp(['sub:',num2str(i_sub),', trial:',num2str(i_trial),', syn:',num2str(i_syn),', ' ])
            [W_best, S_best, M_rec] = nnmfEMGKin0720(M_EMG,M_kin,i_syn,option);
            P(i_sub).Synergy.EMGkin.(option.side)(i_trial).syn(i_syn).W_best = W_best;
            P(i_sub).Synergy.EMGkin.(option.side)(i_trial).syn(i_syn).S_best = S_best;
            P(i_sub).Synergy.EMGkin.(option.side)(i_trial).syn(i_syn).M_rec = M_rec;
        end

    end

end
save('P.mat','P');

%% EMG>gof
SYN = 1:8;
option.side = 'Right';  % Right: right side, Left: left side;
option.type = 'kin';
for i_sub = subject
    i_sub
    TrialNum =length(P(i_sub).EMG.Right);
    for i_trial = 1:TrialNum
        M = P(i_sub).(option.type).Right(i_trial).Mn;
        
        for i_syn = SYN
            M_rec = P(i_sub).Synergy.(option.type).(option.side)(i_trial).syn(i_syn).M_rec;
            P(i_sub).Synergy.(option.type).(option.side)(i_trial).gof.VAF(i_syn) = vaf1(M_rec, M, 0);
            P(i_sub).Synergy.(option.type).(option.side)(i_trial).gof.RSQ(i_syn) = rsq1(M_rec, M, 0);
            P(i_sub).Synergy.(option.type).(option.side)(i_trial).gof.vaf(i_syn,:) = vaf1(M_rec, M, 1);
            P(i_sub).Synergy.(option.type).(option.side)(i_trial).gof.rsq(i_syn,:) = rsq1(M_rec, M, 1);
            
            [r2_bootstat_lb, r2_bootstat_ub] = myBootStrap(100,95,'@rsq1',M,M_rec);
            [vaf_bootstat_lb, vaf_bootstat_ub] = myBootStrap(100,95,'@vaf1',M,M_rec);
            P(i_sub).Synergy.(option.type).(option.side)(i_trial).gof.BootStrap.r2_bootstat_lb(i_syn) = r2_bootstat_lb;
            P(i_sub).Synergy.(option.type).(option.side)(i_trial).gof.BootStrap.r2_bootstat_ub(i_syn) = r2_bootstat_ub;
            P(i_sub).Synergy.(option.type).(option.side)(i_trial).gof.BootStrap.vaf_bootstat_lb(i_syn) = vaf_bootstat_lb;
            P(i_sub).Synergy.(option.type).(option.side)(i_trial).gof.BootStrap.vaf_bootstat_ub(i_syn) = vaf_bootstat_ub;
            
            % Factor Analysis
        end
    end
end
save('P.mat','P');

%% EMG>number
SYN = 1:8;
option.side = 'Right';  % Right: right side, Left: left side;
option.type = 'kin';
vaf_th=0.9;    VAF_th = 0.95;
rsq_th=0.6;     RSQ_th = 0.8;
rsq_boot_th = 0.8;     vaf_boot_th=0.95;

for i_sub = subject
    i_sub
    TrialNum =length(P(i_sub).EMG.Right);
    for i_trial = 1:TrialNum
        gof = P(i_sub).Synergy.(option.type).Right(i_trial).gof  ;
        

        P(i_sub).Synergy.(option.type).(option.side)(i_trial).number.VAF = min(find(gof.VAF>=VAF_th));
        P(i_sub).Synergy.(option.type).(option.side)(i_trial).number.RSQ = min(find(gof.RSQ>=RSQ_th));
        P(i_sub).Synergy.(option.type).(option.side)(i_trial).number.vaf = min(find(min(gof.vaf')>=vaf_th));
        P(i_sub).Synergy.(option.type).(option.side)(i_trial).number.rsq = min(find(min(gof.rsq')>=rsq_th));
        if isempty(P(i_sub).Synergy.(option.type).(option.side)(i_trial).number.rsq)
            P(i_sub).Synergy.(option.type).(option.side)(i_trial).number.rsq=8;
        end

        P(i_sub).Synergy.(option.type).(option.side)(i_trial).number.BootStrap_r2 = min(find(gof.BootStrap.r2_bootstat_lb>=rsq_boot_th));
        P(i_sub).Synergy.(option.type).(option.side)(i_trial).number.BootStrap_vaf = min(find(gof.BootStrap.vaf_bootstat_lb>=vaf_boot_th));

        % Factor Analysis
        
    end
end
save('P.mat','P');


