function [M, std_val] = preNNMF0510(h_trial,i_trial, option)

switch option.type
    case 'emg'
        switch option.side
            case 'r'
                M = h_trial.EMG.Right(i_trial).M_R;                 
            case 'l'
                M = h_trial.EMG.Left(i_trial).M_L;
            case 'b'
                M = [h_trial.EMG.Right(i_trial).M_R ;...
                     h_trial.EMG.Left(i_trial).M_L];
        end
        std_val =  std(M,0,2)+eps;  
        M = M ./ std_val;
    case 'kin'
         switch option.side
            case 'r'
                M = kinProcess(h_trial.kin2.Right(i_trial).M_R);
            case 'l'
                M = kinProcess(h_trial.kin2.Right(i_trial).M_L);
            case 'b'
                M = kinProcess([h_trial.kin2.Right(i_trial).M_R ;...
                                h_trial.kin2.Right(i_trial).M_L]);
         end
        std_val =  std(M,0,2)+eps;  
        M = M ./ std_val;
    
    case 'emgkin'
        switch option.side
            case 'r'
                M1 = h_trial.EMG.Right(i_trial).M_R;
                M2 = kinProcess(h_trial.kin2.Right(i_trial).M_R);                
            case 'l'
                M1 =h_trial.EMG.Left(i_trial).M_L;
                M2 = kinProcess(h_trial.kin2.Right(i_trial).M_L);
            case 'b'
                M1 = [h_trial.EMG.Right(i_trial).M_R ;...
                     h_trial.EMG.Left(i_trial).M_L];
                M2 =  kinProcess([h_trial.kin2.Right(i_trial).M_R ;...
                                h_trial.kin2.Right(i_trial).M_L]);
        end
        std_val1 =  std(M1,0,2)+eps;
        std_val2 =  std(M2,0,2)+eps;
        std_val = [std_val1; std_val2];
        
        M = [M1; M2];
        M = M ./ std_val;
end

switch option.condition
    case 'VAF'        
        option.vaf_th = 1;  
        option.corr_th = 1;
        option.rsq_th = 1;
        option.RSQ_th = 1;
        option.rsqCI_th = 1;
        option.vafCI_th = 1;
    case 'vaf'
        option.VAF_th =1;  
        option.corr_th =1;
        option.rsq_th =1;
        option.RSQ_th = 1;
        option.rsqCI_th = 1;
        option.vafCI_th = 1;
    case 'rsq'
        option.VAF_th = 1;  
        option.vaf_th = 1;
        option.corr_th = 1;
        option.RSQ_th = 1;
        option.rsqCI_th = 1;
        option.vafCI_th = 1;
    case 'corr'
        option.VAF_th =1;  
        option.vaf_th =1;
        option.rsq_th =1;
        option.RSQ_th = 1;
        option.rsqCI_th = 1;
        option.vafCI_th = 1;
    case 'vafCI'
        option.vaf_th =1;
        option.rsq_th =1;
        option.RSQ_th = 1;       
        option.corr_th = 1;
        option.rsqCI_th = 1;
        option.VAF_th = 1;
    case 'rsqCI'
        option.vaf_th =1;
        option.rsq_th =1;
        option.VAF_th = 1;       
        option.corr_th = 1;
        option.RSQ_th = 1;
        option.vafCI_th = 1;
end