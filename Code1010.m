clc
clear
close all
%% calculate balance index
% COP is matrix of m*3 where m is number of frames(freq=1200)
% COM is matrix of m*3 where m is number of frames(freq=120)
    
        
BI = calIndex(COP, COM);
        

    
%% function
function BI = calIndex(COP, COM)
COPx = COP(:,1) - mean(COP(:,1));
COPy = COP(:,2) - mean(COP(:,2));
COMx = COM(:,1) - mean(COM(:,1));
COMy = COM(:,2) - mean(COM(:,2));
COMz = COM(:,3) - mean(COM(:,3));

COPr = sqrt(COPx .^2 + COPy .^2);
r90CI = sort(COPr);
BI.r90CI = r90CI(floor(.9 * length(COPr)));
BI.COPx_std = std(COPx);
BI.COPy_std = std(COPy);
BI.COMx_std = std(COMx);
BI.COMy_std = std(COMy);
BI.COMz_std = std(COMz);


BI.COPrv_std = std(diff(COPr)*1200);


end
