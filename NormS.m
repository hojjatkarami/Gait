function [t,Snew] = NormS(S,Event)
Snew=[];
% if i_test==5 || i_test==6    
%     Norm  = 0:100/4:100;
Norm=[0 11.5 50 60 100];
% elseif i_test==4
%     Norm  = [0 100];
% end

for j=1:size(S,1)
    dim=length(S(j,:));
    [t, Snew(j,:)] = rescale111(1:dim,S(j,:),Event,Norm,100);
    
    
end