function plot_open(cb,evendata)
%cb is the handle of the axes that was clicked
%click on the whitespace within and axes and not on the line object
id = round(cb.CurrentPoint(1,1));
if id<1 || id>16
    return
end
figure;
subplot(2,2,1:2)
hold on
plot(M{1}(id,:),'linewidth',1.5,'DisplayName','Original')
title(label{id})

vaf=[];
for i=1:length(M)
    sig = rec{i}(id,:);
%     sig = sig/max(sig);
    plot(sig,'DisplayName',[num2str(i),' Synergy'])
    vaf(i) = vaf1(M{i}(id,:),rec{i}(id,:),1);
    rsq(i) = rsq1(M{i}(id,:),rec{i}(id,:),1);

end
legend
subplot(2,2,3)
plot(vaf)
subplot(2,2,4)
plot(rsq)
