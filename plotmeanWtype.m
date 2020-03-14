function plotmeanWtype(h1,W,Wtype,id,name,label,supTitle)
n=length(W);
N=length(Wtype);
rgb = maxdistcolor(N+3,@srgb_to_Lab);    % create maximum dixtinct color
for i=1:N
    meanOfType{i} =  mean(Wtype{i},2);
    stdOfType{i} = std(Wtype{i},0,2);
%     bad = find((meanOfType{i}-2*stdOfType{i})<0);
%     meanOfType{i}(bad) = 0;
end

% h1=figure(1001);
set(h1,'name',supTitle,'units','normalized','outerposition',[0 0 1 1]);% suptitle(supTitle)

for i=1:N
    
    subplot(N,1,i); hold on; 
    bar(meanOfType{i},'FaceColor',rgb(i,:))
    errorbar( meanOfType{i},  1*stdOfType{i},'LineStyle','none','Color',rgb(end,:))
    plot(Wtype{i},'k.','LineStyle','none')
    

    xticks(1:length(label))
    xtickangle(45)
    xticklabels(label)
%     a = get(gca,'XTickLabel');    set(gca,'XTickLabel',a,'fontsize',6)
end



