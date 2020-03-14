function plotWtype1122(h1,W,Wstd,Wtype,id,name,label,supTitle)
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
    for j=1:n        
        t = find(id{j} == i);
        subplot(N,n+1,(i-1)*(n+1) + j); hold on;
        set(gca,'YTick',[],'Tag',num2str((i-1)*(n+1) + j),'ButtonDownFcn',{@Click_Wtype,label}) 
        if (i==1)
           title(name{j})
           
        end
        if (j==1)
           ylabel(['Type ',num2str(i)])
        end
        if ~isempty(t)
            bar(W{j}(:,t),'FaceColor',rgb(i,:))
            errorbar( W{j}(:,t),  Wstd{j}(:,t),'LineStyle','none','Color',rgb(end,:))

            xlabel(num2str(corr(W{j}(:,t),mean(Wtype{i},2)),'%.2f'))
%             if i==N
                xticks(1:length(label))
                xtickangle(45)
                xticklabels(label)
                a = get(gca,'XTickLabel');    set(gca,'XTickLabel',a,'fontsize',4);
%             end
        end    
    end
    subplot(N,n+1,(i)*(n+1) ); hold on;
    set(gca,'XTick',[],'YTick',[],'Tag',num2str((i)*(n+1)),'ButtonDownFcn',@Click_Wtype)
    title(['Mean'])
        
    bar(meanOfType{i},'FaceColor',rgb(i,:))
    errorbar( meanOfType{i},  2*stdOfType{i},'LineStyle','none','Color',rgb(end,:))
    plot(Wtype{i},'k*','LineStyle','none')
    

    xticks(1:length(label))
    xtickangle(45)
    xticklabels(label)
%     a = get(gca,'XTickLabel');    set(gca,'XTickLabel',a,'fontsize',6)
end



