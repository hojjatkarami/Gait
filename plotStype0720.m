function plotStype0720(h2,S,Stype,id,name,label,supTitle)
% 0720 version: t variable removed
n=length(S);
N=length(Stype);
rgb = maxdistcolor(N+3,@srgb_to_Lab);    % create maximum dixtinct color
% for i=1:N
%     meanOfType{i} =  mean(Wtype{i},2);
%     stdOfType{i} = std(Wtype{i},0,2);
% %     bad = find((meanOfType{i}-2*stdOfType{i})<0);
% %     meanOfType{i}(bad) = 0;
% end
% h2=figure(1002);
set(h2,'name',supTitle,'units','normalized','outerposition',[0 0 1 1]);
% suptitle(supTitle)


for i=1:N
    for j=1:n        
        z = find(id{j} == i);
        subplot(N,n+1,(i-1)*(n+1) + j); hold on;
        set(gca,'XTick',[],'YTick',[],'Tag',num2str((i-1)*(n+1) + j),'ButtonDownFcn',@Click_Stype) 
        ylim([0 6])
        if (i==1)
           title(name{j})           
        end
        if (j==1)
           ylabel(['Type ',num2str(i)])
        end
        if ~isempty(z)
%             dim=length(S{j}(z,:));
%             [t2, y2] = rescale111(1:dim,S{j}(z,:),Event{j},Norm{j},1000);
            plot(S{j}(z,:),'Color',rgb(i,:))
%             xlabel(num2str(corr(W{j}(:,t),mean(Wtype{i},2)),'%.2f'))
        end    
    end
    subplot(N,n+1,(i)*(n+1) ); hold on;
    set(gca,'XTick',[],'YTick',[],'Tag',num2str((i)*(n+1)),'ButtonDownFcn',@Click_Stype2)
    ylim([0 6])
    title(['Mean'])
%     plot([0:.1:99.9],mean(Stype{i},1),'Color',rgb(i,:))
    plot([0:1:99],mean(Stype{i},1),'Color',rgb(i,:))
    x=[[0:1:99],fliplr([0:1:99])];
%     x=[[0:.1:99.9],fliplr([0:.1:99.9])];
    y=[mean(Stype{i},1)-std(Stype{i},0,1),...
    fliplr(mean(Stype{i},1)+std(Stype{i},0,1))];
    s=fill(x,y,'k','EdgeColor','none');
    alpha(s,.1)
%     bar(meanOfType{i},'FaceColor',rgb(i,:))
%     errorbar( meanOfType{i},  2*stdOfType{i},'LineStyle','none','Color',rgb(end,:))
%     plot(Wtype{i},'k*','LineStyle','none')
%     xticks(1:16)
%     xtickangle(45)
%     xticklabels(label)
end
