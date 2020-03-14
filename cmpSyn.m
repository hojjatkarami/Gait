function cmpSyn(data,rep,label,partition)
    n = length(data);   % number of sets
    m=[];
    for i=1:n
        m(i) = size(data{i},2); % number of vectors in each set
    end
    corrMat={};
    for i=2:n
        for j=1:i-1
            corrMat{i,j} = NDP2(data{i},data{j});
            
        end
    end
    type=1;
    id = zeros(n,max(m));
    for i=1:n
        for k=1:m(i)
            
            if id(i,k)~=0
                continue;
            end                
            id(i,k)=type;
            for j=i+1:n     % compare with vectors of next set       
                
                t = find(corrMat{j,i}(:,k) ~= 0);
                if length(t)>1  % if to vector of one set is detected then choose the best one
                   max_corr = max(corrMat{j,i}(:,k));
                   t = find(corrMat{j,i}(:,k) == max_corr);
                   t2 = find(corrMat{j,i}(:,k) ~= max_corr);
                   corrMat{j,i}(t2,k)=0;
                end
                id(j,t)=type;
                if (id(j,t)~=0)
                    if length(find(id(i,:)==id(j,t)))==0
                       id(i,k) = id(j,t);
                       type = type-1;
                       break;
                    end
                end
            end
            type=type+1;
        end
    end
%     Y = tsne(vec,'Algorithm','exact','Distance','euclidean');
%     plot(Y(:,1),Y(:,2),'*','MarkerSize',10,'DisplayName',['Group ']);
    
%         [idx,cent,sumdist] = kmeans(data1,n,'Distance','cosine','Replicates',30);
figure('units','normalized','outerposition',[0 0 1 1]);
type = max(max(id));
rgb = maxdistcolor(type+1,@srgb_to_Lab);
group.type = cell(1,type);
group.id = id;
group.typeNo = type;
group.n = n;
for i=1:type
    
    for j=1:n
        x = find(id(j,:)==i);
        if isempty(x)
            continue;
        end
        subplot(type,n,(i-1)*(n) + j); hold on;
%         if i==1
%             title(['Set',num2str(j)])
%         end
        set(gca,'XTick',[],'YTick',[],'ButtonDownFcn',@createnew_fig)
%         if j==1
%             ylabel(['Type ',num2str(i)])
%         end
        group.type{i}=[group.type{i} , data{j}(:,x)];
        bar(data{j}(:,x),'FaceColor',rgb(i,:))
        for ii = 1:length(partition)
            plot(partition(ii)*[1 1]+0.5 , [0 1],'color','black') 
        end
    end
    
    
%         h=subplot(type,n+1,(i)*(n+1) ); hold on;
%         if i==1
%             title(['Mean'])
%         end
%         bar(mean(group.type{i},2),'FaceColor',rgb(i,:))
%         errorbar(mean(group.type{i},2), std(group.type{i},0,2),'LineStyle','none','Color',rgb(end,:))
%         set(gca,'XTick',[],'YTick',[],'ButtonDownFcn',@createnew_fig)
%         for ii = 1:length(partition)
%                 plot(partition(ii)*[1 1]+0.5 , [0 1],'color','black') 
%         end
%         
%     
%     if n==1
%        delete(h); 
%     end
%     xticks(1:16)
%     xtickangle(45)
%     xticklabels(label(1:16))
    
end
        
        


    for i=1:type

        for j=1:n
            x = find(id(j,:)==i);
            if isempty(x)
                continue;
            end
            subplot(type,n,(i-1)*(n) + j); hold on;
%             title(['rep:',num2str(rep(j,x)),' Corr:',num2str(corr(data{j}(:,x),mean(group.type{i},2)),'%.2f'),' NDP:',num2str(NDP(data{j}(:,x),mean(group.type{i},2)),'%.2f')])
            title([' Corr:',num2str(corr(data{j}(:,x),mean(group.type{i},2)),'%.2f'),' NDP:',num2str(NDP(data{j}(:,x),mean(group.type{i},2)),'%.2f')])
%             xticks(1:16)
%             xtickangle(45)
%             xticklabels(label(1:16))


        end    


    end

        