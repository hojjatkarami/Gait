function group = kmean1(N,data,data1,label,partition,name)

n = length(data);   % number of set
 m=[];   % number of synergies in each set
 W=[];
for i=1:n
    m(i) = size(data{i},2); % number of vectors in each set
    W=[W, data{i}];
end
 
[idx,cent,sumdist] = kmeans(W',N,'Distance','correlation','Replicates',30);
% plot W
h1 = figure('units','normalized','outerposition',[0 0 1 1]);

type = N;
rgb = maxdistcolor(type+1,@srgb_to_Lab);    % create maximum dixtinct color
group.rgb = rgb;
group.Wtype = cell(1,type);
group.Stype = cell(1,type);
for i=1:n
   id(i,1:m(i))=idx(1:m(i));
   idx(1:m(i))=[];
end
group.id = id;
group.WtypeNo = type;
group.n = n;

for i=1:type
    
    for j=1:n
        x = find(id(j,:)==i);
        if isempty(x)
            continue;
        end
        group.Wtype{i}=[group.Wtype{i} , data{j}(:,x)];
        
    end
end

for k=1:n
for i=1:N

M=mean(group.Wtype{i},2);
t=find(id(k,:)==i);
if length(t)<=1
    continue;
end
corr_p=0;
h=0;
for j=1:length(t)
    ccorr=corr(data{k}(:,t(j)), M);
    if ccorr>corr_p
        corr_p=ccorr;
        h=t(j);
    end
end

id(k,t)=0;
id(k,h)=i;
end
end

LP5=designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',4,'SampleRate',120,'DesignMethod','butter');
for i=1:type
    
    for j=1:n
        x = find(id(j,:)==i);
        if isempty(x)
            continue;
        end
        group.Wtype{i}=[group.Wtype{i} , data{j}(:,x)];
        group.Stype{i}=[group.Stype{i} ; filtfilt(LP5, data1{j}.Y(x,:))];
    end
end





for i=1:type
    
    for j=1:n
        x = find(id(j,:)==i);
        if isempty(x)
            continue;
        end
        group.Wtype{i}=[group.Wtype{i} , data{j}(:,x)];
        group.Stype{i}=[group.Stype{i} ; data1{j}.Y(x,:)];
        
        subplot(type,n+1,(i-1)*(n+1) + j); hold on;
        set(gca,'XTick',[],'YTick',[],'ButtonDownFcn',@createnew_fig)        
        bar(data{j}(:,x),'FaceColor',rgb(i,:))
        for ii = 1:length(partition)    % draw partition line
            plot(partition(ii)*[1 1]+0.5 , [0 1],'color','black') 
        end
    end
    
        % draw mean of type
        h=subplot(type,n+1,(i)*(n+1) ); hold on;
        if i==1
            title(['Mean'])
        end
        bar(mean(group.Wtype{i},2),'FaceColor',rgb(i,:))
        errorbar(mean(group.Wtype{i},2), std(group.Wtype{i},0,2),'LineStyle','none','Color',rgb(end,:))
        plot(group.Wtype{i},'k*','LineStyle','none')
        set(gca,'XTick',[],'YTick',[],'ButtonDownFcn',@createnew_fig)
        for ii = 1:length(partition)
                plot(partition(ii)*[1 1]+0.5 , [0 1],'color','black') 
        end
        
    
    if n==1
       delete(h); 
    end
    xticks(1:16)
    xtickangle(45)
    xticklabels(label)
    
end 


if n~=1
    for i=1:type

        for j=1:n
            x = find(id(j,:)==i);
            if isempty(x)
                continue;
            end
            subplot(type,n+1,(i-1)*(n+1) + j); hold on;
            title(['C:',num2str(corr(data{j}(:,x),mean(group.Wtype{i},2)),'%.2f'),' N:',num2str(NDP(data{j}(:,x),mean(group.Wtype{i},2)),'%.2f')])


        end    


    end
end
suptitle(name)
    


% plot S
h2 = figure('units','normalized','outerposition',[0 0 1 1]);

for i=1:type
    
    for j=1:n
        x = find(id(j,:)==i);
        if isempty(x)
            continue;
        end
%         group.Wtype{i}=[group.Wtype{i} , data{j}(:,x)];
        
        
        subplot(type,n+1,(i-1)*(n+1) + j); hold on;
        set(gca,'XTick',[],'YTick',[],'ButtonDownFcn',@createnew_fig)        
        plot(data1{j}.X(x,:),data1{j}.Y(x,:),'Color',rgb(i,:))
         
        l=[0:100/9:100];
            for ii=1:length(l)
                line([1 1]*l(ii),[0 1],'linestyle','--');
                
            end
    end

%         draw mean of type
        h=subplot(type,n+1,(i)*(n+1) ); hold on;
        if i==1
            title(['Mean'])
        end
%         plot(data1{1}.X(x,:),mean(group.Stype{i},1),'Color',rgb(i,:))
             plot([0:.1:99.9],mean(group.Stype{i},1),'Color',rgb(i,:))
%             plot(mean(group.Stype{i},1),'Color',rgb(i,:))
%         plot(mean(group.Stype{i},1)-std(group.Stype{i},0,1),'Color',rgb(i,:))
%         plot(mean(group.Stype{i},1)+std(group.Stype{i},0,1),'Color',rgb(i,:))
        x=[data1{1}.X(x,:),fliplr(data1{1}.X(x,:))];
        x=[[0:.1:99.9],fliplr([0:.1:99.9])];

        y=[mean(group.Stype{i},1)-std(group.Stype{i},0,1),...
            fliplr(mean(group.Stype{i},1)+std(group.Stype{i},0,1))];
        s=fill(x,y,'k','EdgeColor','none');
        alpha(s,.1)
        
        set(gca,'XTick',[],'YTick',[],'ButtonDownFcn',@createnew_fig)
        l=[0:100/9:100];
            for ii=1:length(l)
                line([1 1]*l(ii),[0 1],'linestyle','--');
                
            end
        
    
%     if n==1
%        delete(h); 
%     end
%     xticks(1:16)
%     xtickangle(45)
%     xticklabels(label)
    
end 
suptitle(name)

% plot W and S

h3 = figure('units','normalized','outerposition',[0 0 1 1]);

for i=1:N
    
    % plot S
    subplot(N,2,2*i); hold on;
                plot([0:.1:99.9],mean(group.Stype{i},1),'Color',rgb(i,:))
%         plot(mean(group.Stype{i},1)-std(group.Stype{i},0,1),'Color',rgb(i,:))
%         plot(mean(group.Stype{i},1)+std(group.Stype{i},0,1),'Color',rgb(i,:))
%         x=[data1{1}.X(x,:),fliplr(data1{1}.X(x,:))];
        x=[[0:.1:99.9],fliplr([0:.1:99.9])];

        y=[mean(group.Stype{i},1)-std(group.Stype{i},0,1),...
            fliplr(mean(group.Stype{i},1)+std(group.Stype{i},0,1))];
        s=fill(x,y,'k','EdgeColor','none');
        alpha(s,.1)
        
        set(gca,'XTick',[],'YTick',[],'ButtonDownFcn',@createnew_fig)
        l=[0:100/9:100];
            for ii=1:length(l)
                line([1 1]*l(ii),[0 1],'linestyle','--');
                
            end
    
    % plot W
    subplot(N,2,2*i-1); hold on;
    bar(mean(group.Wtype{i},2),'FaceColor',rgb(i,:))
        errorbar(mean(group.Wtype{i},2), std(group.Wtype{i},0,2),'LineStyle','none','Color',rgb(end,:))
        plot(group.Wtype{i},'k*','LineStyle','none')
        set(gca,'XTick',[],'YTick',[],'ButtonDownFcn',@createnew_fig)
        
    
        
    
end

    xticks(1:28)
    xtickangle(45)
    xticklabels(label)

  