function tSNEPlot(Y, TYPE, h)

label=[TYPE.(h)];
label = label(:);
id = [TYPE.id];
id = id(:);
N = max(id);
%%
distance = {'euclidean'  'seuclidean'  'cityblock'  'chebychev'...
    'minkowski'  'mahalanobis'  'cosine'  'correlation'...
    'spearman'  'hamming'  'jaccard'};
distance = {'euclidean'   'cityblock' ...
    'minkowski' 'cosine'  'correlation'...
    'spearman'};
distance = {'seuclidean'};
%% 13) plot t-SNE
allColor = 'mcrgby';
% X=[];
% N=6;
rgb = maxdistcolor(N+3,@srgb_to_Lab);    % create maximum dixtinct color
plotColor = rgb(1:N,:);


for i_dis = 1:length(distance)
    if length(distance)>1
        subplot(3,5,i_dis)
    end
    
    
    figure
    hold on
    axis equal
    % gscatter(Y(:,1),Y(:,2),TYPE)
    U = unique(label);
    
    %     for i=1:length(U)
    i=1;
    % plot omitted
    % bad corr
    %         indx1 = find(id==100+i);
    %         plot(Y(indx1,1),Y(indx1,2),'color','k','linestyle','none',...
    %             'marker','x','DisplayName','Omitted','MarkerSize',10,...
    %             'linewidth',1,'Handlevisibility','off')
    % repeatetive
    %         indx2 = find(id==200+i);
    %         plot(Y(indx2,1),Y(indx2,2),'color','k','linestyle','none',...
    %             'marker','x','DisplayName','Omitted','MarkerSize',10,...
    %             'linewidth',1,'Handlevisibility','off')
    
    indx3 = find(id==0);
    %         gscatter(Y(indx3,1), Y(indx3,2), label(indx3))
    length(indx3)/length(id)
    Y(indx3,:) = [];
    label(indx3) = [];
    % figure
    gscatter(Y(:,1),Y(:,2),label)
    
    
    % plot group
    %         indx3 = find(id~=100+i && id~=200+i );
    %         plot(Y(:,1),Y(:,2),'color',plotColor(i,:),'linestyle','none',...
    %             'marker','o','DisplayName',['G' num2str(i)],'linewidth',1,...
    %             'MarkerEdgeColor','k','MarkerFaceColor',plotColor(i,:),...
    %             'MarkerSize',6)
    %
    %         C = [mean(Y(indx,1)) mean(Y(indx,2))];
    %         r = 2* std(sqrt(Y(indx,1).^2 + Y(indx,2).^2));
    %
    %         th=0:pi/50:2*pi;
    %         xunit=r*cos(th)+C(1);
    %         yunit=r*sin(th)+C(2);
    %         plot(xunit,yunit,'--k','LineWidth',1,'Handlevisibility','off')
    
    
    
    %     end
    
    title([h,' - ',distance{i_dis}])
end
legend

% TYPE1 = Tag_sub(Tag_sub==[5 6 7]);
% X = X(:,find(Tag_sub==[5 6 7]));
% Y = tsne(X','Algorithm','exact','Distance',distance{i_dis});
%     gscatter(Y(:,1),Y(:,2),TYPE1)
%     title(distance{i_dis})

