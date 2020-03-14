function [id,Wtype,Wmean, Wstd] = kmeanW1122(data, N,minCorr)


n = length(data);   % number of set
 m=[];   % number of synergies in each set
 W=[];
for i=1:n
    m(i) = size(data{i},2); % number of vectors in each set
    W=[W, data{i}];
end
[idx,cent,sumdist] = kmeans(W',N,'Distance','correlation','Replicates',100,'Display','off');
for i=1:n
   id{i}=idx(1:m(i));
   idx(1:m(i))=[];
end


Wtype = cell(1,N);
for i=1:N    
    for j=1:n
        x = find(id{j}==i);
        if isempty(x)
            continue;
        end
        Wtype{i}=[Wtype{i} , data{j}(:,x)];
        
    end
end



%% remove repeatitive items
for i=1:N
    M=mean(Wtype{i},2);
    for k=1:n
        
        t=find(id{k}==i);
        if length(t)<=1
            continue;
        end
        corr_p=0;
        best=0;
        for x=t'
            ccorr=corr(data{k}(:,x), M);
            if ccorr>corr_p
                corr_p=ccorr;
                best=x;
            end
        end
        
        bad = t(t~=best);        
        id{k}(bad)=0;
        data{k}(:,bad) = 0;  
%         Wtype{i}(:,bad)=[];

    end
end

Wtype = cell(1,N);
for i=1:N    
    for j=1:n
        x = find(id{j}==i);
        if isempty(x) || id{j}(x)==0
            continue;
        end
        Wtype{i}=[Wtype{i} , data{j}(:,x)];
        
    end
end

%% remove bad corr under #minCorr
corrMat = {};
for i=1:N
    M=mean(Wtype{i},2);
    for j=1:n
        t=find(id{j}==i);
        if isempty(t)
            continue;
        end
        corrMat{j}(t) = corr(data{j}(:,t),M)';
        if corrMat{j}(t)<minCorr
            id{j}(t)=0;
            data{j}(:,t) = 0;  
%             Wtype{i}(:,t)=[];
            corrMat{j}(t)=[];
        end
        
    end
end


Wtype = cell(1,N);
for i=1:N    
    for j=1:n
        x = find(id{j}==i);
        if isempty(x) || id{j}(x)==0
            continue;
        end
        Wtype{i}=[Wtype{i} , data{j}(:,x)];
        
    end
end
%% calculate mean of W
for i=1:N
    Wmean(:,i) = mean(Wtype{i},2);
    Wstd(:,i) = std(Wtype{i},0,2);
    
    
end
