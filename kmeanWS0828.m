function [id,Wtype,dataW,Stype,dataS] = kmeanWS0828(dataW,dataS,type, N,minCorr)


n = length(dataW);   % number of set
 m=[];   % number of synergies in each set
 W=[];
for i=1:n
    m(i) = size(dataW{i},2); % number of vectors in each set
    W=[W, dataW{i}];
end
[idx,cent,sumdist] = kmeans(W',N,'Distance','correlation','Replicates',1,'Display','final');
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
        Wtype{i}=[Wtype{i} , dataW{j}(:,x)];
        
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
            ccorr=corr(dataW{k}(:,x), M);
            if ccorr>corr_p
                corr_p=ccorr;
                best=x;
            end
        end
        
        bad = t(t~=best);        
        id{k}(bad)=0;
        dataW{k}(:,bad) = 0;  
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
        Wtype{i}=[Wtype{i} , dataW{j}(:,x)];
        
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
        corrMat{j}(t) = corr(dataW{j}(:,t),M)';
        if corrMat{j}(t)<minCorr
            id{j}(t)=0;
            dataW{j}(:,t) = 0;  
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
        Wtype{i}=[Wtype{i} , dataW{j}(:,x)];
        
    end
end


% return Stype
Stype = cell(1,N);
%     for i=1:N    
%         for j=1:n
%             x = find(id{j}==i);
%             if isempty(x) || id{j}(x)==0
%                 continue;
%             end
%             Stype{i}=[Stype{i} ; dataS{j}(x,:)];
%         end
%     end
