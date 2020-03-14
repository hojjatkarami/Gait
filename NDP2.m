function out = NDP2(data1, data2)
    n1=size(data1,2);
    n2=size(data2,2);
    out=[];
    for i=1:n1
        for j=1:n2
            v1=data1(:,i);
            v2=data2(:,j);
            out(i,j) = dot(v1,v2) / (norm(v1) * norm(v2));
        end
    end
    out(out<0.75)=0;
end