function out = vaf1(x, ref, dim)

m=size(x,1);
n=size(x,2);

if dim==1
    out = zeros(m,1);
    for i =1:m
    out(i,1) = 1-(norm(x(i,:)-ref(i,:),2)/norm(ref(i,:),2))^2;
    end
elseif dim==2
    out = zeros(1,n);
    for i =1:n
    out(1,n) = 1-(norm(x(:,i)-ref(:,i),2)/norm(ref(:,i),2))^2;
    end
else
    out = 1-(norm(x - ref,2)/norm(ref,2))^2;
end
end