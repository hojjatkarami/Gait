function out = InvkinProcess(M)
out=[];
 len = size(M,1);
for i=1:len/2
    out(i,:) = M(2*i-1,:)-M(2*i,:);
end

