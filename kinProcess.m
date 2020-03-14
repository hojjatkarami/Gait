function out = kinProcess(M)

 len = size(M,1);
for i=1:len
temp = M(i,:);       
M1(i,:) = temp;
temp_pos = temp;
temp_neg = -temp;
temp_pos(find(temp_pos<0))=0;
temp_neg(find(temp_neg<0))=0;
%                 M_pos(i,:) = temp_pos;
%                 M_neg(i,:) = temp_neg;
M2(2*i-1,:) = temp_pos;
M2(2*i,:) = temp_neg;
end
out = M2;
%             M2 = [M_pos;M_neg];
%              m=max(M2,[],2);
%              m(find(m==0)) = 1;
%              M2 = M2./m;
%              std_val=std(M2,0,2);
%               std_valM(find(std_valM==0)) =1;
%              M2=M2./std_valM;
%              MM=[M;M2];