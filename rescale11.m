function [t2, y2] = rescale11(x1,y1,s1,s2,len)
x2=[];
n=length(s1)-1;
s1=double(s1);
s2=double(s2);
for i=1:n-1
   r(i)=floor((s1(i+1)-s1(i))/(s1(end)-s1(1))*len); 
end
if n==1
   r=[1000];
else
  r(i+1) = len-sum(r);   
end
t2=[];
t1=[];
y2=[];
temp=[];
for i=1:n
    % a*s1+b=s2
    
    
    step=(s2(i+1)-s2(i))/(r(i));
    
    for j=0:r(i)-1
       temp=[temp, s2(i)+j*step];
       
        
    end
    s1=double(s1);
    s2=double(s2);
    coeff = inv([s1(i),1;s1(i+1),1])*[s2(i);s2(i+1)];
    t1 = [t1,  floor((temp-coeff(2))/coeff(1))+1  ];
%     x2 = [x2 , coeff(1)*[s1(i):s1(i+1)]+coeff(2)];
    
%     if i~=n
%         x2(end)=[];
%     end
t2=[t2, temp];
temp=[];
end
y2 = y1(t1-s1(1)+1);


