function [Mat2,t] = rescale0720(Mat, PurtFrame, PurtFrameNorm)

% each row of matrix is a vector
t = zeros(1,100);

for i = 1:(length(PurtFrameNorm)-1)
    a1 = PurtFrameNorm(i);  a2 = PurtFrameNorm(i+1);    b1 = PurtFrame(i);  b2 = PurtFrame(i+1);      
    c=a1;
    while c<=a2
        t(c) = round( (c-a1)/(a2-a1)*(b2-b1) ) + b1;
        c = c+1;        
    end   
end

Mat2(:,1:100) = Mat(:,t);




