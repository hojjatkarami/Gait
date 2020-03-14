function Mat2 = myspline(Mat,PurtFrame)
PurtFrame = double(PurtFrame);
n=size(Mat,1);
Norm=[1 12 50 60 101];
m=length(PurtFrame)-1;
for i=1:n
    for j=1:m
        x = PurtFrame(j):PurtFrame(j+1);
        y = Mat(i,x);
        xq = Norm(j):(Norm(j+1)-1);
        Mat2(i,xq) =  spline(x,y,xq);
    
    end
end




