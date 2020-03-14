function out = calFeature(data,windowSize,fea)

% dataG is t by 8 matrix
% dataF calculates desired feature "fea" for each frame of dataG
% based on windowSize
if length(data)==0
    return;
end
allFeatures = {'MAV','RMS','WL','ZC','SSC','LD','MAX','VAR'};
T=windowSize;
n = size(data,1);   % length of signal
m = size(data,2);   % dimension of signal
for i=1:windowSize:n
    for j = 1:m
        eval(['out((i-1)/windowSize+1,j) = j' fea '(data(i:min(n,i+windowSize-1),j));'])
    end
end




end

function MAV=jMAV(X)
MAV=mean(abs(X));
end
function RMS=jRMS(X)
RMS=sqrt(mean(X.^2));
end
function WL=jWL(X)
N=length(X); WL=0;
for i=2:N
    WL=WL+abs(X(i)-X(i-1));
end
end
function ZC=jZC(X)
thres = 0;
N=length(X); ZC=0;
for i=1:N-1
    if ((X(i) > 0 && X(i+1) < 0) || (X(i) < 0 && X(i+1) > 0)) ...
            && (abs(X(i)-X(i+1)) >= thres)
        ZC=ZC+1;
    end
end
end
function SSC=jSSC(X)
thres=0;
N=length(X); SSC=0;
for i=2:N-1
    if ((X(i) > X(i-1) && X(i) > X(i+1)) || (X(i) < X(i-1) && X(i) < X(i+1))) ...
            && ((abs(X(i)-X(i+1)) >= thres) || (abs(X(i)-X(i-1)) >= thres))
        SSC=SSC+1;
    end
end
end
function LD=jLD(X)
N=length(X); Y=0;
for k=1:N
    Y=Y+log(abs(X(k)));
end
LD=exp(Y/N);
end
function MFL=jMFL(X)

MFL = log10(sqrt(sum(diff(X).^2)));

end

function MAX=jMAX(X)
N=length(X);

MAX = max(X);

end

function VAR=jVAR(X)
N=length(X);

VAR = var(X);

end


function MSR=jMSR(X)
N=length(X);

MSR = sum(sqrt(X)) / N;

end


