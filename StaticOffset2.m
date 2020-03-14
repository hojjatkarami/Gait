clc
clear all
close all

LR={'Left' 'Right'};
Co={'PelCo' 'ThiCo' 'TibCo' 'UNTibCo' 'FootCo'};
xyz={'x' 'y' 'z'};

load Static%load Trajectory file
load Balance %load Trajectory file
sub=1;  % subject Number

% Segments Coordinate Defenition___________________________________________
Static=Coordinate(Static);
Balance=Coordinate(Balance);

% Rotation Matrix for Static Offset _______________________________________
for i=1:length(Co)
    for k=1:length(xyz)
        if i==1
            Static0.(Co{i}).(xyz{k})=mean(Static.(Co{i}).(xyz{k}),1);  
        else
            for h=1:length(LR)     %Left,Right
                Static0.(Co{i}).(LR{h}).(xyz{k})=mean(Static.(Co{i}).(LR{h}).(xyz{k}),1);
            end
        end
    end
end

% Euler Angle Calculation__________________________________________________
Static2=EulerAngle(Static,Static0);
Balance2=EulerAngle(Balance,Static0);

T(sub).Static=Static;
T(sub).Static2=Static2;
T(sub).Balance=Balance;
T(sub).Balance2=Balance2;
save('T','T')

