function T = ISB(T, subject)

LR={'Left' 'Right'};
Co={'PelCo' 'ThiCo' 'TibCo' 'UNTibCo' 'FootCo'};
xyz={'x' 'y' 'z'};
% Kinematic angle extraction
for i_sub=subject
    disp(['subject:',num2str(i_sub)])
    for side=[2]
        Static.Trajectory=T(i_sub).Static.(LR{side});
%         out = checkISBDimensions(Static.Trajectory);
%         if out==1
%             disp(['sub:' num2str(i_sub) ' static' ' side:' LR{side}])
%             continue;
%         end
        Static=Coordinate(Static);  % Segments Coordinate Defenition
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
        Static2=EulerAngle(Static,Static0);     % Euler Angle Calculation
%         T(i_sub).Static.Static0.(LR{side})=Static0;
%         T(i_sub).Static.Static2.(LR{side})=Static2;
        
        for i_trial=1:length(T(i_sub).Trajectory.(LR{side}))
            kin.Trajectory = T(i_sub).Trajectory.(LR{side})(i_trial);
            out = checkISBDimensions(kin.Trajectory);
%             if out==1
%                 disp(['sub:' num2str(i_sub) ' kin' ' trial:' num2str(i_trial) ' side:' LR{side}])
%                 continue;
%             end
            kin=Coordinate(kin);  % Segments Coordinate Defenition
            kin2=EulerAngle(kin,Static0);   % Euler Angle Calculation
            T(i_sub).KIN.(LR{side})(i_trial).value = kin2;
            
            clear kin kin2
            
        end
        
        clear Static Static0 Static2
        
    end
end

disp('ISB done')
end

function out = checkISBDimensions(traj)
% this function check if dimension of all markers in a trajectory set is
% consistent or not
out=0;
allFields= fieldnames(traj);
t = size(traj.(allFields{1}),1);
for i=2:length(allFields)
    if size(traj.(allFields{i}),1) ~= t
        warning('wrong marker dimension in traj -> %s',allFields{i})
        out=1;
    end
end

end

