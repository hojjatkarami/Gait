function T = PhasicNorm(T, subject)

LR={'Left' 'Right'};
% phasic time normalization
for i_sub=subject
    disp(['subject:',num2str(i_sub)])
    %     %% Remove Static Offset and Calculate Kinematic angles
    %
    
    % Calculate angles of plate and find purt frames
    % add all and purtFrame field to Event
    for i_trial= 1:length(T(i_sub).Trajectory.Left)
        T(i_sub).Events.Left(i_trial).all=double([T(i_sub).Events.Left(i_trial).LFC1 T(i_sub).Events.Left(i_trial).RTO...
            T(i_sub).Events.Left(i_trial).RFC T(i_sub).Events.Left(i_trial).LTO T(i_sub).Events.Left(i_trial).LFC2]);
        T(i_sub).Events.Left(i_trial).PurtFrame = T(i_sub).Events.Left(i_trial).all - T(i_sub).Events.Left(i_trial).all(1)+1;
    end
    
    for i_trial= 1:length(T(i_sub).Trajectory.Right)
        T(i_sub).Events.Right(i_trial).all=double([T(i_sub).Events.Right(i_trial).RFC1 T(i_sub).Events.Right(i_trial).LTO...
            T(i_sub).Events.Right(i_trial).LFC T(i_sub).Events.Right(i_trial).RTO T(i_sub).Events.Right(i_trial).RFC2]);
        T(i_sub).Events.Right(i_trial).PurtFrame = T(i_sub).Events.Right(i_trial).all - T(i_sub).Events.Right(i_trial).all(1)+1;
        
    end
    
    % phasic time normalization for EMG and kinematics
    for side=[1 2]
        for i_trial= 1:length(T(i_sub).Trajectory.(LR{side}))
            
            PurtFrame = T(i_sub).Events.(LR{side})(i_trial).PurtFrame;
            PurtFrameNorm = round(PurtFrame/PurtFrame(end)*100);
            if PurtFrameNorm(1)==0
                PurtFrameNorm(1)=1;
            end
            T(i_sub).Events.(LR{side})(i_trial).PurtFrameNorm = PurtFrameNorm;
            
            Mat1 = T(i_sub).EMG2.(LR{side})(i_trial).M;
            [Mat2,t] = rescale0720(Mat1, PurtFrame, PurtFrameNorm);
            T(i_sub).EMG2.(LR{side})(i_trial).Mn = Mat2;
            T(i_sub).Events.(LR{side})(i_trial).t = t;
            T(i_sub).Events.(LR{side})(i_trial).tNorm = (t/t(end)*100);
            Mat1 = T(i_sub).KIN.(LR{side})(i_trial).M;
            [Mat2,t] = rescale0720(Mat1, PurtFrame, PurtFrameNorm);
            T(i_sub).KIN.(LR{side})(i_trial).Mn = Mat2;
        end
    end
    
end