function Mat2 = Coordinate(Mat)
    Frames=size(Mat.Trajectory.LASI,1);
    
    % Pelvic Coordination______________________________________________________
    for i=1:Frames %Number Of Frames
        Mat.PelCo.Org(i,1:3)=(Mat.Trajectory.RASI(i,1:3)+Mat.Trajectory.LASI(i,1:3))/2;
        MP(i,1:3)=(Mat.Trajectory.RPSI(i,1:3)+Mat.Trajectory.LPSI(i,1:3))/2;
        y=Mat.Trajectory.LASI(i,1:3)-Mat.Trajectory.RASI(i,1:3);
        z=cross(Mat.Trajectory.LASI(i,1:3)-MP(i,1:3),y);
        x=cross(y,z);

        Mat.PelCo.x(i,1:3)=x/norm(x);        % Obliquity
        Mat.PelCo.y(i,1:3)=y/norm(y);        % Tilt
        Mat.PelCo.z(i,1:3)=z/norm(z);        % Rotation

        clear x y z
    end


    % Hip Joint Center_________________________________________________________
    mm=14/2;      %Marker Radius
    theta=0.5;
    beta=0.314;

    for i=1:Frames %Number Of Frames
        LegL(i,1:3)=norm(Mat.Trajectory.LASI(i,1:3)-Mat.Trajectory.LANK_M(i,1:3));
        LegR(i,1:3)=norm(Mat.Trajectory.RASI(i,1:3)-Mat.Trajectory.RANK_M(i,1:3));
        AsisD(i,1:3)=norm(Mat.Trajectory.LASI(i,1:3)-Mat.Trajectory.RASI(i,1:3));

        AsisTrocDistL(i,1:3)=0.1288*LegL(i,1:3)-48.56;
        AsisTrocDistR(i,1:3)=0.1288*LegR(i,1:3)-48.56;
        C(i,:)=(LegL(i,1:3)+LegR(i,1:3))/2*0.115-15.3;
        aa(i,:)=AsisD(i,1:3)/2;

        Lhjc(i,:)=[C(i,1:3)*cos(theta)*sin(beta)-(AsisTrocDistL(i,1:3)+mm)*cos(beta) ...
            -(C(i,1:3)*sin(theta)-aa(i,1:3)) ...
            -C(i,1:3)*cos(theta)*cos(beta)-(AsisTrocDistL(i,1:3)+mm)*sin(beta)];

        Rhjc(i,:)=[C(i,1:3)*cos(theta)*sin(beta)-(AsisTrocDistR(i,1:3)+mm)*cos(beta) ...
            (C(i,1:3)*sin(theta)-aa(i,1:3)) ...
            -C(i,1:3)*cos(theta)*cos(beta)-(AsisTrocDistR(i,1:3)+mm)*sin(beta)];

        Rot_Pel=[Mat.PelCo.x(i,1:3)' Mat.PelCo.y(i,1:3)' Mat.PelCo.z(i,1:3)'];
        Mat.HJC_GL.Left(i,1:3)=(Mat.PelCo.Org(i,1:3)'+Rot_Pel*Lhjc(i,1:3)')';
        Mat.HJC_GL.Right(i,1:3)=(Mat.PelCo.Org(i,1:3)'+Rot_Pel*Rhjc(i,1:3)')';
    end


    % Thigh Coordinate_________________________________________________________
    for i=1:Frames %Number Of Frames
        z=Mat.HJC_GL.Left(i,1:3)-(Mat.Trajectory.LKNE(i,1:3)+Mat.Trajectory.LKNE_M(i,1:3))/2;
        x=cross(Mat.Trajectory.LKNE(i,1:3)-Mat.Trajectory.LKNE_M(i,1:3),z);
        y=cross(z,x);

        Mat.ThiCo.Left.Org(i,1:3)=(Mat.Trajectory.LKNE(i,1:3)+Mat.Trajectory.LKNE_M(i,1:3))/2;
        Mat.ThiCo.Left.x(i,1:3)=x/norm(x);
        Mat.ThiCo.Left.y(i,1:3)=y/norm(y);
        Mat.ThiCo.Left.z(i,1:3)=z/norm(z);
        clear x y z

        z=Mat.HJC_GL.Right(i,1:3)-(Mat.Trajectory.RKNE(i,1:3)+Mat.Trajectory.RKNE_M(i,1:3))/2;
        x=cross(z,Mat.Trajectory.RKNE(i,1:3)-Mat.Trajectory.RKNE_M(i,1:3));
        y=cross(z,x);

        Mat.ThiCo.Right.Org(i,1:3)=(Mat.Trajectory.RKNE(i,1:3)+Mat.Trajectory.RKNE_M(i,1:3))/2;
        Mat.ThiCo.Right.x(i,1:3)=x/norm(x);
        Mat.ThiCo.Right.y(i,1:3)=y/norm(y);
        Mat.ThiCo.Right.z(i,1:3)=z/norm(z);
        clear x y z
    end


    % Shank Coordination_______________________________________________________
    for i=1:Frames%Number Of Frames
        % Tortioned Tibia for ankle angles
        y=Mat.Trajectory.LANK(i,1:3)-Mat.Trajectory.LANK_M(i,1:3);
        x=cross(y,(Mat.Trajectory.LKNE(i,1:3)+Mat.Trajectory.LKNE_M(i,1:3))/2-Mat.Trajectory.LANK_M(i,1:3));
        z=cross(x,y);

        Mat.TibCo.Left.Org(i,1:3)=(Mat.Trajectory.LANK(i,1:3)+Mat.Trajectory.LANK_M(i,1:3))/2;
        Mat.TibCo.Left.x(i,1:3)=x/norm(x);
        Mat.TibCo.Left.y(i,1:3)=y/norm(y);
        Mat.TibCo.Left.z(i,1:3)=z/norm(z);
        clear x y z

        y=Mat.Trajectory.RANK_M(i,1:3)-Mat.Trajectory.RANK(i,1:3);
        x=cross(y,(Mat.Trajectory.RKNE(i,1:3)+Mat.Trajectory.RKNE_M(i,1:3))/2-Mat.Trajectory.RANK_M(i,1:3));
        z=cross(x,y);

        Mat.TibCo.Right.Org(i,1:3)=(Mat.Trajectory.RANK(i,1:3)+Mat.Trajectory.RANK_M(i,1:3))/2;
        Mat.TibCo.Right.x(i,1:3)=x/norm(x);
        Mat.TibCo.Right.y(i,1:3)=y/norm(y);
        Mat.TibCo.Right.z(i,1:3)=z/norm(z);
        clear x y z


        % UnTortioned Tibia for knee angles
        z=(Mat.Trajectory.LKNE(i,1:3)+Mat.Trajectory.LKNE_M(i,1:3))/2-(Mat.Trajectory.LANK(i,1:3)+Mat.Trajectory.LANK_M(i,1:3))/2;
        x=cross(Mat.ThiCo.Left.y(i,1:3),z);
        y=cross(z,x);

        Mat.UNTibCo.Left.Org(i,1:3)=(Mat.Trajectory.LANK(i,1:3)+Mat.Trajectory.LANK_M(i,1:3))/2;
        Mat.UNTibCo.Left.x(i,1:3)=x/norm(x);
        Mat.UNTibCo.Left.y(i,1:3)=y/norm(y);
        Mat.UNTibCo.Left.z(i,1:3)=z/norm(z);
        clear x y z

        z=(Mat.Trajectory.RKNE(i,1:3)+Mat.Trajectory.RKNE_M(i,1:3))/2-(Mat.Trajectory.RANK(i,1:3)+Mat.Trajectory.RANK_M(i,1:3))/2;
        x=cross(Mat.ThiCo.Right.y(i,1:3),z);
        y=cross(z,x);

        Mat.UNTibCo.Right.Org(i,1:3)=(Mat.Trajectory.RANK(i,1:3)+Mat.Trajectory.RANK_M(i,1:3))/2;
        Mat.UNTibCo.Right.x(i,1:3)=x/norm(x);
        Mat.UNTibCo.Right.y(i,1:3)=y/norm(y);
        Mat.UNTibCo.Right.z(i,1:3)=z/norm(z);
        clear x y z
    end


    % Foot Coordination________________________________________________________
    for i=1:Frames%Number Of Frames
        % The Main Foot Segment
        z=Mat.Trajectory.LTOE(i,1:3)-(Mat.Trajectory.LANK(i,1:3)+Mat.Trajectory.LANK_M(i,1:3))/2;
        y=cross((Mat.Trajectory.LKNE(i,1:3)+Mat.Trajectory.LKNE_M(i,1:3))/2-(Mat.Trajectory.LANK(i,1:3)+Mat.Trajectory.LANK_M(i,1:3))/2,z);
        x=cross(y,z);

        Mat.FootCo.Left.Org(i,1:3)=(Mat.Trajectory.LANK(i,1:3)+Mat.Trajectory.LANK_M(i,1:3))/2;
        Mat.FootCo.Left.x(i,1:3)=x/norm(x);
        Mat.FootCo.Left.y(i,1:3)=y/norm(y);
        Mat.FootCo.Left.z(i,1:3)=z/norm(z);
        clear x y z

        z=Mat.Trajectory.RTOE(i,1:3)-(Mat.Trajectory.RANK(i,1:3)+Mat.Trajectory.RANK_M(i,1:3))/2;
        y=cross((Mat.Trajectory.RKNE(i,1:3)+Mat.Trajectory.RKNE_M(i,1:3))/2-(Mat.Trajectory.RANK(i,1:3)+Mat.Trajectory.RANK_M(i,1:3))/2,z);
        x=cross(y,z);

        Mat.FootCo.Right.Org(i,1:3)=(Mat.Trajectory.RANK(i,1:3)+Mat.Trajectory.RANK_M(i,1:3))/2;
        Mat.FootCo.Right.x(i,1:3)=x/norm(x);
        Mat.FootCo.Right.y(i,1:3)=y/norm(y);
        Mat.FootCo.Right.z(i,1:3)=z/norm(z);
        clear x y z
    end
    Mat2=Mat;
end

