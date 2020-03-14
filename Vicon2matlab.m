clc
% clear all
close all

% load P
sub=4; %SubjectNumber::::::::::::::::::::::::::::::::::::::::::::::::::

vicon=ViconNexus();
P(sub).Name=vicon.GetSubjectNames;
P(1).Name='Hojat Karami';
P(2).Name='Sina Esmaeili';
P(3).Name='M.H Mashaghi';
P(4).Name='M.R Fereydoni';
P(5).Name='Mohsen Bagheri';
P(6).Name='Ehsan Abdoulahi';
P(7).Name='Parsa Riazi';
P(8).Name='Mohamad Azimi';
Name(sub)=vicon.GetSubjectNames;


LAng={'LAbsAnkleAngle' 'LAnkleAngles' 'LFootProgressAngles' 'LHipAngles' 'LKneeAngles' 'LPelvisAngles' 'LSpineAngles' 'LThoraxAngles'};
RAng={'RAbsAnkleAngle' 'RAnkleAngles' 'RFootProgressAngles' 'RHipAngles' 'RKneeAngles' 'RPelvisAngles' 'RSpineAngles' 'RThoraxAngles'};
LKin={'LGroundReactionForce' 'LGroundReactionMoment' 'LNormalisedGRF' 'LAnkleForce' 'LAnkleMoment' 'LAnklePower' 'LHipForce' 'LHipMoment' 'LHipPower' 'LKneeForce' 'LKneeMoment' 'LKneePower' 'LWaistForce' 'LWaistMoment' 'LWaistPower'};
RKin={'RGroundReactionForce' 'RGroundReactionMoment' 'RNormalisedGRF' 'RAnkleForce' 'RAnkleMoment' 'RAnklePower' 'RHipForce' 'RHipMoment' 'RHipPower' 'RKneeForce' 'RKneeMoment' 'RKneePower' 'RWaistForce' 'RWaistMoment' 'RWaistPower'};
Lmrk={'LASI' 'RASI' 'LPSI' 'RPSI' 'LTHI' 'LKNE' 'LTIB' 'LANK' 'LHEE' 'LTOE' 'LTRO' 'LANK-M' 'LKNE-M' 'LSHO' 'CLAV' 'STRN' 'T10' 'C7' 'LMET5' 'LMET1' 'RTHI' 'RKNE' 'RTIB' 'RANK' 'RHEE' 'RTOE' 'RTRO' 'RANK-M' 'RKNE-M'  'RSHO'  'RMET5' 'RMET1' 'LTHI-B' 'LTHI-F' 'LTIB-B' 'LTIB-F' 'LPEL' 'RTHI-B' 'RTHI-F' 'RTIB-B' 'RTIB-F'  'RPEL'  }; %
Rmrk={'LASI' 'RASI' 'LPSI' 'RPSI' 'RTHI' 'RKNE' 'RTIB' 'RANK' 'RHEE' 'RTOE' 'RTRO' 'RANK-M' 'RKNE-M'  'RSHO' 'CLAV' 'STRN' 'T10' 'C7' 'RMET5' 'RMET1' 'LTHI' 'LKNE' 'LTIB' 'LANK' 'LHEE' 'LTOE' 'LTRO' 'LANK-M' 'LKNE-M' 'LSHO'  'LMET5' 'LMET1' 'RTHI-B' 'RTHI-F' 'RTIB-B' 'RTIB-F'  'RPEL'  'LTHI-B' 'LTHI-F' 'LTIB-B' 'LTIB-F'  'LPEL'}; % 
Lmrk2={'LASI' 'RASI' 'LPSI' 'RPSI' 'LTHI' 'LKNE' 'LTIB' 'LANK' 'LHEE' 'LTOE' 'LTRO'  'LANK_M' 'LKNE_M' 'LSHO' 'CLAV' 'STRN' 'T10' 'C7' 'LMET5' 'LMET1' 'RTHI' 'RKNE' 'RTIB' 'RANK' 'RHEE' 'RTOE' 'RTRO' 'RANK_M' 'RKNE_M' 'RSHO' 'RMET5' 'RMET1' 'LTHI_B' 'LTHI_F' 'LTIB_B' 'LTIB_F'  'LPEL' 'RTHI_B' 'RTHI_F' 'RTIB_B' 'RTIB_F'  'RPEL' };%  
Rmrk2={'LASI' 'RASI' 'LPSI' 'RPSI' 'RTHI' 'RKNE' 'RTIB' 'RANK' 'RHEE' 'RTOE' 'RTRO'  'RANK_M' 'RKNE_M' 'RSHO' 'CLAV' 'STRN' 'T10'  'C7' 'RMET5' 'RMET1' 'LTHI' 'LKNE' 'LTIB' 'LANK' 'LHEE' 'LTOE' 'LTRO'  'LANK_M' 'LKNE_M' 'LSHO' 'LMET5' 'LMET1' 'RTHI_B' 'RTHI_F' 'RTIB_B' 'RTIB_F' 'RPEL' 'LTHI_B' 'LTHI_F' 'LTIB_B' 'LTIB_F' 'LPEL'};% 
Ls={'LAJC' 'LCL' 'LFE' 'LFO' 'LHJC' 'LKJC' 'LSJC' 'LTI' 'LTO' 'PEL' 'TRX'};
Rs={'RAJC' 'RCL' 'RFE' 'RFO' 'RHJC' 'RKJC' 'RSJC' 'RTI' 'RTO' 'PEL' 'TRX'};
LEMG={'LTA' 'LPL' 'LSOL' 'LGC' 'LRF' 'LHM'...
       'LVL' 'LIP' 'LGMAX' 'LGMED' 'LAD' 'LTFL' 'LIC' 'LLG' 'LRA' 'LEO'};
LEMG2={'LTA' 'LPL' 'LSOL' 'LGC' 'LRF' 'LHM'...
        'LVL' 'LIP' 'LGMAX' 'LGMED' 'LAD' 'LTFL' 'LIC' 'LLG' 'LRA' 'LEO'};
% LEMG2={ 'LIPS' 'LAD' 'LRF' 'LTA' 'LMH' 'LLH' 'LGC' 'LGMAX'};%'LSOL' 'LGLUT_M'};
REMG={'RTA' 'RPL' 'RSOL' 'RGC' 'RRF' 'RHM' 'RVL' 'RIP' 'RGMAX' 'RGMED' 'RAD' 'RTFL' 'RIC' 'RLG' 'RRA' 'REO'} ;
REMG2={'RTA' 'RPL' 'RSOL' 'RGC' 'RRF' 'RHM' 'RVL' 'RIP' 'RGMAX' 'RGMED' 'RAD' 'RTFL' 'RIC' 'RLG' 'RRA' 'REO'} ;
% REMG={'RIPS' 'RAD' 'RRF' 'RTA' 'RST' 'RBF' 'RGC' 'RGMAX'};%'RSOL' 'RGLUT-M'};
% REMG2={'RIPS' 'RAD' 'RRF' 'RTA' 'RMH' 'RLH' 'RGC' 'RGMAX'};%'RSOL' 'RGLUT_M'};
% REMG={ 'RTA' 'RPL' 'RSOL'};

l=0; r=0;
flag=0;

P(sub).ForcePlate.Left=[];
P(sub).ForcePlate.Right=[];
P(sub).Kinetics.Left=[];
P(sub).Kinetics.Right=[];
P(sub).Kinematics=[];
P(sub).JC=[];
P(sub).Trajectory=[];
P(sub).Events=[];
P(sub).EMG=[];
P(sub).Static=[];
%% Static
clear fp fp1

GAIT='static'
if flag~=0
    if l~=0
        l=size(P(sub).Trajectory.Left,2);
    end
    if r~=0
        r=size(P(sub).Trajectory.Right,2);
    end
end

%Devices--------------------------
DeviceNames=vicon.GetDeviceNames;
DeviceID=cell(length(DeviceNames),1);
DeviceType=cell(length(DeviceNames),1);
DeviceRate=cell(length(DeviceNames),1);
DeviceOutputID=cell(length(DeviceNames),1);
ChannelID=cell(length(DeviceNames),1);
for i = 1:length(DeviceNames)
    DeviceID{i}=vicon.GetDeviceIDFromName(DeviceNames{i});
    [~,DeviceType{i}, DeviceRate{i}, DeviceOutputID{i}]=vicon.GetDeviceDetails(DeviceID{i});
end
u=1;
for j=1:length(DeviceType)
    if strcmp(DeviceType{j},'ForcePlate')
        fp(u)=DeviceID(j);
        u=u+1;
    end
end
fp=double(cell2mat(fp));
fp1=min(fp);

% markers --------------------------------
NxsMrk=vicon.GetMarkerNames(Name{sub});
Iml=0; Imr=0;
for i=1:length(Lmrk)
    lm=vicon.HasTrajectory((Name{sub}),(Lmrk{i}));
    if lm==0
        Iml=[Iml i];
    end
    rm=vicon.HasTrajectory((Name{sub}),(Rmrk{i}));
    if rm==0
        Imr=[Imr i];
    end
    clear lm lr
end
maxframe=vicon.GetFrameCount;

    for j=1:length(Lmrk2)
            ch=find (Iml==j);
            if size(ch,2)==0
                [temp1 temp2 temp3 temp4]=vicon.GetTrajectory((Name{sub}),(Lmrk{j}));
                temp=[temp1;temp2;temp3; temp4];
                temp=temp';
                P(sub).Static.Left.(Lmrk2{j})=temp(488:1000,:);
%                 50:maxframe/2,:
% 200
                clear temp temp1 temp2 temp3 temp4
            else
                P(sub).Static.Left.(Lmrk2{j})=[];
            end
            clear ch
    end
    
      for j=1:length(Rmrk2)
            ch=find (Imr==j);
            if size(ch,2)==0
                [temp1 temp2 temp3 temp4]=vicon.GetTrajectory((Name{sub}),(Rmrk{j}));
                temp=[temp1;temp2;temp3; temp4];
                temp=temp';
                P(sub).Static.Right(1).(Rmrk2{j})=temp(488:1000,:);
                clear temp temp1 temp2 temp3 temp4
            else
                P(sub).Static.Right(1).(Rmrk2{j})=[];
            end
            clear ch
       end



clear maxframe

%% Kinematics & L & Trajectory & Events
   clear fp fp1
GAIT='GAIT19';  %Trial Name :::::::::::::::::::::::::::::::::::::::::::::::

if flag~=0
    if l~=0
        l=size(P(sub).Trajectory.Left,2);
    end
    if r~=0
        r=size(P(sub).Trajectory.Right,2);
    end
end

%Devices--------------------------
DeviceNames=vicon.GetDeviceNames;
DeviceID=cell(length(DeviceNames),1);
DeviceType=cell(length(DeviceNames),1);
DeviceRate=cell(length(DeviceNames),1);
DeviceOutputID=cell(length(DeviceNames),1);
ChannelID=cell(length(DeviceNames),1);
for i = 1:length(DeviceNames)
    DeviceID{i}=vicon.GetDeviceIDFromName(DeviceNames{i});
    [~,DeviceType{i}, DeviceRate{i}, DeviceOutputID{i}]=vicon.GetDeviceDetails(DeviceID{i});
end
u=1;
for j=1:length(DeviceType)
    if strcmp(DeviceType{j},'ForcePlate')
        fp(u)=DeviceID(j);
        u=u+1;
    end
end
fp=double(cell2mat(fp));
fp1=min(fp);

% markers --------------------------------
NxsMrk=vicon.GetMarkerNames(Name{sub});
Iml=0; Imr=0;
for i=1:length(Lmrk)
    lm=vicon.HasTrajectory((Name{sub}),(Lmrk{i}));
    if lm==0
        Iml=[Iml i];
    end
    rm=vicon.HasTrajectory((Name{sub}),(Rmrk{i}));
    if rm==0
        Imr=[Imr i];
    end
    clear lm lr
end

NxsOutput=vicon.GetModelOutputNames(Name{sub});
% kinematics angles------------------------------
Ial=0; Iar=0;
for i=1:length(LAng)
    lm=strcmp((LAng{i}),NxsOutput);
    if max(lm)==0
        Ial=[Ial i];
    end
    rm=strcmp((RAng{i}),NxsOutput);
    if max(rm)==0
        Iar=[Iar i];
    end
    clear lm lr
end

% JointCenter------------------------------
Ijl=0; Ijr=0;
for i=1:length(Ls)
    lm=strcmp((Ls{i}),NxsOutput);
    if max(lm)==0
        Ijl=[Ijl i];
    end
    rm=strcmp((Rs{i}),NxsOutput);
    if max(rm)==0
        Ijr=[Ijr i];
    end
    clear lm lr
end

clear LEvent REvent S h LFootOff RFootOff LGeneral RGeneraL
LEvent=vicon.GetEvents((Name{sub}),'Left','Foot Strike');
REvent=vicon.GetEvents((Name{sub}),'Right','Foot Strike');
LFootOff=vicon.GetEvents((Name{sub}),'Left','Foot Off');
RFootOff=vicon.GetEvents((Name{sub}),'Right','Foot Off');
LGeneral=vicon.GetEvents((Name{sub}),'Left','General');
RGeneral=vicon.GetEvents((Name{sub}),'Right','General');
if size(LEvent,2)>1
    if LFootOff(1)<LEvent(1)
        LFootOff=LFootOff(2:end);
    end
    for i=1:length(LEvent)
        lfs(i)=double(LEvent(i));
    end
    for i=1:length(LFootOff)
        lfo(i)=double(LFootOff(i));
    end
    clear LEvent; LEvent=lfs; clear lfs;
    clear LFootOff; LFootOff=lfo; clear lfo;
end
if size(REvent,2)>1
    if RFootOff(1)<REvent(1)
        RFootOff=RFootOff(2:end);
    end
    for i=1:length(REvent)
        rfs(i)=double(REvent(i));
    end
    for i=1:length(RFootOff)
        rfo(i)=double(RFootOff(i));
    end
    clear REvent; REvent=rfs; clear rfs;
    clear RFootOff; RFootOff=rfo; clear rfo;
end
S=vicon.GetTrialRegionOfInterest;
if size(LEvent,2)>1
    for i=l+1:l+length(LEvent)-1
        if i>10
            break
        end
        for j=1:length(LAng)
            ch=find(Ial==j);
            if size(ch,2)==0
                temp=vicon.GetModelOutput((Name{sub}),(LAng{j}));
                temp=temp';
                P(sub).Kinematics.Left(i).(LAng{j})=temp(LEvent(i-l):LEvent(i-l+1),:);
                clear temp
            else
                P(sub).Kinematics.Left(i).(LAng{j})=[];
            end
            clear ch
        end
        for j=1:length(Ls)
            ch=find(Ijl==j);
            if size(ch,2)==0
                temp=vicon.GetModelOutput((Name{sub}),(Ls{j}));
                temp=temp';
                P(sub).JC.Left(i).(Ls{j})=temp(LEvent(i-l):LEvent(i-l+1),:);
                clear temp
            else
                P(sub).JC.Left(i).(Ls{j})=[];
            end
            clear ch
        end
        for j=1:length(Lmrk2)
            ch=find (Iml==j);
            if size(ch,2)==0
                [temp1 temp2 temp3 temp4]=vicon.GetTrajectory((Name{sub}),(Lmrk{j}));
                temp=[temp1;temp2;temp3; temp4];
                temp=temp';
                P(sub).Trajectory.Left(i).(Lmrk2{j})=temp(LEvent(i-l):LEvent(i-l+1),:);
                clear temp temp1 temp2 temp3 temp4
            else
                P(sub).Trajectory.Left(i).(Lmrk2{j})=[];
            end
            clear ch
        end
        P(sub).Events.Left(i).StancePhase=(LFootOff(i-l)-LEvent(i-l))/(LEvent(i-l+1)-LEvent(i-l));
        P(sub).Events.Left(i).LFC1=LEvent(i-l);
        P(sub).Events.Left(i).RTO=LGeneral(2*(i-l)-1);
        P(sub).Events.Left(i).RFC=LGeneral(2*(i-l));
        P(sub).Events.Left(i).LTO=LFootOff(i-l);
        P(sub).Events.Left(i).LFC2=LEvent(i-l+1);
        P(sub).Events.Left(i).FileName=GAIT;   
    end
end

if size(REvent,2)>1
    for i=r+1:r+length(REvent)-1
        if i>10
            break
        end
        for j=1:length(RAng)
            ch=find(Iar==j);
            if size(ch,2)==0
                temp=vicon.GetModelOutput((Name{sub}),(RAng{j}));
                temp=temp';
                P(sub).Kinematics.Right(i).(RAng{j})=temp(REvent(i-r):REvent(i-r+1),:);
                clear temp
            else
                P(sub).Kinematics.Right(i).(RAng{j})=[];
            end
            clear ch
        end
        for j=1:length(Rs)
            ch=find(Ijr==j);
            if size(ch,2)==0
                temp=vicon.GetModelOutput((Name{sub}),(Rs{j}));
                temp=temp';
                P(sub).JC.Right(i).(Rs{j})=temp(REvent(i-r):REvent(i-r+1),:);
                clear temp
            else
                P(sub).JC.Right(i).(Rs{j})=[];
            end
            clear ch
        end
        for j=1:length(Rmrk2)
            ch=find (Imr==j);
            if size(ch,2)==0
                [temp1 temp2 temp3 temp4]=vicon.GetTrajectory((Name{sub}),(Rmrk{j}));
                temp=[temp1;temp2;temp3; temp4];
                temp=temp';
                P(sub).Trajectory.Right(i).(Rmrk2{j})=temp(REvent(i-r):REvent(i-r+1),:);
                clear temp temp1 temp2 temp3 temp4
            else
                P(sub).Trajectory.Right(i).(Rmrk2{j})=[];
            end
            clear ch
        end
        P(sub).Events.Right(i).StancePhase=(RFootOff(i-r)-REvent(i-r))/(REvent(i-r+1)-REvent(i-r));
        P(sub).Events.Right(i).RFC1=REvent(i-r);
        P(sub).Events.Right(i).LTO=RGeneral(2*(i-r)-1);
        P(sub).Events.Right(i).LFC=RGeneral(2*(i-r));
        P(sub).Events.Right(i).RTO=RFootOff(i-r);
        P(sub).Events.Right(i).RFC2=REvent(i-r+1);
        P(sub).Events.Right(i).FileName=GAIT;  
    end
end




%EMG_______________________________________________________________________
if size(LEvent,2)>1
    
    Channels=LEMG;
    Channels2=LEMG2;
    for i = 1:length(DeviceNames)
        if strcmp(DeviceType{i},'Other')
            for k=l+1:l+length(LEvent)-1
                if k<11
                    for j = 1:length(Channels)
                        ChannelID{j} = vicon.GetDeviceChannelIDFromName(DeviceID{i}, ceil(double(j)/j), Channels{j});
                        temp=vicon.GetDeviceChannel(DeviceID{i}, ceil(double(j)/j), ChannelID{j});
                        temp=temp';
                        P(sub).EMG.Left(k).(Channels2{j})=temp(LEvent(k-l)*10-10:LEvent(k-l+1)*10);
                        clear temp
                    end
                end
            end
        end
    end
    clear Channels ChannelID temp
end

if size(REvent,2)>1
    Channels=REMG;
    Channels2=REMG2;
    for i = 1:length(DeviceNames)
        if strcmp(DeviceType{i},'Other')
            for k=r+1:r+length(REvent)-1
                if k<11
                    for j = 1:length(Channels)
                        ChannelID{j} = vicon.GetDeviceChannelIDFromName(DeviceID{i}, ceil(double(j)/j), Channels{j});
                        temp=vicon.GetDeviceChannel(DeviceID{i}, ceil(double(j)/j), ChannelID{j});
                        temp=temp';
                        P(sub).EMG.Right(k).(Channels2{j})=temp(REvent(k-r)*10-10:REvent(k-r+1)*10);
                        clear temp
                    end
                end
            end
        end
    end
    clear Channels ChannelID temp
end



% Detect Correct ForcePlate Cycle__________________________________________
LFoot={'LHEE'  'LTOE' 'LMET5' 'LMET1'};
RFoot={'RHEE'  'RTOE' 'RMET5' 'RMET1'};
cycleL=[0 0];
cycleR=[0 0];
if size(LEvent,2)>1
    for i=l+1:l+length(LEvent)-1
        %check fp1--------------------------------
        flag2=0;
        for j=1:length(LFoot)
            if isempty(P(sub).Trajectory.Left(i).(LFoot{j}))==0
                tempL=P(sub).Trajectory.Left(i).(LFoot{j})(1:LFootOff(i-l)-LEvent(i-l),:);
                if isempty(tempL)==0
                    if min(tempL(:,1))<-10 || max(tempL(:,1))>510 || min(tempL(:,2))<-10 || max(tempL(:,2))>310
                        flag2=1;
                        break
                    end
                end
                clear tempL
            end
        end
        if flag2==0
            for j=1:length(RFoot)
                [temp1 temp2 temp3 temp4]=vicon.GetTrajectory((Name{sub}),(RFoot{j}));
                temp=[temp1;temp2;temp3; temp4];
                temp=temp';
                tempR1=temp(LGeneral(2*(i-l)):LFootOff(i-l),:);
                tempR2=temp(LEvent(i-l):LGeneral(2*(i-l)-1),:);
                if (min(tempR1(:,1))>10 && max(tempR1(:,1))<490 && min(tempR1(:,2))>10 && max(tempR1 (:,2))<290) || ...
                        (min(tempR2(:,1))>10 && max(tempR2(:,1))<490 && min(tempR2(:,2))>10 && max(tempR2 (:,2))<290)
                    flag2=1;
                    break
                end
                clear tempR1 tempR2 temp temp1 temp2 temp3 temp4
            end
        end
        if flag2==0
            cycleL=[cycleL; i 1];
        end
        
        %check fp2--------------------------------
        flag2=0;
        for j=1:length(LFoot)
            if isempty(P(sub).Trajectory.Left(i).(LFoot{j}))==0
                tempL=P(sub).Trajectory.Left(i).(LFoot{j})(1:LFootOff(i-l)-LEvent(i-l),:);
                if isempty(tempL)==0
                    if min(tempL(:,1))<490 || max(tempL(:,1))>1040 || min(tempL(:,2))<-10 || max(tempL(:,2))>610
                        flag2=1;
                        break
                    end
                end
                clear tempL
            end
        end
        if flag2==0
            for j=1:length(RFoot)
                [temp1 temp2 temp3 temp4]=vicon.GetTrajectory((Name{sub}),(RFoot{j}));
                temp=[temp1;temp2;temp3; temp4];
                temp=temp';
                tempR1=temp(LGeneral(2*(i-l)):LFootOff(i-l),:);
                tempR2=temp(LEvent(i-l):LGeneral(2*(i-l)-1),:);
                if (min(tempR1(:,1))>510 && max(tempR1(:,1))<990 && min(tempR1(:,2))>10 && max(tempR1(:,2))<590) || ...
                        (min(tempR2(:,1))>510 && max(tempR2(:,1))<990 && min(tempR2(:,2))>10 && max(tempR2(:,2))<590)
                    flag2=1;
                    break
                end
                clear tempR1 tempR2 temp temp1 temp2 temp3 temp4
            end
        end
        if flag2==0
            cycleL=[cycleL; i 2];
        end
    end
end

if size(REvent,2)>1
    for i=r+1:r+length(REvent)-1
        %check fp1--------------------------------
        flag2=0;
        for j=1:length(RFoot)
            if isempty(P(sub).Trajectory.Right(i).(RFoot{j}))==0
                tempR=P(sub).Trajectory.Right(i).(RFoot{j})(1:RFootOff(i-r)-REvent(i-r),:);
                if isempty(tempR)==0
                    if min(tempR(:,1))<-10 || max(tempR(:,1))>510 || min(tempR(:,2))<-10 || max(tempR(:,2))>310
                        flag2=1;
                        break
                    end
                end
                clear tempR
            end
        end
        if flag2==0
            for j=1:length(RFoot)
                [temp1 temp2 temp3 temp4]=vicon.GetTrajectory((Name{sub}),(LFoot{j}));
                temp=[temp1;temp2;temp3; temp4];
                temp=temp';
                tempL1=temp(RGeneral(2*(i-r)):RFootOff(i-r),:);
                tempL2=temp(REvent(i-r):RGeneral(2*(i-r)-1),:);
                if (min(tempL1(:,1))>10 && max(tempL1(:,1))<490 && min(tempL1(:,2))>10 && max(tempL1(:,2))<290) || ...
                        (min(tempL2(:,1))>10 && max(tempL2(:,1))<490 && min(tempL2(:,2))>10 && max(tempL2(:,2))<290)
                    flag2=1;
                    break
                end
                clear tempL temp temp1 temp2 temp3 temp4
            end
        end
        if flag2==0
            cycleR=[cycleR; i 1];
        end
        
        %check fp2--------------------------------
        flag2=0;
        for j=1:length(RFoot)
            if isempty(P(sub).Trajectory.Right(i).(RFoot{j}))==0
                tempR=P(sub).Trajectory.Right(i).(RFoot{j})(1:RFootOff(i-r)-REvent(i-r),:);
                if isempty(tempR)==0
                    if min(tempR(:,1))<490 || max(tempR(:,1))>1010 || min(tempR(:,2))<-10 || max(tempR(:,2))>610
                        flag2=1;
                        break
                    end
                end
                clear tempR
            end
        end
        if flag2==0
            for j=1:length(RFoot)
                [temp1 temp2 temp3 temp4]=vicon.GetTrajectory((Name{sub}),(LFoot{j}));
                temp=[temp1;temp2;temp3; temp4];
                temp=temp';
                tempL1=temp(RGeneral(2*(i-r)):RFootOff(i-r),:);
                tempL2=temp(REvent(i-r):RGeneral(2*(i-r)-1),:);
                if (min(tempL1(:,1))>510 && max(tempL1(:,1))<990 && min(tempL1(:,2))>10 && max(tempL1(:,2))<590) || ...
                        (min(tempL2(:,1))>510 && max(tempL2(:,1))<990 && min(tempL2(:,2))>10 && max(tempL2(:,2))<590)
                    flag2=1;
                    break
                end
                clear tempL
            end
        end
        if flag2==0
            cycleR=[cycleR; i 2];
        end
    end
end


% kinetics_________________________________________________________________
Ikl=0; Ikr=0;
for i=1:length(LKin)
    lm=strcmp((LKin{i}),NxsOutput);
    if max(lm)==0
        Ikl=[Ikl i];
    end
    rm=strcmp((RKin{i}),NxsOutput);
    if max(rm)==0
        Ikr=[Ikr i];
    end
    clear lm lr
end

if size(LEvent,2)>1
    for i=l+1:l+length(LEvent)-1
        if l>10
            break
        end
        ch0=find(cycleL(:,1)==i);
        if size(ch0,1)~=0
            for j=1:length(LKin)
                ch=find(Ikl==j);
                if size(ch,2)==0
                    temp=vicon.GetModelOutput((Name{sub}),(LKin{j}));
                    temp=temp';
                    P(sub).Kinetics.Left(i).(LKin{j})=temp(LEvent(i-l):LEvent(i-l+1),:);
                    clear temp temp2 Max Ind
                else
                    P(sub).Kinetics.Left(i).(LKin{j})=[];
                end
                clear ch
            end
        end
        clear ch0
    end
end

if size(REvent,2)>1
    for i=r+1:r+length(REvent)-1
        if r>10
            break
        end
        ch0=find(cycleR(:,1)==i);
        if size(ch0,1)~=0
            for j=1:length(RKin)
                ch=find(Ikr==j);
                if  size(ch,2)==0
                    temp=vicon.GetModelOutput((Name{sub}),(RKin{j}));
                    temp=temp';
                    P(sub).Kinetics.Right(i).(RKin{j})=temp(REvent(i-r):REvent(i-r+1),:);
                    clear temp temp2 Max Ind
                else
                    P(sub).Kinetics.Right(i).(RKin{j})=[];
                end
                clear ch
            end
        end
        clear ch0
    end
end


% GRF && COP_______________________________________________________________
Channels={'Fx' 'Fy' 'Fz' 'Mx' 'My' 'Mz' 'Cx' 'Cy' 'Cz'};

if size(LEvent,2)>1
    for i=l+1:l+length(LEvent)-1
        ch0=find(cycleL(:,1)==i);
        if size(ch0,1)~=0
            if i<=10
                for j=1:length(DeviceNames)
                    if strcmp(DeviceType{j},'ForcePlate')
                        for k=1:length(Channels)
                            ChannelID{k} = vicon.GetDeviceChannelIDFromName(DeviceID{j}, ceil(double(k)/3), Channels{k});
                            temp3=vicon.GetDeviceChannel(DeviceID{j}, ceil(double(k)/3), ChannelID{k});
                            temp3=temp3';
                            if DeviceID{j}==fp1%5 || DeviceID{j}==2
                                Left.FP1.(Channels{k})=temp3(LEvent(i-l)*10-10:LEvent(i-l+1)*10);
                            else
                                Left.FP2.(Channels{k})=temp3(LEvent(i-l)*10-10:LEvent(i-l+1)*10);
                            end
                            clear temp3
                        end
                    end
                end
                if isempty(ch0)==0
                    if cycleL(ch0,2)==1
                        for j=1:3
                            P(sub).ForcePlate.Left(i).Force(:,j)=Left.FP1.(Channels{j});
                            P(sub).ForcePlate.Left(i).Moment(:,j)=Left.FP1.(Channels{j+3});
                            P(sub).ForcePlate.Left(i).COP(:,j)=Left.FP1.(Channels{j+6});
                        end
                    else
                        for j=1:3
                            P(sub).ForcePlate.Left(i).Force(:,j)=Left.FP2.(Channels{j});
                            P(sub).ForcePlate.Left(i).Moment(:,j)=Left.FP2.(Channels{j+3});
                            P(sub).ForcePlate.Left(i).COP(:,j)=Left.FP2.(Channels{j+6});
                        end
                    end
                end
            end
        end
        clear ch0 Left
    end   
end

if size(REvent,2)>1
    for i=r+1:r+length(REvent)-1
        ch0=find(cycleR(:,1)==i);
        if size(ch0,1)~=0
            if i<=10
                for j=1:length(DeviceNames)
                    if strcmp(DeviceType{j},'ForcePlate')
                        for k=1:length(Channels)
                            ChannelID{k} = vicon.GetDeviceChannelIDFromName(DeviceID{j}, ceil(double(k)/3), Channels{k});
                            temp3=vicon.GetDeviceChannel(DeviceID{j}, ceil(double(k)/3), ChannelID{k});
                            temp3=temp3';
                            %                         Left(i-1).(Channels{j})=temp3(LEvent(k-l)*10:LEvent(k-l+1)*10);
                            if DeviceID{j}==fp1%5 || DeviceID{j}==2
                                Right.FP1.(Channels{k})=temp3(REvent(i-r)*10-10:REvent(i-r+1)*10);
                            else
                                Right.FP2.(Channels{k})=temp3(REvent(i-r)*10-10:REvent(i-r+1)*10);
                            end
                            clear temp3
                        end
                    end
                end
                if isempty(ch0)==0
                    if cycleR(ch0,2)==1
                        for j=1:3
                            P(sub).ForcePlate.Right(i).Force(:,j)=Right.FP1.(Channels{j});
                            P(sub).ForcePlate.Right(i).Moment(:,j)=Right.FP1.(Channels{j+3});
                            P(sub).ForcePlate.Right(i).COP(:,j)=Right.FP1.(Channels{j+6});
                        end
                    else
                        for j=1:3
                            P(sub).ForcePlate.Right(i).Force(:,j)=Right.FP2.(Channels{j});
                            P(sub).ForcePlate.Right(i).Moment(:,j)=Right.FP2.(Channels{j+3});
                            P(sub).ForcePlate.Right(i).COP(:,j)=Right.FP2.(Channels{j+6});
                        end
                    end
                end
            end
        end
        clear ch0 Right
    end    
end

lr=fieldnames(P(sub).Kinematics);
if length(lr)==1
    if strcmp(lr,'Left')
        l=size(P(sub).Kinematics.Left,2);
    else
        r=size(P(sub).Kinematics.Right,2);
    end
else
l=size(P(sub).Kinematics.Left,2);
r=size(P(sub).Kinematics.Right,2);
end

if isempty(P(sub).Kinetics.Left)==0
    l2=size(P(sub).Kinetics.Left,2);
    if l2<l
        P(sub).Kinetics.Left(l).(LKin{1})=[];
        P(sub).ForcePlate.Left(l).Force=[];
    end
end
if isempty(P(sub).Kinetics.Right)==0
    r2=size(P(sub).Kinetics.Right,2);
    if r2<r
        P(sub).Kinetics.Right(r).(RKin{1})=[];
        P(sub).ForcePlate.Right(r).Force=[];
    end
end

if l>=10 && r>=10
    flag=0;
else
    flag=1;
end
clear Left Right fp fp1 ChannelID DeviceID DeviceNames DeviceOutpotID DeviceRate DeviceType
% save('P','P')
disp('done')