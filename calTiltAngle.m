function [tiltAP, tiltML] = calTiltAngle(tilt1, tilt2, tilt3, tilt4)

d = (tilt2(:,1:3) - tilt4(:,1:3));
d = sqrt(d(:,1).^2 + d(:,2).^2 + d(:,3).^2);
dz = tilt2(:,3) - tilt4(:,3);
tiltAP = asind(dz ./ d);

d = (tilt1(:,1:3) - tilt3(:,1:3));
d = sqrt(d(:,1).^2 + d(:,2).^2 + d(:,3).^2);
dz = tilt1(:,3) - tilt3(:,3);
tiltML = asind(dz ./ d);




