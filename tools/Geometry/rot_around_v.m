function R = rot_around_v(e_rot,phi)
% ROT_AROUND_V - matrix for rotation PHI radians around E_ROT
%   
% Calling:
% R = rot_around_v(e_rot,phi)



%   Copyright © 2002 Bjorn Gustavsson, <bjorn.gustavsson@irf.se>
%   This is free software, licensed under GNU GPL version 2 or later

r = [cos(phi) -sin(phi) 0;
     sin(phi)  cos(phi) 0;
     0         0        1];
e_Rot = e_rot./sum(e_rot.^2)^(1/2);
az = atan2(e_Rot(1),e_Rot(2));
ze = acos(e_Rot(3));
T = [1 0        0;
     0 cos(ze)  -sin(ze);
     0 sin(ze) cos(ze)]*...
    [cos(az)  -sin(az) 0;
     sin(az) cos(az) 0;
     0         0        1];

R = T'*r*T;
