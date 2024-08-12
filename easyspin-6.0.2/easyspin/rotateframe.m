% rotateframe   Rotate a frame around a given axis
%
%    ang = rotateframe(ang0,nRot,rho)
%    [ang,R] = rotateframe(ang0,nRot,rho)
%
%    Rotates a frame described by the three Euler angles in ang0 around the axis
%    given in nRot by the rotation angles listed in rho.
%
%    Input:
%     ang0   Euler angles describing the frame orientation, in radians
%     nRot   letter or vector specifying the rotation axis
%              nRot = [1;0;0]
%              nRot = 'x' is the x axis
%              vectors do not need to be normalized 
%     rho    rotation angle, or array of rotation angles, for the rotation
%               around nRot, in radians
%
%    Output:
%     ang    list of Euler angle sets for the rotated frames, one per row.
%     R      list of rotation matrices for the rotated frames, one per row.
%
%    Example:
%       ang0 = [0 45 0]*pi/180;
%       nRot = [1;0;0];
%       rho = (0:30:180)*pi/180;
%       ang = rotateframe(ang0,nRot,rho);
