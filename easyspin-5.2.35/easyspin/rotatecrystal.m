% rotatecrystal   Rotate a crystal around a given axis
%
%    angles_rot = rotatecrystal(CryOri,nRot,rho)
%
%    Rotates the crystal frame descibed by the three Euler angles
%    in CryOri around the axis given in nL by the rotation angles
%    listed in rho.
%
%    Input
%     CryOri: Euler angles for the crystal->lab frame transformation
%             (same as Exp.CrystalOrientation)
%     nRot:   rotation axis, in lab coordinates
%             E.g, nRot = [1;0;0] is the lab x axis, xL
%     rho:    rotation angle, or list of rotation angles, for the rotation
%             around nRot
%
%    Output:
%     angles_rot: list of Euler angle sets for the rotated frames, one per row.
%                 (can be directly used as Exp.CrystalOrientation in
%                 simulations).
%
%    Example:
%       COri0 = [0 45 0]*pi/180;
%       nRot = [1;0;0];
%       rho = (0:30:180)*pi/180;
%       COri_list = rotatecrystal(nRot,rho);
%       Exp.CrystalOrientation = COri_list;
