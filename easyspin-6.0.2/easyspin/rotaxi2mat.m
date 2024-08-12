% rotaxi2mat   Compute rotation matrix from rotation axis and rotation angle
%
%    R = rotaxi2mat(n,rho)
%
%    Generates the matrix representing the rotation around the axis given by
%    the 3-element vector n by the angle rho (in radians). n does not have to
%    be normalized.
%
%    There are shortcuts for a few specific directions of the rotation axis:
%      n='x' implies n=[1;0;0]
%      n='y'         n=[0;1;0]
%      n='z'         n=[0;0;1]
%      n='xy'        n=[1;1;0]
%      n='xz'        n=[1;0;1]
%      n='yz'        n=[0;1;1]
%      n='xyz'       n=[1;1;1]
%
%    Example:
%      % A rotation by 2*pi/3 (120 degrees) around the axis [1;1;1]
%      % permutes the x, y and z axes: x->y, y->z, z->x
%      R = rotaxi2mat('xyz',2*pi/3)
%
%      % To apply the rotation, use:
%      v = [1;0;0];    % vector to rotate
%      v_rot = R.'*v;  % rotated vector; note transpose!
