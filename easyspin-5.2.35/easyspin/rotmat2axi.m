% rotmat2axi   Rotation axis and angle form rotation matrix
%
%    [n,rho] = rotmat2axi(R)
%
%    Determines the axis n (a 3-element vector) and the angle
%    rho (in  radians) of the rotation described by the 3x3
%    rotation matrix R.
%
%    Example:
%      % A rotation matrix describing a rotation around y
%      R = erot(0,pi/3,0);
%      % rotmax2axi yields axis y and angle pi/3
%      [n,rho] = rotmat2axi(R)
