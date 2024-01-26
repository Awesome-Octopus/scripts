% rotaxi2mat   Compute rotation matrix from axis and angle
%
%    R = rotaxi2mat(n,rho)
%
%    Generates the matrix representing the rotation
%    around the axis given by the 3-element vector n
%    by the angle rho (in radians).
%
%    There are three shortcuts for x, y, and z axis:
%      n=1 implies n=[1;0;0]
%      n=2 means   n=[0;1;0]
%      n=3 gives   n=[0;0;1].
%
%    Example:
%      % A rotation by 2*pi/3 around the axis [1;1;1]
%      % permutes the x, y and z axes.
%      R = rotaxi2mat([1;1;1],2*pi/3)
