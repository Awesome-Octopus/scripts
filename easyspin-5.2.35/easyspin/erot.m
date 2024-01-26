% erot  Rotation/transformation matrix from Euler angles
%
%   Rp = erot(Angles)
%   Rp = erot(alpha,beta,gamma)
%
%   [xc,yc,zc] = erot(...,'cols')
%   [xr,yr,zr] = erot(...,'rows')
%
%   Computes a 3x3 rotation/transformation matrix Rp from
%   a vector of 3 Euler angles.
%
%   Input
%   - Angles: vector containing the three
%     Euler angles (in radians) that define
%     the rotation. [alpha, beta, gamma] rotate
%     the coordinate system around [z,y',z'']
%     counterclockwise in in that order.
%   - alpha, beta, gamma: the three Euler angles
%     as defined above.
%   - 'cols', 'rows': tells erot to return either
%     the three rows or the three columns of the
%     rotation matrix separately
%
%   Output
%   - Rp: matrix for the passive rotation/coordinate transformation
%       vec2 = Rp*vec1 or
%       mat2 = Rp*mat1*Rp.'
%   - xc,yc,zc: three columns of the rotation matrix
%       such that Rp = [xc yc zc]
%   - xr,yr,zr: three rows of the rotation matrix, as
%       column vectors, such that Rp = [xr.'; yr.'; zr.']
