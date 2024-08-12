% vec2ang  Convert cartesian vectors to polar angles 
%
%   [phi,theta] = vec2ang(v)
%   ang = vec2ang(v)
%   ... = vec2ang(x,y,z)
%
%   Converts cartesian vectors in v or x, y, z to polar angles phi and theta.
%
%   Inputs:
%     v can be a 3xN array (list of column vectors) or Nx3 array (list of row
%        vectors). If 3x3, column vectors are assumed.
%     v can be a vector shorthand such as 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz'.
%     x,y,z are arrays with x, y, and z coordinates of the vectors.
%
%   Outputs:
%     theta is the angle down from the z axis to the vector v
%     phi is the angle between the x axis and the projection of the vector v
%        onto the xy plane, anticlockwise.
%
%   If only one output is requested, then an array with the phi values in the
%   first row and the theta values in the secon row is returned (ang = [phi; theta]).
%
%   Angles are in radians. The input vectors don't have to be normalized.
