% vec2ang  Convert cartesian vectors to polar angles 
%
%   [phi,theta] = vec2ang(x)
%   ang = vec2ang(x)
%
%   Converts cartesian vectors in x to polar angles.
%   theta is the angle down from the z axis to the
%   vector, phi is the angle between the x axis and
%   the projection of the vector onto the xy plane,
%   anticlockwise.
%
%   If only one output is requested, then an array with
%   the phi values in the first row and the theta values
%   in the secon row is returned (ang = [phi; theta]).
%
%   x must be a 3xN array (list of column vectors) or
%   Nx3 array (list of row vectors). If 3x3, column vectors
%   are assumed.
%
%   Angles are in radians. The vectors in x don't have
%   to be normalized.
%
