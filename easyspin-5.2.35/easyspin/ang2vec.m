% ang2vec  Polar angles to cartesian coordinates 
%
%   V = ang2vec(phi,theta)
%   [x,y,z] = ang2vec(phi,theta)
%
%   Converts angles to cartesian unit vector.
%   Angles have to be in radians. V is a 3xN
%   matrix containing the cartesian vectors
%   along columns. theta is the angle down
%   from the z axis, phi is the angle in the xy
%   plane between the x axis and the vectors xy
%   projection (counterclockwise). phi and
%   theta must be arrays with the same number of
%   elements.
