% eulang  Euler angles from rotation matrix
%
%   Angles = eulang(R)
%   [alpha,beta,gamma] = eulang(R)
%
%   Returns the three Euler angles [alpha, beta, gamma]
%   (in radians) of the rotation matrix R, which must be
%   a 3x3 real matrix with determinant very close to +1.
%
%   For a definition of the Euler angles, see erot().
%
%   [alpha,beta,gamma] and [alpha+-pi,-beta,gamma+-pi]
%   give the same rotation matrix. eulang() returns the
%   set with beta>=0.
