% eulang  Euler angles from rotation matrix
%
%   angles = eulang(R)
%   [alpha,beta,gamma] = eulang(R)
%
%   Returns the three Euler angles alpha, beta and gamma (in radians) of the
%   rotation matrix R, which must be a 3x3 real matrix with determinant very
%   close to +1.
%
%   [alpha,beta,gamma] and [alpha+-pi,-beta,gamma+-pi]
%   give the same rotation matrix. eulang() returns the
%   set with beta>=0.
%
%   If the matrix is close to orthogonal, a neighoring orthogonal matrix is
%   calculated using singular-value decomposition.
