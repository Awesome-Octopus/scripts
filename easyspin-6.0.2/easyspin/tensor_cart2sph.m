% tensor_cart2sph  Irreducbile spherical tensors from rank-2 cartesian tensor
%
%    [T0,T1,T2] = tensor_cart2sph(Tc)
%
% tensor_cart2sph converts a cartesian tensor (3x3 matrix) into
% its three irreducible spherical tensors (ranks 0, 1, and 2).
%
% Input:
%    Tc  ... cartesian tensor (real-valued); 3x3 matrix, or 3x1 array of
%            principal values, or single number if isotropic
%
% Output:
%    T0  ... rank-0 irreducbile spherical tensor (scalar)
%    T1  ... rank-1 irreducible spherical tensor (3x1 array)
%    T2  ... rank-2 irreducible spherical tensor (5x1 array)
