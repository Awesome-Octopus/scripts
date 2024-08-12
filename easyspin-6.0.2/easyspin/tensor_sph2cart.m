% tensor_sph2cart  Rank-2 cartesian tensor from irreducible spherical tensors
%
%    Tc = tensor_sph2cart(T0,T1,T2)
%
% tensor_sph2cart converts the three irreducible spherical tensors (rank 0, 1,
% and 2) to the corresponding rank-2 cartesian tensor.
%
% Input:
%    T0 ... rank-0 spherical tensor (1 element, real-valued)
%    T1 ... rank-1 spherical tensor (3 elements)
%    T2 ... rank-2 spherical tensor (5 elements)
%
% Output:
%    Tc ... rank-2 cartesian tensor (3x3 array, real-valued)
