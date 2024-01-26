% spherharm    Spherical and tesseral harmonics
%
%   y = spherharm(L,M,theta,phi)
%   y = spherharm(L,M,theta,phi,'c')
%   y = spherharm(L,M,theta,phi,'s')
%
%   Computes spherical harmonic Y_L,M(theta,phi)
%   with L>=0 and |M|<L. theta and phi can be
%   scalars or arrays (of the same size). If they
%   are arrays, the spherical harmonic is computed
%   element-by-element.
%
%   If the options 'c' or 's' are included, real spherical
%   harmonics called tesseral harmonics (corresponding to
%   the commonly used orbitals in chemistry) are returned.
%   'c' specifies the function containing cos(m*phi), and
%   's' specifies the function containing sin(m*phi).
%   For the tesseral harmonics, M  is restricted
%   to non-negative values. For M=0, the 'c' and the 's'
%   tesseral harmonic and the spherical harmonic are identical.
%
%   The spherical harmonics include the Condon-Shortley
%   phase convention. The tesseral harmonics ('c' and 's')
%   do not include it.
