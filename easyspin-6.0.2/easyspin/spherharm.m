% spherharm    Complex- and real-valued spherical harmonics
%
%   y = spherharm(L,M,theta,phi)
%   y = spherharm(L,M,theta,phi,'r')
%
%   Computes spherical harmonic Y_L,M(theta,phi) with L>=0 and -L<=M<=L.
%
%   theta is the angle down from the z axis (colatitude), and phi is the
%   counterclockwise angle off the x axis in the xy plane (longitude).
%   theta and phi can be scalars or arrays (of the same size). Both angles
%   are assumed to be in units of radians.
%
%   If the option 'r' is included, the real-valued spherical harmonics
%   are returned, using cos(M*phi) for M>=0, and sin(abs(M)*phi) for M<0.
%
%   The complex-valued spherical harmonics evaluated by spherharm() include the
%   Condon-Shortley phase (-1)^M.
%
%   The sign of the real-valued harmonics is defined such that they give
%   nonnegative values near theta=0 and phi=0 for all L and M.
