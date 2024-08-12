% plegendre    Legendre polynomials and associated Legendre polynomials
%
%    y = plegendre(L,z)
%    y = plegendre(L,M,z)
%    y = plegendre(L,M,z,CSphase)
%
%    Computes the Legendre polynomial of degree L, or the associated Legendre
%    polynomial of degree L and order M, evaluated at z.
%
%    L and M are integers with L >= 0  and |M| <= L.
%    z can be a scalar or an array of real numbers with -1 <= z <= 1.
%
%    The optional input CSphase specifies whether the Condon-Shortley phase
%    (-1)^M should be included (true) or not (false). The default is true.
