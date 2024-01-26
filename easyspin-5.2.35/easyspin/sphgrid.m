% sphgrid  Spherical triangular grids
%
%   vecs = sphgrid(Symmetry,nKnots)
%   [vecs,Weights] = sphgrid(Symmetry,nKnots,'c')
%   [phi,theta] = sphgrid(Symmetry,nKnots)
%   [phi,theta,Weights] = sphgrid(Symmetry,nKnots)
%
%   Returns a set of unique orientations together with
%   the covered solid angle by each, for a given
%   symmetry.
%
%   Input:
%   - Symmetry: the point group in Schoenflies notation
%     ('Ci', 'Dinfh', 'C2h', etc.)
%   - nKnots: number of knots along a quarter of a
%     meridian, must be at least 2
%
%   Output:
%   - vecs: an 3xm array of orientations (unit column vectors)
%   - Weights: associated solid angles, 1xm vector, sum is 4*pi
%   - phi,theta: polar angles in radians
%
%   Only centrosymmetric point groups and C1 are supported.
