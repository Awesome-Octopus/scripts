% sphgrid  Spherical grid
%
%   grid = sphgrid(GridSymmetry,GridSize)
%   [grid,tri] = sphgrid(GridSymmetry,GridSize)
%
%   Returns a set of unique orientations together with
%   the covered solid angle by each, for a given symmetry.
%
%   Input:
%   - Symmetry:   the desired point group of the grid n Schoenflies notation
%                  ('Ci', 'Dinfh', 'C2h', etc.)
%                 only centrosymmetric point groups and C1 are supported.
%   - GridSize:   number of knots along a quarter of a meridian, at least 2
%
%   Output:
%   - grid: a structure with fields
%      .GridSymmetry input grid symmetry
%      .GridSize     input grid size
%      .phi,.theta:  m-element arrays of polar angles, in radians
%      .vecs:        3xm array of orientations (unit column vectors)
%      .weights:     m-element array of associated weights, sum is 4*pi
%   - tri: a structure with fields
%      .idx          triangulation array, one triangle per row
%      .areas        triangle areas (solid angles), sum is 4*pi
%
