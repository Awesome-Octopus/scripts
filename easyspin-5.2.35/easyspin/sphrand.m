% sphrand  Generate random points over unit sphere
%
%   vecs = sphrand(N)
%   vecs = sphrand(N,k)
%   [phi,theta] = sphrand(...)
%   [x,y,z] = sphrand(...)
%
%   Generates N randomly distributed point over k
%   octants of the unit sphere. The underlying
%   distribution is uniform.
%
%   Input
%   - N: number of points
%   - k: number of octants, can be 1, 2, 4 or 8.
%          1   x>0,y>0,z>0
%          2   y>0,z>0
%          4   z>0 (upper hemisphere)
%          8   full sphere
%        4 is the default.
%
%    Output
%   - vecs: 3xN array with column vectors
%   - phi, theta: polar angles in radians
