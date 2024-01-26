% sphtri  Triangulation of standard grid
%
%   tri = sphtri(Symmetry,n)
%
%   Computes the triangulation of the standard
%   EasySpin grid with n knots along a quarter
%   of a meridian for point group Symmetry.
%   tri is a Nx3 matrix containing indices
%   into the grid vector.
%   The triangulation of point groups with open
%   phi intervals (C4h,C6h,S6,C2h) is done after
%   closing the phi interval.
