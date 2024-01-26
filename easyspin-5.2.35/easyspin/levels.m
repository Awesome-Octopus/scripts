% levels  Energy levels of a spin system 
%
%   En = levels(SpinSystem,Ori,B)
%   En = levels(SpinSystem,phi,theta,B)
%   [En,Ve] = levels(...);
%
%   Calculates energy levels of a spin system.
%
%   Input:
%   - SpinSystem: spin system specification structure
%   - phi,theta: vectors of orientation angles of
%     magnetic field (in radians)
%   - Ori: orientations of the field in the molecular frame
%       a) nx2 array of Euler angles (phi,theta), or
%       b) nx3 array of Euler angles (phi,theta,chi), or
%       c) 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz' for
%          special directions
%   - B: vector of magnetic field magnitudes (mT)
%
%   Output:
%   - En: array containing all energy eigenvalues (in MHz), sorted,
%     for all possible (phi,theta,B) or (Ori,B) combinations. Depending
%     on the dimensions of Ori, phi,theta and B, out can be up
%     to 4-dimensional. The dimensions are in the order phi,
%     theta (or Ori), field, level number.
