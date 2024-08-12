% levels  Energy levels of a spin system 
%
%   E = levels(SpinSystem,Ori,B)
%   E = levels(SpinSystem,phi,theta,B)
%   [E,V] = levels(...)
%
%   Calculates energy levels of a spin system.
%
%   Input:
%   - SpinSystem: spin system structure
%   - phi,theta: vectors of polar angles (in radians) specifying orientation
%                of magnetic field in molecular frame
%   - Ori: orientations of the field in the molecular frame
%       a) nx2 array of Euler angles (phi,theta), or
%       b) nx3 array of Euler angles (phi,theta,chi), or
%       c) 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz' for special directions
%       phi and theta are the polar angles; the field direction is
%       independent of chi
%   - B: array of magnetic field magnitudes (mT)
%
%   Output:
%   - E: array containing all energy eigenvalues (in MHz), sorted,
%        for all possible (phi,theta,B) or (Ori,B) combinations. Depending
%        on the dimensions of Ori, phi, theta and B, E can be up
%        to 4-dimensional. The dimensions are in the order phi,
%        theta (or Ori), field, level number.
%   - V: array of eigenvectors
%
%   Example:
%
%    Sys = struct('S',5/2,'D',1000);
%    B = linspace(0,600,601);  % mT
%    E = levels(Sys,'z',B);
%    plot(B,E);
