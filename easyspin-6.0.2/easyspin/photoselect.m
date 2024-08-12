% photoselect Calculate orientation-dependent photo-selection weight
%
%   weight = photoselect(tdm,ori,k,alpha)
%
% This function calculates a selection weight between 0 and 1 for photo-excited
% spin centers, 1 indicating full excitation and 0 indicating no
% excitation. It needs information about the electric transition dipole moment
% orientation (tdm), the spin center orientation (ori), the light
% excitation direction (k) and the light polarization (alpha).
%
% Inputs:
%   tdm         Orientation of the transition dipole moment vector in the
%               molecular frame. There are three ways to input this:
%               - a letter or letter combination, e.g. 'x', 'z', 'xz', '-y', etc.
%               - a three-element vector [mx my mz]
%               - two spherical angles [phim thetam] (in radians)
%   ori         [phi theta chi] or [phi theta] angles of lab frame in molecular frame.
%               phi and theta determine the direction of the lab z axis (zL, aligned
%               with B0), chi additionally determines the direction of lab x and y
%               (xL and yL). If chi is omitted, the integral over chi from
%               0 to 2*pi is calculated.
%   k           Orientation of the propagation direction of the light excitation beam
%               in the lab frame.
%               There are three ways to input this:
%               - a letter or letter combination, e.g. 'x', 'z', 'xz', '-y', etc.
%               - a three-element vector [kx ky kz]
%               - two spherical angles [phik thetak] (in radians)
%               The most common situation is k = 'y' (perpendicular to
%               B0) or k = [pi/2 pi/2] or k = [0; 1; 0]
%   alpha       Polarization angle (in radians) that indicates the orientation
%               of the E-field vector in the plane perpendicular to the propagation
%               direction k. Use alpha=NaN to indicate depolarized light.
%               For k = 'y', alpha = 0 or pi puts the E-field along the lab z
%               axis (zL), parallel to B0; alpha = pi/2 or -pi/2 puts the E-field
%               along the lab x axis (perpendicular to zL and B0)
%
% Outputs:
%   weight      photoselection weight, between 0 and 1
%               0 means no excitation at all, 1 indicates complete excitation
%
% Example:
%
%   k = 'y';  % propagation direction along lab y axis
%   photoselect(tdm,ori,k,pi/2)  % E-field along lab x axis
%   photoselect(tdm,ori,k,0)     % E-field along lab z axis
