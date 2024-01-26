% rotplane   Set of vectors/orientations in a plane
%
%     v = rotplane(n,chi)
%     v = rotplane(n,chi12,nChi)
%     [phi,theta] = ...
%
%     Computes vectors in a plane perpendicular to
%     a given direction n.
%
%     n      3-element vector specifying the rotation axis
%     chi    array of plane angles, in radians
%     chi12  [chi1 chi2], start and end angle
%            in plane, in radians
%            rotplane(n,[chi1 chi2],nChi) corresponds to
%            rotplane(n,linspace(chi1,chi2,nChi))
%    
%     v      array of vectors giving the orientations
%            in the plane; one vector per column
%     phi,theta   angles describing the orientations
