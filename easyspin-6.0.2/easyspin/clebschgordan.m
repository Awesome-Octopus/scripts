% clebschgordan  Clebsch-Gordan coefficient
%
%   v = clebschgordan(j1,j2,j,m1,m2,m)
%   v = clebschgordan(jm1,jm2,jm)
%   v = clebschgordan(j,m)
%
%   Returns the Clebsch-Gordan coefficient, also called vector coupling
%   coefficient
%
%      (j1,j2,m1,m2|j1,j2,j,m)
%
%   involved in the coupling of two angular momenta j1 and j2 to
%   resultant angular momemtum j.
%
%   j1,j2,m1,m2 are the quantum number for the uncoupled representation,
%   and j1,j2,j,m are the quantum numbers for the coupled representation.
%
%   Definitons for alternative input forms
%   a)  jm1 = [j1 m1], jm2 = [j2 m2], jm = [j m]
%   b)  j = [j1 j2 j], m = [m1 m2 m]
