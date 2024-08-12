% diptensor  Calculate dipolar coupling tensor between two spins
%
%   T = diptensor(g1,g2,rvec)
%   T = diptensor(nuc1,nuc2,rvec)
%   T = diptensor(g1,nuc2,rvec)
%   T = diptensor(nuc1,g2,rvec)
%
% Calculates the dipolar coupling tensor T (3x3 matrix, in MHz) between two
% spins using the inter-spin distance vector rvec (in nm). For electrons,
% provide g tensor information in g1 and g2. For nuclei, provide the
% isotope in nuc1 and nuc2.
%
% The tensor T is for the Hamiltonian H = S1*T*S2, where S1 are the spin operators
% of the two spins, in units of hbar.
%
% The tensor T is defined in the same frame as the vector rvec and the g
% tensors g1 and g2, which is typically the molecular frame.
%
% Inputs:
%   g1     g tensor of electron spin 1: isotropic (1 number) or full tensor
%          (3x3 matrix) or principal values plus Euler angles {gpv,gFrame}
%          Euler angles are in radians and defined as in Sys.gFrame
%   g2     same for spin 2
%   nuc1   nuclear isotope: '1H', '14N', etc.
%   nuc2   same for spin 2
%   rvec   3-element vector pointing from spin1 to spin2, in nm
%
% Outputs:
%   T      dipolar coupling tensor (3x3 matrix), in MHz
