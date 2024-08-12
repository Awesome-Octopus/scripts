% wignerd  Wigner rotation matrices
%
%   D = wignerd(J,alpha,beta,gamma);
%   D = wignerd(J,alpha,beta,gamma,phase);
%   Dm1m2 = wignerd(Jm1m2,alpha,beta,gamma);
%   Dm1m2 = wignerd(Jm1m2,alpha,beta,gamma,phase);
%   d = wignerd(J,beta);
%   d = wignerd(J,beta,phase);
%   dm1m2 = wignerd(Jm1m2,beta);
%   dm1m2 = wignerd(Jm1m2,beta,phase);
%
%   This function computes the Wigner rotation matrix with elements D^J_{m1,m2},
%   where m1 and m2 run from J to -J. If m1 and m2 are given in Jm1m2 = [J m1 m2],
%   only the single matrix element D^J_{m1,m2} is calculated.
%
%   Input:
%     J      ... rank J; J = 0, 1/2, 1, 3/2, 2, 5/2, etc.
%     Jm1m2  ... rank and two indices, [J m1 m2]; with m1,m2 = J, J-1, ..., -J
%     alpha, beta, gamma ... Euler angles, in radians
%             If only beta is given, alpha and gamma are assumed zero.
%     phase  ... sign convention for the rotation operator, '+' or '-'
%            '-': expm(-1i*alpha*Jz)*expm(-1i*beta*Jy)*expm(-1i*gamma*Jz)
%                (as used in Brink/Satcher, Zare, Sakurai, Varshalovich, Biedenharn/Louck, Mehring).
%            '+': expm(+1i*alpha*Jz)*expm(+1i*beta*Jy)*expm(+1i*gamma*Jz)
%                (as used in Edmonds, Mathematica).
%            The default sign convention is '-'.
%
%  Output: 
%     D      ... Wigner rotation matrix D^J, size (2J+1)x(2J+1)
%     Dm1m2  ... Wigner rotation matrix element D^J_{m1,m2}
%     d      ... reduced Wigner rotation matrix d^J, size (2J+1)x(2J+1)
%     dm1m2  ... reduced Wigner rotation matrix element d^J_{m1,m2}
%
%   The basis is ordered +J..-J from left to right and from top to bottom,
%   so that e.g. the output matrix element D(1,2) corresponds to D^J_{J,J-1}.
