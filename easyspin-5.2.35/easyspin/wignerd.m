% wignerd  Compute Wigner D-matrix
%
%   D = wignerd(J,angles);
%   D = wignerd(J,angles,'-');
%   D = wignerd(J,angles,'+');
%
%   This function computes the Wigner rotation matrix D^J_(m1,m2)
%   where m1 and m2 run from J to -J.
%
%   angles ... [alpha beta gamma], Euler angles in radians, or
%              beta, second Euler angle in radians
%   J      ... rank
%   D      ... Wigner rotation matrix, size (2J+1)x(2J+1)
%
%   If angles contains only one number, the Wigner small-d matrix
%   d^J_(m1,m2) is computed.
%
%   If '+' is given (or the third argument is omitted), wignerd computes
%   the (2J+1)x(2J+1) matrix
%     <Jm1|expm(1i*alpha*Jz)*expm(1i*beta*Jy)*expm(1i*gamma*Jz)|Jm2>
%     (sign convention as in Edmonds, Mathematica)
%   If '-' is given, wignerd returns the (2J+1)x(2J+1) matrix representation of
%     <Jm1|expm(-1i*alpha*Jz)*expm(-1i*beta*Jy)*expm(-1i*gamma*Jz)|Jm2>
%     (sign convention as in Brink/Satcher, Zare, Sakurai, Varshalovich)
%
%   The basis is ordered +J..-J from left to right and from top to bottom,
%   so that e.g. the output matrix element D(1,2) corresponds to D^J_(J,J-1).
