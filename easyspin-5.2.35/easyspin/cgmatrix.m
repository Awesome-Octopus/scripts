% cgmatrix     Transformation matrix between coupled and uncoupled representations
%
%   U2C = cgmatrix(S1,S2)
%
%   This function calculates the matrix representing the transformation from
%   uncoupled to coupled representation for a spin system consisting of two
%   spins with quantum numbers S1 and S2.
%
%   To transform vectors and matrices from the uncoupled to the coupled
%   representation, use
%      U2C = cgmatrix(S1,S2);
%      psi_c = U2C*psi_u;
%      H_c = U2C*H_u*U2C';
%
%   Input:
%      S1, S2: spin quantum numbers (1/2, 1, 3/2, 2, 5/2, etc)
%
%   Output:
%      U2C:    transformation matrix
%
%  Basis ordering convention - uncoupled basis:
%     top level: decreasing mS1 = S1,...,-S1
%     next level: decreasing mS2 = S2,...,-S2
%     |mS1,mS2> = 
%       |S1,S2>, |S1,S2-1>, ..., |S1-1,S2>, |S1-1,S2-1>, ..., |-S1,-S2>
%  Basis ordering convention - coupled basis
%     top level:  decreasing total spin St = S1+S2...|S1-S2|
%     next level: decreasing mSt = St...-St
%     |Stot,mStot> = 
%       |S1+S2,S1+S2>,|S1+S2,S1+S2-1>,...|S1+S2-1,S1+S2-1>,...|abs(S1-S2),-abs(S1-S2)>
%
%  cgmatrix() uses clebschgordan() to calculate the matrix elements of U2C.
