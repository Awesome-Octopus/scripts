% sham  Spin Hamiltonian 
%
%   [F,Gx,Gy,Gz] = sham(sys)
%   H = sham(sys, B)
%   ... = sham(sys, B,'sparse')
%
%   Constructs a spin Hamiltonian.
%
%   Input:
%   - 'sys': spin system specification structure
%   - 'B': 1x3 vector specifying the magnetic field [mT]
%   - 'sparse': if present, sparse instead of full matrices are returned
%
%   Output:
%   - 'H': complete spin Hamiltonian [MHz]
%   - 'F',Gx','Gy','Gz' ([MHz] and [MHz/mT])
%      H = F + B(1)*Gx + B(2)*Gy + B(3)*Gz
