% sigeq  Thermal equilibrium spin density 
%
%   sigma = sigeq(Ham, Temp)
%   sigma = sigeq(Ham, Temp, 'pol')
%
%   Calculates the equilibrium density matrix
%   at temperature Temp of the system described by
%   the Hamiltonian Ham.
%
%   Input:
%   - Ham: Hamiltonian of the spin system (in MHz)
%   - Temp: temperature of the system (in K)
%   - 'pol': if present, only the polarization part
%       will be returned.
%
%   Output:
%   - sigma: density matrix at thermal equilibrium
%       in the same basis as Ham
