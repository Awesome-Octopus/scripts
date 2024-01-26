% symm  Determines spin Hamiltonian symmetry 
%
%   Group = symm(Sys)
%   [Group,R] = symm(Sys)
%
%   Determines the point group of the Hamiltonian
%   of a spin sytem together with its symmetry frame.
%
%   Input:
%   - Sys: Spin system specification structure
%
%   Output:
%   - Group: Schoenfliess point group symbol, one of
%     'Ci','C2h','D2h','C4h','D4h','S6','D3d',
%     'C6h','D6h','Th','Oh',Dinfh','O3'.
%   - R: Rotation matrix containing the axes of the
%     symmetry frame along columns.
