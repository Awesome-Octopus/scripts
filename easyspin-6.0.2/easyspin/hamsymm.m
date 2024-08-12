% hamsymm  Determine spin Hamiltonian symmetry 
%
%   PGroup = hamsymm(Sys)
%   [PGroup,R] = hamsymm(Sys)
%
%   Determines the point group of the Hamiltonian
%   of a spin sytem and the associated symmetry frame.
%
%   Input:
%   - Sys: Spin system structure
%
%   Output:
%   - PGroup: Schoenfliess point group symbol, one of
%     'Ci','C2h','D2h','C4h','D4h','S6','D3d','C6h','D6h','Th','Oh',Dinfh','O3'.
%   - R: Rotation matrix with the axes of the symmetry frame along columns.
