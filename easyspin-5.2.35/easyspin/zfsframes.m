% zfsframes   D tensor in various frame conventions
%
%    zfsframes(D1,D2,D3)
%    zfsframes([D1 D2 D3])
%
%    Determines and displays the zero-field splitting
%    parameters D and E for all possible choices of the
%    eigenframe of the D tensor, with the given principal
%    values D1, D2, and D3. Conventions are listed.
%
%    Example:
%      zfsframes(-1,-1,2)   % axial tensor
%      zfsframes(-1,0,1)    % rhombic tensor
%      zfsframes(-1,0,2)    % general tensor
