% resfreqs_matrix Compute resonance frequencies for frequency-swept EPR
%
%   ... = resfreqs_matrix(Sys,Exp)
%   ... = resfreqs_matrix(Sys,Exp,Opt)
%   [Pos,Int] = resfreqs_matrix(...)
%   [Pos,Int,Wid] = resfreqs_matrix(...)
%   [Pos,Int,Wid,Trans] = resfreqs_matrix(...)
%
%   Computes frequency-domain EPR line positions, intensities and widths using
%   matrix diagonalization.
%
%   Input:
%    Sys: spin system structure
%    Exp: experimental parameter settings
%      Field               static field, in mT
%      Range               frequency sweep range, [numin numax], in GHz
%      CenterField         frequency sweep range, [center sweep], in GHz
%      Temperature         temperature, in K
%      CrystalOrientation  nx3 array of Euler angles (in radians) for crystal orientations
%      CrystalSymmetry     crystal symmetry (space group etc.)
%      MolFrame            Euler angles (in radians) for molecular frame orientation
%      mwPolarization      'linear', 'circular+', 'circular-', 'unpolarized'
%      Mode                excitation mode: 'perpendicular', 'parallel', [k_tilt alpha_pol]
%    Opt: additional computational options
%      Verbosity           level of detail of printing; 0, 1, 2
%      Transitions         nx2 array of level pairs
%      Threshold           cut-off for transition intensity, between 0 and 1
%      Hybrid              0 or 1, switches hybrid mode off or on
%      HybridCoreNuclei    for hybrid mode, nuclei to include in exact core
%      Sites               list of crystal sites to include (default []: all)
%
%   Output:
%    Pos     line positions (in mT)
%    Int     line intensities
%    Wid     Gaussian line widths, full width half maximum (FWHM)
%    Trans   list of transitions included in the computation
