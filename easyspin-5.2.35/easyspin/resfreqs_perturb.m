% resfreqs_perturb Compute resonance frequencies for frequency-swept EPR
%
%   ... = resfreqs_perturb(Sys,Exp)
%   ... = resfreqs_perturb(Sys,Exp,Opt)
%   [Pos,Int] = resfreqs_perturb(...)
%   [Pos,Int,Wid] = resfreqs_perturb(...)
%   [Pos,Int,Wid,Trans] = resfreqs_perturb(...)
%
%   Computes frequency-domain EPR line positions, intensities and widths using
%   perturbation theory.
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
%      PerturbOrder        perturbation order; 1 or 2
%      Sites               list of crystal sites to include (default []: all)
%
%   Output:
%    Pos     line positions (in mT)
%    Int     line intensities
%    Wid     Gaussian line widths, full width half maximum (FWHM)
%    Trans   list of transitions included in the computation
