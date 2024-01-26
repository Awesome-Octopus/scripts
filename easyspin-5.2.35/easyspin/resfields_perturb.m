% resfields_perturb  Compute resonance fields for cw EPR 
%
%   ... = resfields_perturb(Sys,Exp)
%   ... = resfields_perturb(Sys,Exp,Opt)
%   [Pos,Int] = resfields_perturb(...)
%   [Pos,Int,Wid] = resfields_perturb(...)
%   [Pos,Int,Wid,Trans] = resfields_perturb(...)
%
%   Computes cw EPR line positions, intensities and widths using
%   perturbation theory.
%
%   Input:
%    Sys: spin system structure
%    Exp: experimental parameter settings
%      mwFreq              microwave frequency, in GHz
%      Range               field sweep range, [Bmin Bmax], in mT
%      CenterField         field sweep range, [center sweep], in mT
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
