% pepper  Computation of powder cw EPR spectra 
%
%   pepper(Sys,Exp)
%   pepper(Sys,Exp,Opt)
%   spec = pepper(...)
%   [x,spec] = pepper(...)
%   [x,spec,Trans] = pepper(...)
%
%   Calculates field-swept and frequency-swept cw EPR spectra.
%
%   Input:
%    Sys: parameters of the paramagnetic system
%      S, g, Nucs, A, Q, D, ee,
%      gFrame, AFrame, QFrame, DFrame, eeFrame
%      lw, lwpp
%      HStrain, gStrain, AStrain, DStrain
%      B2, B4, B6 etc.
%    Exp: experimental parameter settings
%      mwFreq              microwave frequency, in GHz (for field sweeps)
%      Range               sweep range, [sweepmin sweepmax], in mT (for field sweep)
%      CenterSweep         sweep range, [center sweep], in mT (for field sweeps
%      Field               static field, in mT (for frequency sweeps)
%      mwRange             sweep range, [sweepmin sweepmax], in GHz (for freq. sweeps)
%      mwCenterSweep       sweep range, [center sweep], in GHz (for freq. sweeps)
%      nPoints             number of points
%      Harmonic            detection harmonic: 0, 1 (default), 2
%      ModAmp              peak-to-peak modulation amplitude, in mT (field sweeps only)
%      mwPhase             detection phase (0 = absorption, pi/2 = dispersion)
%      Temperature         temperature, in K
%      CrystalOrientation  nx3 array of Euler angles (in radians) for crystal orientations
%      CrystalSymmetry     crystal symmetry (space group etc.)
%      MolFrame            Euler angles (in radians) for molecular frame orientation
%      mwPolarization      'linear', 'circular+', 'circular-', 'unpolarized'
%      Mode                excitation mode: 'perpendicular', 'parallel', k_tilt, [k_tilt alpha_pol]
%      Ordering            coefficient for non-isotropic orientational distribution
%    Opt: computational options
%      Method              'matrix', 'perturb1', 'perturb2'='perturb'
%      Output              'summed', 'separate'
%      Verbosity           0, 1, 2
%      nKnots              N1, [N1 Ninterp]
%      Transitions, Threshold
%      Symmetry, SymmFrame,
%      Intensity, Freq2Field, Sites
%
%   Output:
%    x       field axis (in mT) or frequency axis (in GHz)
%    spec    spectrum
%    Trans   transitions included in the calculation
%
%   If no output argument is given, the simulated spectrum is plotted.
