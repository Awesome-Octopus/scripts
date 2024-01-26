% endorfrq  Compute ENDOR frequencies and intensities 
%
%    Pos = endorfrq(Sys,Exp)
%    Pos = endorfrq(Sys,Exp,Opt)
%    [Pos,Int] = endorfrq(...)
%    [Pos,Int,Tra] = endorfrq(...)
%
%    Sys:  spin system structure
%
%    Exp: experimental parameter structure
%      mwFreq              spectrometer frequency [GHz]
%      Field               magnetic field [mT]
%      Temperature         in K (optional, by default off (NaN))
%      ExciteWidth         FWHM of excitation [MHz] (optional, default inf)
%      Range               frequency range [MHz] (optional, default [])
%      CrystalOrientation  nx3 array of Euler angles (in radians) for crystal orientations
%      CrystalSymmetry     crystal symmetry (space group etc.)
%      MolFrame            Euler angles (in radians) for molecular frame orientation
%
%    Opt: computational options structure
%      Verbosity          log level (0, 1 or 2)
%      Transitions        [mx2 array] of level pairs
%      Threshold          for transition selection.
%      Nuclei             index for nuclei (1 to #nuclei)
%      Intensity          'on','off'
%      Enhancement        'on','off'
%      Sites              list of crystal sites to include (default []: all)
%
%    See the documentation for a detailed description of the arguments.
