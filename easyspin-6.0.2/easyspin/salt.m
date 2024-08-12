% salt  ENDOR spectra simulation
%
%   salt(Sys,Exp)
%   salt(Sys,Exp,Opt)
%   spec = salt(...)
%   [rf,spec] = salt(...)
%   [rf,spec,info] = salt(...)
%
%   Input:
%   - Sys: spin system structure specification
%       S, g, gFrame, Nucs, A, AFrame, Q, QFrame etc.
%       lwEndor     FWHM Endor line width [Gaussian,Lorentzian]
%   - Exp: experiment specification
%       mwFreq        spectrometer frequency, in GHz
%       Field         magnetic field, in mT
%       Range         radiofrequency range [low,high], in MHz
%       nPoints       number of points
%       Temperature   temperature of the sample, by default off (NaN)
%       ExciteWidth   ENDOR excitation width, FWHM, in MHz
%   - Opt: simulation options
%       Verbosity           0, 1, 2
%       Method        'matrix', 'perturb1', 'perturb2'='perturb'
%       separate      '', 'components', 'transitions', 'sites', 'orientations'
%       GridSize      grid size;  N1, [N1 Ninterp]
%       Transitions, Threshold, GridSymmetry
%       Intensity, Enhancement, Sites
%
%   Output:
%   - rf:     the radiofrequency axis, in MHz
%   - spec:   the spectrum or spectra
%   - info:    structure with details of the calculation
%
%   If no output argument is given, the simulated spectrum is plotted.
