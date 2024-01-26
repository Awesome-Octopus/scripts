% orisel  Orientation selection 
%
%   orisel(Sys,Exp)
%   orisel(Sys,Exp,Opt)
%   Weights = ...
%   [Weights,Trans] = ...
%
%   Calculates weights for orientation selectivity.
%
%   Input:
%   - Sys    spin system
%   - Exp    experimental parameters
%      Field         magnetic field [mT]
%      mwFreq        spectrometer frequency [GHz]
%      ExciteWidth   FWHM excitation bandwidth [MHz]
%      Orientations  [phi;theta]
%   - Opt    options structure
%      nKnots       number of knots
%      Symmetry     point group: C1, Ci, D2h, or Dinfh
%
%   Output:
%   - Weights   vector of selectivity weights, one per orientation
%   - Trans     list of transitions included
%
%   If the function is called with no output arguments,
%   the selectivity weights are plotted.
