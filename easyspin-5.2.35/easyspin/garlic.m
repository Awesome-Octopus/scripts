% garlic    Simulates isotropic and fast-motion cw EPR spectra 
%
%   garlic(Sys,Exp)
%   garlic(Sys,Exp,Opt)
%   spec = ...
%   [B,spec] = ...
%
%   Computes the solution cw EPR spectrum of systems with
%   an unpaired electron and arbitrary numbers of nuclear spins.
%
%   Sys: spin system structure
%     g            isotropic g factor or 3 principal values of g
%     Nucs         string with comma-separated list of isotopes
%     n            vector of number of equivalent nuclei (default all 1)
%     A            vector of hyperfine couplings [MHz]
%     lw           vector with FWHM line widths [mT]
%                   1 element:  GaussianFWHM
%                   2 elements: [GaussianFWHM LorentzianFWHM]
%     lwpp         peak-to-peak line widths [mT], same format as Sys.lw
%
%     tcorr        correlation time for fast-motion linewidths [s]
%                  If omitted or zero, the isotropic spectrum is computed.
%     logtcorr     log10 of the correlation time for fast-motion linewidths [s]
%                  If logtcorr is given, tcorr is ignored.
%
%   Exp:  experimental parameter settings
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
%
%  Opt:  simulation parameters
%      Verbosity    log level (0 none, 1 normal, 2 very verbose)
%      Method       method used to calculate line positions:
%                     'exact' - Breit-Rabi solver
%                     'perturb1','perturb2','perturb3','perturb4','perturb5'
%      AccumMethod  method used to construct spectrum
%                     'binning' - stick spectrum with line shape convolution
%                     'template' - interpolative construction using line shape template
%                     'explicit' - explicit evaluation of line shape for each line
%      IsoCutoff    relative isotopologue abundance cutoff threshold
%                     between 0 and 1, default 1e-6
%
%   Output
%     B                magnetic field axis [mT]
%     spec             spectrum [arbitrary units]
%
%     If no output parameter is specified, the simulated spectrum
%     is plotted.
