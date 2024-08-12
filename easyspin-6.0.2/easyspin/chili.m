% chili    Simulation of cw EPR spectra in the slow motional regime
%
%   chili(Sys,Exp,Opt)
%   spc = chili(...)
%   [B,spc] = chili(...)
%   [nu,spc] = chili(...)
%   [B,spc,info] = chili(...)
%   [nu,spc,info] = chili(...)
%
%   Computes a slow-motion cw EPR spectrum.
%
%   Sys: spin system structure
%
%     Sys.tcorr       rotational correlation time (in seconds)
%     Sys.logtcorr    log10 of rotational correlation time (in seconds)
%     Sys.Diff        diffusion rate (s^-1)
%     Sys.logDiff     log10 of diffusion rate (s^-1)
%
%         All fields can have 1 (isotropic), 2 (axial) or 3 (rhombic) elements.
%         Precedence: logtcorr > tcorr > logDiff > Diff.
%
%     Sys.DiffFrame   Euler angles describing the orientation of the
%                     diffusion tensor in the molecular frame (default [0 0 0])
%     Sys.lw          vector with FWHM residual broadenings
%                     1 element:  GaussianFWHM
%                     2 elements: [GaussianFWHM LorentzianFWHM]
%                     units: mT for field sweeps, MHz for frequency sweeps
%     Sys.lwpp        peak-to-peak line widths, same format as Sys.lw
%     Sys.Exchange    spin exchange rate (microsecond^-1)
%     Sys.Potential   orientational potential coefficients
%                       [L1 M1 K1 lambda1; L2 M2 K2 lambda2; ...]
%
%    Exp: experimental parameter settings
%      mwFreq         microwave frequency, in GHz (for field sweeps)
%      Range          sweep range, [sweepmin sweepmax], in mT (for field sweeps)
%      CenterSweep    sweep range, [center sweep], in mT (for field sweeps)
%      Field          static field, in mT (for frequency sweeps)
%      mwRange        sweep range, [sweepmin sweepmax], in GHz (for freq. sweeps)
%      mwCenterSweep  sweep range, [center sweep], in GHz (for freq. sweeps)
%      nPoints        number of points
%      Harmonic       detection harmonic: 0, 1, 2
%      ModAmp         peak-to-peak modulation amplitude, in mT (field sweeps only)
%      mwPhase        detection phase (0 = absorption, pi/2 = dispersion)
%      Temperature    temperature, in K
%
%   Opt: simulation options
%      LLMK           basis set parameters, [evenLmax oddLmax Mmax Kmax]
%      evenK          whether to use only even K values (true/false)
%      highField      whether to use the high-field approximation (true/false)  
%      pImax          maximum nuclear coherence order for basis
%      GridSize       grid size for powder simulation
%      PostConvNucs   nuclei to include perturbationally via post-convolution
%      Verbosity      0: no display, 1: show info
%      GridSymmetry   grid symmetry to use for powder simulation
%      separate       subspectra output, '' (default) or 'components'
%
%   Output:
%     B               magnetic field axis vector, in mT (for field sweeps)
%     nu              frequency axis vector, in GHz (for frequency sweeps)
%     spc             simulated spectrum, arbitrary units
%     info            structure containing details about the calculation
%
%     If no output arguments are specified, chili plots the simulated spectrum.
