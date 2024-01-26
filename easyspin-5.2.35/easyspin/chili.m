% chili    Simulation of cw EPR spectra in the slow motional regime
%
%   chili(Sys,Exp,Opt)
%   spc = chili(...)
%   [B,spc] = chili(...)
%
%   Computes the slow-motion cw EPR spectrum of systems with
%   one electron and one nuclear spin.
%
%   Sys: spin system structure
%
%     Sys.tcorr           rotational correlation time (in seconds)
%     Sys.logtcorr        log10 of rotational correlation time (in seconds)
%     Sys.Diff            diffusion rate (s^-1)
%     Sys.logDiff         log10 of diffusion rate (s^-1)
%
%         All fields can have 1 (isotropic), 2 (axial) or 3 (rhombic) elements.
%         Precedence: logtcorr > tcorr > logDiff > Diff.
%
%     Sys.DiffFrame       Euler angles of the diffusion tensor (default [0 0 0])
%     Sys.lw              vector with FWHM residual broadenings
%                         1 element:  GaussianFWHM
%                         2 elements: [GaussianFWHM LorentzianFWHM]
%                         field sweep: mT, frequency sweep: MHz
%     Sys.lwpp            peak-to-peak line widths, same format as Sys.lw
%     Sys.lambda          ordering potential coefficients
%                         [lambda20 lambda22 lambda40 lambda42 lambda44]
%     Sys.Exchange        Heisenberg exchange frequency (MHz)
%
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
%
%   Opt: simulation options
%      LLKM            basis size: [evenLmax oddLmax Kmax Mmax]
%      PostConvNucs    nuclei to include perturbationally via post-convolution
%      Verbosity       0: no display, 1: show info
%      nKnots          number of knots for powder simulation
%      Symmetry        symmetry to use for powder simulation
%
%   Output:
%     B      magnetic field axis vector, in mT
%     spc    simulated spectrum, arbitrary units
%
%     If no output arguments are specified, chili plots the
%     simulated spectrum.
