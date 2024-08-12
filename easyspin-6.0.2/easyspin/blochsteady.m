% blochsteady   Steady-state solutions of Bloch equations
%
%  blochsteady(g,T1,T2,DeltaB0,B1,modAmp,modFreq)
%  blochsteady(g,T1,T2,DeltaB0,B1,modAmp,modFreq,Options)
%  [t,My] = blochsteady(...)
%  [t,Mx,My,Mz] = blochsteady(...)
%
%  Computes periodic steady-state solution of the Bloch equations
%  for a single spin-1/2 in the presence of a sinusoidal field modulation.
%
%  Inputs:
%    g        g value of the electron spin (S = 1/2)
%    T1       longitudinal relaxation time constant, µs
%    T2       transverse relaxation time constant, µs
%
%    DeltaB0  offset from resonance field, mT
%    B1       microwave field amplitude, mT
%    modAmp   peak-to-peak modulation field amplitude, mT
%    modFreq  modulation frequency, kHz
%
%    Options  calculation options
%      .Verbosity   whether to print information (0 or 1; 0 default)
%      .nPoints     number of points, chosen automatically by default
%      .kmax        highest Fourier order, chosen automatically by default
%      .Method      calculation method for time-domain signal
%                   'td'   explicit evolution in time-domain
%                   'fft'  using inverse Fourier transform (default)
%
%  Outputs:
%    t        time axis, µs
%    Mx       in-phase signal (dispersion)
%    My       quadrature signal (absorption)
%    Mz       longitudinal magnetization
