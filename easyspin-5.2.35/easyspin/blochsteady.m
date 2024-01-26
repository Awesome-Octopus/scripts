% blochsteady   Steady-state solutions of Bloch equations
%
%  blochsteady(g,T1,T2,deltaB0,B1,ModAmp,ModFreq)
%  blochsteady(g,T1,T2,deltaB0,B1,ModAmp,ModFreq,Options)
%  [t,My] = blochsteady(...)
%  [t,Mx,My,Mz] = blochsteady(...)
%
%  Computes periodic steady-state solution of the Bloch equations
%  for a single spin-1/2 in the presence of field modulation.
%
%  Inputs:
%    g        g value of the electron spin (S = 1/2)
%    T1       longitudinal relaxation time, us
%    T2       transverse relaxation time, us
%
%    deltaB0  offset from resonance field, mT
%    B1       microwave field amplitude, mT
%    ModAmp   peak-to-peak field modulation amplitude, mT
%    ModFreq  modulation frequency, kHz
%
%  Outputs:
%    t        time axis, us
%    Mx       in-phase signal (dispersion)
%    My       quadrature signal (absorption)
%    Mz       longitudinal magnetization
