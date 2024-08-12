% propint  Compute propagation operator 
%
%   U = propint(H0,H1,t,freq)
%   U = propint(H0,H1,t,freq,phase)
%   U = propint(H0,H1,t,freq,phase,n)
%   [U,SaveU] = propint(...)
%
%   U = propint(SaveU,t,freq)
%   U = propint(SaveU,t,freq,phase)
%
%   Computes a propagator by integrating
%
%      expm(-i*2*pi*(H0+H1*cos(2*pi*freq*t+phase))*t)
%
%   over the interval t(1) to t(2) using steps of
%   dt = 1/freq/n. Defaults: phase = 0, n = 256.
%   t = t2 is equivalent to t = [0 t2]. H0 and
%   H1 should be in frequency unit, t in complementary
%   time units, e.g. GHz and ns or MHz and �s.
%
%   SaveU is a cell array containing intermediate
%   results that can be reused by supplying SaveU
%   instead of H0 and H1 in the next call. The
%   speed-up is huge.
%
%   Example: An resonant pi/2 pulse for a spin-1 would be
%
%      Sz = sop(1,'z'); Sy = sop(1,'y'); % S=1 system
%      mwFreq = 10e3; tp = 0.010; % 10 GHz, 10 ns
%      U = propint(mwFreq*Sz,1/2/tp*Sy,tp,mwFreq)
%
%   Example: A resonance pi pulse on an two spin system
%
%      spins = [1/2 1/2]
%      Sz = sop(spins,'ze'); Sy = sop(spins,'ye');
%      mwFreq = 10e3; tp = 0.010; % 10 GHz, 10 ns
%      U = propint(mwFreq*Sz,1/tp*Sy,tp,mwFreq);
%      Density = -Sz; U*Density*U'
