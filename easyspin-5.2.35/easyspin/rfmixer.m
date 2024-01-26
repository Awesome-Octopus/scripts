% rfmixer  Digital up- or downconversion
%
%   [tOut,signalOut] = rfmixer(tIn,signalIn,mwFreq,type) 
%   [tOut,signalOut] = rfmixer(tIn,signalIn,mwFreq,type,Opt) 
%
%   Mixes the input signal with the LO frequency mwFreq. Depending on the 
%   selected mixer type, the function acts as a double-sideband (DSB) mixer,
%   single-sideband (SSB) mixer or performs IQ modulation, demodulation or 
%   frequency shifting.
%
%   Resampling of the input signal before up/downconversion can be
%   performed if the new sampling time step is given in Opt.dt.
%
%   If no outputs are requested, the results are plotted.
%
%   Input:
%   - tIn:        time axis, in �s
%   - signalIn:   input signal vector, in-phase part only for DSB,
%                 SSB or IQ demodulation mode, in-phase and quadrature
%                 part for IQ mixer and IQ frequency shift operation
%   - mwFreq:     LO frequency, in GHz
%                 (for IQ frequency shifts, the direction of the shift
%                 has to be included)
%   - type :      mode of operation, the options are:
%                 'DSB' = double sideband mixer
%                 'USB' = single sideband mixer, upper sideband selected
%                 'LSB' = single sideband mixer, lower sideband selected
%                 'IQmod' = IQ modulator
%                 'IQdemod' = IQ demodulator
%                 'IQshift' = IQ frequency shifter (up- or downconversion)
%   - Options:
%     Opt.dt =     time step for the output signal, in �s.
%                  If no time step for resampling is given and the input
%                  time axis is too large, a new time step is computed as
%                  1/(2*Opt.OverSampleFactor*maxfreq).
%     Opt.OverSampleFactor     = factor for determining new time step in
%                                signal resampling (default = 1.25).
%     Opt.InterpolationMethod  = interpolation method for signal resampling
%                                (default = 'spline').
%     Opt.BandwidthThreshold   = threshold for input bandwidth determination
%                                (default = 0.1).
%     Opt.HilbertThreshold     = threshold for amplitude and phase FT overlap,
%                                the Hilbert transform corresponds to the
%                                quadrature signal only if amplitude and phase
%                                FT do not overlap (Ref. 1). If the maximum of
%                                the product of the normalized amplitude and
%                                phase FT (extracted from real input signal
%                                by Hilbert transform) is larger than the
%                                specified threshold, recovery of the I and
%                                Q data from the real input signal is not
%                                possible (default = 0.05).
%     Opt.NoiseCutoffThreshold = threshold for selection of part of the input
%                                signal used to compute amplitude and cos(phase)
%                                FT overlap (the evaluation of the validity
%                                of the Hilbert transform is affected by
%                                the noise level in the signal and baseline).
%
%   Output:
%   - tOut:          time axis for the output signal, in �s.
%   - signalOut:     output signal vector, in-phase component only for DSB,
%                    SSB or IQmod, in-phase and quadrature component for
%                    IQdemod and IQshift
%
% References:
% 1. Bedrosian's product theorem is explained in:
%    Boashash, B., Estimating and interpreting the instantaneous frequency 
%    of a signal. I. Fundamentals
%    Proc. IEEE 80, 520 (1992). (DOI: 10.1109/5.135376)
%
