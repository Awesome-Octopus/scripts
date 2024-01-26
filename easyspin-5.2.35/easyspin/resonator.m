% resonator      Simulation of/compensation for the effect of the resonator
%                on a pulse
%
%  [tOut,signalOut] = resonator(tIn,signalIn,mwFreq,nu,TransferFunction,'simulate')
%  [tIn,signalIn] = resonator(tOut,signalOut,mwFreq,nu,TransferFunction,'compensate')
%  [tOut,signalOut] = resonator(tIn,signalIn,mwFreq,nu0,QL,'simulate')
%  [tIn,signalIn] = resonator(tOut,signalOut,mwFreq,nu0,QL,'compensate')
%  ... = resonator(...,Opt)
%
%  If the option 'simulate' is selected, resonator() simulates the effect
%  of the resonator on the input signal signal0.
%  If the option 'compensate' is selected, resonator() adapts the pulse 
%  shape to compensate for the resonator transfer function.
%
%  The resonator transfer function can be given directly by providing
%  a frequency axis and the corresponding transfer function as the second
%  and third input argument. If the input is real, it is assumed to be the
%  magnitude response and the phase response is estimated as described in
%  reference 1. The transfer function is expanded over the required frequency
%  range by exponentially approaching an ideal transfer function with
%  parameters estimated by fit to the provided transfer function.
%
%  Alternatively, the resonator center frequency and the loaded Q-value
%  can be provided as the second and third input argument and the transfer
%  function is calculated based on the ideal transfer function of an RLC
%  series circuit (for details see reference 1).
%
% References:
% 1. Doll, A., Pribitzer, S., Tschaggelar, R., Jeschke, G., Adiabatic and
%    fast passage ultra-wideband inversion in pulsed EPR
%    J. Magn. Reson. 230, 27-39 (2013). (DOI: 10.1016/j.jmr.2013.01.002)
%
%  Input:
%   - t0                   = time axis for the input signal (in �s)
%   - signal0              = input signal vector
%   - mwFreq               = microwave frequency for the input signal
%                            in GHz
%   - nu/nu0               = frequency axis for the resonator transfer 
%                            function (in GHz) or resonator center 
%                            frequency (in GHz)
%   - TransferFunction/QL = resonator transfer function or magnitude
%                            response or loaded Q-value
%   - 'simulate'/'compensate'
%   - Options structure with the following fields:
%        Opt.CutoffFactor     = cutoff factor for truncation of the impulse 
%                               response function (used for resonator simulation) 
%                               (default:1/1000)
%        Opt.TimeStep         = time step in �s (if it is not provided the
%                               ideal time step is estimated based on the 
%                               Nyquist condition with an oversampling factor)
%        Opt.OverSampleFactor = oversampling factor for the determination of the 
%                               time step (default: 10)
%
%  Output:
%   - t         = time axis for the output signal (in �s)
%   - signal    = signal modified by the resonator transfer function or
%                 compensated for the resonator transfer function
%
