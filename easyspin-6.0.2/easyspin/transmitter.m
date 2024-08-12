% transmitter      Simulation of/compensation for the effect of transmitter
%                  nonlinearity on a pulse
%
%  signal = transmitter(signal0,Ain,Aout,'simulate')
%  signal = transmitter(signal0,Ain,Aout,'compensate')
%
%  If the option 'simulate' is selected, transmitter() simulates the effect
%  of transmitter nonlinearity provided as output amplitude (Aout) as a
%  function of input amplitude (Ain) on the input signal signal0.
%
%  If the option 'compensate' is selected, transmitter() adapts the signal
%  provided in signal0 to compensate for the defined transmitter nonlinearity.
%
%  Input:
%   - signal0   = input signal vector, for 'simulate' the amplitude of signal0
%                 refers to the amplitude scale provided in Ain, for 'compensate'
%                 it refers to the amplitude scale Aout
%   - Ain       = input amplitudes
%   - Aout      = output amplitudes
%   - 'simulate'/'compensate'
%
%  Output:
%   - signal    = signal with simulated amplitude compression or signal compensated
%                 for the transmitter nonlinearity
%
