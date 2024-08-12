% pulse      Pulse definition
%
% [t,IQ] = pulse(Par)
% [t,IQ] = pulse(Par,Opt)
%
% [t,IQ,modulation] = pulse(Par,Opt)
%
% Input:
%   Par = structure containing the following fields:
%     Par.tp          = pulse length, in µs
%     Par.TimeStep    = time step for waveform definition, in µs (default:
%                       determined automatically based on pulse parameters)
%     Par.Flip        = pulse flip angle, in radians (see Ref. 1)
%                       (default: pi), ignored if Par.Amplitude or
%                       Par.Qcrit are given
%     Par.Amplitude   = pulse amplitude, in MHz; ignored if Par.Qcrit is
%                       given
%     Par.Qcrit       = critical adiabaticity, used to calculate pulse
%                       amplitude for frequency-swept pulses [1]; if given
%                       takes precedence over Par.Amplitude and Par.Flip
%     Par.Frequency   = pulse frequency; center frequency for amplitude
%                       modulated pulses, [start-frequency end-frequency]
%                       for frequency swept pulses; (default: 0)
%     Par.Phase       = phase for the pulse in radians (default: 0 = +x)
%     Par.Type        = pulse shape name in a string with structure 'AM/FM'
%                       (or just 'AM'), where AM refers to the amplitude
%                       modulation function and FM to the frequency
%                       modulation function. The available options are
%                       listed below.
%                       If only a single keyword is given, the FM function
%                       is set to 'none'. For a pulse with constant
%                       amplitude, the AM needs to be specified as
%                       'rectangular'.
%                       Different AM functions can be multiplied by
%                       concatenating the keywords as 'AM1*AM2/FM'.
%                       (default: 'rectangular')
%     Par.*           = value for the pulse parameters defining the
%                       specified pulse shape. The pulse parameters
%                       required for each of the available modulation
%                       functions are listed below.
%     Par.I, .Q, .IQ  = I and Q data describing an arbitrary pulse.
%                       The time axis is reconstructed based on Par.tp
%                       and the length of the I and Q vectors, all other
%                       input parameters (Amplitude, Flip, Frequency,
%                       Phase, etc.) are ignored.
%   To compensate for the resonator bandwidth to get uniform adiabaticity
%   (see Ref. 2), define:
%     Par.FrequencyResponse  = frequency axis, in GHz, and resonator frequency 
%                              response (ideal or experimental, real-valued
%                              input is interpreted as magnitude response)
%     Par.mwFreq             = microwave frequency for the experiment, in GHz
%   or
%     Par.ResonatorFrequency = resonator center frequency, in GHz
%     Par.ResonatorQL        = loaded resonator Q-value
%     Par.mwFreq             = microwave frequency for the experiment, in GHz
%
%   Opt = optional structure with the following fields
%    Opt.OverSampleFactor = oversampling factor for the determination of the 
%                           time step (default: 10)
%
% Available pulse modulation functions:
%   - Amplitude modulation: rectangular, gaussian, sinc, halfsin, quartersin,
%                           tanh2, sech, WURST, Gaussian pulse cascades (G3, G4, Q3, Q5
%                           custom coefficients using 'GaussianCascade', see
%                           private/GaussianCascadeCoefficients.txt for 
%                           details), Fourier-series pulses (I-BURP 1/2,
%                           E-BURP 1/2, U-BURP, RE-BURP, SNOB i2/i3, 
%                           custom coefficients using 'FourierSeries', see 
%                           private/FourierSeriesCoefficients.txt for details)
%   - Frequency modulation: none, linear, tanh, uniformQ
%
% The parameters required for the different modulation functions are:
