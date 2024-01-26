% pulse      Pulse definition
%
% [t,IQ] = pulse(Par)
% [t,IQ] = pulse(Par,Opt)
%
% [t,IQ,modulation] = pulse(Par,Opt)
%
% Input:
%   Par = structure containing the following fields:
%     Par.tp          = pulse length, in �s
%     Par.TimeStep    = time step for waveform definition, in �s (default:
%                       determined automatically based on pulse parameters)
%     Par.Flip        = pulse flip angle, in radians (see Ref. 1)
%                       (default: pi)
%     Par.Amplitude   = pulse amplitude, in MHz; ignored if Par.Flip or
%                       Par.Qcrit are given
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
%     Par.FrequencyResponse  = frequency axis in GHz and resonator frequency 
%                              response (ideal or experimental, real-valued
%                              input is interpreted as magnitude response)
%     Par.mwFreq             = microwave frequency for the experiment in GHz
%   or
%     Par.ResonatorFrequency = resonator center frequency in GHz
%     Par.ResonatorQL        = loaded resonator Q-value
%     Par.mwFreq             = microwave frequency for the experiment in GHz
%
%   Opt = optional structure with the following fields
%    Opt.OverSampleFactor = oversampling factor for the determination of the 
%                           time step (default: 10)
%
% Available pulse modulation functions:
%   - Amplitude modulation: rectangular, gaussian, sinc, halfsin, quartersin,
%                           sech, WURST, Gaussian pulse cascades (G3, Q3, 
%                           custom coefficients using 'GaussianCascade', see
%                           private/GaussianCascadeCoefficients.txt for 
%                           details), Fourier-series pulses (I-BURP 1/2,
%                           SNOB i2/i3, custom coefficients using 
%                           'FourierSeries',see 
%                           private/FourierSeriesCoefficients.txt for details)
%   - Frequency modulation: none, linear, tanh, uniformQ
%
% The parameters required for the different modulation functions are:
% Amplitude modulation:
% 'rectangular'         - none
% 'gaussian'            - tFWHM     = FWHM in �s
%                       Alternatively:
%                       - trunc     = truncation parameter (0 to 1)
% 'sinc'                - zerocross = width between the first zero-
%                                     crossing points in �s
% 'sech'                - beta      = dimensionless truncation parameter
%                       - n         = order of the secant function argument
%                                     (default: 1);
%                                     asymmetric sech pulses can be
%                                     obtained by specifying two values,
%                                     e.g. [6 1]
% 'WURST'               - nwurst    = WURST n parameter (determining the
%                                     steepness of the amplitude function)
% 'halfsin'             - none
% 'quartersin'          - trise     = rise time in �s for quarter sine
%                                     weighting at the pulse edges
% 'GaussianCascade'     - A0        = list of relative amplitudes
%                       - x0        = list of positions (in fractions of
%                                     tp)
%                       - FWHM     = list of FWHM (in fractions of tp)
% 'FourierSeries'       - A0        = initial amplitude coefficient
%                       - An        = list of Fourier coefficients for cos
%                       - Bn        = list of Fourier coefficients for sin
%
% Frequency modulation:
% For all frequency-modulated pulses, the frequency sweep range needs to be
% defined in Par.Frequency (e.g. Par.Frequency = [-50 50])
% 'linear'              no additional parameters
% 'tanh'                - beta      = dimensionless truncation parameter
% 'uniformQ'            The frequency modulation is calculated as the
%                       integral of the squared amplitude modulation
%                       function (for nth order sech/tanh pulses or in
%                       general to obtain offset-independent adiabaticity
%                       pulses, see Ref. 3.
%
% Output:   t          = time axis for defined waveform in �s
%           IQ         = real and imaginary part of the pulse function
%           modulation = structure with amplitude (modulation.A, in MHz),
%                        frequency (modulation.freq, in MHz) and phase
%                        (modulation.phase, in rad) modulation functions
%
% References:
% 1. The conversion from flip angles or critical adiabaticities to amplitudes 
%    is performed using the approximations described in:
%    Jeschke, G., Pribitzer, S., Doll, A. Coherence Transfer by Passage
%    Pulses in Electron Paramagnetic Resonance Spectroscopy.
%    J. Phys. Chem. B 119, 13570�13582 (2015). (DOI: 10.1021/acs.jpcb.5b02964)
% 2. Chirps with variable rate to compensate for the resonator bandwidth.
%    The bandwidth compensation is implemented as described in:
%    Doll, A., Pribitzer, S., Tschaggelar, R., Jeschke, G., Adiabatic and
%    fast passage ultra-wideband inversion in pulsed EPR
%    J. Magn. Reson. 230, 27-39 (2013). (DOI: 10.1016/j.jmr.2013.01.002)
%    and
%    Pribitzer, S., Doll, A. & Jeschke, G. SPIDYAN, a MATLAB library for
%    simulating pulse EPR experiments with arbitrary waveform excitation.
%    J. Magn. Reson. 263, 45�54 (2016). (DOI: 10.1016/j.jmr.2015.12.014)
% 3. Shaped pulses with offset-independent adiabaticity and determination
%    of their frequency modulation functions is described in:
%    Garwood, M., DelaBarre, L., The return of the frequency sweep: 
%    designing adiabatic pulses for contemporary NMR.
%    J. Magn. Reson. 153, 155-177 (2001). (DOI: 10.1006/jmre.2001.2340)
%
