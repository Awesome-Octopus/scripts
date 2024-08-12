% rapidscan2spc  Conversion of sinusoidal rapid-scan signal to EPR spectrum
%
%   [B,spc] = rapidscan2spc(rsSignal,rsAmp,rsFreq,g)
%
% Converts a sinusoidal rapid-scan time-domain EPR signal to the corresponding
% field-swept EPR spectrum.
%
% Inputs:
%   rsSignal time-domain rapid-scan signal (My+1i*Mx) over one period of
%            the field modulation (vector); My is the absorption, and Mx
%            is the dispersion
%   rsAmp    peak-to-peak modulation amplitude (in mT)
%   rsFreq   modulation frequency (in kHz)
%   g        (optional) g value, needed for field-domain rapid scan.
%            If omitted, g = 2.0023193 is assumed.
%
% Outputs:
%   dB       field offset axis (in mT)
%   spc      spectrum (both absorption and dispersion)
%
% For this deconvolution to work, the rapid-scan signals of the up- and
% the down-sweeps should not overlap, and the data should not be saturated.
%
% Additionally, the input signal needs to be correctly quadrature phased, time
% shifted to align with a cosine field modulation, and background corrected.
